#include <cmath>
#include <algorithm>
#include <random>
#include <iostream>

#include "vec3.hpp"

const double EPS = 1e-5;
const double INF = 1e20;
const int LEASTDEPTH = 5;
const int MAXDEPTH = 20;

double rand01() { return rand() / double(RAND_MAX); }

struct Ray {
  Vec3 org, dir;
  Ray(const Vec3 &org, const Vec3 &dir) : org(org), dir(dir) {}
};

enum ReflectionType {
  DIFFUSE,
  SPECULAR,
  REFRACTION,
};

struct Sphere {
  Vec3 pos, em, col;  // position, emission, color
  double r;           // radius
  ReflectionType ref_type;
  Sphere(const double r, const Vec3 &pos, const Vec3 &em, const Vec3 &col,
         const ReflectionType ref_type)
      : r(r), pos(pos), em(em), col(col), ref_type(ref_type) {}
  // returns distance to this sphere, or 0 if nohit
  double intersect(const Ray &ray) const {
    Vec3 op = pos - ray.org;
    const double b = dot(op, ray.dir);
    const double det = b * b - dot(op, op) + r * r;
    if (det >= 0.0) {
      const double sqrt_det = sqrt(det);
      const double t1 = b - sqrt_det;
      if (t1 > EPS) return t1;
      const double t2 = b + sqrt_det;
      if (t2 > EPS) return t2;
    }
    return 0.0;
  }
};

const int N = 9;
Sphere spheres[N] = {
    Sphere(1e5, Vec3(1e5 + 1, 40.8, 81.6), Vec3(), Vec3(.75, .25, .25), DIFFUSE),    // left
    Sphere(1e5, Vec3(-1e5 + 99, 40.8, 81.6), Vec3(), Vec3(.25, .25, .75), DIFFUSE),  // right
    Sphere(1e5, Vec3(50, 40.8, 1e5), Vec3(), Vec3(.75, .75, .75), DIFFUSE),          // back
    Sphere(1e5, Vec3(50, 40.8, -1e5 + 170), Vec3(), Vec3(), DIFFUSE),                // front
    Sphere(1e5, Vec3(50, 1e5, 81.6), Vec3(), Vec3(.75, .75, .75), DIFFUSE),          // bottom
    Sphere(1e5, Vec3(50, -1e5 + 81.6, 81.6), Vec3(), Vec3(.75, .75, .75), DIFFUSE),  // top
    Sphere(16.5, Vec3(27, 16.5, 47), Vec3(), Vec3(1, 1, 1) * .999, SPECULAR),        // mirror
    Sphere(16.5, Vec3(73, 16.5, 78), Vec3(), Vec3(1, 1, 1) * .999, REFRACTION),      // glass
    Sphere(600, Vec3(50, 681.6 - .27, 81.6), Vec3(12, 12, 12), Vec3(), DIFFUSE)};    // light

inline bool intersect_scene(const Ray &ray, double *t, int *id) {
  *t = INF;
  *id = -1;
  for (int i = 0; i < N; i++) {
    const double d = spheres[i].intersect(ray);
    if (d > 0.0 && d < *t) {
      *t = d;
      *id = i;
    }
  }
  return *t < INF;
}

void tangent_space(const Vec3 &n, Vec3 &u, Vec3 &v) {
  u = fabs(n.x) > 0.1 ? normalize(cross(Vec3(0, 1, 0), n)) : normalize(cross(Vec3(1, 0, 0), n));
  v = cross(n, u);
}

Vec3 radiance(const Ray &ray, const int depth) {
  double t;
  int id;
  if (!intersect_scene(ray, &t, &id)) return Vec3(0, 0, 0);  // background color

  const Sphere &obj = spheres[id];
  const Vec3 hit_pos = ray.org + t * ray.dir;
  const Vec3 norm2out = normalize(hit_pos - obj.pos);
  // ray toward inside : normal is directed toward outside
  // ray toward outside : normal is directed toward inside
  bool into = dot(norm2out, ray.dir) < 0.0;
  const Vec3 n = into ? norm2out : -norm2out;

  // Russian Roulette
  Vec3 ref_coef = obj.col;
  const double p1 = std::max(obj.col.x, std::max(obj.col.y, obj.col.z));
  if (depth > MAXDEPTH) return obj.em;
  if (depth > LEASTDEPTH) {
    if (rand01() > p1) return obj.em;
    ref_coef = ref_coef / p1;
  }
  switch (obj.ref_type) {
    case DIFFUSE: {
      Vec3 u, v;
      tangent_space(n, u, v);
      const double t = 2 * M_PI * rand01(), r2 = rand01();
      const double r = sqrt(r2), r_ = sqrt(1 - r2);
      Vec3 dir = normalize(u * cos(t) * r + v * sin(t) * r + n * r_);
      return obj.em + ref_coef * radiance(Ray(hit_pos, dir), depth + 1);
    }
    case SPECULAR: {
      Vec3 dir = ray.dir - 2 * n * dot(n, ray.dir);
      return obj.em + ref_coef * radiance(Ray(hit_pos, dir), depth + 1);
    }
    case REFRACTION: {
      const double n_air = 1.0, n_glass = 1.5;
      const double nnt = into ? n_air / n_glass : n_glass / n_air;
      const double ddn = dot(ray.dir, n);
      const double cos2t = 1 - nnt * nnt * (1 - ddn * ddn);
      const Vec3 tr_dir = ray.dir - 2 * n * dot(n, ray.dir);  // total reflection
      if (cos2t < 0.0) {
        return obj.em + ref_coef * radiance(Ray(hit_pos, tr_dir), depth + 1);
      } else {
        const Vec3 rf_dir = normalize(ray.dir * nnt - n * (ddn * nnt + sqrt(cos2t)));  // refraction
        const double a = n_glass - n_air, b = n_glass + n_air;
        const double c = 1 - (into ? -ddn : dot(rf_dir, norm2out));
        const double r0 = (a * a) / (b * b);
        const double re = r0 + (1 - r0) * pow(c, 5);
        const double tr = 1 - re;
        if (depth > 2) {
          const double p2 = 0.25 + 0.5 * re;
          if (rand01() < p2) {
            ref_coef = ref_coef / p2;
            return obj.em + ref_coef * radiance(Ray(hit_pos, tr_dir), depth + 1) * re;
          } else {
            ref_coef = ref_coef / (1 - p2);
            return obj.em + ref_coef * radiance(Ray(hit_pos, rf_dir), depth + 1) * tr;
          }
        } else {
          return obj.em + ref_coef * (radiance(Ray(hit_pos, tr_dir), depth + 1) * re +
                                      radiance(Ray(hit_pos, rf_dir), depth + 1) * tr);
        }
      }
    }
  }
}

inline double clamp(double x) { return x < 0 ? 0 : x > 1 ? 1 : x; }
inline int toInt(double x) { return int(pow(clamp(x), 1 / 2.2) * 255 + .5); }
int main() {
  int w = 1024, h = 768, samps = 100;
  Ray cam(Vec3(50, 52, 295.6), normalize(Vec3(0, -0.042612, -1)));
  Vec3 cx = Vec3(w * .5135 / h), cy = normalize(cross(cx, cam.dir)) * .5135, r,
       *img = new Vec3[w * h];
  for (int y = 0; y < h; y++) {
    fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps * 4, 100. * y / (h - 1));
    for (unsigned short x = 0; x < w; x++) {
      for (int sy = 0, i = (h - y - 1) * w + x; sy < 2; sy++) {
        for (int sx = 0; sx < 2; sx++) {
          r = Vec3();
          for (int s = 0; s < samps; s++) {
            const double r1 = 2 * rand01(), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
            const double r2 = 2 * rand01(), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
            const Vec3 d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
                           cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.dir;
            r = r + radiance(Ray(cam.org + d * 140, normalize(d)), 0) * (1. / samps);
          }
          img[i] = img[i] + Vec3(clamp(r.x), clamp(r.y), clamp(r.z)) * .25;
        }
      }
    }
  }
  FILE *f = fopen("image.ppm", "w");  // Write image to PPM file.
  fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
  for (int i = 0; i < w * h; i++)
    fprintf(f, "%d %d %d ", toInt(img[i].x), toInt(img[i].y), toInt(img[i].z));
}
