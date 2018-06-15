extern crate image;
extern crate nalgebra as na;

use image::{Pixel, RgbImage};
use na::Vector3;

struct Ray {
    o: Vector3<f64>,
    d: Vector3<f64>,
}

struct Sphere {
    p: Vector3<f64>,
    r: f64,
}

impl Sphere {
    fn intersect(&self, ray: &Ray, tmin: f64, tmax: f64) -> Option<Hit> {
        let op = self.p - ray.o;
        let b = Vector3::dot(&op, &ray.d);
        let det = b * b - Vector3::dot(&op, &op) + self.r * self.r;
        if det < 0.0 {
            return None;
        }
        let t1 = b - det.sqrt();
        if tmin < t1 && t1 < tmax {
            return Some(Hit {
                t: t1,
                sphere: Some(self),
            });
        }
        let t2 = b + det.sqrt();
        if tmin < t2 && t2 < tmax {
            return Some(Hit {
                t: t2,
                sphere: Some(self),
            });
        }
        None
    }
}

struct Hit<'a> {
    t: f64,
    sphere: Option<&'a Sphere>,
}

struct Scene {
    spheres: Vec<Sphere>,
}

impl Scene {
    fn intersect(&self, ray: &Ray, tmin: f64, tmax: f64) -> Option<SurfInfo> {
        let mut t = tmax;
        let mut sphere: Option<&Sphere> = None;
        for s in &self.spheres {
            let h = s.intersect(ray, tmin, t);
            match h {
                Some(val) => {
                    t = val.t;
                    sphere = val.sphere;
                }
                None => continue,
            }
        }
        match sphere {
            Some(s) => {
                let p = ray.o + ray.d * t;
                let n = (p - s.p) / s.r;
                Some(SurfInfo { t, p, n, sphere })
            }
            None => None,
        }
    }
}

struct SurfInfo<'a> {
    t: f64,
    p: Vector3<f64>,
    n: Vector3<f64>,
    sphere: Option<&'a Sphere>,
}

fn main() {
    const W: u32 = 1200;
    const H: u32 = 800;

    let spheres = vec![Sphere {
        p: Vector3::new(0.0, 0.0, 0.0),
        r: 1.0,
    }];
    let scene = Scene { spheres };

    let mut img = RgbImage::new(W, H);
    for x in 0..W {
        for y in 0..H {
            let w = W as f64;
            let h = H as f64;
            let ray = Ray {
                o: Vector3::new(2.0 * x as f64 / w - 1.0, 2.0 * y as f64 / h - 1.0, 5.0),
                d: Vector3::new(0.0, 0.0, -1.0),
            };

            let h = scene.intersect(&ray, 0.0, 1e10);
            match h {
                Some(_) => img.put_pixel(x, y, Pixel::from_channels(255, 0, 255, 255)),
                None => img.put_pixel(x, y, Pixel::from_channels(0, 0, 0, 0)),
            }
        }
    }
    img.save("sample.png").unwrap();
}
