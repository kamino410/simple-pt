#include <cmath>
#include <string>

// vector like GLSL
struct Vec3 {
  double x, y, z;
  Vec3(const double x = 0, const double y = 0, const double z = 0) : x(x), y(y), z(z) {}
  inline Vec3 operator+(const Vec3 &b) const { return Vec3(x + b.x, y + b.y, z + b.z); }
  inline Vec3 operator-(const Vec3 &b) const { return Vec3(x - b.x, y - b.y, z - b.z); }
  inline Vec3 operator*(double c) const { return Vec3(c * x, c * y, c * z); }
  inline Vec3 operator*(const Vec3 &b) const { return Vec3(x * b.x, y * b.y, z * b.z); }
  inline Vec3 operator/(double c) const { return Vec3(x / c, y / c, z / c); }
  inline std::string to_string() const {
    return "(" + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ")";
  }
};
inline Vec3 operator*(const double c, const Vec3 &v) { return v * c; }
inline Vec3 operator+(const Vec3 &v) { return v; }
inline Vec3 operator-(const Vec3 &v) { return Vec3(-v.x, -v.y, -v.z); }
inline const double dot(const Vec3 &a, const Vec3 &b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline const Vec3 cross(const Vec3 &a, const Vec3 &b) {
  return Vec3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}
inline const double length(const Vec3 &v) { return sqrt(dot(v, v)); }
inline const Vec3 normalize(const Vec3 &v) { return v / length(v); }
