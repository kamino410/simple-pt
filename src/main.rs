extern crate image;
extern crate nalgebra as na;

use image::{Pixel, RgbImage};
use na::{Perspective3, Vector3};
use std::cmp::{max, min};

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

fn tonemap(v: f64) -> u8 {
    min(max((v.powf(1.0 / 2.2) * 255.0) as u32, 0), 255) as u8
}

fn main() {
    const W: u32 = 1200;
    const H: u32 = 800;

    let eye = Vector3::new(50.0, 52.0, 295.6);
    let center = eye + Vector3::new(0.0, -0.042612, -1.0);
    let up = Vector3::new(0.0, 1.0, 0.0);
    let fov = 30.0 * 3.141593 / 180.0;
    let aspect = W as f64 / H as f64;
    let wE = Vector3::normalize(&(eye - center));
    let uE = Vector3::normalize(&Vector3::cross(&up, &wE));
    let vE = Vector3::cross(&wE, &uE);

    let spheres = vec![
        Sphere {
            p: Vector3::new(1e5 + 1.0, 40.8, 81.6),
            r: 1e5,
        },
        Sphere {
            p: Vector3::new(-1e5 + 99.0, 40.8, 81.6),
            r: 1e5,
        },
        Sphere {
            p: Vector3::new(50.0, 40.8, 1e5),
            r: 1e5,
        },
        Sphere {
            p: Vector3::new(50.0, 1e5, 81.6),
            r: 1e5,
        },
        Sphere {
            p: Vector3::new(50.0, -1e5 + 81.6, 81.6),
            r: 1e5,
        },
        Sphere {
            p: Vector3::new(27.0, 16.5, 47.0),
            r: 16.5,
        },
        Sphere {
            p: Vector3::new(73.0, 16.5, 78.0),
            r: 16.5,
        },
        Sphere {
            p: Vector3::new(50.0, 681.6 - 0.27, 81.6),
            r: 600.0,
        },
    ];
    let scene = Scene { spheres };

    let mut img = RgbImage::new(W, H);
    for x in 0..W {
        for y in 0..H {
            let tf = (fov * 0.5f64).tan();
            let rpx = 2.0 * x as f64 / W as f64 - 1.0;
            let rpy = 2.0 * y as f64 / H as f64 - 1.0;
            let w = Vector3::normalize(&Vector3::new(aspect * tf * rpx, tf * rpy, -1.0));
            let ray = Ray {
                o: eye,
                d: uE * w.x + vE * w.y + wE * w.z,
            };

            let h = scene.intersect(&ray, 0.0, 1e10);
            match h {
                Some(h) => img.put_pixel(
                    x,
                    H - y - 1,
                    Pixel::from_channels(
                        tonemap(h.n.x.abs()),
                        tonemap(h.n.y.abs()),
                        tonemap(h.n.z.abs()),
                        255,
                    ),
                ),
                None => img.put_pixel(x, H - y - 1, Pixel::from_channels(0, 0, 0, 0)),
            }
        }
    }
    img.save("sample.png").unwrap();
}
