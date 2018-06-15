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
        if det < 0f64 {
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
    fn intersect(&self, ray: &Ray, tmin: f64, tmax: f64) -> Option<Hit> {
        let mut tmax = tmax;
        let mut res: Option<Hit> = None;
        for s in &self.spheres {
            let h = s.intersect(ray, tmin, tmax);
            match h {
                Some(val) => {
                    tmax = val.t;
                    res = Some(val);
                }
                None => continue,
            }
        }
        res
    }
}

fn main() {
    const W: u32 = 1200;
    const H: u32 = 800;

    let spheres = vec![Sphere {
        p: Vector3::new(0f64, 0f64, 0f64),
        r: 1f64,
    }];
    let scene = Scene { spheres };

    let mut img = RgbImage::new(W, H);
    for x in 0..W {
        for y in 0..H {
            let w = W as f64;
            let h = H as f64;
            let ray = Ray {
                o: Vector3::new(2f64 * x as f64 / w - 1f64, 2f64 * y as f64 / h - 1f64, 5f64),
                d: Vector3::new(0f64, 0f64, -1f64),
            };

            let h = scene.intersect(&ray, 0f64, 1e10);
            match h {
                Some(_) => img.put_pixel(x, y, Pixel::from_channels(255, 0, 255, 255)),
                None => img.put_pixel(x, y, Pixel::from_channels(0, 0, 0, 0)),
            }
        }
    }
    img.save("sample.png").unwrap();
}
