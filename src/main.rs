extern crate image;
extern crate nalgebra;
extern crate rand;
extern crate rayon;

use image::RgbImage;
use nalgebra::Vector3;
use rand::prelude::*;
use rayon::prelude::*;
use std::cmp::{max, min};

struct Ray {
    o: Vector3<f64>,
    d: Vector3<f64>,
}

struct Sphere {
    p: Vector3<f64>,
    r: f64,
    rf: Vector3<f64>,
    le: Vector3<f64>,
    spc: bool,
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
                Some(SurfInfo { p, n, sphere: s })
            }
            None => None,
        }
    }
}

struct SurfInfo<'a> {
    p: Vector3<f64>,
    n: Vector3<f64>,
    sphere: &'a Sphere,
}

fn tangent_space(n: &Vector3<f64>) -> (Vector3<f64>, Vector3<f64>) {
    let s = if n.z > 0.0 { 1.0 } else { -1.0 };
    let a = -1.0 / (s + n.z);
    let b = n.x * n.y * a;
    (
        Vector3::new(1.0 + s * n.x * n.x * a, s * b, -s * n.x),
        Vector3::new(b, s + n.y * n.y * a, -n.y),
    )
}

fn tonemap(v: f64) -> u8 {
    min(max((v.powf(1.0 / 2.2) * 255.0) as u32, 0), 255) as u8
}

fn main() {
    const W: u32 = 1200;
    const H: u32 = 800;
    const SPP: u32 = 1000;
    const MAXDEPTH: u32 = 10;

    let eye = Vector3::new(50.0, 52.0, 295.6);
    let center = eye + Vector3::new(0.0, -0.042612, -1.0);
    let up = Vector3::new(0.0, 1.0, 0.0);
    let fov = 30.0 * 3.141593 / 180.0;
    let aspect = W as f64 / H as f64;
    let we = Vector3::normalize(&(eye - center));
    let ue = Vector3::normalize(&Vector3::cross(&up, &we));
    let ve = Vector3::cross(&we, &ue);

    let spheres = vec![
        Sphere {
            p: Vector3::new(1e5 + 1.0, 40.8, 81.6),
            r: 1e5,
            rf: Vector3::new(0.75, 0.25, 0.25),
            le: Vector3::new(0.0, 0.0, 0.0),
            spc: false,
        },
        Sphere {
            p: Vector3::new(-1e5 + 99.0, 40.8, 81.6),
            r: 1e5,
            rf: Vector3::new(0.25, 0.25, 0.75),
            le: Vector3::new(0.0, 0.0, 0.0),
            spc: false,
        },
        Sphere {
            p: Vector3::new(50.0, 40.8, 1e5),
            r: 1e5,
            rf: Vector3::new(0.75, 0.75, 0.75),
            le: Vector3::new(0.0, 0.0, 0.0),
            spc: false,
        },
        Sphere {
            p: Vector3::new(50.0, 1e5, 81.6),
            r: 1e5,
            rf: Vector3::new(0.75, 0.75, 0.75),
            le: Vector3::new(0.0, 0.0, 0.0),
            spc: false,
        },
        Sphere {
            p: Vector3::new(50.0, -1e5 + 81.6, 81.6),
            r: 1e5,
            rf: Vector3::new(0.75, 0.75, 0.75),
            le: Vector3::new(0.0, 0.0, 0.0),
            spc: false,
        },
        Sphere {
            p: Vector3::new(27.0, 16.5, 47.0),
            r: 16.5,
            rf: Vector3::new(0.999, 0.999, 0.999),
            le: Vector3::new(0.0, 0.0, 0.0),
            spc: true,
        },
        Sphere {
            p: Vector3::new(73.0, 16.5, 78.0),
            r: 16.5,
            rf: Vector3::new(0.999, 0.999, 0.999),
            le: Vector3::new(0.0, 0.0, 0.0),
            spc: false,
        },
        Sphere {
            p: Vector3::new(50.0, 681.6 - 0.27, 81.6),
            r: 600.0,
            rf: Vector3::new(0.0, 0.0, 0.0),
            le: Vector3::new(12.0, 12.0, 12.0),
            spc: false,
        },
    ];
    let scene = Scene { spheres };

    let mut img = RgbImage::new(W, H);
    let tf = (fov * 0.5f64).tan();

    img.enumerate_pixels_mut()
        .collect::<Vec<(u32, u32, &mut image::Rgb<u8>)>>()
        .par_iter_mut()
        .for_each(|(x, y, pixel)| {
            let mut rng = rand::thread_rng();
            let mut il = Vector3::new(0.0, 0.0, 0.0);
            for _ in 0..SPP {
                let rpx = 2.0 * (*x as f64 + rng.gen::<f64>()) / W as f64 - 1.0;
                let rpy = 2.0 * ((H - *y - 1) as f64 + rng.gen::<f64>()) / H as f64 - 1.0;
                let w = Vector3::normalize(&Vector3::new(aspect * tf * rpx, tf * rpy, -1.0));
                let mut ray = Ray {
                    o: eye,
                    d: ue * w.x + ve * w.y + we * w.z,
                };

                let mut l = Vector3::new(0.0, 0.0, 0.0);
                let mut th = Vector3::new(1.0, 1.0, 1.0);
                for _ in 0..MAXDEPTH {
                    let h = scene.intersect(&ray, 1e-4, 1e10);
                    if let Some(h) = h {
                        l.x += th.x * h.sphere.le.x;
                        l.y += th.y * h.sphere.le.y;
                        l.z += th.z * h.sphere.le.z;

                        ray.o = h.p;
                        ray.d = {
                            let n = if Vector3::dot(&h.n, &-ray.d) > 0.0 {
                                h.n
                            } else {
                                -h.n
                            };
                            if h.sphere.spc {
                                ray.d + 2.0 * Vector3::dot(&n, &-ray.d) * n
                            } else {
                                let (u, v) = tangent_space(&n);
                                let d = {
                                    let r = rng.gen::<f64>().sqrt();
                                    let t = 2.0 * 3.141593 * rng.gen::<f64>();
                                    let x = r * t.cos();
                                    let y = r * t.sin();
                                    let z = 1.0 - x * x - y * y;
                                    let z = if z < 0.0 {
                                        0.0
                                    } else {
                                        z.sqrt()
                                    };
                                    Vector3::new(x, y, z)
                                };
                                u * d.x + v * d.y + n * d.z
                            }
                        };
                        th.x *= h.sphere.rf.x;
                        th.y *= h.sphere.rf.y;
                        th.z *= h.sphere.rf.z;
                        if th.x < 0.0 && th.y < 0.0 && th.z < 0.0 {
                            break;
                        }
                    } else {
                        break;
                    }
                }
                il += l / SPP as f64;
            }
            pixel[0] = tonemap(il.x);
            pixel[1] = tonemap(il.y);
            pixel[2] = tonemap(il.z);
        });

    img.save("sample.png").unwrap();
}
