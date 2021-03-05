use std::fs;
use std::io::Write;
use derive_more::{Add,Sub, Mul, Div};

// Define struct for all time-evolving quantities
#[derive(Debug, Add, Sub, Mul, Div)]
struct Variables {
    x:  f64,
    y:  f64,
    vx: f64,
    vy: f64,
}

// Main RK4 algorithm
fn rk4(fxx: &dyn Fn(f64, f64, f64, f64, f64) -> f64, 
       fyy: &dyn Fn(f64, f64, f64, f64, f64) -> f64,
       fvx: &dyn Fn(f64, f64, f64, f64, f64) -> f64,
       fvy: &dyn Fn(f64, f64, f64, f64, f64) -> f64,
       t: f64, x: f64, y: f64, vx: f64, vy: f64, dt: f64) -> Variables {

    let kx1 =  fxx(t, x, y, vx, vy);
    let ky1 =  fyy(t, x, y, vx, vy);
    let kvx1 = fvx(t, x, y, vx, vy);
    let kvy1 = fvy(t, x, y, vx, vy);

    let kx2 =  fxx(t + 0.5 * dt, x + dt * 0.5 * kx1, y + dt * 0.5 * ky1, vx + dt * 0.5 * kvx1, vy + dt * 0.5 * kvy1);
    let ky2 =  fyy(t + 0.5 * dt, x + dt * 0.5 * kx1, y + dt * 0.5 * ky1, vx + dt * 0.5 * kvx1, vy + dt * 0.5 * kvy1);
    let kvx2 = fvx(t + 0.5 * dt, x + dt * 0.5 * kx1, y + dt * 0.5 * ky1, vx + dt * 0.5 * kvx1, vy + dt * 0.5 * kvy1);
    let kvy2 = fvy(t + 0.5 * dt, x + dt * 0.5 * kx1, y + dt * 0.5 * ky1, vx + dt * 0.5 * kvx1, vy + dt * 0.5 * kvy1);

    let kx3 =  fxx(t + 0.5 * dt, x + dt * 0.5 * kx2, y + dt * 0.5 * ky2, vx + dt * 0.5 * kvx2, vy + dt * 0.5 * kvy2);
    let ky3 =  fyy(t + 0.5 * dt, x + dt * 0.5 * kx2, y + dt * 0.5 * ky2, vx + dt * 0.5 * kvx2, vy + dt * 0.5 * kvy2);
    let kvx3 = fvx(t + 0.5 * dt, x + dt * 0.5 * kx2, y + dt * 0.5 * ky2, vx + dt * 0.5 * kvx2, vy + dt * 0.5 * kvy2);
    let kvy3 = fvy(t + 0.5 * dt, x + dt * 0.5 * kx2, y + dt * 0.5 * ky2, vx + dt * 0.5 * kvx2, vy + dt * 0.5 * kvy2);

    let kx4 =  fxx(t + dt, x + dt * kx3, y + dt * ky3, vx + dt * kvx3, vy + dt * kvy3);
    let ky4 =  fyy(t + dt, x + dt * kx3, y + dt * ky3, vx + dt * kvx3, vy + dt * kvy3);
    let kvx4 = fvx(t + dt, x + dt * kx3, y + dt * ky3, vx + dt * kvx3, vy + dt * kvy3);
    let kvy4 = fvy(t + dt, x + dt * kx3, y + dt * ky3, vx + dt * kvx3, vy + dt * kvy3);

    return Variables {
        x:  x  + (1.0 / 6.0) * dt * (kx1 + 2.0 * kx2 + 2.0 * kx3 + kx4),
        y:  y  + (1.0 / 6.0) * dt * (ky1 + 2.0 * ky2 + 2.0 * ky3 + ky4),
        vx: vx + (1.0 / 6.0) * dt * (kvx1 + 2.0 * kvx2 + 2.0 * kvx3 + kvx4),
        vy: vy + (1.0 / 6.0) * dt * (kvy1 + 2.0 * kvy2 + 2.0 * kvy3 + kvy4),
    };
}

// Functions of evolving quantities pulled from pertinent differential equations
fn f_x(_t: f64, _x: f64, _y: f64, vx: f64, _vy: f64) -> f64 {
    vx
}
fn f_y(_t: f64, _x: f64, _y: f64, _vx: f64, vy: f64) -> f64 {
    vy
}
fn f_vx(_t: f64, _x: f64, _y: f64, _vx: f64, _vy: f64) -> f64 {
    0.0
}
fn f_vy(_t: f64, _x: f64, _y: f64, _vx: f64, _vy: f64) -> f64 {
    -9.8
}

fn main() {
    let mut t  = 0.0;
    let mut x  = 0.0;
    let mut y  = 0.0;
    let mut vx = 4.0;
    let mut vy = 4.0;
    let step   = 0.001;
    let mut counter = 0;
    let mut vecx: Vec<f64> = Vec::new();
    let mut vecy: Vec<f64> = Vec::new();
 
    while y >= 0.0 {
 
        let va = rk4(&f_x, &f_y, &f_vx, &f_vy, t, x, y, vx, vy, step);

        x  = va.x;
        y  = va.y;
        vx = va.vx;
        vy = va.vy;

        vecx.push(x);
        vecy.push(y);

        t += step;

        counter += 1;

    }

    // Write data file
    let file = fs::File::create("solution.dat").unwrap();
    for i in 0..(counter-1) {

        writeln!(&file, "{:.6} {:.6}", vecx[i as usize], vecy[i as usize]).unwrap();

    }

}
