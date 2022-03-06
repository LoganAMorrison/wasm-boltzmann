use cyphus_diffeq::ode::*;
use cyphus_diffeq::radau::*;
use haliax_thermal_functions::standard_model::{sm_energy_density, sm_sqrt_gstar};
use haliax_thermal_functions::thermal_functions::neq;
use ndarray::prelude::*;
use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

const MPLANK: f64 = 1.220910e19;

// When the `wee_alloc` feature is enabled, use `wee_alloc` as the global
// allocator.
#[cfg(feature = "wee_alloc")]
#[global_allocator]
static ALLOC: wee_alloc::WeeAlloc = wee_alloc::WeeAlloc::INIT;

#[wasm_bindgen]
#[derive(Serialize, Deserialize)]
pub struct Config {
    pub x_start: f64,
    pub x_end: f64,
    pub mass: f64,
    pub n: i32,
    pub sigma: f64,
}

#[wasm_bindgen]
impl Config {
    pub fn new() -> Config {
        Config {
            x_start: 1.0,
            x_end: 500_f64,
            mass: 100.0,
            n: 0,
            sigma: 1e-9,
        }
    }

    fn verify(&self) -> Result<(), &'static str> {
        if self.x_start < 0.0 {
            return Err("Invalid x-start.");
        }
        if self.x_end < 0.0 {
            return Err("Invalid x-end.");
        }
        if self.sigma < 0.0 {
            return Err("Invalid sigma.");
        }
        if self.mass < 0.0 {
            return Err("Invalid mass.");
        }
        if self.n < 0 {
            return Err("Invalid n.");
        }
        if self.x_end < self.x_start {
            return Err("Invalid x-end or x-start.");
        }
        Ok(())
    }
}

struct Model {
    mass: f64,
    n: i32,
    sigma: f64,
}

#[wasm_bindgen]
#[derive(Serialize, Deserialize)]
pub struct Solution {
    logxs: Vec<f64>,
    ws: Vec<f64>,
    status: String,
}

#[wasm_bindgen]
impl Solution {
    pub fn size(&self) -> usize {
        return self.logxs.len();
    }
    pub fn get_logx(&self, i: usize) -> f64 {
        return self.logxs[i];
    }
    pub fn get_w(&self, i: usize) -> f64 {
        return self.ws[i];
    }
    pub fn get_status(&self) -> String {
        return self.status.clone();
    }
}

fn solve_boltzmann(config: Config) -> Solution {
    match config.verify() {
        Err(e) => {
            return Solution {
                logxs: vec![],
                ws: vec![],
                status: format!("Errored with: {}", e).to_owned(),
            };
        }
        _ => {}
    };

    let mass = config.mass;
    let model = Model {
        mass,
        n: config.n,
        sigma: config.sigma,
    };

    let logxspan = (config.x_start.ln(), config.x_end.ln());
    let pre = -(std::f64::consts::PI / 45.0).sqrt() * model.sigma * model.mass * MPLANK;
    let t0 = mass * (-logxspan.0).exp();
    let w0 = neq(t0, model.mass, 2.0, 2).ln() - sm_energy_density(t0).ln();
    let winit = array![w0];

    let dwdlx = |mut dw: ArrayViewMut1<f64>, w: ArrayView1<f64>, logx: f64, model: &Model| {
        let x = logx.exp();
        let temp = model.mass / x;
        let pf = pre * sm_sqrt_gstar(temp) / x.powi(model.n + 1);
        let weq = neq(temp, model.mass, 2.0, 2).ln() - sm_energy_density(temp).ln();
        dw[0] = pf * (w[0].exp() - (2.0 * weq - w[0]).exp());
    };
    let jac = |mut j: ArrayViewMut2<f64>, w: ArrayView1<f64>, logx: f64, model: &Model| {
        let x = logx.exp();
        let temp = model.mass / x;
        let pf = pre * sm_sqrt_gstar(temp) / x.powi(model.n + 1);
        let weq = neq(temp, model.mass, 2.0, 2).ln() - sm_energy_density(temp).ln();
        j[[0, 0]] = pf * (w[0].exp() + (2.0 * weq - w[0]).exp());
    };
    let mut integrator = OdeIntegratorBuilder::default(&dwdlx, winit, logxspan, Radau5, model)
        .dfdu(&jac)
        .reltol(1e-7)
        .abstol(1e-7)
        .build();
    integrator.integrate();

    let mut status = String::new();

    match integrator.sol.retcode {
        OdeRetCode::Success => {
            status = "Success".to_owned();
        }
        OdeRetCode::Stiff => {
            status = "Too stiff".to_owned();
        }
        OdeRetCode::Failure => {
            status = "Failure".to_owned();
        }
        OdeRetCode::MaxIters => {
            status = "Maximum iterations exceeded".to_owned();
        }
        OdeRetCode::DtLessThanMin => {
            status = "Step size too small".to_owned();
        }
        OdeRetCode::SingularMatrix => {
            status = "Singular jacobian".to_owned();
        }
        _ => {}
    }

    let mut sol = Solution {
        logxs: vec![],
        ws: vec![],
        status,
    };
    sol.logxs.reserve(integrator.sol.ts.len());
    sol.ws.reserve(integrator.sol.ts.len());
    for (t, u) in (&mut integrator.sol).into_iter() {
        sol.logxs.push(*t);
        sol.ws.push(u[0]);
        println!("{}, {}", t, u);
    }
    return sol;
}

#[wasm_bindgen]
pub fn solve(config: Config) -> Solution {
    // let config = config
    //     .into_serde::<Config>()
    //     .map_err(|_| JsValue::from_str("Failed to deserialize configuration."));
    let solution = solve_boltzmann(config);
    solution
}
