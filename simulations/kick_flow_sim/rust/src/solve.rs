use differential_equations::{prelude::*};
use nalgebra::SVector;

/* 
Functions herein are hard-coded for 50 possible alpha values
Replace all instances of 50 when changing that.
*/

// Analytical solution for a single strain, used to determine when a patch is at equilibrium
pub fn get_biomass(alpha: f64) -> f64 {
    let s_fe = 5f64;
    let s_apo = 1f64;
    let dr = 1f64;
    let dm = 0.1f64;
    let kf = 1f64;
    let gamma = 3f64;
    let vh = 10f64;
    let kh = 1f64;
    let epsilon = 1f64;

    let g = (1.0 - alpha) * gamma * vh;
    let a = kf * dm.powi(2) * vh.powi(2) / g.powi(2)
                - alpha * epsilon * kf * dm * vh / g;
    let b = 2.0 * kf * dm.powi(2) * vh * kh * dr / g.powi(2)
                - kh * alpha * epsilon * kf * dm * dr / g
                - dm * (dr.powi(2) + kf * (s_apo + s_fe)) * vh / g
                + alpha * epsilon * kf * s_fe;
    let c = kf * dm.powi(2) * kh.powi(2) * dr.powi(2) / g.powi(2)
                - dm * kh * dr * (dr.powi(2) + kf * (s_apo + s_fe)) / g
                + kf * s_apo * s_fe;

    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 { return 0.0; }

    let sqrt_disc = disc.sqrt();
    let m1 = (-b + sqrt_disc) / (2.0 * a);
    let m2 = (-b - sqrt_disc) / (2.0 * a);

    match (m1 > 0.0, m2 > 0.0) {
        (true, false) => m1,
        (false, true) => m2,
        (true, true)  => m1.min(m2),
        _             => 0.0,
    }
}


pub type PatchState = SVector<f64, 50>;

struct SiderophoreModel {
    alphas: Vec<f64>,
    equilibria: Vec<f64>
}

impl ODE<f64, PatchState> for SiderophoreModel {
    fn diff(&self, _t: f64, y: &PatchState, dy: &mut PatchState) {
        let s_fe: f64 = 5.0;
        let s_apo: f64 = 1.0;
        let d_r: f64 = 1.0;
        let d_m: f64 = 0.1;
        let k_f: f64 = 1.0;
        let gamma: f64 = 3.0;
        let v_h: f64 = 10.0;
        let k_h: f64 = 1.0;
        let epsilon: f64 = 1.0;

        let size = 50;

        // Calculate total biomass (sum of M_i)
        let total_biomass: f64 = y.sum();

        // Calculate sum of alpha_i * M_i
        let alpha_biomass_sum: f64 = self.alphas.iter()
            .zip(y.iter())
            .map(|(alpha, m)| alpha * m)
            .sum();

        // Calculate A
        let a = s_apo + epsilon * alpha_biomass_sum;

        // Calculate F*
        let kf_term = k_f * (s_fe + a) + d_r.powi(2);
        let sqrt_term = (kf_term.powi(2) - 4.0 * k_f.powi(2) * s_fe * a).sqrt();
        let f_star = (kf_term - sqrt_term) / (2.0 * k_f);

        // Denominator for the growth rate equation
        let denominator = v_h * total_biomass + k_h * d_r;

        // Compute derivatives for each strain
        for i in 0..size {
            let growth_rate = ((1.0 - self.alphas[i]) * gamma * v_h * f_star) / denominator;
            dy[i] = (growth_rate - d_m) * y[i];
        }
    }

    fn event(&self, _t: f64, y: &PatchState) -> ControlFlag<f64, PatchState> {
        // Determine who's alive (biomass > 1e-3)
        let living: Vec<bool> = y.iter()
                                 .map(|x| *x > 1.0e-3)
                                 .collect();

        let living_nr = living.iter().filter(|b| **b).count();
        
        // Terminate if the patch is empty
        if living_nr == 0 {
            return ControlFlag::Terminate("Patch is empty".to_string());
        }
        
        // Continue if multiple strains are still competing
        if living_nr > 1 {
            return ControlFlag::Continue;
        }
        
        // Check if the single living strain is near equilibrium
        let ind = living.iter().position(|x| *x).unwrap();
        if (y[ind] - self.equilibria[ind]).abs() < 1.0e-3 {
            ControlFlag::Terminate("Equilibrium was reached".to_string())
        } else {
            ControlFlag::Continue
        }
    }
}

pub fn solve(alphas: Vec<f64>, y0: PatchState, tend: f64) -> (PatchState, bool) {
    let equilibria: Vec<f64> = alphas.iter().map(|x| get_biomass(*x)).collect();
    
    // Prepare model
    let ode = SiderophoreModel { alphas, equilibria };
    let model = ODEProblem::new(ode, 0.0, tend, y0);
    
    // Solve
    let mut method = DOPRI5::new().rtol(1e-6).atol(1e-9);
    
    // Set initial step size for small end times
    let h0 = if tend < 1e-4 { 
        tend 
    } else { 
        tend.min(0.001) // Use the smaller of 0.001 or the interval
    };
    
    method = method.h0(h0);
    
    match model.solve(&mut method) {
        Ok(solution) => {
            // Check if termination condition (equilibrium) was met
            if let Status::Interrupted(ref _reason) = solution.status {
                (*solution.last().unwrap().1, true)
            } else {
                (*solution.last().unwrap().1, false)
            }
        }
        Err(e) => panic!("Error: {:?}", e),
    }
}