use differential_equations::{prelude::*};
use nalgebra::SVector;

/* 
Functions herein are hard-coded for 50 possible alpha values
plus three resource state variables. Replace all instances of 50 and 53 when changing that.
*/

// Analytical solution for a single strain, used to determine when a patch is at equilibrium
pub fn get_biomass(alpha: f64) -> f64 {
    let (fein, d,   kf,   gamma,  vh,   kh,     epsilon, rho) = (
        1f64, 1f64, 1f64, 10f64,  1f64, 0.1f64, 2f64,    1f64);   
    // stable positive solution from mathematica 
    let sol = (gamma * (-1.0 + alpha) * (d.powi(4) + d.powi(3) * (-kf * kh * (-2.0 + rho) + 
           vh * gamma * (-1.0 + alpha)) + fein * kf * vh * gamma.powi(2) * epsilon * 
        (-1.0 + alpha).powi(2) * alpha + d * fein * kf * gamma * (-1.0 + alpha) * (vh - 
           vh * rho + epsilon * alpha) + d.powi(2) * kf * (fein - fein * rho + 
           kh * gamma * epsilon * (-1.0 + alpha) * alpha) + (4.0 * d.powi(5) * kf * kh * (d + 
              vh * gamma * (-1.0 + alpha)) * (d * rho - gamma * epsilon * (-1.0 + alpha) * 
              alpha) + (d.powi(4) + d.powi(3) * (-kf * kh * rho + vh * gamma * (-1.0 + alpha)) +
              fein * kf * vh * gamma.powi(2) * epsilon * (-1.0 + alpha).powi(2) * 
              alpha + d * fein * kf * gamma * (-1.0 + alpha) * (vh - vh * rho + 
                epsilon * alpha) + d.powi(2) * kf * (fein - fein * rho + 
                kh * gamma * epsilon * (-1.0 + alpha) * alpha)).powi(2)).sqrt())) / 
                (2.0 * d * kf * (d + vh * gamma * (-1.0 + alpha)) * (d * (-1.0 + rho) - 
        gamma * epsilon * (-1.0 + alpha) * alpha));
    if sol.is_nan() || sol.is_sign_negative() {
        return 0f64;
    }
    return sol;        
}


pub type PatchState = SVector<f64, 53>;

struct SiderophoreModel {
    alphas: Vec<f64>,
    equilibria: Vec<f64>
}

impl ODE<f64, PatchState> for SiderophoreModel {
    fn diff(&self, _t: f64, y: &PatchState, dy: &mut PatchState) {
        let d = 1.0;
        let vh = 1.0;
        let gamma = 10.0;
        let epsilon = 2.0;
        let kh = 0.1;
        let kf = 1.0;
        let rho = 1.0;
        let sfe = 1.0;

        
        let size = 50;

        let uptake = vh * y[size + 1] / (kh + y[size + 1]);
        let chelation = kf * y[size + 2] * y[size];

        let biomass = y.rows(0, 50);

        for i in 0..(size) {
            dy[i] = (1.0 - self.alphas[i]) * gamma * uptake * y[i] - d * y[i];
        }

        // Apo
        dy[size] = self.alphas.iter()
                              .zip(biomass)
                              .map(|(alpha, mass)| 
                                   mass * (epsilon * alpha + rho * uptake))
                              .sum::<f64>() 
                    - chelation - d * y[size];
        // Holo 
        dy[size + 1] = chelation - uptake * biomass.sum() - d * y[size + 1];
        // Fe 
        dy[size + 2] = sfe - chelation - d * y[size + 2];
    }

    // Stop when the patch reaches the equilibrium biomass or death
    fn event(&self, _t: f64, y: &PatchState) -> ControlFlag {
        // Determine who's alive (biomass > 1e-3)
        let living: Vec<bool> = y.rows_range(0..50)
                                 .iter()
                                 .map(|x| *x > 1.0e-3)
                                 .collect();

        let living_nr = living.iter().filter(|b| **b).count();
        // If no one is alive, stop
        if living_nr == 0 {
            return ControlFlag::Terminate("Patch is empty".to_string())
        }
        // If multiple strains exist, we must not be at equilibrium
        if living_nr > 1 {
            return ControlFlag::Continue
        }
        // Check if the one living strain is near equilibrium
        let ind = living.iter().position(|x| *x ).unwrap();
        if (y[ind] - self.equilibria[ind]).abs() < 1.0e-3 {
            ControlFlag::Terminate("Equilibrium was reached".to_string())
        }
        else {
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
    // Help it set an initial step size for very small end times
    if tend < 0.5 {method = method.h0(0.0001)}
    match model
        .solve(&mut method) // Solve the system
    {
        // Return final state
        Ok(solution) => {
            // Check if the solver stopped due to the event command
            // This indicates the patch is at equilibrium
            if let Status::Interrupted(ref _reason) = solution.status {
                return (*solution.last().unwrap().1, true);
            }
            else {
                return (*solution.last().unwrap().1, false);
            }
        }
        Err(e) => panic!("Error: {:?}", e),
    };
}

