use rand::prelude::*;
use rand_distr::Normal;
use std::cmp::Ordering;

#[derive(Debug, Clone, PartialEq)]
pub struct Strain {
    pub alpha: f64,
    pub biomass: f64
}
impl Strain {
    pub fn new(alpha: f64) -> Strain {
        Strain {
            alpha,
            biomass: get_biomass(alpha)
        }
    }

    pub fn default() -> Strain {
        Strain {alpha: std::f64::NAN, biomass: 0f64}
    }
}

impl Ord for Strain {
    fn cmp(&self, other: &Self) -> Ordering {
        self.alpha.partial_cmp(&other.alpha).unwrap().reverse()
    }
}

impl PartialOrd for Strain {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Eq for Strain {}



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


pub fn get_toxin_biomass(alpha: f64) -> f64 {
    let (d, foodIn, toxinIn) = (0.1f64, 1.0f64, 0.9f64);
    let sol: f64 = (2.0 * d.powi(3) + d * (toxinIn + foodIn * (alpha - 1.0)) + d.powi(2) * alpha + foodIn * (alpha - 1.0) * alpha
            - (-4.0 * d.powi(3) * toxinIn * alpha + (d * (toxinIn + foodIn * (alpha - 1.0)) + d.powi(2) * alpha + foodIn * (alpha - 1.0) * alpha).powi(2)).sqrt())
            / (2.0 * d * (d + alpha));
    if sol.is_nan() || sol.is_sign_negative() {
        return 0f64
    }
    return sol
}


pub fn get_yeast_biomass(alpha: f64) -> f64 {
    let (d, supply, omega) = (0.2f64, 1f64, 0.2f64);
    let sol: f64 = -(1.0 / (2.0 * alpha * (omega - 1.0) - 2.0 * omega)) * 
        (supply * (alpha - 1.0) * (alpha * (omega - 2.0) - omega) + d * (-1.0 + (alpha - 1.0) * omega) + 
        (
            ( supply * (alpha - 1.0) * (alpha * (omega - 2.0) - omega) + d * (1.0 + alpha * (omega - 2.0)) ).powi(2) -
            4.0 * d * supply * (alpha - 1.0).powi(2) * (omega - omega.powi(2) + alpha * (2.0 + (omega - 2.0) * omega))
        ).sqrt()
    );


    if sol.is_nan() || sol.is_sign_negative() {
        return 0f64
    }
    return sol
}


// Truncated normal distribution by rejection
pub fn trunc_norm<R: Rng + ?Sized>(mean: f64, std_dev: f64, min: f64, max: f64, rng: &mut R) -> f64 {
    loop {
        let x: f64 = Normal::new(mean, std_dev).unwrap().sample(rng);
        if min <= x && x <= max {
            return x;
        }
    }
}
