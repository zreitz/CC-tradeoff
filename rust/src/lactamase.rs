use rand::prelude::*;
use rand_chacha::ChaCha8Rng;
use rand_distr::{Normal, WeightedIndex};
use serde::{Serialize, Deserialize};
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::env;
use std::fs;
use std::fs::OpenOptions;
use std::io::prelude::*;
use std::path::PathBuf;
use zzz::ProgressBarIterExt as _;


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
    let (d, foodIn, toxinIn) = (0.1f64, 1.0f64, 0.9f64);
    let sol: f64 = - (2.0 * d.powi(3) + d * (toxinIn + foodIn * (alpha - 1.0)) + d.powi(2) * alpha + foodIn * (alpha - 1.0) * alpha
            - (-4.0 * d.powi(3) * toxinIn * alpha + (d * (toxinIn + foodIn * (alpha - 1.0)) + d.powi(2) * alpha + foodIn * (alpha - 1.0) * alpha).powi(2)).sqrt())
            / (2.0 * d * (d + alpha));
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


#[derive(Serialize, Deserialize, Debug)]
struct Session {
    description: String,
    rounds: usize,
    write_every: usize,
    mortality: f64,
    max_biomass: f64,
    mut_stdev: f64,
    mut_uniform: f64,
    mut_range: [f64; 2],
    seed: u64,
    initial: Vec<f64>
}

impl Session {
    fn from_json(path: &String) -> Session {
        // Get parameters
        let mut file = fs::File::open(path)
                .expect("Couldn't open session file");
        let mut data = String::new();
        file.read_to_string(&mut data).unwrap();
        serde_json::from_str(&data)
                .expect("session file was misformatted")
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
struct Resident {
    time: Option<usize>,
    alpha: f64,
    frequency: f64,
    pressure: f64
}

// Use iterative niche shadows to determine survivors
fn get_survivors(mut strains: BinaryHeap<Strain>, mortality: &f64) 
            -> BinaryHeap<Strain> {
    let mut survivors: BinaryHeap<Strain> = BinaryHeap::new();
    let mut threshold: f64 = *mortality;

    while !strains.is_empty() {
        let this: Strain = strains.pop().unwrap();
        if this.biomass > threshold {
            threshold = this.biomass.powi(2) / threshold;
            survivors.push(this);
        }
    }
    return survivors
}

// Use iterative formula to calculate frequencies
//  Returns a vector of alpha, frequency, pressure
//  Consumes the strain list; give it a clone
fn get_freqs(mut strains: BinaryHeap<Strain>, mortality: &f64) 
            -> Vec<Resident> {
    let mut results: Vec<Resident> = vec![];
    let mut freq_sums: f64 = 0.0;
    let mut mass_sums: f64 = 0.0;

    while !strains.is_empty() {
        let this: Strain = strains.pop().unwrap();
        let freq: f64 = 1.0 - freq_sums - (mortality + mass_sums) / this.biomass;
        assert!(freq.signum() == 1.0, 
                "ERROR: found negative frequency; run get_survivors() first!");
        freq_sums += freq;
        mass_sums += freq * this.biomass;
        results.push(Resident{
            time: None,
            alpha: this.alpha, 
            frequency: freq,
            pressure: freq * this.biomass });
    }

    return results
}

// Takes a heap of strains and a vec of their pressures
// Generates a mutant and returns the survivors
fn mutate_and_resolve<R: Rng + ?Sized>(
                    mut strains: BinaryHeap<Strain>,
                    residents: Vec<Resident>, 
                    session: &Session, 
                    rng: &mut R
                    ) -> BinaryHeap<Strain>{

    let (alphas, weights): (Vec<f64>, Vec<f64>) = residents
                        .into_iter().map(|x| (x.alpha, x.pressure)).unzip();

    // Draw if nothing happens
    let mut_prob: f64 = &weights.iter().sum() / session.max_biomass;
    let fate: f64 = rng.gen::<f64>();
    if fate >= session.mut_uniform + mut_prob {
        return strains
    }

    let new_alpha: f64; 
    // Uniform distribution
    if fate < session.mut_uniform {
        new_alpha = rng.gen_range(session.mut_range[0]..session.mut_range[1]);
    }
    
    // Truncated normal distribution from a parent strain
    else {
        // Choose parent strain biased by abundance and mutate it
        let dist = WeightedIndex::new(weights).unwrap();
        let parent = alphas[dist.sample(rng)];
        new_alpha = trunc_norm(parent, session.mut_stdev,
            session.mut_range[0], session.mut_range[1], rng)
    }
    strains.push(Strain::new(new_alpha));
    return get_survivors(strains, &session.mortality)
}


fn simulate_pop(mut strains: BinaryHeap<Strain>, session: &Session) 
                -> Vec<Resident> {
    // initialize RNG
    let mut rng = ChaCha8Rng::seed_from_u64(session.seed);

    // Results storage: vec<vec<(alpha, freq)>>
    let mut census: Vec<Resident> = vec![];

    // Initialize population
    let mut residents: Vec<Resident>;

    residents = get_freqs(strains.clone(), &session.mortality);
    for r in &mut residents {
        r.time = Some(0);
    }

    // Save initial freqs
    census.extend_from_slice(&residents);

    // Main loop
    for round in (0..(session.rounds)).into_iter().progress() {
        strains = mutate_and_resolve(strains, residents, &session, &mut rng);
        residents = get_freqs(strains.clone(), &session.mortality);
        
        if round % session.write_every == 0 || round == session.rounds {
            for r in &mut residents {
                r.time = Some(round);
            }
            census.extend_from_slice(&residents);
        }
    }

    return census
}


fn main() {
    /*
    let session = Session {
            description: "test".to_string(),
            rounds: 10,
            write_every: 1,
            mortality: 1.5,
            max_biomass: 0.0,
            mut_stdev: 1e-4,
            mut_uniform: 1.0,
            mut_range: [0.05, 0.85],
            seed: 93117,
            initial: vec!(0.075),
            results: vec!()
    };
     */
    let args: Vec<String> = env::args().collect();
    // Get parameters
    assert!(args.len() > 1, "ERROR. Example usage: $ cargo run -- in_config [outpath]");
    let session = Session::from_json(&args[1]);

    // Prepare outpath
    let mut outpath: PathBuf;
    if args.len() == 2 {
        outpath = std::env::current_dir()
        .expect("could not get directory")
        .join("results");
        fs::create_dir_all(&outpath)
                .expect("could not create output directory");
        outpath.push(chrono::offset::Local::now()
                .format("%Y-%m-%d_%H-%M-%S.csv").to_string());
    }
    else {
        outpath = PathBuf::from(&args[2]);
        fs::create_dir_all(&outpath.parent().unwrap())
                .expect("could not create output directory");
    }
    assert!(! outpath.is_file(), "File exists. Exiting.");

    // Initialize metapopulation
    let strains: BinaryHeap<Strain> = BinaryHeap::from(
                session.initial
                .iter().map(|alpha| Strain::new(*alpha))
                .collect::<Vec<Strain>>());
    
    // MAIN SIMULATION and add to session
    let results = simulate_pop(strains, &session);
    
    // Write the main results as CSV
    let mut wtr = csv::Writer::from_path(&outpath).unwrap();
    for row in results {
        wtr.serialize(row).unwrap();
    }
    wtr.flush().unwrap();

    // Serialize the session to JSON, this will be a comment at the bottom
    let mut serialized_stats = serde_json::to_string(&session).unwrap();
    serialized_stats.insert_str(0, "# ");
    let mut file = OpenOptions::new().append(true).open(outpath).unwrap();
    file.write_all(serialized_stats.as_bytes()).unwrap();
}
