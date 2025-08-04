extern crate nalgebra as na;
use na::{*};

use rand::prelude::*;
use rand::seq::index::sample;
use rand_chacha::ChaCha8Rng;
use rand_distr::weighted::WeightedIndex;
use rand_distr::{Binomial, Uniform};
use rayon::prelude::*;
use serde::{Serialize, Deserialize};
use std::fs::OpenOptions;
use std::path::PathBuf;
use std::{env, fs};
use std::io::prelude::*;

use solve::solve;


/* 
Functions herein are hard-coded for 50 possible alpha values
plus three resource state variables. Replace all instances of 50 to 53 when changing that.
*/

static INIT_HOLO: f64 = 0.2;

#[derive(Serialize, Deserialize, Debug)]
struct Session {
    // "meta" parameters
    description: String,
    replicates: usize,
    rounds: usize,
    seed: u64,
    write_every: usize,
    print_every: usize,
    
    // metapopulation setup
    pop_size: usize,
    inoc_size: f64,

    mort_to_disp: f64,
    patch_lifetime: f64,
    uniform_fraction: f64
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

// Matrix has 53 rows (state variables) and [pop_size] columns (patches)
type PatchMatrix = Matrix<f64, U53, Dyn, VecStorage<f64, U53, Dyn>>;

#[derive(Debug)]
struct Population {
    // Represents patch concentrations. Each col is a patch.
    states: PatchMatrix,
    alphas: Vec<f64>,
    at_eq: Vec<bool>,
    time: f64
}

impl Population {

    fn initialize<R: Rng + ?Sized>(session: &Session, rng: &mut R) -> Population {
        //let alphas: Vec<f64> = (0..50).map(|x| x as f64 / 100.0).collect();
        let alphas: Vec<f64> = draw_strains(
            50, 0, 4000, rng);
        
        // Make a small population to calculate equilibria for each strain
        let mut eq_patches: PatchMatrix = PatchMatrix::repeat(50, 0.0);
        eq_patches.set_row(52, &RowDVector::repeat(50, 3.0));

        let mut eq_pop = Population{
            states: eq_patches, alphas: alphas.clone(), at_eq: vec![false; 50], time: 0.0};
        for (i, alpha) in (0..50).enumerate() {
            eq_pop.inoculate(i, alpha, 2.0);
        }
        
        // Run to equilibrium
        eq_pop.step_forward(50.0);

        // Generate a new patch matrix
        let mut patches: PatchMatrix = PatchMatrix::repeat(session.pop_size, 0.0);
        // Last row is iron - initialize to equilibrium
        patches.set_row(52, &RowDVector::repeat(session.pop_size, 1.0));
        // Add holo siderophore to prevent death with small inocula
        patches.set_row(51, &RowDVector::repeat(session.pop_size, INIT_HOLO));

        Population{
            states: patches,
            alphas: alphas,
            at_eq: vec![true; session.pop_size],
            time: 0.0
        }
    }

    /*
    Increase the value of a given (strain, patch) coordinate by inoc_size
     */
    fn inoculate(&mut self, patch: usize, invader: usize, inoc_size: f64) {
        self.states[(invader, patch)] += inoc_size;
        self.at_eq[patch] = false;
    }

    fn step_forward(&mut self, stepsize: f64) {
        self.states
            .par_column_iter_mut()
            .zip(self.at_eq.par_iter_mut())
            .for_each(|(mut col, at_eq)| {
                // Do nothing if the patch is at equilibrium
                if *at_eq {return};
    
                let patch: SVector<f64, 53> = col.clone_owned();
    
                // Step forward with the numerical ODE solver
                let (mut result, new_at_eq) = solve(
                    self.alphas.clone(),
                    patch,
                    stepsize
                );
    
                // Cull dead strains (<1e-3)
                let mut living = 0;
                for strain in 0..50 {
                    if result[strain] < 1.0e-3 {
                        result[strain] = 0.0;
                    }
                    else {
                        living += 1
                    }
                }
                if living == 0 {    // Reset patch
                    result[50] = 0.0;
                    result[51] = INIT_HOLO;
                    result[52] = 1.0;
                }
    
                // Update the states matrix and residents list
                col.copy_from(&result);
                *at_eq = new_at_eq;
            });
    }

    fn turnover(&mut self, patch: usize) {
        self.states.set_column(patch, &SVector::repeat( 0.0));
        self.states[(51, patch)] = INIT_HOLO;
        self.states[(52, patch)] = 1.0;
        self.at_eq[patch] = true;
    }
}


// Draw n strains from a uniform distribution
fn draw_strains<R: Rng + ?Sized>(n: usize, min: u64, max: u64, rng: &mut R) -> Vec<f64> {
    let dist = Uniform::new(min, max).unwrap();
    let mut samples = Vec::with_capacity(n);

    while samples.len() < n {
        let x = dist.sample(rng);
        if !samples.contains(&x){
            samples.push(x);
        }
    }
    let mut strains = samples.iter()
            .map(|x| *x as f64 / 10000.).collect::<Vec<f64>>();
    strains.sort_by(f64::total_cmp);
    return strains
}


/*
Perform one round of kick / flow
    */
fn simulate_round<R: Rng + ?Sized>(mut population: Population, session: &Session, rng: &mut R, print_stats: bool) -> Population {
    // Calculate colonization pressures for each strain
    let col_pool:SVector<f64, 50> = population.states.column_sum()
                            .remove_fixed_rows::<3>(50);
    let col_pressure:f64 = col_pool.iter().sum();

    // Uniform distribution for the immigrants
    let uniform_dist = Uniform::new(0,50).unwrap();

    // If the population is dead, inoculate patch 0 with a random immigrant, step forward, and return
    if col_pressure < 1e-3 {
        population.inoculate(0, uniform_dist.sample(rng), session.inoc_size);
        population.step_forward(50.0);
        population.time += 50.0;
        return population
    }

    // Calculate rates
    let prob_wiped = 1.0 / session.patch_lifetime; // mu
    let prob_col = prob_wiped / session.mort_to_disp * col_pressure / session.pop_size as f64; 
    let prob_event = prob_wiped + prob_col;

    // Dynamically set time step to scale probability so that 5% of patches are kicked
    // This is just for housekeeping, shouldn't affect dynamics
    let stepsize = - 0.95.ln() / (prob_event);

    if print_stats {
        println!("Time step: {}", stepsize);
    }
    
    // 5% of patches have an event happen
    let num_events = Binomial::new(session.pop_size as u64, 0.05)
            .unwrap().sample(rng);

    // Draw the number of colonizations (the rest are wiped)
    let num_col = Binomial::new(num_events, prob_col / prob_event)
            .unwrap().sample(rng) as usize;

    // Draw number of colonizations that are from uniform
    let num_unif = Binomial::new(num_col as u64, session.uniform_fraction)
            .unwrap().sample(rng) as usize;

    // Draw who is being kicked (without replacement)
    let kicked:Vec<usize> = sample(rng, session.pop_size, num_events as usize).into_vec();

    // Prepare distributions to draw colonizers from
    let weighted_dist = WeightedIndex::new(&col_pool)
            .expect(&format!("Error making WeightedIndex for \n{:?}\n", col_pool));
    
    for (i, victim) in kicked.into_iter().enumerate() {
        if i < num_col {
            let colonizer = 
                if i < num_unif { uniform_dist.sample(rng) }
                else { weighted_dist.sample(rng) };
            population.inoculate(victim, colonizer, session.inoc_size);
        }
        else {
            population.turnover(victim);
        }
    }

    if print_stats {
        println!("Flowing {} patches", population.at_eq.iter().filter(|x| !*x).count());
    }

    // Flow
    population.step_forward(stepsize);
    population.time += stepsize;

    return population
}


fn prepare_outpath(args: Vec<String>) -> PathBuf {
    let mut outpath: PathBuf;
    // If output file is not provided, make a timestamped file
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
    return outpath;
    }

fn main() {
    env::set_var("RUST_BACKTRACE", "1");

    // Get parameters from config file
    let args: Vec<String> = env::args().collect();
    assert!(args.len() > 1, "ERROR. Example usage: $ cargo run -- params.json");
    let session = Session::from_json(&args[1]);

    // Prepare output file
    let outpath: PathBuf = prepare_outpath(args);

    // Serialize the session to JSON, this will be a comment
    let mut serialized_stats = serde_json::to_string(&session).unwrap();
    serialized_stats.insert_str(0, "# ");
    serialized_stats.push_str("\n");

    let mut handle = OpenOptions::new()
        .create(true)
        .write(true)
        .open(outpath)
        .unwrap();
    handle.write_all(serialized_stats.as_bytes()).unwrap();


    // Prepare to write results as csv
    let mut wtr = csv::Writer::from_writer(handle);

    // MAIN LOOP
    let mut seed = session.seed;
    for replicate in 0..session.replicates {
        // initialize RNG
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        
        // Initialize matrix representing patches. Each col is a patch
        let mut population = Population::initialize(&session, &mut rng);
        
        // header row gives the alphas
        let mut header = population.alphas.iter()
                .map(|x| {format!("{:.4}", x)})
                .collect::<Vec<String>>();
        header.push("time".to_string());
        header.push(seed.to_string());
        wtr.serialize(header).unwrap();

        for i in 0..=session.rounds {
            //eprintln!("{}", population.states.column(0));
            if i % session.write_every == 0 {
                let mut total_biomasses:Vec<f64> = population.states
                    .column_sum()
                    .data.0[0][0..50]
                    .to_vec();
                total_biomasses.push(population.time);
                total_biomasses.push(seed as f64);
                wtr.serialize(total_biomasses).unwrap();
            }
            if i % session.print_every == 0 {
                eprintln!("{replicate}..{i}: {:.3?}", population.states.map(|x| if x > 0.0 {1} else {0}).column_sum().map(|x| x as f64 / session.pop_size as f64).data);
                population = simulate_round(population, &session, &mut rng, true);
            }
            else {
                population = simulate_round(population, &session, &mut rng, false);
            }
        }
        wtr.flush().unwrap();

        seed += 1;
    }
}