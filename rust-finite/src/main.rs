use rand::prelude::*;
use rand::seq::index::sample;
use rand_chacha::ChaCha8Rng;
use rand_distr::{Binomial, Normal};
use rand_distr::weighted::WeightedIndex;
use serde::{Serialize, Deserialize};
use std::collections::BTreeMap;
use std::env;
use std::fs;
use std::fs::OpenOptions;
use std::io::prelude::*;
use std::path::PathBuf;

use rayon::prelude::*;
use std::sync::mpsc;
use std::thread;


// ============================================================
//  Strain
// ============================================================

#[derive(Debug, Clone, PartialEq)]
pub struct Strain {
    pub alpha:   f64,
    pub biomass: f64,
}

impl Strain {
    pub fn new(alpha: f64) -> Strain {
        Strain { alpha, biomass: get_biomass(alpha) }
    }

    pub fn is_viable(&self) -> bool {
        self.biomass > 0.0
    }
}


// ============================================================
//  Analytical biomass model
// ============================================================

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


// ============================================================
//  Within-patch invasion
// ============================================================

fn resolve_invasion(resident: &Strain, invader: &Strain) -> Option<Strain> {
    let winner = if invader.alpha < resident.alpha { invader } else { resident };
    if winner.is_viable() { Some(winner.clone()) } else { None }
}


// ============================================================
//  Mutation helper
// ============================================================

// Truncated normal distribution by rejection
pub fn trunc_norm<R: Rng + ?Sized>(mean: f64, std_dev: f64, min: f64, max: f64, rng: &mut R) -> f64 {
    loop {
        let x: f64 = Normal::new(mean, std_dev).unwrap().sample(rng);
        if min <= x && x <= max {
            return x;
        }
    }
}

// ============================================================
//  Session (JSON config)
// ============================================================

// One entry in the `initial` array, e.g. {"alpha": 0.1, "patches": 50}
#[derive(Serialize, Deserialize, Debug)]
struct InitStrain {
    alpha:   f64,
    patches: usize,
}

#[derive(Serialize, Deserialize, Debug)]
struct Session {
    description: String,
    replicates:  usize,
    rounds:      usize,
    write_every: usize,
    seed:        u64,

    // Metapopulation geometry
    pop_size:         usize,
    patch_lifetime:   f64,  // mean patch lifespan (time units)
    mort_to_disp:     f64,  // ratio: extinction rate / colonization rate
    uniform_fraction: f64,  // fraction of colonizations drawn uniformly
    decay_factor:     f64,  // fraction of seed bank that survives (0.0 for no memory)

    // Mutation
    max_biomass: f64,  // normaliser for biomass-scaled mutation probability
    mut_prob:    f64,  // baseline per-dispersal mutation probability
    mut_stdev:   f64,  // std-dev of truncated-normal mutant step
    mut_range:   [f64; 2],

    // Initial condition: list of {alpha, patches} pairs
    // Example:
    //   "initial": [
    //     {"alpha": 0.10, "patches": 100},
    //     {"alpha": 0.30, "patches":  50}
    //   ]
    initial: Vec<InitStrain>,
}

impl Session {
    fn from_json(path: &String) -> Session {
        let mut file = fs::File::open(path).expect("Couldn't open session file");
        let mut data = String::new();
        file.read_to_string(&mut data).unwrap();
        serde_json::from_str(&data).expect("session file was misformatted")
    }
}


// ============================================================
//  Population
// ============================================================

struct Population {
    patches: Vec<Option<Strain>>,
    time:    f64,
    seed_bank: BTreeMap<u64, f64>
}

impl Population {
    fn initialize(session: &Session) -> Population {
        let mut patches: Vec<Option<Strain>> =
            (0..session.pop_size).map(|_| None).collect();

        let mut slot = 0;
        for init in session.initial.iter() {
            let s = Strain::new(init.alpha);
            if s.is_viable() {
                for _ in 0..init.patches {
                    if slot < session.pop_size {
                        patches[slot] = Some(s.clone());
                        slot += 1;
                    }
                }
            }
        }

        Population { patches, time: 0.0, seed_bank: BTreeMap::new() }
    }

    fn col_pressure(&self) -> f64 {
        self.patches.iter()
            .filter_map(|p| p.as_ref())
            .map(|s| s.biomass)
            .sum()
    }
}


// ============================================================
//  Census output
// ============================================================

#[derive(Serialize, Deserialize, Debug)]
struct Record {
    replicate:   usize,
    round:       usize,
    time:        f64,
    alpha:       f64,
    patch_count: usize,
    biomass:     f64
}

fn take_census(population: &Population, replicate: usize, time: f64, round: usize) -> Vec<Record> {
    let mut map: BTreeMap<u64, (usize, f64)> = BTreeMap::new();
    for s in population.patches.iter().filter_map(|p| p.as_ref()) {
        let entry = map.entry(s.alpha.to_bits()).or_insert((0, 0.0));
        entry.0 += 1;
        entry.1 += s.biomass;
    }
    map.into_iter()
        .map(|(key, (patch_count, biomass))| Record {
            replicate,
            round,
            time,
            alpha: f64::from_bits(key),
            patch_count,
            biomass,
        })
        .collect()
}


// ============================================================
//  Single simulation round
//  Returns None if the metapopulation has gone extinct.
// ============================================================

fn simulate_round<R: Rng + ?Sized>(
        mut population: Population,
        session: &Session,
        rng: &mut R) -> Option<Population> {

        
    // -- Calculate colonization pressures, adding seed bank if relevant

    // Live pressure from current patches
    let live_pressure: f64 = population.col_pressure();

    // Decay seed bank
    for val in population.seed_bank.values_mut() {
        *val *= session.decay_factor;
    }
    population.seed_bank.retain(|_, v| *v > 1e-3);

    let bank_pressure: f64 = population.seed_bank.values().sum();
    let col_pressure: f64 = live_pressure + bank_pressure;

    if col_pressure < 1e-3 {
        return None;
    }

    // --- Event rates ---
    let prob_wiped = 1.0 / session.patch_lifetime;
    let prob_col   = prob_wiped / session.mort_to_disp
                   * col_pressure / session.pop_size as f64;
    let prob_event = prob_wiped + prob_col;

    let stepsize = -0.99_f64.ln() / prob_event;
    population.time += stepsize;

    // --- Draw events ---
    let num_events = (Binomial::new(session.pop_size as u64, 0.01)
        .unwrap().sample(rng) as usize)
        .min(session.pop_size);
    let num_col  = Binomial::new(num_events as u64, prob_col / prob_event)
        .unwrap().sample(rng) as usize;
    let num_unif = Binomial::new(num_col as u64, session.uniform_fraction)
        .unwrap().sample(rng) as usize;

    let victims: Vec<usize> = sample(rng, session.pop_size, num_events).into_vec();

    // Make distribution of possible strains weighted by biomass
    let (mut alphas, mut biomasses): (Vec<f64>, Vec<f64>) = population.patches.iter()
        .filter_map(|p| p.as_ref())
        .map(|s| (s.alpha, s.biomass))
        .unzip();

    for (&bits, &pressure) in population.seed_bank.iter() {
        alphas.push(f64::from_bits(bits));
        biomasses.push(pressure);
    }
    let weighted_dist = WeightedIndex::new(&biomasses).unwrap();

    for (i, victim) in victims.into_iter().enumerate() {
        if i < num_col {
            let colonizer_alpha = if i < num_unif {
                rng.random_range(session.mut_range[0]..session.mut_range[1])
            } else {
                let parent = alphas[weighted_dist.sample(rng)];
                if rng.random::<f64>() < session.mut_prob {
                    trunc_norm(parent, session.mut_stdev,
                               session.mut_range[0], session.mut_range[1], rng)
                } else {
                    parent
                }
            };

            let colonizer = Strain::new(colonizer_alpha);
            population.patches[victim] = match &population.patches[victim] {
                None           => colonizer.is_viable().then_some(colonizer),
                Some(resident) => resolve_invasion(resident, &colonizer),
            };
        } else {
            population.patches[victim] = None;
        }
    }

    // Inject current round's live pressures into bank
    if session.decay_factor > 0.0 {
        for s in population.patches.iter().filter_map(|p| p.as_ref()) {
            *population.seed_bank.entry(s.alpha.to_bits()).or_insert(0.0) += s.biomass;
        }
    }

    Some(population)
}


// ============================================================
//  Single replicate
// ============================================================

fn simulate_replicate(session: &Session, replicate: usize, seed: u64) -> Vec<Record> {
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let mut population = Population::initialize(session);
    let mut census: Vec<Record> = vec![];

    census.extend(take_census(&population, replicate, population.time, 0));

    for round in 1..session.rounds {
        match simulate_round(population, session, &mut rng) {
            // Record the time of extinction
            None => {
                census.push(Record {
                    replicate,
                    round,
                    time: 0.0,
                    alpha: f64::NAN,
                    patch_count: 0,
                    biomass: 0.0,
                });
                break;
            },
            Some(pop) => {
                population = pop;
                if round % session.write_every == 0 || round == session.rounds - 1 {
                    census.extend(take_census(&population, replicate, population.time, round));
                    if round % 10000 == 0 {
                        eprintln!("replicate {replicate} | round {round} | t={:.2} | occupied={}/{}",
                            population.time,
                            population.patches.iter().filter(|p| p.is_some()).count(),
                            session.pop_size);
                    }
                }
            }
        }
    }

    census
}


// ============================================================
//  Output path helper
// ============================================================

fn prepare_outpath(args: &[String]) -> PathBuf {
    let mut outpath: PathBuf;
    if args.len() == 2 {
        outpath = std::env::current_dir()
            .expect("could not get directory")
            .join("results");
        fs::create_dir_all(&outpath).expect("could not create output directory");
        outpath.push(chrono::offset::Local::now()
            .format("%Y-%m-%d_%H-%M-%S.csv").to_string());
    } else {
        outpath = PathBuf::from(&args[2]);
        fs::create_dir_all(outpath.parent().unwrap())
            .expect("could not create output directory");
    }
    assert!(!outpath.is_file(), "File exists. Exiting.");
    outpath
}


// ============================================================
//  main
// ============================================================

fn main() {
    let args: Vec<String> = env::args().collect();
    assert!(
        args.len() > 1,
        "ERROR. Example usage: $ cargo run -- in_config.json [outpath.csv]"
    );

    let session = Session::from_json(&args[1]);
    let outpath = prepare_outpath(&args);

    // Initialize an MPSC channel to send data to the writer thread
    let (tx, rx) = mpsc::channel::<Vec<Record>>();

    // Spawn a dedicated thread for CSV writing to prevent I/O blocking in the worker pool
    let outpath_clone = outpath.clone();
    let writer_thread = thread::spawn(move || {
        let mut wtr = csv::Writer::from_path(&outpath_clone).unwrap();
        
        // Loop runs until all `tx` senders are dropped
        for results in rx {
            for row in results {
                wtr.serialize(row).unwrap();
            }
            wtr.flush().unwrap();
        }
    });

    // Execute replicates in parallel
    (0..session.replicates)
        .into_par_iter()
        .for_each_with(tx, |tx, replicate| {
            let seed = session.seed + replicate as u64;
            
            let results = simulate_replicate(&session, replicate, seed);
            
            // Transmit results to the writer thread
            tx.send(results).expect("Failed to send results to writer thread");
        });

    writer_thread.join().expect("Writer thread panicked");

    // Append trailing JSON comment for provenance sequentially
    let mut serialized = serde_json::to_string(&session).unwrap();
    serialized.insert_str(0, "# ");
    let mut file = OpenOptions::new().append(true).open(&outpath).unwrap();
    file.write_all(serialized.as_bytes()).unwrap();
}