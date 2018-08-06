#[macro_use]
extern crate structopt;
#[macro_use]
extern crate lazy_static;
extern crate nalgebra as na;
extern crate rayon;
extern crate rug;

use std::collections::VecDeque;
use std::path::PathBuf;
use std::io::prelude::*;
use std::fs::File;
use na::Vector3;
use rayon::prelude::*;
use std::f64::consts::PI;
use rug::Float;

// Constant Declarations for whole crate
const THETA_BIN_SIZE: f64 = 0.01;

// Constants for doing lattice sums. Non-obvious
// numerical values computed through Mathematica 11.2.0
const Z_OFFSET: f64 = 0.816496580927726;

const X_OFFSET_B: f64 = 0.5;
const Y_OFFSET_B: f64 = 0.2886751345948129;

const X_OFFSET_C: f64 = -0.5;
const Y_OFFSET_C: f64 = -0.2886751345948129;

const X_HEX_1: f64 = 1.0;
const Y_HEX_1: f64 = 0.0;
const X_HEX_2: f64 = 0.5;
const Y_HEX_2: f64 = 0.866025403784439;

// lazy_static declarations of useful vectors
lazy_static! {
    static ref OFFSET_ALL: Vector3<f64> = Vector3::new(0.0, 0.0, Z_OFFSET);
    static ref OFFSET_B: Vector3<f64> = Vector3::new(X_OFFSET_B, Y_OFFSET_B, 0.0);
    static ref OFFSET_C: Vector3<f64> = Vector3::new(X_OFFSET_C, Y_OFFSET_C, 0.0);
    static ref HEX_1: Vector3<f64> = Vector3::new(X_HEX_1, Y_HEX_1, 0.0);
    static ref HEX_2: Vector3<f64> = Vector3::new(X_HEX_2, Y_HEX_2, 0.0);
}

// StructOpt derivation is based on Github sample code for structopt
///A program to calculate order metrics for periodic Barlow packings.
#[derive(StructOpt, Debug)]
#[structopt(name = "Periodic Barlow Order Metric Calculator")]
pub struct Opt {
    ///Choose order metric.
    ///Panics on invalid request.
    ///Choices:
    ///theta (mostly for debugging, preprogrammed bin size),
    ///T,
    ///T*,
    ///tau,
    ///hyperuniformity
    #[structopt(short = "m", long = "metric")]
    pub metric: String,

    ///Distance cutoff, in units of particle diameter
    #[structopt(short = "c", long = "cutoff")]
    pub cutoff: f64,

    ///Input file. Format:
    ///ABCABC\nABABAB
    #[structopt(short = "i", long = "input", parse(from_os_str))]
    pub input: PathBuf,

    ///Output file (optional)
    #[structopt(short = "o", long = "output", parse(from_os_str))]
    pub output: Option<PathBuf>,

    ///Convergence parameter for hyperuniformity B coefficient
    #[structopt(short = "a", long = "alpha")]
    pub alpha: Option<f64>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Layer {
    A,
    B,
    C,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Spec ( VecDeque<Layer> );

// Spec invariants (may panic if not upheld):
// Non-empty
// Starts with A
// Ends with B or C
impl Spec {
    fn len(&self) -> usize {
        self.0.len()
    }

    fn get(&self, n: usize) -> Layer {
        self.0[n]
    }

    pub fn new_with_check(raw_spec: &str) -> Spec {
        use Layer::*;

        // Parse the string read in from file
        let spec: Spec = Spec ( raw_spec.chars()
            .map(|c| match c {
                'A' => A,
                'B' => B,
                'C' => C,
                _ => panic!("Invalid spec string"), }
            )
            .collect() );
        
        // Check invariants
        assert!(*spec.0.front().expect("Spec must be non-empty") == A);
        assert!(*spec.0.back().expect("Spec must be non-empty") == B 
                || *spec.0.back().expect("Spec must be non-empty") == C);
        
        spec
    }
    
    fn shift(&self, n: usize) -> Spec {
        use Layer::*;
        
        // Rotate list of layers
        let mut new_spec = self.clone();
        for _ in 0..n {
            let front = new_spec.0.pop_front().expect("Spec can't be empty.");
            new_spec.0.push_back(front);
        }
        
        // Rename layers so beginning is still A
        match *new_spec.0.front().expect("Spec must be non-empty") {
            A => (),
            B => {
                for layer in &mut new_spec.0 {
                    *layer = match *layer {
                        A => C,
                        B => A,
                        C => B,
                    };
                }
            },
            C => {
                for layer in &mut new_spec.0 {
                    *layer = match *layer {
                        A => B,
                        B => C,
                        C => A,
                    };
                }
            },
        }
        new_spec
    }
}

pub fn theta_hex(
    layer: Layer, // A B or C
    l_idx: i64, //How many layers above or below center are we?
    max_lat_idx: i64, //How far are we going out in lattice sum?
    theta: &mut Vec<f64>) //Holder for theta series
{
    use Layer::*;
    let l_f64 = l_idx as f64;
    
    let offset = match layer {
        A => l_f64 * *OFFSET_ALL,
        B => l_f64 * *OFFSET_ALL + *OFFSET_B,
        C => l_f64 * *OFFSET_ALL + *OFFSET_C,
    };

    for i in (-max_lat_idx)..(max_lat_idx + 1) {
        for j in (-max_lat_idx)..(max_lat_idx + 1) {
            let vec = offset + (i as f64) * *HEX_1 + (j as f64) * *HEX_2;
            let r: f64 = na::norm(&vec);
            let idx = (r/THETA_BIN_SIZE) as usize;
            if idx < theta.len() {
                theta[idx] += 1.0;
            }
        }
    }

}

pub fn compute_theta(spec: Spec, opt: &Opt) {
    // turn radial cutoff into integer lattice cutoff
    // probably overkill
    let max_lat_idx = (2.0*opt.cutoff) as i64;
    
    let spec_len = spec.len();
    let spec_i64 = spec_len as i64;

    // set up theta series bins
    let n_bins = (opt.cutoff/THETA_BIN_SIZE) as usize;
    let domain: Vec<f64> = (0..n_bins)
        .map(|x| (x as f64)*THETA_BIN_SIZE)
        .collect();
    let mut theta: Vec<f64> = vec![0.0; n_bins];
    
    // iterate over particle types
    for i in 0..spec_len {
        println!("Working on particle type: {}", i);
        let curr_spec = spec.shift(i);
        // iterate over layers
        for j in (-max_lat_idx)..(max_lat_idx + 1) {
            //println!("Working on layer: {}", j);
            // standard modulo hack, see for example
            // https://stackoverflow.com/questions/31210357/is-there-a-modulus-not-remainder-function-operation
            let layer = curr_spec.get((((j%spec_i64)+spec_i64)%spec_i64) as usize);
            theta_hex(layer, j, max_lat_idx, &mut theta);
        }
    }

    // Output theta series, optionally to file
    let mut output_file;
    if let Some(ref output_path) = opt.output {
        output_file = Some(File::create(output_path)
                           .expect("Error creating output file"));
    } else {
        output_file = None;
    }
    for (i,t) in theta.iter_mut().enumerate() {
        *t /= spec_len as f64;
        println!("{}, {}", domain[i], *t);
        match output_file {
            Some(ref mut buffer) => {
                write!(buffer, "{}, {}\n", domain[i], *t)
                    .expect("Error writing to output file");
            },
            None => (),
        }
    }
}

pub fn hyper_hex(
    layer_pos: Layer, // A B or C
    layer_neg: Layer,
    l_idx: i64, //How many layers above or below center are we?
    max_lat_idx: i64, //How far are we going out in lattice sum?
    alpha: f64 //exponential convergence parameter
    ) -> f64
{
    use Layer::*;
    let l_f64 = l_idx as f64;
    
    let offset_pos = match layer_pos {
        A => l_f64 * *OFFSET_ALL,
        B => l_f64 * *OFFSET_ALL + *OFFSET_B,
        C => l_f64 * *OFFSET_ALL + *OFFSET_C,
    };
    
    let offset_neg = match layer_neg {
        A => -l_f64 * *OFFSET_ALL,
        B => -l_f64 * *OFFSET_ALL + *OFFSET_B,
        C => -l_f64 * *OFFSET_ALL + *OFFSET_C,
    };

    let mut b: f64 = 0.0;

    for i in (0..(max_lat_idx + 1)) {
        for j in (0..(max_lat_idx + 1)) {
            // +++
            let vec = offset_pos + (i as f64) * *HEX_1 + (j as f64) * *HEX_2;
            let r: f64 = na::norm(&vec);
            b += r*((-alpha*r*r).exp());
            // ++-
            if j != 0 {
                let vec = offset_pos + (i as f64) * *HEX_1 - (j as f64) * *HEX_2;
                let r: f64 = na::norm(&vec);
                b += r*((-alpha*r*r).exp());
            }
            // +-+
            if i != 0 {
                let vec = offset_pos - (i as f64) * *HEX_1 + (j as f64) * *HEX_2;
                let r: f64 = na::norm(&vec);
                b += r*((-alpha*r*r).exp());
            }
            // +--
            if (i != 0) && (j != 0) {
                let vec = offset_pos - (i as f64) * *HEX_1 - (j as f64) * *HEX_2;
                let r: f64 = na::norm(&vec);
                b += r*((-alpha*r*r).exp());
            }
            // -++
            if l_idx != 0 {
                let vec = offset_neg + (i as f64) * *HEX_1 + (j as f64) * *HEX_2;
                let r: f64 = na::norm(&vec);
                b += r*((-alpha*r*r).exp());
            }
            // -+-
            if (l_idx != 0) && (j != 0) {
                let vec = offset_neg + (i as f64) * *HEX_1 - (j as f64) * *HEX_2;
                let r: f64 = na::norm(&vec);
                b += r*((-alpha*r*r).exp());
            }
            // --+
            if (l_idx != 0) && (i != 0) {
                let vec = offset_neg - (i as f64) * *HEX_1 + (j as f64) * *HEX_2;
                let r: f64 = na::norm(&vec);
                b += r*((-alpha*r*r).exp());
            }
            // ---
            if (l_idx != 0) && (i != 0) && (j != 0) {
                let vec = offset_neg - (i as f64) * *HEX_1 - (j as f64) * *HEX_2;
                let r: f64 = na::norm(&vec);
                b += r*((-alpha*r*r).exp());
            }
        }
    }
    b
}


pub fn compute_hyperuniformity(spec: Spec, opt: &Opt) {
    let alpha = opt.alpha
        .expect("Must supply convergence parameter alpha for B calculation.");

    let max_lat_idx = (2.0*opt.cutoff) as i64;

    let spec_len = spec.len();
    let spec_i64 = spec_len as i64;

    let mut b: f64 = 0.0; //B surficial variance coefficient
    
    for i in 0..spec_len {
        println!("Working on particle type: {}", i);
        let curr_spec = spec.shift(i);
        let temp_vec: Vec<f64> = (0..(max_lat_idx + 1))
            .collect::<Vec<i64>>()
            .par_iter()
            .map(|&j| hyper_hex(
                curr_spec.get((((j%spec_i64)+spec_i64)%spec_i64) as usize),
                curr_spec.get((((-j%spec_i64)+spec_i64)%spec_i64) as usize),
                j,
                max_lat_idx,
                alpha
                )
            )
            .collect();
            //println!("{:?}", temp_vec);
            //.sum::<f64>();
            b += temp_vec.iter().sum::<f64>();
    }

    println!("b is {}", b);
    
    // normalize by # of particle types
    b /= spec_len as f64;
    
    // account for the -1 in the definition of h
    b -= 4.0*PI*(2.0f64.sqrt())/(2.0*alpha*alpha);

    println!("B Coefficient (non-normalized) = {}", b);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn valid_parse() {
        use Layer::*;
        let mut deque = VecDeque::new();
        deque.push_back(A);
        deque.push_back(B);
        deque.push_back(C);
        deque.push_back(B);
        let output = Spec ( deque );
        assert!(Spec::new_with_check("ABCB") == output);
    }
    
    #[test]
    #[should_panic]
    fn invalid_parse1() {
        Spec::new_with_check("ABCA");
    }

    #[test]
    #[should_panic]
    fn invalid_parse2() {
        Spec::new_with_check("BCA");
    }

    #[test]
    fn test_shift_2() {
        assert!(Spec::new_with_check("ABCAB").shift(2) == Spec::new_with_check("ABCBC"));
    }
}


