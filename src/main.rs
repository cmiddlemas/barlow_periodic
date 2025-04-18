extern crate barlow_periodic;
extern crate structopt;
extern crate chrono;

use barlow_periodic::*;

use std::fs::File;
use std::io::prelude::*;
use structopt::StructOpt;
use chrono::prelude::*;

fn main() {

    // From the chrono docs, heard about chrono from
    // rust reddit and stackoverflow
    // was originally using std::time, with
    // associated docs as reference
    println!("Starting at time: {}", Utc::now());
    
    let opt = Opt::from_args();

    println!("Got argument struct : {:?}", opt);
    
    println!("Parsing Barlow strings from {:?}", opt.input);
    
    let mut input_file = File::open(&opt.input)
        .expect("Couldn't open input file");
    let mut barlow_buffer = String::new();
    input_file.read_to_string(&mut barlow_buffer)
        .expect("Error reading input file");
    
    for line in barlow_buffer.lines() {
        println!("Working on packing: {}", line);
        let layer_spec = Spec::new_with_check(line);
        match opt.metric.as_str() {
            "theta" => compute_theta(layer_spec, &opt),
            "hyperuniformity" => compute_hyperuniformity(layer_spec, &opt),
            "integral" => compute_h(layer_spec, &opt),
            "ap-integral" => ap_compute_h(layer_spec, &opt),
            "ap-hyperuniformity" => ap_compute_hyper(layer_spec, &opt),
            _ => panic!("Invalid request, or not implemented."),
        }
    }
    
    // same as above
    println!("Exiting at time: {}", Utc::now());

}
