// extern crate csv;

use std::collections::HashMap;
use std::vec::Vec;
use std::f64::consts::PI;
use std::error;
// use csv::Writer;

const BETA: f64 = 0.5;
const MIN_R: u64 = 5695325;
const MAX_R: u64 = 5695325;
const VALID_RESIDUALS: [u64; 5] = [1,4,5,6,9];

fn main() -> Result<(), Box<dyn error::Error>>{
   let res = get_max_lattice_points(MIN_R, MAX_R);

    let max_count = res.values().max().unwrap();
    let mut max_rads = find_keys_for_value(&res,max_count);
    max_rads.sort();

    println!("From radius = {} to {}", MIN_R, MAX_R);
    println!("The maximum number of lattice points for beta = {} is {:?}", BETA, max_count);
    println!("The corresponding squared radii are {:?}", max_rads);

    // let mut wtr = Writer::from_path("jarnik1to10.csv")?;
    // wtr.write_record(&["Radius", "Max Lattice Points"])?;
    // for (k,v) in res.iter() {
    //     wtr.serialize(&[k, v])?;
    // }
    // wtr.flush()?;
    Ok(())
}

fn find_keys_for_value<'a>(map: &'a HashMap<u64,u64>, value: &u64) -> Vec<&'a u64> {
    map.iter()
        .filter_map(|(key, val)| if val == value { Some(key) } else { None })
        .collect()
}

fn get_max_lattice_points<'a>(min_r: u64, max_r: u64) -> HashMap<u64, u64> {
    let mut rad_counts: HashMap<u64,u64> = HashMap::new();

    for R in min_r..=max_r {
        let radius = (R as f64).sqrt();

        let mut lattice_points: Vec<(f64, f64)> = Vec::new();
        let mut angles: Vec<f64> = Vec::new();

        // lattice_points.push((0.0, radius));

        for i in (((radius / (2 as f64).sqrt()).ceil() as u64)..=(radius.floor() as u64)).rev() {
            let i = i as f64;
            let cand = radius.powf(2.0) - i.powf(2.0);
            if cand == 0.0 {
                lattice_points.push((0.0, i));
            } else if VALID_RESIDUALS.contains(&((cand as u64) % 10)) {
                let temp = cand.sqrt();
                if temp.fract() == 0.0 {
                    lattice_points.push((temp, i));
                }
            }
        }
        // println!{"{:?}",lattice_points};
        let lp_len = lattice_points.len();
        if lp_len > 0 {
            angles.push(2.0 * lattice_points[0].1.atan2(0.0));

            for i in 0..(lp_len-1) {
                let (x1, y1) = lattice_points[i];
                let (x2, y2) = lattice_points[i + 1];
                let a1 = y1.atan2(x1);
                let a2 = y2.atan2(x2);
                angles.push(a1 - a2)
            }
        } else {
            rad_counts.insert(R, 0);
            continue
        }

        angles.push( 2.0 * lattice_points[lp_len-1].1.atan2(lattice_points[lp_len-1].0) - PI/2.0 );
        let mut rev_angles:Vec<f64> = Vec::new();
        rev_angles.extend(angles.clone().iter().rev());
        angles.extend_from_slice(&rev_angles[1..rev_angles.len()]);
        angles.extend(&mut angles.clone().iter());

        let mut count: u64 = 0;

        let arc_length: f64 = radius.powf(BETA);
        let max_ang: f64 = arc_length / radius;

        for i in 0..angles.len() {
            let mut running_ang_sum: f64 = 0.0;
            let mut running_count = 1;

            for ang in &angles[i..] {
                running_ang_sum += *ang;
                if running_ang_sum >= max_ang {
                    break
                }
                running_count += 1;
            }
            if running_count > count {
                count = running_count;
            }
        }
        rad_counts.insert(R, count);
    }
    rad_counts
}
