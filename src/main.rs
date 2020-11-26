// extern crate csv;

use std::collections::HashMap;
use std::vec::Vec;
use std::f64::consts::PI;
use std::error;
use std::cmp::min;
// use csv::Writer;

const BETA: f64 = 1.0/3.0;
const N : usize = 4;
const MIN_R: u64 = 65;
const MAX_R: u64 = 65;
const VALID_RESIDUALS: [u64; 5] = [1,4,5,6,9];

fn main() -> Result<(), Box<dyn error::Error>>{
    // let lp_counts = get_latticepoint_counts(MIN_R, MAX_R);
    let lp_coeffs = get_arclength_coeffs(MIN_R, MAX_R);

    // let max_count = lp_counts.values().max().unwrap();
    let min_coeffs = lp_coeffs.values().fold(f64::INFINITY, |a, &b| a.min(b));
    // let mut max_rads = find_keys_for_value(&res,max_count);
    // max_rads.sort();

    println!("From radius = {} to {}", MIN_R, MAX_R);
    // println!("The maximum number of lattice points for beta = {:.4} is {:?}", BETA, max_count);
    println!("The minimum ratio of arc_length containing {} lattice points to R^beta is {:.4} for beta = {:.4}", N, min_coeffs, BETA);
    // println!("The corresponding squared radii are {:?}", max_rads);

    /* Uncomment if you want to output data to a CSV */
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

fn get_lattice_points(radius_squared: u64) -> Option<Vec<(f64,f64)>> {
    let radius = (radius_squared as f64).sqrt();
    let mut lattice_points: Vec<(f64, f64)> = Vec::new();

    for i in (((radius / (2 as f64).sqrt()).ceil() as u64)..=(radius.floor() as u64)).rev() {
        let i = i as f64;
        let cand = (radius.powi(2) - i.powi(2)).round();
        if cand == 0.0 {
            lattice_points.push((0.0, i));
        } else if VALID_RESIDUALS.contains(&((cand as u64) % 10)) {
            let temp = cand.sqrt().round();
            if temp.powi(2) == cand {
                lattice_points.push((temp, i));
            }
        }
    }
    let lp_len = lattice_points.len();
    if lp_len > 0 {
        for pt in lattice_points.clone().iter().rev() {
            lattice_points.push((pt.1,pt.0));
        }
        let mut fourth_quad_lp = Vec::new();

        for pt in lattice_points.iter().rev() {
            fourth_quad_lp.push((pt.0, -pt.1));
        }
        if radius.round().powi(2) == radius_squared as f64 {
            lattice_points.extend(fourth_quad_lp.split_first().unwrap().1);
        } else {
            lattice_points.extend(fourth_quad_lp)
        }
        // dbg!(lattice_points.clone());

        Some(lattice_points)
    } else {
        None
    }
}

fn get_arclength_coeffs(min_r: u64, max_r: u64) -> HashMap<u64, f64> {
    let mut coeffs: HashMap<u64,f64> = HashMap::new();

    for rad_squared in min_r..=max_r {
        let lattice_points = if let Some(res) = get_lattice_points(rad_squared) {
            res
        } else {
            // no lattice points, skip to next r^2 value
            coeffs.insert(rad_squared, f64::INFINITY);
            continue
        };

        let mut min_coeff: f64 = 0.0;
        for arc_v in lattice_points.windows(N) {
            let (x1, y1) = arc_v[0];
            let (x2, y2) = arc_v[N-1];
            let a1 = y1.atan2(x1);
            let a2 = y2.atan2(x2);
            let coeff = (a1 - a2) * (rad_squared as f64).powf(0.5*(1.0-BETA));
            if min_coeff == 0.0 || min_coeff > coeff {
                min_coeff = coeff;
            } else {
                continue
            }
        }
        coeffs.insert(rad_squared, min_coeff);
    }
    coeffs
}

fn get_latticepoint_counts<'a>(min_r: u64, max_r: u64) -> HashMap<u64, u64> {
    let mut rad_counts: HashMap<u64,u64> = HashMap::new();

    for radius_squared in min_r..=max_r {
        let radius = (radius_squared as f64).sqrt();

        let lattice_points = if let Some(res) = get_lattice_points(radius_squared) {
            res
        } else {
            // no lattice points, skip to next r^2 value
            rad_counts.insert(radius_squared, 0);
            continue
        };

        // compute angles between lattice points
        let mut angles = Vec::new();
        for i in 0..(lattice_points.len()-1) {
            let (x1, y1) = lattice_points[i];
            let (x2, y2) = lattice_points[i + 1];
            let a1 = y1.atan2(x1);
            let a2 = y2.atan2(x2);
            angles.push(a1 - a2)
        }

        let mut count: u64 = 0;

        let arc_length: f64 = 2.0_f64.sqrt() * radius.powf(BETA);
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
        rad_counts.insert(radius_squared, count);
    }
    rad_counts
}
