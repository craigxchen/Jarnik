use std::collections::HashMap;
use std::vec::Vec;
use std::f64::consts::PI;

const BETA: f64 = 1.0/3.0;
const MIN_R: u64 = 10;
const MAX_R: u64 = 10;
const VALID_RESIDUALS: [u64; 5] = [1,4,5,6,9];

fn main() {
   let res = get_max_lattice_points(MIN_R, MAX_R);

    println!("The maximum number of lattice points for beta = {} is {}", &BETA, &res)
}

fn get_max_lattice_points<'a>(min_r: u64, max_r: u64) -> u64 {
    let mut rad_counts: HashMap<u64,u64> = HashMap::new();

    for radius in min_r..max_r+1 {
        let radius = radius as f64;

        let mut lattice_points: Vec<(f64, f64)> = Vec::new();
        let mut angles: Vec<f64> = Vec::new();

        lattice_points.push((0.0, radius));

        for i in (((radius / (2 as f64).sqrt()).ceil() as u64)..(radius as u64)).rev() {
            let i = i as f64;
            let cand = radius.powf(2.0) - i.powf(2.0);
            if VALID_RESIDUALS.contains(&((cand as u64) % 10)) {
                let temp = cand.sqrt();
                if temp.fract() == 0.0 {
                    lattice_points.push((temp, i));
                }
            }
        }
        let lp_len = lattice_points.len();
        // println!("{:?}",lattice_points);
        for i in 0..(lp_len-1) {
            let (x1, y1) = lattice_points[i];
            let (x2, y2) = lattice_points[i+1];
            let a1 = y1.atan2(x1);
            let a2 = y2.atan2(x2);
            angles.push(a1 - a2)
        }
        angles.push( 2.0 * lattice_points[lp_len-1].1.atan2(lattice_points[lp_len-1].0) - PI/2.0 );
        let mut rev_angles:Vec<f64> = Vec::new();
        rev_angles.extend(angles.clone().iter().rev());
        angles.extend_from_slice(&rev_angles[1..rev_angles.len()]);
        angles.extend(&mut angles.clone().iter());

        let mut count: u64 = 0;

        let arc_length: f64 = (2 as f64).sqrt() * radius.powf(BETA);
        let max_ang: f64 = arc_length / radius;

        for i in 0..angles.len() {
            let mut running_ang_sum: f64 = 0.0;
            let mut running_count = 1;

            for ang in &angles[i..] {
                running_ang_sum += *ang;
                if running_ang_sum > max_ang {
                    break
                }
                running_count += 1;
            }
            if running_count > count {
                count = running_count;
            }
        }
        rad_counts.insert(radius as u64, count);
    }
    let max_count = rad_counts.values().max().cloned();
    match max_count {
        Some(x) => x,
        None => unimplemented!(),
    }
}
