use std::collections::HashMap;
use std::vec::Vec;

pub const BETA: f64 = 1.0 / 3.0;
pub const N: usize = 4;
pub const MIN_R: u64 = 5;
pub const MAX_R: u64 = 5;

pub fn find_keys_for_value<'a, K: std::cmp::PartialEq>(
    map: &'a HashMap<u64, K>,
    value: &K,
) -> Vec<&'a u64> {
    map.iter()
        .filter_map(|(key, val)| if val == value { Some(key) } else { None })
        .collect()
}

fn get_lattice_points(radius_squared: u64) -> Option<Vec<(f64, f64)>> {
    if radius_squared % 4 == 3 {
        return None;
    }
    let radius = (radius_squared as f64).sqrt();
    let mut lattice_points: Vec<(f64, f64)> = Vec::new();

    for i in (((radius / (2 as f64).sqrt()).ceil() as u64)..=(radius.floor() as u64)).rev() {
        let i = i as f64;
        let cand = (radius.powi(2) - i.powi(2)).round();
        if cand == 0.0 {
            lattice_points.push((0.0, i));
        } else if (cand as u64) % 10 == 1
            || (cand as u64) % 10 == 4
            || (cand as u64) % 10 == 5
            || (cand as u64) % 10 == 6
            || (cand as u64) % 10 == 9
        {
            let temp = cand.sqrt().round();
            if temp.powi(2) == cand {
                lattice_points.push((temp, i));
            }
        }
    }
    let lp_len = lattice_points.len();
    if lp_len > 0 {
        for pt in lattice_points.clone().iter().rev() {
            lattice_points.push((pt.1, pt.0));
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

pub fn get_arclength_coeffs(min_r: u64, max_r: u64) -> HashMap<u64, f64> {
    let mut coeffs: HashMap<u64, f64> = HashMap::new();

    for rad_squared in min_r..=max_r {
        let lattice_points = if let Some(res) = get_lattice_points(rad_squared) {
            res
        } else {
            // no lattice points, skip to next r^2 value
            continue;
        };
        // dbg!(lattice_points.clone());

        let mut min_coeff: f64 = f64::INFINITY;
        for arc_v in lattice_points.windows(N) {
            // dbg!(arc_v);
            let (x1, y1) = arc_v[0];
            let (x2, y2) = arc_v[N - 1];
            let a1 = y1.atan2(x1);
            let a2 = y2.atan2(x2);
            let coeff = (a1 - a2) * (rad_squared as f64).powf(0.5 * (1.0 - BETA));
            if min_coeff > coeff {
                min_coeff = coeff;
            } else {
                continue;
            }
        }
        coeffs.insert(rad_squared, min_coeff);
    }
    coeffs
}

pub fn get_latticepoint_counts<'a>(min_r: u64, max_r: u64) -> HashMap<u64, u64> {
    let mut rad_counts: HashMap<u64, u64> = HashMap::new();

    for radius_squared in min_r..=max_r {
        let radius = (radius_squared as f64).sqrt();

        let lattice_points = if let Some(res) = get_lattice_points(radius_squared) {
            res
        } else {
            // no lattice points, skip to next r^2 value
            rad_counts.insert(radius_squared, 0);
            continue;
        };

        // compute angles between lattice points
        let mut angles = Vec::new();
        for i in 0..(lattice_points.len() - 1) {
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
                    break;
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

#[cfg(test)]
mod tests {
    use crate::get_arclength_coeffs;

    #[test]
    fn arc_ratio() {
        // comparing values for degenerate 4-tuples with the paper
        // "Close Lattice Points on Circles" by Cilleruelo and Granville
        let temp = get_arclength_coeffs(5, 5);
        let ratio = temp.get(&5_u64).unwrap();
        let abs_dif = ratio - 3.786395353643682;
        assert!(abs_dif < 1e-10);

        let temp = get_arclength_coeffs(65, 65);
        let ratio = temp.get(&65_u64).unwrap();
        let abs_dif = ratio - 4.174688308044826;
        assert!(abs_dif < 1e-10);
    }
}
