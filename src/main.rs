use jarnik::{
    get_arclength_coeffs,
    find_keys_for_value,
    // get_latticepoint_counts,
    BETA, MIN_R, MAX_R, N,
};
use std::error;
// use csv::Writer;

fn main() -> Result<(), Box<dyn error::Error>>{
    // let lp_counts = get_latticepoint_counts(MIN_R, MAX_R);
    let lp_coeffs = get_arclength_coeffs(MIN_R, MAX_R);

    // let max_count = lp_counts.values().max().unwrap();
    let min_coeffs = lp_coeffs.values().fold(f64::INFINITY, |a, &b| a.min(b));


    println!("From radius = {} to {}", MIN_R, MAX_R);
    // println!("The maximum number of lattice points for beta = {:.4} is {:?}", BETA, max_count);
    // let mut max_rads = find_keys_for_value(&res,max_count);
    // max_rads.sort();
    // println!("The corresponding squared radii are {:?}", max_rads);

    println!("The minimum ratio of arc_length containing {} lattice points to R^beta is {:} for beta = {:.4}", N, min_coeffs, BETA);
    let mut max_rads = find_keys_for_value(&lp_coeffs,&min_coeffs);
    max_rads.sort();
    println!("The corresponding squared radii are {:?}", max_rads);

    /* Uncomment if you want to output data to a CSV */
    // let mut wtr = Writer::from_path("jarnik_arc_coeffs_5lp_1to1e8.csv")?;
    // wtr.write_record(&["Radius", "Ratio"])?;
    // for (k,v) in lp_coeffs.iter() {
    //     // if not enough lattice points no need to save data
    //     // the "inf" comes from the default value we assign
    //     // to the arc-length ratio
    //     if !v.is_infinite() {
    //         wtr.serialize(&(k, v))?;
    //     }
    // }
    // wtr.flush()?;
    Ok(())
}


