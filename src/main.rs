use jarnik::numerics::clusters::{
    cluster_points, find_max_clusters, geometry_for_cluster, ArcCluster, ClusterGeometry,
};
use jarnik::numerics::cuts::{analyze_cuts, enumerate_cuts, CutDiagnostics};
use jarnik::numerics::exact::{isqrt_u128, parse_i128, parse_u64_list};
use jarnik::numerics::factor::{
    analyze_factorization, analyze_factorization_cached, factor_u64, pow_u128, primes_up_to,
    represent_prime_as_sum_of_squares, FactorizationData, PrimeFactor, SplitPrime,
};
use jarnik::numerics::half_angle::{analyze_half_angles, HalfAngleDiagnostics};
use jarnik::numerics::legacy_arc::{
    find_keys_for_value, get_arclength_coeffs, BETA, MAX_R, MIN_R, N,
};
use jarnik::numerics::points::{generate_lattice_points_from_factorization, LatticePointRecord};
use jarnik::numerics::relation_height::{
    analyze_full_relation_height, CutTableRecord, FullRelationHeightDiagnostics,
};
use jarnik::numerics::transition_patterns::{
    analyze_transition_patterns, TransitionPatternDiagnostics,
};
use jarnik::{
    analyze_circle, find_extensions, search_sign_product_counterexamples, Certificate,
    CertificateSearch, ExtensionSearchConfig,
};
use serde::Serialize;
use std::collections::{BTreeMap, HashMap};
use std::convert::TryFrom;
use std::env;
use std::error;
use std::fs::File;
use std::io::{self, BufWriter, Write};

fn main() -> Result<(), Box<dyn error::Error>> {
    let args: Vec<String> = env::args().skip(1).collect();
    if args.is_empty() {
        print_help();
        return Ok(());
    }

    match args[0].as_str() {
        "arc" => run_arc(&args[1..])?,
        "search" => run_search(&args[1..])?,
        "analyze-n" => run_analyze_n(&args[1..])?,
        "cg-family" => run_cg_family(&args[1..])?,
        "random-squarefree" => run_random_squarefree(&args[1..])?,
        "known-cert" => run_known_cert(),
        "generate" => run_generate(&args[1..])?,
        "extend-known" => run_extend_known(&args[1..])?,
        "verify-circle" => run_verify_circle(&args[1..])?,
        "sign-search" => run_sign_search(&args[1..])?,
        "help" | "--help" | "-h" => print_help(),
        other => {
            eprintln!("unknown command `{}`\n", other);
            print_help();
            std::process::exit(2);
        }
    }

    Ok(())
}

#[derive(Clone, Debug, Serialize)]
struct ClusterAnalysisRecord {
    #[serde(flatten)]
    geometry: ClusterGeometry,
    factors: Vec<PrimeFactor>,
    split_primes: Vec<SplitPrime>,
    cluster_exponents: Vec<Vec<u32>>,
    cut_diagnostics: CutDiagnostics,
    full_relation_height: FullRelationHeightDiagnostics,
    full_cut_table: Vec<CutTableRecord>,
    transition_patterns: TransitionPatternDiagnostics,
    half_angle: HalfAngleDiagnostics,
}

#[derive(Clone, Debug, Serialize)]
struct PointExponentRecord {
    point: [i128; 2],
    angle: f64,
    exponents: Vec<u32>,
}

#[derive(Clone, Debug, Serialize)]
struct AnalyzeNOutput {
    n: u64,
    arc_scale: f64,
    point_count: usize,
    max_cluster_size: usize,
    factorization: FactorizationData,
    lattice_points: Vec<PointExponentRecord>,
    maximum_clusters: Vec<ClusterAnalysisRecord>,
}

#[derive(Clone, Debug, Serialize)]
struct TopExample {
    n: u64,
    cluster_size: usize,
    angular_width: f64,
    points: Vec<[i128; 2]>,
}

#[derive(Clone, Debug, Serialize)]
struct SearchSummary {
    record_type: String,
    processed_n: u64,
    num_circles_with_points: u64,
    max_cluster_seen: usize,
    top_examples: Vec<TopExample>,
}

#[derive(Clone, Debug, Serialize)]
struct CgFamilyRecord {
    family: String,
    family_index: u64,
    requested_index_start: u64,
    requested_index_end: u64,
    radius_squared: Option<u128>,
    radius: Option<f64>,
    generated_points: Vec<[i128; 2]>,
    analyzed: bool,
    skipped_reason: Option<String>,
    analysis: Option<AnalyzeNOutput>,
}

#[derive(Clone, Debug, Serialize)]
struct RandomSquarefreeRecord {
    trial: usize,
    seed: u64,
    selected_primes: Vec<u64>,
    n: u64,
    point_count: usize,
    max_cluster_size: usize,
    maximum_clusters: Vec<ClusterAnalysisRecord>,
}

fn run_search(args: &[String]) -> Result<(), Box<dyn error::Error>> {
    let n_max = required_flag_u64(args, "--n-max")?;
    let arc_scale = optional_flag_f64(args, "--arc-scale", 1.0)?;
    let min_cluster = optional_flag_usize(args, "--min-cluster", 4)?;
    let summary_every = optional_flag_u64(args, "--summary-every", 100_000)?;
    let out = optional_flag(args, "--out")?;
    let mut writer = output_writer(out)?;

    let primes = primes_up_to(isqrt_u128(n_max as u128) as u64 + 1);
    let mut representation_cache = HashMap::new();
    let mut num_circles_with_points = 0_u64;
    let mut max_cluster_seen = 0_usize;
    let mut top_examples = Vec::new();

    for n in 1..=n_max {
        let factorization = analyze_factorization_cached(n, &primes, &mut representation_cache);
        if !factorization.is_sum_of_two_squares {
            maybe_emit_summary(
                n,
                summary_every,
                num_circles_with_points,
                max_cluster_seen,
                &top_examples,
            )?;
            continue;
        }

        let points = generate_lattice_points_from_factorization(&factorization);
        if points.is_empty() {
            maybe_emit_summary(
                n,
                summary_every,
                num_circles_with_points,
                max_cluster_seen,
                &top_examples,
            )?;
            continue;
        }
        num_circles_with_points += 1;

        let clusters = find_max_clusters(&points, n, arc_scale);
        let best_size = clusters
            .iter()
            .map(|cluster| cluster.size())
            .max()
            .unwrap_or(0);
        max_cluster_seen = max_cluster_seen.max(best_size);
        if let Some(first_cluster) = clusters.first() {
            update_top_examples(
                &mut top_examples,
                top_example(n, first_cluster, &points),
                10,
            );
        }

        if best_size >= min_cluster {
            for cluster in clusters
                .iter()
                .filter(|cluster| cluster.size() == best_size)
            {
                let record = build_cluster_analysis(n, arc_scale, &factorization, &points, cluster);
                write_json_line(&mut writer, &record)?;
            }
        }

        maybe_emit_summary(
            n,
            summary_every,
            num_circles_with_points,
            max_cluster_seen,
            &top_examples,
        )?;
    }

    if summary_every == 0 || n_max % summary_every != 0 {
        emit_summary(
            n_max,
            num_circles_with_points,
            max_cluster_seen,
            &top_examples,
        )?;
    }
    writer.flush()?;
    Ok(())
}

fn run_analyze_n(args: &[String]) -> Result<(), Box<dyn error::Error>> {
    let n = required_flag_u64(args, "--n")?;
    let arc_scale = optional_flag_f64(args, "--arc-scale", 1.0)?;
    let out = optional_flag(args, "--out")?;
    let output = analyze_n_value(n, arc_scale)?;
    let mut writer = output_writer(out)?;
    serde_json::to_writer_pretty(&mut writer, &output)?;
    writeln!(writer)?;
    writer.flush()?;
    Ok(())
}

fn run_cg_family(args: &[String]) -> Result<(), Box<dyn error::Error>> {
    let n_start = optional_flag_u64(args, "--n-start", 2)?;
    let n_end = optional_flag_u64(args, "--n-end", 15)?;
    let arc_scale = optional_flag_f64(args, "--arc-scale", 1.0)?;
    let out = optional_flag(args, "--out")?;
    let mut writer = output_writer(out)?;

    if n_start < 2 {
        return Err("cg-family indices start at 2 because the formula uses F_{n-2}".into());
    }
    if n_end < n_start {
        return Err("--n-end must be at least --n-start".into());
    }

    let max_fib_index = n_end
        .checked_mul(3)
        .and_then(|value| value.checked_add(3))
        .ok_or("requested Fibonacci index overflowed u64")?;
    let fibs = fibonacci_numbers_u128(max_fib_index as usize)?;

    for family_index in n_start..=n_end {
        let (radius_squared, generated_points) =
            match fibonacci_family_circle(family_index as usize, &fibs) {
                Ok(value) => value,
                Err(reason) => {
                    let record = CgFamilyRecord {
                        family: "cilleruelo-granville-fibonacci".to_string(),
                        family_index,
                        requested_index_start: n_start,
                        requested_index_end: n_end,
                        radius_squared: None,
                        radius: None,
                        generated_points: Vec::new(),
                        analyzed: false,
                        skipped_reason: Some(reason),
                        analysis: None,
                    };
                    write_json_line(&mut writer, &record)?;
                    continue;
                }
            };

        let (analysis, skipped_reason) = if radius_squared > u64::MAX as u128 {
            (
                None,
                Some("radius_squared exceeds u64 analysis limit".to_string()),
            )
        } else {
            let factorization = fibonacci_family_factorization(
                family_index as usize,
                radius_squared as u64,
                &fibs,
            )?;
            (Some(analyze_factorized_n(factorization, arc_scale)), None)
        };

        let record = CgFamilyRecord {
            family: "cilleruelo-granville-fibonacci".to_string(),
            family_index,
            requested_index_start: n_start,
            requested_index_end: n_end,
            radius_squared: Some(radius_squared),
            radius: Some((radius_squared as f64).sqrt()),
            generated_points,
            analyzed: analysis.is_some(),
            skipped_reason,
            analysis,
        };
        write_json_line(&mut writer, &record)?;
    }
    writer.flush()?;
    Ok(())
}

fn run_random_squarefree(args: &[String]) -> Result<(), Box<dyn error::Error>> {
    let num_primes = required_flag_usize(args, "--num-primes")?;
    let trials = required_flag_usize(args, "--trials")?;
    let arc_scale = optional_flag_f64(args, "--arc-scale", 1.0)?;
    let seed = optional_flag_u64(args, "--seed", 0x4a52_4e49_4b_u64)?;
    let prime_limit = optional_flag_u64(args, "--prime-limit", 100)?;
    let out = optional_flag(args, "--out")?;
    let mut writer = output_writer(out)?;

    let prime_pool = primes_up_to(prime_limit)
        .into_iter()
        .filter(|prime| prime % 4 == 1)
        .collect::<Vec<_>>();
    if prime_pool.len() < num_primes {
        return Err(format!(
            "--prime-limit {} only gives {} primes congruent to 1 mod 4, need {}",
            prime_limit,
            prime_pool.len(),
            num_primes
        )
        .into());
    }

    let mut rng = SplitMix64::new(seed);
    let mut completed = 0_usize;
    let mut attempts = 0_usize;
    let max_attempts = trials.saturating_mul(10_000).max(10_000);
    while completed < trials {
        attempts += 1;
        if attempts > max_attempts {
            return Err(format!(
                "could not produce {} non-overflowing squarefree products after {} attempts; lower --num-primes or --prime-limit",
                trials, attempts
            )
            .into());
        }

        let selected_primes = sample_distinct(&prime_pool, num_primes, &mut rng);
        let factorization = match squarefree_factorization_from_primes(&selected_primes) {
            Some(factorization) => factorization,
            None => continue,
        };
        let analysis = analyze_factorized_n(factorization, arc_scale);
        completed += 1;
        let record = RandomSquarefreeRecord {
            trial: completed,
            seed,
            selected_primes,
            n: analysis.n,
            point_count: analysis.point_count,
            max_cluster_size: analysis.max_cluster_size,
            maximum_clusters: analysis.maximum_clusters,
        };
        write_json_line(&mut writer, &record)?;
    }

    writer.flush()?;
    Ok(())
}

fn analyze_n_value(n: u64, arc_scale: f64) -> Result<AnalyzeNOutput, Box<dyn error::Error>> {
    let primes = primes_up_to(isqrt_u128(n as u128) as u64 + 1);
    let factorization = analyze_factorization(n, &primes);
    Ok(analyze_factorized_n(factorization, arc_scale))
}

fn analyze_factorized_n(factorization: FactorizationData, arc_scale: f64) -> AnalyzeNOutput {
    let n = factorization.n;
    let points = generate_lattice_points_from_factorization(&factorization);
    let clusters = find_max_clusters(&points, n, arc_scale);
    let max_cluster_size = clusters
        .iter()
        .map(|cluster| cluster.size())
        .max()
        .unwrap_or(0);
    let maximum_clusters = clusters
        .iter()
        .filter(|cluster| cluster.size() == max_cluster_size)
        .map(|cluster| build_cluster_analysis(n, arc_scale, &factorization, &points, cluster))
        .collect::<Vec<_>>();
    let lattice_points = points
        .iter()
        .map(|point| PointExponentRecord {
            point: point.point_array(),
            angle: point.angle,
            exponents: point.exponents.clone(),
        })
        .collect::<Vec<_>>();

    AnalyzeNOutput {
        n,
        arc_scale,
        point_count: points.len(),
        max_cluster_size,
        factorization,
        lattice_points,
        maximum_clusters,
    }
}

fn build_cluster_analysis(
    n: u64,
    arc_scale: f64,
    factorization: &FactorizationData,
    points: &[LatticePointRecord],
    cluster: &ArcCluster,
) -> ClusterAnalysisRecord {
    let cluster_refs = cluster_points(cluster, points);
    let cuts = enumerate_cuts(&factorization.split_primes, &cluster_refs);
    let cut_diagnostics = analyze_cuts(&cuts, cluster.size());
    let log_r = 0.5 * (n as f64).ln();
    let full_relation_height = analyze_full_relation_height(
        &cuts,
        cluster.size(),
        cut_diagnostics.total_cut_weight,
        log_r,
    );
    let transition_patterns =
        analyze_transition_patterns(&cuts, cluster.size(), cut_diagnostics.total_cut_weight);
    let half_angle = analyze_half_angles(n, &cluster_refs);
    let cluster_exponents = cluster_refs
        .iter()
        .map(|point| point.exponents.clone())
        .collect::<Vec<_>>();

    ClusterAnalysisRecord {
        geometry: geometry_for_cluster(n, arc_scale, cluster, points),
        factors: factorization.factors.clone(),
        split_primes: factorization.split_primes.clone(),
        cluster_exponents,
        cut_diagnostics,
        full_relation_height,
        full_cut_table: if cluster.size() >= 4 {
            cuts.iter()
                .map(|cut| CutTableRecord {
                    prime: cut.prime,
                    threshold: cut.threshold,
                    weight: cut.weight,
                    pattern: cut.pattern.clone(),
                    variation: cut.variation,
                })
                .collect()
        } else {
            Vec::new()
        },
        transition_patterns,
        half_angle,
    }
}

fn top_example(n: u64, cluster: &ArcCluster, points: &[LatticePointRecord]) -> TopExample {
    TopExample {
        n,
        cluster_size: cluster.size(),
        angular_width: cluster.angular_width,
        points: cluster
            .indices
            .iter()
            .map(|&idx| points[idx].point_array())
            .collect(),
    }
}

fn update_top_examples(top_examples: &mut Vec<TopExample>, example: TopExample, limit: usize) {
    top_examples.push(example);
    top_examples.sort_by(|a, b| {
        b.cluster_size
            .cmp(&a.cluster_size)
            .then_with(|| a.angular_width.partial_cmp(&b.angular_width).unwrap())
            .then_with(|| a.n.cmp(&b.n))
    });
    top_examples.truncate(limit);
}

fn maybe_emit_summary(
    processed_n: u64,
    summary_every: u64,
    num_circles_with_points: u64,
    max_cluster_seen: usize,
    top_examples: &[TopExample],
) -> Result<(), Box<dyn error::Error>> {
    if summary_every != 0 && processed_n % summary_every == 0 {
        emit_summary(
            processed_n,
            num_circles_with_points,
            max_cluster_seen,
            top_examples,
        )?;
    }
    Ok(())
}

fn emit_summary(
    processed_n: u64,
    num_circles_with_points: u64,
    max_cluster_seen: usize,
    top_examples: &[TopExample],
) -> Result<(), Box<dyn error::Error>> {
    let summary = SearchSummary {
        record_type: "summary".to_string(),
        processed_n,
        num_circles_with_points,
        max_cluster_seen,
        top_examples: top_examples.to_vec(),
    };
    eprintln!("{}", serde_json::to_string(&summary)?);
    Ok(())
}

fn squarefree_factorization_from_primes(selected_primes: &[u64]) -> Option<FactorizationData> {
    let mut sorted = selected_primes.to_vec();
    sorted.sort_unstable();
    let mut n = 1_u64;
    let mut factors = Vec::new();
    let mut split_primes = Vec::new();
    for prime in sorted {
        n = n.checked_mul(prime)?;
        let (a, b) = represent_prime_as_sum_of_squares(prime)?;
        factors.push(PrimeFactor { prime, exponent: 1 });
        split_primes.push(SplitPrime {
            prime,
            exponent: 1,
            a: a as i128,
            b: b as i128,
        });
    }
    Some(FactorizationData {
        n,
        factors,
        split_primes,
        two_exponent: 0,
        inert_scale: 1,
        is_sum_of_two_squares: true,
    })
}

fn fibonacci_numbers_u128(max_index: usize) -> Result<Vec<u128>, Box<dyn error::Error>> {
    let mut fibs = vec![0_u128; max_index + 1];
    if max_index >= 1 {
        fibs[1] = 1;
    }
    for index in 2..=max_index {
        fibs[index] = fibs[index - 1]
            .checked_add(fibs[index - 2])
            .ok_or_else(|| format!("Fibonacci number F_{} exceeds u128", index))?;
    }
    Ok(fibs)
}

fn fibonacci_family_circle(
    family_index: usize,
    fibs: &[u128],
) -> Result<(u128, Vec<[i128; 2]>), String> {
    if family_index < 2 {
        return Err("family index must be at least 2".to_string());
    }
    let f = |idx: usize| -> Result<u128, String> {
        fibs.get(idx)
            .copied()
            .ok_or_else(|| format!("missing Fibonacci number F_{}", idx))
    };

    let product = f(2 * family_index - 1)?
        .checked_mul(f(2 * family_index + 1)?)
        .ok_or_else(|| "radius_squared overflowed u128".to_string())?
        .checked_mul(f(2 * family_index + 3)?)
        .ok_or_else(|| "radius_squared overflowed u128".to_string())?
        .checked_mul(5)
        .ok_or_else(|| "radius_squared overflowed u128".to_string())?;
    if product % 2 != 0 {
        return Err("unexpected odd Fibonacci family radius numerator".to_string());
    }
    let radius_squared = product / 2;

    let half_x = f(3 * family_index + 3)?;
    let half_y = f(3 * family_index)?;
    if half_x % 2 != 0 || half_y % 2 != 0 {
        return Err("unexpected odd base point coordinate numerator".to_string());
    }
    let base_x = u128_to_i128(half_x / 2)?;
    let base_y = u128_to_i128(half_y / 2)?;
    let sign = if family_index % 2 == 0 {
        1_i128
    } else {
        -1_i128
    };

    let shifts = [
        (
            -2 * u128_to_i128(f(family_index - 1)?)?,
            2 * u128_to_i128(f(family_index + 2)?)?,
        ),
        (
            -u128_to_i128(f(family_index - 2)?)?,
            u128_to_i128(f(family_index + 1)?)?,
        ),
        (
            u128_to_i128(f(family_index - 1)?)?,
            -u128_to_i128(f(family_index + 2)?)?,
        ),
        (
            u128_to_i128(f(family_index)?)?,
            -u128_to_i128(f(family_index + 3)?)?,
        ),
    ];

    let mut points = Vec::new();
    for (dx, dy) in shifts {
        let x = base_x
            .checked_add(sign * dx)
            .ok_or_else(|| "x coordinate overflowed i128".to_string())?;
        let y = base_y
            .checked_add(sign * dy)
            .ok_or_else(|| "y coordinate overflowed i128".to_string())?;
        let x_abs = x.unsigned_abs();
        let y_abs = y.unsigned_abs();
        let norm = x_abs
            .checked_mul(x_abs)
            .and_then(|value| value.checked_add(y_abs.checked_mul(y_abs)?))
            .ok_or_else(|| "point norm overflowed u128".to_string())?;
        if norm != radius_squared {
            return Err(format!(
                "generated point ({},{}) has norm {}, expected {}",
                x, y, norm, radius_squared
            ));
        }
        points.push([x, y]);
    }

    Ok((radius_squared, points))
}

fn fibonacci_family_factorization(
    family_index: usize,
    radius_squared: u64,
    fibs: &[u128],
) -> Result<FactorizationData, String> {
    let f = |idx: usize| -> Result<u128, String> {
        fibs.get(idx)
            .copied()
            .ok_or_else(|| format!("missing Fibonacci number F_{}", idx))
    };

    let mut components = vec![
        5_u128,
        f(2 * family_index - 1)?,
        f(2 * family_index + 1)?,
        f(2 * family_index + 3)?,
    ];
    let mut divided_by_two = false;
    for component in &mut components {
        if *component % 2 == 0 {
            *component /= 2;
            divided_by_two = true;
            break;
        }
    }
    if !divided_by_two {
        return Err("could not divide Fibonacci family product by 2".to_string());
    }

    let mut check_product = 1_u128;
    let mut factor_map: BTreeMap<u64, u32> = BTreeMap::new();
    for component in components {
        check_product = check_product
            .checked_mul(component)
            .ok_or_else(|| "component product overflowed u128".to_string())?;
        let component_u64 = u64::try_from(component)
            .map_err(|_| "Fibonacci component exceeds u64 factorization limit".to_string())?;
        for factor in factor_u64(component_u64) {
            *factor_map.entry(factor.prime).or_insert(0) += factor.exponent;
        }
    }
    if check_product != radius_squared as u128 {
        return Err(format!(
            "component product {} did not match radius_squared {}",
            check_product, radius_squared
        ));
    }

    factorization_from_factor_map(radius_squared, factor_map)
}

fn factorization_from_factor_map(
    n: u64,
    factor_map: BTreeMap<u64, u32>,
) -> Result<FactorizationData, String> {
    let mut factors = Vec::new();
    let mut split_primes = Vec::new();
    let mut two_exponent = 0_u32;
    let mut inert_scale = 1_u128;
    let mut is_sum_of_two_squares = true;

    for (prime, exponent) in factor_map {
        factors.push(PrimeFactor { prime, exponent });
        match prime % 4 {
            1 => {
                if let Some((a, b)) = represent_prime_as_sum_of_squares(prime) {
                    split_primes.push(SplitPrime {
                        prime,
                        exponent,
                        a: a as i128,
                        b: b as i128,
                    });
                } else {
                    is_sum_of_two_squares = false;
                }
            }
            2 => {
                two_exponent = exponent;
            }
            3 => {
                if exponent % 2 == 1 {
                    is_sum_of_two_squares = false;
                } else {
                    inert_scale *= pow_u128(prime as u128, exponent / 2);
                }
            }
            _ => unreachable!(),
        }
    }

    Ok(FactorizationData {
        n,
        factors,
        split_primes,
        two_exponent,
        inert_scale,
        is_sum_of_two_squares,
    })
}

fn u128_to_i128(value: u128) -> Result<i128, String> {
    if value > i128::MAX as u128 {
        Err("value exceeds i128".to_string())
    } else {
        Ok(value as i128)
    }
}

fn sample_distinct(pool: &[u64], count: usize, rng: &mut SplitMix64) -> Vec<u64> {
    let mut values = pool.to_vec();
    for i in 0..count {
        let j = i + rng.next_usize(values.len() - i);
        values.swap(i, j);
    }
    let mut selected = values[..count].to_vec();
    selected.sort_unstable();
    selected
}

struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    fn next_usize(&mut self, upper: usize) -> usize {
        if upper <= 1 {
            return 0;
        }
        (self.next_u64() as usize) % upper
    }
}

fn output_writer(path: Option<&str>) -> Result<Box<dyn Write>, Box<dyn error::Error>> {
    match path {
        Some(path) if path != "-" => Ok(Box::new(BufWriter::new(File::create(path)?))),
        _ => Ok(Box::new(BufWriter::new(io::stdout()))),
    }
}

fn write_json_line<W: Write + ?Sized, T: Serialize>(
    writer: &mut W,
    value: &T,
) -> Result<(), Box<dyn error::Error>> {
    serde_json::to_writer(&mut *writer, value)?;
    writeln!(writer)?;
    Ok(())
}

fn optional_flag<'a>(args: &'a [String], flag: &str) -> Result<Option<&'a str>, String> {
    let equals_prefix = format!("{}=", flag);
    for (index, arg) in args.iter().enumerate() {
        if arg == flag {
            let value = args
                .get(index + 1)
                .ok_or_else(|| format!("{} requires a value", flag))?;
            if value.starts_with("--") {
                return Err(format!("{} requires a value", flag));
            }
            return Ok(Some(value.as_str()));
        }
        if arg.starts_with(&equals_prefix) {
            return Ok(Some(&arg[equals_prefix.len()..]));
        }
    }
    Ok(None)
}

fn required_flag<'a>(args: &'a [String], flag: &str) -> Result<&'a str, String> {
    optional_flag(args, flag)?.ok_or_else(|| format!("missing required flag {}", flag))
}

fn required_flag_u64(args: &[String], flag: &str) -> Result<u64, String> {
    required_flag(args, flag)?
        .parse::<u64>()
        .map_err(|err| format!("could not parse {} as u64: {}", flag, err))
}

fn required_flag_usize(args: &[String], flag: &str) -> Result<usize, String> {
    required_flag(args, flag)?
        .parse::<usize>()
        .map_err(|err| format!("could not parse {} as usize: {}", flag, err))
}

fn optional_flag_u64(args: &[String], flag: &str, default: u64) -> Result<u64, String> {
    optional_flag(args, flag)?
        .map(|value| {
            value
                .parse::<u64>()
                .map_err(|err| format!("could not parse {} as u64: {}", flag, err))
        })
        .unwrap_or(Ok(default))
}

fn optional_flag_usize(args: &[String], flag: &str, default: usize) -> Result<usize, String> {
    optional_flag(args, flag)?
        .map(|value| {
            value
                .parse::<usize>()
                .map_err(|err| format!("could not parse {} as usize: {}", flag, err))
        })
        .unwrap_or(Ok(default))
}

fn optional_flag_f64(args: &[String], flag: &str, default: f64) -> Result<f64, String> {
    optional_flag(args, flag)?
        .map(|value| {
            value
                .parse::<f64>()
                .map_err(|err| format!("could not parse {} as f64: {}", flag, err))
        })
        .unwrap_or(Ok(default))
}

fn run_arc(args: &[String]) -> Result<(), Box<dyn error::Error>> {
    let min_r = optional_i128(args, 0)?.unwrap_or(MIN_R as i128) as u64;
    let max_r = optional_i128(args, 1)?.unwrap_or(MAX_R as i128) as u64;
    let lp_coeffs = get_arclength_coeffs(min_r, max_r);
    let min_coeffs = lp_coeffs.values().fold(f64::INFINITY, |a, &b| a.min(b));

    println!("From radius-squared = {} to {}", min_r, max_r);
    println!(
        "The minimum ratio of arc_length containing {} lattice points to R^beta is {} for beta = {:.4}",
        N, min_coeffs, BETA
    );
    let mut max_rads = find_keys_for_value(&lp_coeffs, &min_coeffs);
    max_rads.sort();
    println!("The corresponding squared radii are {:?}", max_rads);
    Ok(())
}

fn run_known_cert() {
    let cert = Certificate::known_first_endpoint_certificate();
    println!("{}", cert.to_json_line());
    println!("radius_squared={}", cert.radius_squared());
    println!("n2={} n3={}", cert.n2(), cert.n3());
}

fn run_generate(args: &[String]) -> Result<(), Box<dyn error::Error>> {
    if args.len() < 2 {
        return Err("usage: cargo run -- generate <max_z> <max_det> [limit]".into());
    }

    let max_z = parse_i128(&args[0])?;
    let max_det = parse_i128(&args[1])?;
    let limit = optional_i128(args, 2)?.map(|value| value as usize);
    let report = CertificateSearch {
        min_z: 3,
        max_z,
        max_det,
        limit,
    }
    .run();

    eprintln!(
        "checked_blocks={} congruence_hits={} certificates={}",
        report.checked_blocks,
        report.congruence_hits,
        report.certificates.len()
    );
    for cert in report.certificates {
        println!("{}", cert.to_json_line());
    }
    Ok(())
}

fn run_extend_known(args: &[String]) -> Result<(), Box<dyn error::Error>> {
    let max_w = optional_i128(args, 0)?;
    let primes = if args.len() > 1 {
        parse_u64_list(&args[1])?
    } else {
        ExtensionSearchConfig::default().primes
    };
    let config = ExtensionSearchConfig {
        max_w,
        primes,
        ..ExtensionSearchConfig::default()
    };
    let cert = Certificate::known_first_endpoint_certificate();
    let report = find_extensions(&cert, &config)?;

    println!("upper_w={}", report.upper_w);
    println!("checked_w={}", report.checked_w);
    if let Some(sieve) = report.sieve {
        println!(
            "sieve_modulus={} sieve_residues={} selectivity={:.8}",
            sieve.modulus,
            sieve.residues.len(),
            sieve.selectivity()
        );
    }
    println!("witnesses={}", report.witnesses.len());
    for witness in report.witnesses {
        println!("{}", witness.to_json_line());
    }
    Ok(())
}

fn run_verify_circle(args: &[String]) -> Result<(), Box<dyn error::Error>> {
    if args.is_empty() {
        return Err("usage: cargo run -- verify-circle <radius_squared>".into());
    }

    let radius_squared = parse_i128(&args[0])? as u128;
    let report = analyze_circle(radius_squared);
    println!("radius_squared={}", report.radius_squared);
    println!("lattice_points={}", report.lattice_points.len());
    println!(
        "max_endpoint_chord_cluster={} start_index={}",
        report.max_chord_cluster.count, report.max_chord_cluster.start_index
    );
    println!("max_arc_count_float={}", report.max_arc_count_float);
    for point in report.max_chord_cluster.points {
        println!("({}, {})", point.x, point.y);
    }
    Ok(())
}

fn run_sign_search(args: &[String]) -> Result<(), Box<dyn error::Error>> {
    if args.len() < 4 {
        return Err(
            "usage: cargo run -- sign-search <a> <n> <max_r> <max_t> [min_vectors] [max_results]"
                .into(),
        );
    }

    let a = parse_i128(&args[0])?;
    let n = parse_i128(&args[1])? as usize;
    let max_r = parse_i128(&args[2])?;
    let max_t = parse_i128(&args[3])?;
    let min_vectors = optional_i128(args, 4)?.unwrap_or(5) as usize;
    let max_results = optional_i128(args, 5)?.unwrap_or(10) as usize;

    let results = search_sign_product_counterexamples(a, n, max_r, max_t, min_vectors, max_results);
    println!("results={}", results.len());
    for (factors, collision) in results {
        let factors_text = factors
            .iter()
            .map(|factor| format!("({}, {})", factor.r, factor.t))
            .collect::<Vec<_>>()
            .join(", ");
        println!(
            "factors=[{}] coeffs={:?} masks={:?} products={:?}",
            factors_text, collision.coeffs, collision.masks, collision.products
        );
    }
    Ok(())
}

fn optional_i128(args: &[String], index: usize) -> Result<Option<i128>, String> {
    if let Some(value) = args.get(index) {
        parse_i128(value).map(Some)
    } else {
        Ok(None)
    }
}

fn print_help() {
    println!(
        "Jarnik Gaussian-chain tools

Commands:
  search --n-max N [--arc-scale C] [--min-cluster M] [--out path.jsonl]
  analyze-n --n N [--arc-scale C] [--out path.json]
  cg-family [--n-start A] [--n-end B] [--arc-scale C] [--out path.jsonl]
  random-squarefree --num-primes K --trials T [--arc-scale C] [--seed S] [--prime-limit P] [--out path.jsonl]
  arc [min_radius_squared] [max_radius_squared]
  known-cert
  generate <max_z> <max_det> [limit]
  extend-known [max_w] [comma_primes]
  verify-circle <radius_squared>
  sign-search <a> <n> <max_r> <max_t> [min_vectors] [max_results]

Examples:
  cargo run --release -- search --n-max 100000000 --arc-scale 1.0 --min-cluster 4 --out results.jsonl
  cargo run --release -- analyze-n --n 567454025 --arc-scale 1.0 --out example.json
  cargo run --release -- cg-family --n-start 7 --n-end 15 --arc-scale 1.0 --out cg.jsonl
  cargo run --release -- random-squarefree --num-primes 10 --trials 1000 --arc-scale 1.0 --out random.jsonl
  cargo run -- known-cert
  cargo run -- generate 1597 6
  cargo run -- extend-known
  cargo run -- verify-circle 567454025"
    );
}
