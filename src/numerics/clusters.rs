use crate::numerics::points::LatticePointRecord;
use serde::Serialize;

#[derive(Clone, Debug)]
pub struct ArcCluster {
    pub start_index: usize,
    pub indices: Vec<usize>,
    pub angles: Vec<f64>,
    pub angular_width: f64,
}

impl ArcCluster {
    pub fn size(&self) -> usize {
        self.indices.len()
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct ClusterGeometry {
    pub n: u64,
    pub r: f64,
    pub arc_scale: f64,
    pub angular_window: f64,
    pub cluster_size: usize,
    pub arc_length: f64,
    pub arc_length_over_r_sqrt: f64,
    pub angular_width: f64,
    pub start_index: usize,
    pub points: Vec<[i128; 2]>,
    pub angles: Vec<f64>,
}

pub fn angular_window(n: u64, arc_scale: f64) -> f64 {
    arc_scale * (n as f64).powf(-0.25)
}

pub fn arc_length(n: u64, arc_scale: f64) -> f64 {
    arc_scale * (n as f64).powf(0.25)
}

pub fn find_max_clusters(points: &[LatticePointRecord], n: u64, arc_scale: f64) -> Vec<ArcCluster> {
    if points.is_empty() {
        return Vec::new();
    }

    let window = angular_window(n, arc_scale);
    let len = points.len();
    let mut doubled_angles = Vec::with_capacity(2 * len);
    doubled_angles.extend(points.iter().map(|point| point.angle));
    doubled_angles.extend(
        points
            .iter()
            .map(|point| point.angle + std::f64::consts::TAU),
    );

    let mut clusters = Vec::new();
    let mut best = 0_usize;
    let mut end = 0_usize;
    for start in 0..len {
        if end < start {
            end = start;
        }
        while end + 1 < start + len
            && doubled_angles[end + 1] - doubled_angles[start] <= window + 1e-14
        {
            end += 1;
        }
        let size = end - start + 1;
        if size < best {
            continue;
        }
        let indices = (start..=end).map(|idx| idx % len).collect::<Vec<_>>();
        let angles = doubled_angles[start..=end].to_vec();
        let cluster = ArcCluster {
            start_index: start,
            indices,
            angles,
            angular_width: doubled_angles[end] - doubled_angles[start],
        };
        if size > best {
            best = size;
            clusters.clear();
        }
        clusters.push(cluster);
    }
    clusters
}

pub fn geometry_for_cluster(
    n: u64,
    arc_scale: f64,
    cluster: &ArcCluster,
    points: &[LatticePointRecord],
) -> ClusterGeometry {
    let r = (n as f64).sqrt();
    let length = arc_length(n, arc_scale);
    ClusterGeometry {
        n,
        r,
        arc_scale,
        angular_window: angular_window(n, arc_scale),
        cluster_size: cluster.size(),
        arc_length: length,
        arc_length_over_r_sqrt: length / r.sqrt(),
        angular_width: cluster.angular_width,
        start_index: cluster.start_index,
        points: cluster
            .indices
            .iter()
            .map(|&idx| points[idx].point_array())
            .collect(),
        angles: cluster.angles.clone(),
    }
}

pub fn cluster_points<'a>(
    cluster: &ArcCluster,
    points: &'a [LatticePointRecord],
) -> Vec<&'a LatticePointRecord> {
    cluster.indices.iter().map(|&idx| &points[idx]).collect()
}

#[cfg(test)]
mod tests {
    use crate::numerics::factor::{analyze_factorization, primes_up_to};
    use crate::numerics::points::generate_lattice_points_from_factorization;

    use super::find_max_clusters;

    #[test]
    fn finds_known_cilleruelo_granville_cluster() {
        let n = 567_454_025;
        let primes = primes_up_to((n as f64).sqrt() as u64 + 1);
        let factorization = analyze_factorization(n, &primes);
        let points = generate_lattice_points_from_factorization(&factorization);
        let clusters = find_max_clusters(&points, n, 1.0);
        let best = clusters
            .iter()
            .map(|cluster| cluster.size())
            .max()
            .unwrap_or(0);
        assert_eq!(best, 4);

        let known = [
            (23200_i128, 5405_i128),
            (23189_i128, 5452_i128),
            (23176_i128, 5507_i128),
            (23171_i128, 5528_i128),
        ];
        assert!(clusters.iter().any(|cluster| {
            let coords = cluster
                .indices
                .iter()
                .map(|&idx| (points[idx].x, points[idx].y))
                .collect::<Vec<_>>();
            coords == known
        }));
    }

    #[test]
    fn finds_known_three_point_square_radius_cluster() {
        let n = 142_129;
        let primes = primes_up_to((n as f64).sqrt() as u64 + 1);
        let factorization = analyze_factorization(n, &primes);
        let points = generate_lattice_points_from_factorization(&factorization);
        let clusters = find_max_clusters(&points, n, 1.0);
        let best = clusters
            .iter()
            .map(|cluster| cluster.size())
            .max()
            .unwrap_or(0);
        assert_eq!(best, 3);

        let known = [
            (352_i128, 135_i128),
            (348_i128, 145_i128),
            (345_i128, 152_i128),
        ];
        assert!(clusters.iter().any(|cluster| {
            let coords = cluster
                .indices
                .iter()
                .map(|&idx| (points[idx].x, points[idx].y))
                .collect::<Vec<_>>();
            coords == known
        }));
    }
}
