use crate::numerics::exact::{cmp_angle, isqrt_u128};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Point {
    pub x: i128,
    pub y: i128,
}

impl Point {
    pub fn norm_squared(&self) -> i128 {
        self.x * self.x + self.y * self.y
    }
}

#[derive(Clone, Debug)]
pub struct CircleCluster {
    pub count: usize,
    pub start_index: usize,
    pub points: Vec<Point>,
}

#[derive(Clone, Debug)]
pub struct CircleReport {
    pub radius_squared: u128,
    pub lattice_points: Vec<Point>,
    pub max_chord_cluster: CircleCluster,
    pub max_arc_count_float: usize,
}

pub fn lattice_points_on_circle(radius_squared: u128) -> Vec<Point> {
    let radius_floor = isqrt_u128(radius_squared) as i128;
    let mut points = Vec::new();

    for x in -radius_floor..=radius_floor {
        let x_squared = (x * x) as u128;
        if x_squared > radius_squared {
            continue;
        }
        let y_squared = radius_squared - x_squared;
        let y = isqrt_u128(y_squared);
        if y * y != y_squared {
            continue;
        }

        points.push(Point { x, y: y as i128 });
        if y != 0 {
            points.push(Point { x, y: -(y as i128) });
        }
    }

    points.sort_by(|a, b| cmp_angle(&(a.x, a.y), &(b.x, b.y)));
    points
}

pub fn analyze_circle(radius_squared: u128) -> CircleReport {
    let lattice_points = lattice_points_on_circle(radius_squared);
    let max_chord_cluster = max_endpoint_chord_cluster(radius_squared, &lattice_points);
    let max_arc_count_float = max_arc_cluster_float(radius_squared, &lattice_points);
    CircleReport {
        radius_squared,
        lattice_points,
        max_chord_cluster,
        max_arc_count_float,
    }
}

pub fn max_endpoint_chord_cluster(radius_squared: u128, points: &[Point]) -> CircleCluster {
    if points.is_empty() {
        return CircleCluster {
            count: 0,
            start_index: 0,
            points: Vec::new(),
        };
    }

    let chord_sq_threshold = isqrt_u128(radius_squared) as i128;
    let mut best_count = 1_usize;
    let mut best_start = 0_usize;
    let mut best_points = vec![points[0]];

    for start in 0..points.len() {
        let origin = points[start];
        let mut cluster = vec![origin];
        for offset in 1..points.len() {
            let candidate = points[(start + offset) % points.len()];
            if chord_sq(origin, candidate) <= chord_sq_threshold {
                cluster.push(candidate);
            } else {
                break;
            }
        }
        if cluster.len() > best_count {
            best_count = cluster.len();
            best_start = start;
            best_points = cluster;
        }
    }

    CircleCluster {
        count: best_count,
        start_index: best_start,
        points: best_points,
    }
}

fn max_arc_cluster_float(radius_squared: u128, points: &[Point]) -> usize {
    if points.is_empty() {
        return 0;
    }
    let angle_width = (radius_squared as f64).powf(-0.25);
    let mut angles: Vec<f64> = points
        .iter()
        .map(|point| {
            let mut angle = (point.y as f64).atan2(point.x as f64);
            if angle < 0.0 {
                angle += std::f64::consts::TAU;
            }
            angle
        })
        .collect();
    angles.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let original_len = angles.len();
    let mut doubled = angles.clone();
    doubled.extend(angles.iter().map(|angle| angle + std::f64::consts::TAU));

    let mut best = 1_usize;
    let mut end = 0_usize;
    for start in 0..original_len {
        if end < start {
            end = start;
        }
        while end + 1 < start + original_len
            && doubled[end + 1] - doubled[start] <= angle_width + 1e-14
        {
            end += 1;
        }
        best = best.max(end - start + 1);
    }
    best
}

fn chord_sq(a: Point, b: Point) -> i128 {
    let dx = a.x - b.x;
    let dy = a.y - b.y;
    dx * dx + dy * dy
}

#[cfg(test)]
mod tests {
    use super::{analyze_circle, lattice_points_on_circle, Point};

    #[test]
    fn known_circle_contains_reported_four_point_cluster() {
        let points = lattice_points_on_circle(567_454_025);
        for point in [
            Point { x: 23200, y: 5405 },
            Point { x: 23189, y: 5452 },
            Point { x: 23176, y: 5507 },
            Point { x: 23171, y: 5528 },
        ] {
            assert!(points.contains(&point));
        }

        let report = analyze_circle(567_454_025);
        assert_eq!(report.max_chord_cluster.count, 4);
        assert_eq!(report.max_arc_count_float, 4);
    }
}
