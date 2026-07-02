use crate::numerics::exact::gcd_i128;
use crate::numerics::points::LatticePointRecord;
use serde::Serialize;
use std::collections::BTreeMap;

#[derive(Clone, Debug, Serialize)]
pub struct HalfAngleDiagnostics {
    pub pairs: Vec<HalfAnglePair>,
    pub root: Option<RootHalfAngleDiagnostics>,
}

#[derive(Clone, Debug, Serialize)]
pub struct HalfAnglePair {
    pub i: usize,
    pub j: usize,
    pub c: i128,
    pub s: i128,
    pub a: i128,
    pub b: i128,
    pub d: Option<u128>,
    pub height: Option<f64>,
    pub angle: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct RootHalfAngleDiagnostics {
    pub root_point: [i128; 2],
    pub primitive_scale_h: i128,
    pub primitive_a: i128,
    pub primitive_b: i128,
    pub successors: Vec<RootSuccessor>,
    pub b_distribution: Vec<BCount>,
    pub fixed_b_groups: Vec<FixedBGroup>,
    pub rooted_normal_indices_available: bool,
}

#[derive(Clone, Debug, Serialize)]
pub struct RootSuccessor {
    pub point_index: usize,
    pub point: [i128; 2],
    pub q_a: i128,
    pub q_b: i128,
    pub d: Option<u128>,
    pub height: Option<f64>,
    pub angle: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct BCount {
    pub b: i128,
    pub count: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct FixedBGroup {
    pub b: i128,
    pub point_indices: Vec<usize>,
}

pub fn analyze_half_angles(n: u64, cluster_points: &[&LatticePointRecord]) -> HalfAngleDiagnostics {
    let mut pairs = Vec::new();
    for i in 0..cluster_points.len() {
        for j in (i + 1)..cluster_points.len() {
            if let Some(pair) = half_angle_pair(n, i, j, cluster_points[i], cluster_points[j]) {
                pairs.push(pair);
            }
        }
    }

    let root = cluster_points.first().map(|root_point| {
        let h = gcd_i128(root_point.x, root_point.y).max(1);
        let primitive_a = root_point.x / h;
        let primitive_b = root_point.y / h;
        let successors = pairs
            .iter()
            .filter(|pair| pair.i == 0)
            .map(|pair| {
                let point = cluster_points[pair.j];
                RootSuccessor {
                    point_index: pair.j,
                    point: point.point_array(),
                    q_a: pair.a,
                    q_b: pair.b,
                    d: pair.d,
                    height: pair.height,
                    angle: pair.angle,
                }
            })
            .collect::<Vec<_>>();

        let mut by_b: BTreeMap<i128, Vec<usize>> = BTreeMap::new();
        for successor in &successors {
            by_b.entry(successor.q_b)
                .or_default()
                .push(successor.point_index);
        }
        let b_distribution = by_b
            .iter()
            .map(|(&b, indices)| BCount {
                b,
                count: indices.len(),
            })
            .collect();
        let fixed_b_groups = by_b
            .into_iter()
            .filter(|(_, indices)| indices.len() >= 2)
            .map(|(b, point_indices)| FixedBGroup { b, point_indices })
            .collect();

        RootHalfAngleDiagnostics {
            root_point: root_point.point_array(),
            primitive_scale_h: h,
            primitive_a,
            primitive_b,
            successors,
            b_distribution,
            fixed_b_groups,
            rooted_normal_indices_available: false,
        }
    });

    HalfAngleDiagnostics { pairs, root }
}

fn half_angle_pair(
    n: u64,
    i: usize,
    j: usize,
    left: &LatticePointRecord,
    right: &LatticePointRecord,
) -> Option<HalfAnglePair> {
    let c = right.x * left.x + right.y * left.y;
    let s = right.y * left.x - right.x * left.y;
    let denominator = n as i128 + c;
    if denominator <= 0 || s == 0 {
        return None;
    }
    let gcd = gcd_i128(denominator, s).max(1);
    let mut a = denominator / gcd;
    let mut b = s / gcd;
    if a < 0 {
        a = -a;
        b = -b;
    }
    if b < 0 {
        b = -b;
    }
    let d = checked_square_sum(a, b);
    let height = d.map(|value| (value as f64).sqrt());
    let angle = 2.0 * (b as f64).atan2(a as f64);
    Some(HalfAnglePair {
        i,
        j,
        c,
        s,
        a,
        b,
        d,
        height,
        angle,
    })
}

fn checked_square_sum(a: i128, b: i128) -> Option<u128> {
    let au = a.unsigned_abs();
    let bu = b.unsigned_abs();
    let aa = au.checked_mul(au)?;
    let bb = bu.checked_mul(bu)?;
    aa.checked_add(bb)
}

#[cfg(test)]
mod tests {
    use crate::numerics::points::LatticePointRecord;

    use super::analyze_half_angles;

    #[test]
    fn computes_basic_half_angle_denominator() {
        let points = vec![
            LatticePointRecord {
                x: 5,
                y: 0,
                angle: 0.0,
                exponents: vec![],
            },
            LatticePointRecord {
                x: 3,
                y: 4,
                angle: (4.0_f64).atan2(3.0),
                exponents: vec![],
            },
        ];
        let refs = points.iter().collect::<Vec<_>>();
        let diagnostics = analyze_half_angles(25, &refs);
        assert_eq!(diagnostics.pairs.len(), 1);
        assert_eq!(diagnostics.pairs[0].a, 2);
        assert_eq!(diagnostics.pairs[0].b, 1);
        assert_eq!(diagnostics.pairs[0].d, Some(5));
    }
}
