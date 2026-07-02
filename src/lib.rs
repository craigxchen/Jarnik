pub mod numerics;

pub use numerics::{
    certificate, circle, clusters, cuts, exact, extension, factor, gaussian, half_angle,
    legacy_arc, points, relation_height, sign_product, transition_patterns,
};

pub use numerics::certificate::{Certificate, CertificateSearch, CertificateSearchReport};
pub use numerics::circle::{analyze_circle, lattice_points_on_circle, CircleReport, Point};
pub use numerics::extension::{
    find_extensions, simplified_leading_coefficient, ExtensionCoefficients, ExtensionSearchConfig,
    ExtensionSearchReport, ExtensionWitness, QuadraticSieve,
};
pub use numerics::legacy_arc::{find_keys_for_value, get_arclength_coeffs, BETA, MAX_R, MIN_R, N};
pub use numerics::sign_product::{
    search_sign_product_counterexamples, sign_coefficient_collisions, SignCollision, SignFactor,
};
