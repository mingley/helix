//! Deterministic regulatory modeling with explicit fold and assay boundaries.
//!
//! Scientific feature/fold choices are fixed by the versioned study protocol.
mod exact;
pub mod math;
pub mod motif;
pub mod study;
pub use study::{run, InputRecord, Protocol, StudyResult};
