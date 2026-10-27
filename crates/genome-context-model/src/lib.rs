//! Outcome-independent DNA contracts and small matched-context ridge models.
//!
//! This crate specifies engineering behavior, not a scientific study protocol.
//! Source selection, grouping, folds, endpoints, alpha, and interpretation must
//! be fixed separately before evaluating real outcomes.

pub mod data;
pub mod evaluation;
pub mod features;
pub mod model;
pub mod resampling;

pub use data::{GroupSplit, MatchedPanel, MeasuredOutcome, OutcomeMatrix, SequenceRecord};
pub use evaluation::{evaluate, ContextEvaluation};
pub use features::{
    canonical_reverse_complement, feature_names, sequence_features, SequenceFeatures,
};
pub use model::{ContextFit, ModelKind, ModelPredictions};
pub use resampling::{paired_group_bootstrap, GroupBootstrap, LossWeighting, ModelComparison};
