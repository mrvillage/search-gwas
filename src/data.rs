use std::{collections::HashSet, hash::Hash};

use chrono::{DateTime, Utc};
use rkyv::{Archive, Deserialize, Serialize};

use crate::consts::THRESHOLD;

#[derive(Debug, Archive, Serialize, Deserialize, PartialEq, PartialOrd)]
pub struct Association {
    // sorted
    pub(crate) traits: Vec<u32>,
    pub(crate) p_value: f64,
    // uppercase, sorted
    pub(crate) mapped_gene: Vec<String>,
    pub(crate) accession_id: u32,
    pub(crate) pubmed: u32,
}

impl Eq for Association {}
#[allow(clippy::derive_ord_xor_partial_ord)]
impl Ord for Association {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.partial_cmp(other).unwrap()
    }
}

impl Hash for Association {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.traits.hash(state);
        self.p_value.to_bits().hash(state);
        self.mapped_gene.hash(state);
        self.accession_id.hash(state);
        self.pubmed.hash(state);
    }
}

impl Association {
    #[inline]
    pub fn is_significant(&self) -> bool {
        self.p_value < THRESHOLD
    }

    #[inline]
    pub fn is_associated_with(&self, efo: u32) -> bool {
        self.traits.contains(&efo)
    }
}

#[derive(Clone, Debug, Archive, Serialize, Deserialize)]
pub struct Efo {
    pub(crate) id: u32,
    // uppercase
    pub(crate) label: String,
    pub(crate) parent: Option<u32>,
    pub(crate) children: HashSet<u32>,
    // uppercase
    pub(crate) synonyms: HashSet<String>,
}

impl Hash for Efo {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.id.hash(state);
    }
}

impl PartialEq for Efo {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

impl Eq for Efo {}

#[derive(Debug, Archive, Serialize, Deserialize)]
pub struct Metadata {
    pub(crate) last_updated: DateTime<Utc>,
}
