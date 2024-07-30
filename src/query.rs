use std::collections::HashSet;

use log::debug;
use prettytable::{row, Table};
use rayon::prelude::*;

use crate::{
    data::{Association, Efo},
    files::AzAssociations,
};

pub fn parse_genes(genes: &[String]) -> Vec<String> {
    genes
        .iter()
        .flat_map(|i| i.split(',').map(|i| i.trim().to_uppercase()))
        .collect()
}

pub fn find_efo<'a>(efos: &'a [Efo], label: &str) -> Option<&'a Efo> {
    match efos.iter().find(|i| i.label == *label) {
        None => efos.iter().find(|i| i.synonyms.contains(label)),
        efo => efo,
    }
}

pub fn query(
    efo: &Efo,
    genes: Vec<String>,
    associations: &[Association],
    with_associations: bool,
    with_pubmed_links: bool,
    csv: bool,
) {
    let results = associations
        .iter()
        .filter(|result| result.is_significant() && result.is_associated_with(efo.id))
        .collect::<Vec<_>>();
    println!("{}:", efo.label);
    if results.is_empty() {
        println!("  No significant associations found");
    } else if genes.is_empty() {
        if with_associations {
            let mut table = Table::new();
            table.set_titles(row!["Genes", "P-value", "Accession ID", "PubMed ID"]);
            for assoc in results {
                let pubmed = if with_pubmed_links {
                    format!("https://pubmed.ncbi.nlm.nih.gov/{}", assoc.pubmed)
                } else {
                    assoc.pubmed.to_string()
                };
                table.add_row(row![
                    assoc.mapped_gene.join(", "),
                    format!("{:e}", assoc.p_value),
                    format!("GCST{}", assoc.accession_id.to_string()),
                    pubmed,
                ]);
            }
            if csv {
                let mut buf = Vec::new();
                table.to_csv(&mut buf).unwrap();
                String::from_utf8(buf)
                    .unwrap()
                    .lines()
                    .for_each(|i| println!("  {}", i));
            } else {
                table.to_string().lines().for_each(|i| println!("  {}", i));
            }
        } else {
            let genes = results
                .iter()
                .flat_map(|result| result.mapped_gene.iter())
                .collect::<HashSet<_>>();
            if csv {
                println!(
                    "{}",
                    genes
                        .into_iter()
                        .map(|i| i.as_str())
                        .collect::<Vec<_>>()
                        .join(",")
                );
            } else {
                for gene in genes {
                    println!("  {gene}");
                }
            }
        }
    } else if with_associations {
        for gene in genes {
            let assocs = results
                .iter()
                .filter(|result| result.mapped_gene.contains(&gene))
                .collect::<Vec<_>>();
            println!("  {gene}:");
            if assocs.is_empty() {
                if !csv {
                    println!("    NONE");
                }
            } else {
                let mut table = Table::new();
                table.set_titles(row!["P-value", "Accession ID", "PubMed ID"]);
                for assoc in assocs {
                    let pubmed = if with_pubmed_links {
                        format!("https://pubmed.ncbi.nlm.nih.gov/{}", assoc.pubmed)
                    } else {
                        assoc.pubmed.to_string()
                    };
                    table.add_row(row![
                        format!("{:e}", assoc.p_value),
                        format!("GCST{}", assoc.accession_id.to_string()),
                        pubmed,
                    ]);
                }
                if csv {
                    let mut buf = Vec::new();
                    table.to_csv(&mut buf).unwrap();
                    String::from_utf8(buf)
                        .unwrap()
                        .lines()
                        .for_each(|i| println!("    {}", i));
                } else {
                    table
                        .to_string()
                        .lines()
                        .for_each(|i| println!("    {}", i));
                }
            }
        }
    } else {
        let mut associated = Vec::with_capacity(genes.len());
        let mut not_associated = Vec::with_capacity(genes.len());
        for gene in genes {
            let assoc = results
                .iter()
                .any(|result| result.mapped_gene.contains(&gene));
            if assoc {
                associated.push(gene);
            } else {
                not_associated.push(gene);
            }
        }
        if !associated.is_empty() {
            println!("  ASSOCIATED:");
            if csv {
                println!(
                    "    {}",
                    associated.into_iter().collect::<Vec<_>>().join(",")
                );
            } else {
                for gene in associated {
                    println!("    {gene}");
                }
            }
        }
        if !not_associated.is_empty() {
            println!("  NOT ASSOCIATED:");
            if csv {
                println!(
                    "    {}",
                    not_associated.into_iter().collect::<Vec<_>>().join(",")
                );
            } else {
                for gene in not_associated {
                    println!("    {gene}");
                }
            }
        }
    }
}

pub fn query_az(term: &str, genes: Vec<String>, with_associations: bool, csv: bool) {
    let associations = AzAssociations::new();
    debug!("Loading AZ associations...");
    let associations = ParallelIterator::collect::<Vec<_>>(associations.into_par_iter());
    debug!("Loaded {} AZ associations", associations.len());
    let results = associations
        .into_par_iter()
        .filter(|result| result.is_significant() && result.is_associated_with(term))
        .collect::<Vec<_>>();
    debug!("Found {} significant associations", results.len());
    println!("{}:", term);
    if results.is_empty() {
        println!("  No significant associations found");
    } else if genes.is_empty() {
        if with_associations {
            let mut table = Table::new();
            table.set_titles(row!["Trait", "Genes", "P-value"]);
            for assoc in results {
                table.add_row(row![
                    assoc.trait_,
                    assoc.mapped_gene,
                    format!("{:e}", assoc.p_value),
                ]);
            }
            if csv {
                let mut buf = Vec::new();
                table.to_csv(&mut buf).unwrap();
                String::from_utf8(buf)
                    .unwrap()
                    .lines()
                    .for_each(|i| println!("  {}", i));
            } else {
                table.to_string().lines().for_each(|i| println!("  {}", i));
            }
        } else {
            let genes = results
                .iter()
                .map(|result| result.mapped_gene.as_str())
                .collect::<HashSet<_>>();
            if csv {
                println!("{}", genes.into_iter().collect::<Vec<_>>().join(","));
            } else {
                for gene in genes {
                    println!("  {gene}");
                }
            }
        }
    } else if with_associations {
        for gene in genes {
            let assocs = results
                .iter()
                .filter(|result| result.mapped_gene.contains(&gene))
                .collect::<Vec<_>>();
            println!("  {gene}:");
            if assocs.is_empty() {
                if !csv {
                    println!("    NONE");
                }
            } else {
                let mut table = Table::new();
                table.set_titles(row!["Trait", "P-value"]);
                for assoc in assocs {
                    table.add_row(row![assoc.trait_, format!("{:e}", assoc.p_value),]);
                }
                if csv {
                    let mut buf = Vec::new();
                    table.to_csv(&mut buf).unwrap();
                    String::from_utf8(buf)
                        .unwrap()
                        .lines()
                        .for_each(|i| println!("    {}", i));
                } else {
                    table
                        .to_string()
                        .lines()
                        .for_each(|i| println!("    {}", i));
                }
            }
        }
    } else {
        let mut associated = Vec::with_capacity(genes.len());
        let mut not_associated = Vec::with_capacity(genes.len());
        for gene in genes {
            let assoc = results
                .iter()
                .any(|result| result.mapped_gene.contains(&gene));
            if assoc {
                associated.push(gene);
            } else {
                not_associated.push(gene);
            }
        }
        if !associated.is_empty() {
            println!("  ASSOCIATED:");
            if csv {
                println!(
                    "    {}",
                    associated.into_iter().collect::<Vec<_>>().join(",")
                );
            } else {
                for gene in associated {
                    println!("    {gene}");
                }
            }
        }
        if !not_associated.is_empty() {
            println!("  NOT ASSOCIATED:");
            if csv {
                println!(
                    "    {}",
                    not_associated.into_iter().collect::<Vec<_>>().join(",")
                );
            } else {
                for gene in not_associated {
                    println!("    {gene}");
                }
            }
        }
    }
}
