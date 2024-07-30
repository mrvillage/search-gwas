use std::{
    collections::{HashMap, HashSet},
    env::{current_dir, temp_dir},
    fs::File,
    io::Write,
    path::{Path, PathBuf},
};

use chrono::{DateTime, NaiveDate, NaiveDateTime, TimeZone, Utc};
use csv::DeserializeRecordsIntoIter;
use flate2::read::GzDecoder;
use rand::{distributions::Alphanumeric, thread_rng, Rng};
use rayon::iter::{plumbing::Folder, ParallelBridge, ParallelIterator};
use reqwest::blocking::{Client, ClientBuilder};
use rkyv::ser::serializers::AllocSerializer;

use crate::{
    consts::{OBO_IN_OWL_NS, OWL_NS, RDFS_NS, RDF_NS},
    data::{Association, AzAssociation, Efo, Metadata},
};

#[inline]
fn last_modified_header(client: &Client, url: &str) -> DateTime<Utc> {
    Utc.from_local_datetime(
        &NaiveDateTime::parse_from_str(
            client
                .head(url)
                .send()
                .unwrap()
                .headers()
                .get("Last-Modified")
                .unwrap()
                .to_str()
                .unwrap(),
            "%a, %d %b %Y %H:%M:%S GMT",
        )
        .unwrap(),
    )
    .unwrap()
}

// TODO: make this only check once every few hours or day
#[inline]
fn latest_gwas_date(client: &Client) -> NaiveDate {
    client
        .head("https://www.ebi.ac.uk/gwas/api/search/downloads/alternative")
        .send()
        .unwrap()
        .headers()
        .get("Content-Disposition")
        .unwrap()
        .to_str()
        .unwrap()
        .split('=')
        .last()
        .unwrap()
        .split('_')
        .last()
        .unwrap()[1..]
        .split('.')
        .next()
        .unwrap()
        .parse::<NaiveDate>()
        .unwrap()
}

#[inline]
fn latest_efo_date(client: &Client) -> DateTime<Utc> {
    last_modified_header(client, "https://www.ebi.ac.uk/efo/efo.owl")
}

struct WriteFile<'a> {
    path: &'a PathBuf,
    tmp: PathBuf,
}

impl<'a> WriteFile<'a> {
    fn new(path: &'a PathBuf) -> Self {
        let name = format!(
            "{}.tmp",
            thread_rng()
                .sample_iter(&Alphanumeric)
                .take(20)
                .map(char::from)
                .collect::<String>()
        );
        let tmp = temp_dir().join(name);
        Self { path, tmp }
    }

    fn write_archive<T>(self, data: &T)
    where
        T: rkyv::Serialize<AllocSerializer<256>>,
    {
        let mut tmp = File::create(&self.tmp).unwrap();
        let bytes = rkyv::to_bytes::<T, 256>(data).unwrap();
        tmp.write_all(&bytes).unwrap();
        std::fs::rename(&self.tmp, &self.path).unwrap();
    }

    fn write_str(self, data: &str) {
        let mut tmp = File::create(&self.tmp).unwrap();
        tmp.write_all(data.as_bytes()).unwrap();
        std::fs::rename(&self.tmp, &self.path).unwrap();
    }
}

impl Drop for WriteFile<'_> {
    fn drop(&mut self) {
        match std::fs::metadata(&self.tmp) {
            Ok(_) => {
                if let Err(e) = std::fs::remove_file(&self.tmp) {
                    panic!("Failed to remove temporary file: {}", e);
                }
            },
            Err(e) => match e.kind() {
                std::io::ErrorKind::NotFound => {},
                _ => panic!("Failed to read temporary file metadata: {}", e),
            },
        }
    }
}

fn write_gwas_file(client: &Client, dir: &Path, local: bool) {
    let tsv = associations_tsv_path(dir);
    let processed = associations_path(dir);
    let file = if local {
        println!("Loading local GWAS file...");
        std::fs::read_to_string(tsv).unwrap()
    } else {
        println!("Downloading new GWAS file...");
        let file = client
            .get("https://www.ebi.ac.uk/gwas/api/search/downloads/alternative")
            .send()
            .unwrap()
            .text()
            .unwrap();
        WriteFile::new(&tsv).write_str(&file);
        file
    };

    println!("Processing GWAS file...");
    let headers = file.lines().next().unwrap().split('\t').collect::<Vec<_>>();
    let disease = get_header_position(&headers, "MAPPED_TRAIT_URI");
    let p_value = get_header_position(&headers, "P-VALUE");
    let mapped_gene = get_header_position(&headers, "MAPPED_GENE");
    let accession_id = get_header_position(&headers, "STUDY ACCESSION");
    let link = headers.iter().position(|&header| header == "LINK").unwrap();
    let mut associations = file
        .lines()
        .skip(1)
        .par_bridge()
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .filter_map(|record| {
            let mut traits = record[disease]
                .split(',')
                .map(|disease| disease.split('/').last().unwrap())
                .filter_map(|disease| {
                    if disease.starts_with("EFO_") {
                        Some(disease.split('_').last().unwrap().parse().unwrap())
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>();
            traits.sort();
            // does the ,/-/;/x indicate a gene combination, or two separate genes?
            if record[mapped_gene].trim().is_empty() {
                return None;
            }
            let mut mapped_gene = vec![record[mapped_gene].trim().to_uppercase()];
            mapped_gene.sort();
            Some(Association {
                traits,
                p_value: record[p_value].parse().unwrap(),
                mapped_gene,
                accession_id: record[accession_id][4..].parse().unwrap(),
                pubmed: record[link].split('/').last().unwrap().parse().unwrap(),
            })
        })
        .collect::<Vec<_>>();
    associations.sort();
    associations.dedup();
    WriteFile::new(&processed).write_archive(&associations);

    println!("Processed GWAS file");
}

fn write_efo_file(client: &Client, dir: &Path, local: bool) {
    let owl = efo_owl_path(dir);
    let processed = efo_path(dir);
    let file = if local {
        println!("Loading local EFO file...");
        std::fs::read_to_string(owl).unwrap()
    } else {
        println!("Downloading new EFO file...");
        let file = client
            .get("https://www.ebi.ac.uk/efo/efo.owl")
            .send()
            .unwrap()
            .text()
            .unwrap();
        WriteFile::new(&owl).write_str(&file);
        file
    };

    println!("Processing EFO file...");
    let efo = roxmltree::Document::parse(&file).unwrap();
    let mut efos = efo
        .descendants()
        .par_bridge()
        .filter_map(|node| {
            if node.has_tag_name((OWL_NS, "Class")) {
                if let Some(id) = node.attribute((RDF_NS, "about")) {
                    if let Some(label) = node
                        .children()
                        .find(|node| node.has_tag_name((RDFS_NS, "label")))
                    {
                        let id = id.split('/').last().unwrap();
                        if id.starts_with("EFO_") {
                            let id = id.split('_').last().unwrap().parse().unwrap();
                            let label = label.text().unwrap().trim().to_uppercase();
                            let parent = node
                                .children()
                                .find(|node| node.has_tag_name((RDFS_NS, "subClassOf")))
                                .map(|node| {
                                    let id = node
                                        .attribute((RDF_NS, "resource"))
                                        .unwrap()
                                        .split('/')
                                        .last()
                                        .unwrap();
                                    if id.starts_with("EFO_") {
                                        id.split('_').last().unwrap().parse().unwrap()
                                    } else {
                                        0
                                    }
                                });
                            let synonyms = node
                                .children()
                                .filter(|node| {
                                    node.has_tag_name((OBO_IN_OWL_NS, "hasExactSynonym"))
                                        && node.is_text()
                                })
                                .map(|node| node.text().unwrap().trim().to_uppercase())
                                .collect::<HashSet<_>>();
                            return Some((
                                id,
                                Efo {
                                    id,
                                    label,
                                    parent,
                                    children: HashSet::new(),
                                    synonyms,
                                },
                            ));
                        }
                    }
                }
            }
            None
        })
        .collect::<HashMap<_, _>>();
    for (id, efo) in efos.clone().into_iter() {
        if let Some(subclass) = efo.parent {
            if let Some(parent) = efos.get_mut(&subclass) {
                parent.children.insert(id);
            }
        }
    }
    WriteFile::new(&processed).write_archive(&efos.into_values().collect::<Vec<_>>());

    println!("Processed EFO file");
}

pub fn check_for_updates(dir: &Path, local: bool, force: u8) {
    let client = ClientBuilder::new().timeout(None).build().unwrap();
    let metadata_path = metadata_path(dir);
    match std::fs::read(&metadata_path) {
        Ok(bytes) => {
            let metadata: Metadata = unsafe { rkyv::from_bytes_unchecked(&bytes).unwrap() };
            let now = Utc::now();
            if now - metadata.last_updated < chrono::Duration::days(1) && force == 0 {
                return;
            }
        },
        Err(e) => {
            if e.kind() != std::io::ErrorKind::NotFound {
                panic!("Failed to read metadata file: {}", e);
            }
        },
    }

    let associations_path = associations_path(dir);
    if let Ok(metadata) = std::fs::metadata(associations_path) {
        let latest = latest_gwas_date(&client);
        if DateTime::<Utc>::from(metadata.modified().unwrap())
            .naive_utc()
            .date()
            < latest
            || force == 2
        {
            write_gwas_file(&client, dir, local);
        }
    } else {
        write_gwas_file(&client, dir, false);
    }
    let efo_path = efo_path(dir);
    if let Ok(metadata) = std::fs::metadata(efo_path) {
        let latest = latest_efo_date(&client);
        if DateTime::<Utc>::from(metadata.modified().unwrap()) < latest || force == 2 {
            write_efo_file(&client, dir, local);
        }
    } else {
        write_efo_file(&client, dir, false);
    }

    let bytes = rkyv::to_bytes::<_, 0>(&Metadata {
        last_updated: Utc::now(),
    })
    .unwrap();
    std::fs::write(metadata_path, bytes).unwrap();
}

pub fn get_data_dir() -> PathBuf {
    dirs::data_dir()
        .unwrap_or_else(|| current_dir().unwrap())
        .join("search-gwas")
}

pub fn get_global_dir() -> PathBuf {
    let dir = "/usr/local/share/search-gwas";
    if Path::new(dir).exists() {
        PathBuf::from(dir)
    } else {
        get_data_dir()
    }
}

pub fn get_az_dir() -> PathBuf {
    get_global_dir().join("az470k-proteomics")
}

pub fn associations_path(dir: &Path) -> PathBuf {
    dir.join("associations.rkyv")
}

pub fn associations_tsv_path(dir: &Path) -> PathBuf {
    dir.join("associations.tsv")
}

pub fn efo_path(dir: &Path) -> PathBuf {
    dir.join("efo.rkyv")
}

pub fn efo_owl_path(dir: &Path) -> PathBuf {
    dir.join("efo.owl")
}

pub fn metadata_path(dir: &Path) -> PathBuf {
    dir.join("metadata.rkyv")
}

#[inline]
fn get_header_position(headers: &[&str], header: &str) -> usize {
    headers.iter().position(|&h| h == header).unwrap()
}

pub fn load_associations(dir: &Path) -> Vec<Association> {
    let file = std::fs::read(associations_path(dir)).unwrap();
    unsafe { rkyv::from_bytes_unchecked::<Vec<Association>>(&file).unwrap() }
}

pub fn load_efo(dir: &Path) -> Vec<Efo> {
    let file = std::fs::read(efo_path(dir)).unwrap();
    unsafe { rkyv::from_bytes_unchecked::<Vec<Efo>>(&file).unwrap() }
}

pub struct AzAssociations {
    binary: Option<DeserializeRecordsIntoIter<GzDecoder<File>, AzAssociation>>,
    proteomics: Option<DeserializeRecordsIntoIter<GzDecoder<File>, AzAssociation>>,
    quantitative: Option<DeserializeRecordsIntoIter<GzDecoder<File>, AzAssociation>>,
}

impl AzAssociations {
    pub fn new() -> Self {
        let file = std::fs::File::open(get_az_dir().join("binary.csv.gz"));
        let binary = if let Ok(file) = file {
            let reader = csv::ReaderBuilder::new()
                .has_headers(true)
                .from_reader(flate2::read::GzDecoder::new(file));
            Some(reader.into_deserialize::<AzAssociation>())
        } else {
            None
        };
        let file = std::fs::File::open(get_az_dir().join("proteomics.csv.gz"));
        let proteomics = if let Ok(file) = file {
            let reader = csv::ReaderBuilder::new()
                .has_headers(true)
                .from_reader(flate2::read::GzDecoder::new(file));
            Some(reader.into_deserialize::<AzAssociation>())
        } else {
            None
        };
        let file = std::fs::File::open(get_az_dir().join("quantitative.csv.gz"));
        let quantitative = if let Ok(file) = file {
            let reader = csv::ReaderBuilder::new()
                .has_headers(true)
                .from_reader(flate2::read::GzDecoder::new(file));
            Some(reader.into_deserialize::<AzAssociation>())
        } else {
            None
        };
        Self {
            binary,
            proteomics,
            quantitative,
        }
    }
}

impl Iterator for AzAssociations {
    type Item = AzAssociation;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(iter) = &mut self.binary {
            if let Some(assoc) = iter.next() {
                return Some(assoc.unwrap());
            }
        }
        if let Some(iter) = &mut self.proteomics {
            if let Some(assoc) = iter.next() {
                return Some(assoc.unwrap());
            }
        }
        if let Some(iter) = &mut self.quantitative {
            if let Some(assoc) = iter.next() {
                return Some(assoc.unwrap());
            }
        }
        None
    }
}

impl ParallelIterator for AzAssociations {
    type Item = AzAssociation;

    fn drive_unindexed<C>(self, consumer: C) -> C::Result
    where
        C: rayon::iter::plumbing::UnindexedConsumer<Self::Item>,
    {
        if let Some(iter) = self.binary {
            println!("BINARY");
            iter.into_iter()
                .par_bridge()
                .filter_map(Result::ok)
                .drive_unindexed(consumer)
        } else if let Some(iter) = self.proteomics {
            println!("PROTEOMICS");
            iter.into_iter()
                .par_bridge()
                .filter_map(Result::ok)
                .drive_unindexed(consumer)
        } else if let Some(iter) = self.quantitative {
            println!("QUANTITATIVE");
            iter.into_iter()
                .par_bridge()
                .filter_map(Result::ok)
                .drive_unindexed(consumer)
        } else {
            consumer.into_folder().complete()
        }
    }
}
