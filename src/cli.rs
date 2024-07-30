use clap::{Args, Parser, Subcommand};

use crate::{
    files::{check_for_updates, get_az_dir, load_associations, load_efo},
    query::{find_efo, parse_genes, query, query_az},
    Context,
};

pub trait Run {
    fn run(self, ctx: Context);
}

#[derive(Parser)]
#[command(version)]
pub struct Cli {
    #[command(subcommand)]
    command: Commands,
}

impl Run for Cli {
    #[inline]
    fn run(self, ctx: Context) {
        self.command.run(ctx);
    }
}

#[derive(Subcommand)]
enum Commands {
    #[command(about = "Download the latest GWAS and EFO data if available")]
    Update(Update),
    #[command(about = "Query the GWAS catalog for a trait")]
    Trait(Trait),
    #[command(about = "Update the AstraZeneca PheWAS catalog", hide = true)]
    AzUpdate(AzUpdate),
    #[command(about = "Query the AstraZeneca PheWAS catalog for a trait")]
    AzTrait(AzTrait),
}

impl Run for Commands {
    #[inline]
    fn run(self, ctx: Context) {
        match self {
            Self::Update(update) => update.run(ctx),
            Self::Trait(query) => query.run(ctx),
            Self::AzUpdate(update) => update.run(ctx),
            Self::AzTrait(query) => query.run(ctx),
        }
    }
}

#[derive(Args)]
#[command(about = "Download the latest GWAS and EFO data if available")]
struct Update {
    #[arg(
        short,
        long,
        help = "Forcibly check for updates even if checked recently, specify twice to forcibly redownload the data",
        action = clap::ArgAction::Count,
    )]
    force: u8,
    #[arg(short, long, help = "Reprocess the local files")]
    reprocess: bool,
}

impl Run for Update {
    fn run(self, ctx: Context) {
        check_for_updates(
            &ctx.dir,
            self.reprocess,
            if self.reprocess { 2 } else { self.force },
        );
        println!("Up to date!");
    }
}

#[derive(Args)]
struct Trait {
    #[arg(help = "The EFO label to query")]
    efo: String,
    #[arg(short, long, action = clap::ArgAction::Append, help = "Gene(s) to query")]
    gene: Vec<String>,
    #[arg(
        short = 'a',
        long = "with-associations",
        help = "Show full association data"
    )]
    with_associations: bool,
    #[arg(
        short = 'l',
        long = "with-pubmed-links",
        help = "Show PubMed links instead of IDs"
    )]
    with_pubmed_links: bool,
    #[arg(short, long, help = "Replace tables with CSV output")]
    csv: bool,
}

impl Run for Trait {
    fn run(self, ctx: Context) {
        check_for_updates(&ctx.dir, false, 0);
        let orig = self.efo.trim();
        let genes = parse_genes(&self.gene);
        let efos = load_efo(&ctx.dir);
        let associations = load_associations(&ctx.dir);
        let efo = match find_efo(&efos, &orig.to_uppercase()) {
            Some(efo) => efo,
            None => {
                eprintln!("\"{orig}\" is not a valid EFO label");
                return;
            },
        };
        query(
            efo,
            genes,
            &associations,
            self.with_associations,
            self.with_pubmed_links,
            self.csv,
        );
    }
}

#[derive(Args)]
struct AzUpdate;

impl Run for AzUpdate {
    fn run(self, ctx: Context) {
        let dir = get_az_dir();
    }
}

#[derive(Args)]
struct AzTrait {
    #[arg(help = "The trait to query")]
    trait_: String,
    #[arg(short, long, action = clap::ArgAction::Append, help = "Gene(s) to query")]
    gene: Vec<String>,
    #[arg(
        short = 'a',
        long = "with-associations",
        help = "Show full association data"
    )]
    with_associations: bool,
    #[arg(short, long, help = "Replace tables with CSV output")]
    csv: bool,
}

impl Run for AzTrait {
    fn run(self, _ctx: Context) {
        let orig = self.trait_.trim().to_lowercase();
        let genes = parse_genes(&self.gene);
        query_az(&orig, genes, self.with_associations, self.csv);
    }
}
