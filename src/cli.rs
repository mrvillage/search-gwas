use clap::{Args, Parser, Subcommand};

use crate::{
    files::{check_for_updates, load_associations, load_efo},
    query::{find_efo, parse_genes, query},
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
    Update(Update),
    Trait(Trait),
}

impl Run for Commands {
    #[inline]
    fn run(self, ctx: Context) {
        match self {
            Self::Update(update) => update.run(ctx),
            Self::Trait(query) => query.run(ctx),
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
    efo: String,
    #[arg(short, long, action = clap::ArgAction::Append)]
    gene: Vec<String>,
    #[arg(short = 'a', long = "with-associations", help = "Show associations")]
    with_associations: bool,
    #[arg(short = 'l', long = "with-pubmed-links", help = "Show PubMed links")]
    with_pubmed_links: bool,
    #[arg(short, long, help = "Replace tables with CSV output")]
    csv: bool,
}

impl Run for Trait {
    fn run(self, ctx: Context) {
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
