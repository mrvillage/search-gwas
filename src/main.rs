mod cli;
mod consts;
mod data;
mod files;
mod query;

use std::path::PathBuf;

use clap::Parser;

use crate::{
    cli::{Cli, Run},
    files::get_data_dir,
};

struct Context {
    dir: PathBuf,
}

fn main() {
    let _ = env_logger::Builder::from_env(
        env_logger::Env::default().filter_or("SEARCH_GWAS_LOG", "warn"),
    )
    .try_init();

    let cli = Cli::parse();

    let dir = get_data_dir();
    if !dir.exists() {
        std::fs::create_dir(&dir).unwrap();
    }

    let ctx = Context { dir };

    cli.run(ctx);
}
