use rusqlite::Connection;
use std::{
    collections::HashMap,
    path::{Path, PathBuf},
};

use anyhow::*;
use clap::Parser;

mod align;

#[derive(Parser, Debug)]
#[clap(author, version, about)]
struct Settings {
    #[clap(required = true)]
    infiles: Vec<String>,
    #[clap(short, long)]
    outdir: Option<String>,
    #[clap(short, long)]
    database: String,
    #[clap(short, long, default_value_t = 15)]
    window: usize,
}

struct Gene {
    gene: String,
    protein: String,
    left_landscape: Vec<usize>,
    right_landscape: Vec<usize>,
}

struct Register {
    genes: HashMap<String, Gene>,
}

fn read_db(filename: &str, window: usize) -> Result<Register> {
    let conn = Connection::open(filename)?;
    let mut query =
        conn.prepare("SELECT gene, protein, left_tail_ids, right_tail_ids FROM genome")?;
    let genes = query.query_map([], |r| {
        let gene: String = r.get(0)?;
        let protein: String = r.get(1)?;
        let left_landscape: String = r.get(2)?;
        let mut left_landscape = left_landscape
            .split('.')
            .map(|x| x.parse::<usize>().unwrap())
            .collect::<Vec<_>>();
        left_landscape.reverse();
        left_landscape.truncate(window);
        left_landscape.reverse();

        let right_landscape: String = r.get(3)?;
        let mut right_landscape = right_landscape
            .split('.')
            .map(|x| x.parse::<usize>().unwrap())
            .collect::<Vec<_>>();
        right_landscape.truncate(window);

        std::result::Result::Ok(Gene {
            gene,
            protein,
            left_landscape,
            right_landscape,
        })
    })?;

    Ok(Register {
        genes: genes
            .into_iter()
            .map(|g| {
                let g = g.unwrap();
                (g.protein.clone(), g)
            })
            .collect(),
    })
}

fn process_file(filename: &str, register: &Register, settings: &Settings) -> Result<PathBuf> {
    let mut outfile = if let Some(ref outdir) = settings.outdir {
        PathBuf::from(outdir)
    } else {
        PathBuf::from(Path::new(filename).parent().unwrap())
    };
    outfile.set_file_name(Path::new(filename).with_extension("mat"));

    Ok(outfile)
}

fn main() -> Result<()> {
    let args = Settings::parse();
    let register = read_db(&args.database, args.window)?;

    println!(
        "{}",
        align::score_landscape(
            [1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3],
            [1, 2, 3, 1, 2, 3, 1, 2, 1, 2, 3, 1, 2, 3, 1, 2, 3],
            &|x, y| x.max(y) as f32
        )
    );

    for f in args.infiles.iter() {
        println!("Processing {}", f);
        let out = process_file(&f, &register, &args)?;
        println!("Result written to {:?}\n", out);
    }

    Ok(())
}
