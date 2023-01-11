use colored::Colorize;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum FileError {
    #[error("{} not found", .0.white().bold())]
    NotFound(String),
    #[error("while creating {filename}")]
    WhileCreating { source: std::io::Error, filename: String },
}

#[derive(Error, Debug)]
pub enum DataError {
    #[error("ID {} not found in the specified database", .0.yellow().bold())]
    UnknownId(String),
}
