use colored::Colorize;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum FileError {
    #[error("filed to open {}", .filename.bright_yellow().bold())]
    CannotOpen { source: std::io::Error, filename: String },

    #[error("{} not found", .0.bright_yellow().bold())]
    NotFound(String),

    #[error("while creating {filename}")]
    WhileCreating { source: std::io::Error, filename: String },

    #[error("invalid filename: {}", .0.yellow().bold())]
    InvalidFilename(String),
}

#[derive(Error, Debug)]
pub enum DataError {
    #[error("ID {} not found in the specified database", .0.yellow().bold())]
    UnknownId(String),

    #[error("failed to connect to d")]
    FailedToConnect { source: rusqlite::Error, filename: String },
}

#[derive(Error, Debug)]
pub enum MatrixParseError {
    #[error("size missing in distance matrix")]
    SizeMissing,

    #[error("erroneous line found in distance matrix")]
    ErroneousLine,
}

#[derive(Error, Debug)]
pub enum RuntimeError {
    #[error("species {} not found in the provided species tree", .0.yellow().bold())]
    SpeciesNotFound(String),

    #[error("ID {} not found in the provided database", .0.yellow().bold())]
    IdNotFound(String),

    #[error("failed to read distance matrix {}", .filename.yellow().bold())]
    FailedToReadMatrix { source: anyhow::Error, filename: String },
}
