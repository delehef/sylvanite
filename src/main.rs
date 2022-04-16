use anyhow::*;

mod align;

fn main() -> Result<()> {
    println!(
        "{}",
        align::score_landscape(
            [1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3],
            [1, 2, 3, 1, 2, 3, 1, 2, 1, 2, 3, 1, 2, 3, 1, 2, 3],
            &|x, y| x.max(y) as f32
        )
    );
    Ok(())
}
