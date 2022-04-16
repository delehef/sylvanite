use anyhow::*;

mod align;

fn main() -> Result<()> {
    println!("{}", align::align_sw("AIOPU", "ABCD", &|x, y| x.max(y) as f32));
    Ok(())
}
