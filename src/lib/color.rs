use ansi_term::Color;

/// Color for matching bases in a FASTQ
pub const COLOR_BASES: Color = Color::Red;
/// Color for matching base qualities in a FASTQ
pub const COLOR_QUALS: Color = Color::Fixed(22);
/// Color for a matching read name (head) of a FASTQ
pub const COLOR_HEAD: Color = Color::Fixed(30);
/// Color for all non-matching text in a FASTQ
pub const COLOR_BACKGROUND: Color = Color::Fixed(240);

/// Colors the text with the given color
pub fn color(text: &[u8], color: &Color) -> Vec<u8> {
    let mut colored: Vec<u8> = Vec::with_capacity(text.len());
    color.paint(text).write_to(&mut colored).unwrap();
    colored
}

/// Color for the read name (head) of a FASTQ
pub fn color_head(text: &[u8]) -> Vec<u8> {
    color(text, &COLOR_HEAD)
}

/// Color for the base qualities of a FASTQ
pub fn color_quals(text: &[u8]) -> Vec<u8> {
    color(text, &COLOR_QUALS)
}

/// Color for non-matching bases, quals, and other text
pub fn color_background(text: &[u8]) -> Vec<u8> {
    color(text, &COLOR_BACKGROUND)
}
