use flate2::read::MultiGzDecoder;
use memmap2::Mmap;
use rayon::prelude::*;
use std::error::Error;
use std::{
    fmt::Debug,
    fs::File,
    io::{BufWriter, Read, Write},
    path::Path,
};

pub fn get_sizes<T: AsRef<Path> + Debug>(fasta: T, accession_only: bool) -> Result<Vec<(String, u64)>, Box<dyn Error>> {
    let path = fasta.as_ref();
    let ext = path.extension().unwrap();
    let file = File::open(path)?;

    let lines = match ext.to_str().unwrap() {
        "gz" => with_gz(&file, accession_only)?,
        "fa" | "fasta" | "fna" => raw(&file, accession_only)?,
        _ => panic!("ERROR: Not a fasta. Wrong file format!"),
    };

    Ok(lines)
}

pub fn raw(file: &File, accession_only: bool) -> Result<Vec<(String, u64)>, Box<dyn Error>> {
    let mmap = unsafe { Mmap::map(file)? };
    let lines = chromsize(&mmap, accession_only)?;

    Ok(lines)
}

fn with_gz(file: &File, accession_only: bool) -> Result<Vec<(String, u64)>, Box<dyn Error>> {
    let mmap = unsafe { Mmap::map(file)? };
    let mut decoder = MultiGzDecoder::new(&mmap[..]);

    let mut buffer = Vec::with_capacity(100 * 1024 * 1024); // 100MB buffer
    decoder.read_to_end(&mut buffer)?;

    let lines = chromsize(&buffer, accession_only)?;

    Ok(lines)
}

fn chromsize(data: &[u8], accession_only: bool) -> Result<Vec<(String, u64)>, Box<dyn Error>> {
    let lines = data
        .par_split(|&c| c == b'>')
        // NOTE this filter will only work when '>' is litterally the first byte of the string. It will not pick up the corner case where there are multiple empty lines before the the first sequence header
        .filter(|chunk| !chunk.is_empty())
        .map(|chunk| {
            let mut totals = 0u64;
            let stop = memchr::memchr(b'\n', chunk).unwrap_or(0);
            let name_stop = match accession_only {
                true => memchr::memchr(b' ', &chunk[..stop]).unwrap_or(stop),
                false => stop,
            };
            let chr = unsafe {
                std::str::from_utf8_unchecked(&chunk[..name_stop])
                    .trim()
                    .to_string()
            };
            let data = &chunk[stop + 1..];
            for line in data.split(|&c| c == b'\n') {
                totals += unsafe { std::str::from_utf8_unchecked(line).trim().len() as u64 };
            }
            (chr, totals)
        })
        .collect::<Vec<(String, u64)>>();

    Ok(lines)
}

pub fn writer<T>(sizes: Vec<(String, u64)>, out: T)
where
    T: AsRef<Path> + Debug,
{
    let o = match File::create(out) {
        Ok(f) => f,
        Err(e) => panic!("Error creating file: {}", e),
    };
    let mut writer = BufWriter::new(o);

    let mut it = sizes.iter();
    // Skip first item if empty. This happens when there are multiple empty lines before the first sequence header in the fasta file
    if 0 == sizes[0].1 && 0 == sizes[0].0.len() {
        it.next();
    }
    for (k, v) in it {
        writeln!(writer, "{}\t{}", k, v).unwrap();
    }
}
