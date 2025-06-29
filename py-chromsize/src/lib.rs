use chromsize;
use pyo3::prelude::*;
use std::path::PathBuf;

#[pyfunction]
#[pyo3(signature = (fasta, accession_only=false))]
fn get_chromsizes(
    py: Python,
    fasta: PyObject,
    accession_only: bool,
) -> PyResult<Vec<(String, u64)>> {
    let fasta = PathBuf::from(fasta.extract::<String>(py)?);

    let sizes = chromsize::get_sizes(&fasta, accession_only);
    sizes.map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))
}

#[pyfunction]
#[pyo3(signature = (fasta, output, accession_only=false))]
fn write_chromsizes(
    py: Python,
    fasta: PyObject,
    output: PyObject,
    accession_only: bool,
) -> PyResult<String> {
    let fasta = PathBuf::from(fasta.extract::<String>(py)?);
    let output = PathBuf::from(output.extract::<String>(py)?);

    let sizes = chromsize::get_sizes(&fasta, accession_only);
    if let Ok(sizes) = sizes {
        chromsize::writer(sizes, &output);
        Ok(format!("Chromosome sizes written to {}", output.display()))
    } else {
        Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "Failed to get chromosome sizes",
        ))
    }
}

#[pymodule]
#[pyo3(name = "chromsize")]
fn py_chromsize(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(get_chromsizes, m)?)?;
    m.add_function(wrap_pyfunction!(write_chromsizes, m)?)?;
    Ok(())
}
