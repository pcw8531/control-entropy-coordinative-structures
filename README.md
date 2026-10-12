# Average information across a network quantifies what a coordinative structure achieves in motor control

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21789488.svg)](https://doi.org/10.5281/zenodo.21789488)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.11](https://img.shields.io/badge/Python-3.11-blue.svg)](requirements.txt)

Code and data for the paper. The values the paper reports are recomputed here from the released files, by the scripts and notebooks listed below.

**Chulwook Park** · Department of Physical Education, Seoul National University · [ORCID 0000-0001-8714-5760](https://orcid.org/0000-0001-8714-5760)

---

## Quick start

```bash
pip install -r requirements.txt
python code/verify_published_values.py
```

That script checks the size and hash of every released file it reads, recomputes the values of Sections 3.1 to 3.3 and of the behavioural check in Section 3.4, and exits non-zero if anything fails. It is the gate every other script runs behind.

> On a clean checkout it reports **75 of 75 values reproducing**.

---

## Reproducing the individual tables

| Script | Produces |
| --- | --- |
| `code/s1_ensembles.py` | Table 1, Figure 3, Supplementary Note 1 and Supplementary Tables 1 to 3 |
| `code/s2_measured.py` | Supplementary Tables 5 and 6, Supplementary Notes 2 and 3, and the upper row of Figure 4 |
| `code/s3_matched_nulls.py` | Supplementary Table 4 |
| `code/s9_resampling.py` | s9 and s10. Supplementary Note 4, Supplementary Table 7, Supplementary Fig. 1 and the lower row of Figure 4 |
| `notebooks/s11_controllability_checks.ipynb` | Section 3.4, Supplementary Note 6 and Supplementary Table 8 |
| `notebooks/fig5_controllability_publication_1.ipynb` | Figure 5 |
| `notebooks/movie_s1_3.ipynb` | Supplementary Video 1 |
| `s4`, `s5`, `s6` | The three analyses specified in advance in Supplementary Note 5 |
| `reproduce.py` | Runs the gate and then all four scripts in order. `--gate` runs the gate alone. |

The notebooks run on their own under Jupyter. Nothing has to be run in a particular order beyond the verification gate.

---

## Repository layout

| Path | What it holds |
| --- | --- |
| `code/core.py` | The equation module. Every number the analysis scripts compute traces to one definition here. Its comments keep the equation numbers of an earlier draft. |
| `code/generators.py` | The four matched topology ensembles |
| `code/loaders.py` | Readers and the provenance record |
| `code/s*.py` | The analyses, named as the Supplementary Information names them |
| `code/verify_published_values.py` | The gate. Hashes every input it reads and recomputes the values of Sections 3.1 to 3.3 and the behavioural check of Section 3.4. |
| `notebooks/` | The s11 checks, and the drawing notebooks of Figure 5 and Supplementary Video 1 |
| `data/derived/` | Ten CSVs. Nine are computed from a released file or transcribed from a published table, and each carries a comment header saying where it came from. `s11_controllability_checks.csv` is the output of the s11 notebook. |
| `data/external/` | Files belonging to other groups. See below. |
| `checksums.sha256` | sha256 of every tracked file in the deposit. Verify with `sha256sum -c checksums.sha256`. |
| [`SOURCES.md`](SOURCES.md) | Full provenance for every input file |

Script numbers follow the Supplementary Information rather than running consecutively. s10 is the part of `s9_resampling.py` that meets the two single published profiles with a null, and s11 is a notebook.

---

## Provenance

**Nothing was collected for this paper.** Every measured value comes from a published table or a released file. The three team sports are transcribed from printed tables, because none of those studies released a file. The cortical and musculoskeletal structures are read from files their authors released.

### The two external files

| File | Source | Size | sha256 begins | Licence |
| --- | --- | --- | --- | --- |
| `fig2e.csv` | Source data behind Figure 2E of Murphy et al. (2018), [Zenodo 1069104](https://doi.org/10.5281/zenodo.1069104) | 135 bytes | `516bfd49deeeae4e` | CC BY-SA 4.0 |
| `cov_MOTOR` | `task_covs` release of Weaver et al. (2026), [BrainCompressibility2025](https://github.com/NicholasJWeaver/BrainCompressibility2025) | 80,162 bytes | `e00c39db9c61fada` | none declared |

`fig2e.csv` is released under CC BY-SA 4.0, so it is included here with attribution, and the files derived from it carry the same licence.

`cov_MOTOR` comes from a repository that declares no licence, so it is **not redistributed here**. Fetch it with:

```bash
python data/external/fetch_external.py
```

The script checks that what arrives matches the size and hash above, which is what the Methods record and what the paper was computed from.

For convenience `data/derived/cov_MOTOR_matrix.csv` holds the same 100 by 100 array as plain text. It agrees with the original to better than 5e-13 in every entry. If the upstream licence position changes, the original can be added alongside it.

---

## Reproducibility notes

- **Seeds are fixed and reported.** The simulated ensembles use seeds 42 to 141 at one hundred realisations, and seeds 42 to 91 at fifty for N = 500. Supplementary Table 4 uses seeds 42 to 101 at sixty.
- **Disconnected realisations.** The eigenvector calculation uses the largest connected component, because the leading eigenvector is not unique otherwise, while the degree statistics, the heterogeneity and the assortativity use the full generated graph. The two are read from different objects within one realisation. This is stated because it affects exact reproduction.
- **Rounding.** Eigenvector shares are rounded to twelve decimal places before use. The solver returns values that differ in the last digits between runs, and the figures draw colours from them, so without the rounding the same script renders two images that differ by a level of grey. Twelve decimals leaves eight significant figures on the smallest share in any figure and moves every reported entropy by less than 1e-12 bits.

---

## Two things a reader should know

**The scale-free seed.** The scale-free ensemble grows from a ring lattice seed so that it holds the link budget of the other three, and the size of the seed moves the value. The seed has twelve elements at N = 100 and N = 270, six in the thirteen-element check of Supplementary Note 1, and seven at N = 11 in Supplementary Table 4. At five and six elements no scale-free ensemble exists at the row's density, so those two cells of Supplementary Table 4 are empty. `code/s3_matched_nulls.py` also prints the superseded form grown from a bare seed beside each row, for the record.

**The cortical threshold sweep.** The compression column of `cortical_threshold_sweep.csv` is computed against the logarithm of the number of regions that remain connected at that threshold, not against log2 of 100. The two agree from 1200 links upward, where every region stays connected. The value the paper reports, at every one of the 4950 pairs and with no threshold imposed, is unaffected.

---

## Not included

The code that draws Figures 1 to 4 and 6 is not part of this deposit. It is available from the author on request.

The three analyses specified in advance in Supplementary Note 5 are provided as s4, s5 and s6 and have been checked on simulated data. Their results will be reported once the releases they name have been analysed.

---

## Licence

Code is MIT. See [`LICENSE`](LICENSE).

| Files | Licence |
| --- | --- |
| `muscle_degree_counts.csv`, `musculoskeletal_coupling_shares.csv`, `musculoskeletal_random_null_shares.csv` | CC BY-SA 4.0, inherited from Murphy et al. (2018) |
| All remaining derived files | CC BY 4.0 |

---

## Citing

See [`CITATION.cff`](CITATION.cff), or cite the deposit directly:

> Park, C. (2026). *Average information across a network quantifies what a coordinative structure achieves in motor control. Code and data.* Zenodo. https://doi.org/10.5281/zenodo.21789488
