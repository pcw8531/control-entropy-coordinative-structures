# Average information across a network measures what a coordinative structure buys in motor control

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21789489.svg)](https://doi.org/10.5281/zenodo.21789489)

Code and data for the paper. Everything printed in the manuscript, the SI
Appendix and the figure captions is recomputed here from the released files.

Chulwook Park, Department of Physical Education, Seoul National University.
ORCID 0000-0001-8714-5760.

## Start here

```
pip install -r requirements.txt
python code/verify_published_values.py
```

That script checks the size and hash of every released file, recomputes every
value the paper prints, and exits non-zero if anything fails. It is the gate
every other script runs behind. On a clean checkout it reports 75 of 75 values
reproducing.

Then, for the individual tables:

```
python code/s1_ensembles.py      # Table 1, Figure 3, SI Tables S1 to S5
python code/s2_measured.py       # SI section 4 and 5, Tables S6 and S7
python code/s3_matched_nulls.py  # SI Table S8
python code/s9_resampling.py     # SI Table S9, and the lower row of Figure 4
```

Every value the paper prints comes out of those five scripts. Nothing in the
deposit has to be run in a particular order beyond the verification gate.

## What is in here

    code/core.py              the equation module, Eq 1 to 13 in the Methods
                              numbering. Every number in the paper traces to one
                              definition here.
    code/generators.py        the four matched topology ensembles
    code/loaders.py           readers and the provenance record
    code/s*.py                the analyses, named as the SI names them
    code/verify_published_values.py
                              the gate. Hashes every input and recomputes every
                              published value.

    data/derived/             nine CSVs, every one computed from a released file
                              or transcribed from a published table. Each carries
                              a comment header saying where it came from.
    data/external/            files belonging to other groups. See below.

## Nothing was collected for this paper

Every measured value comes from a published table or a released file. The three
team sports are transcribed from printed tables, because none of those studies
released a file. The cortical and musculoskeletal structures are read from files
their authors released.

## The two external files

**fig2e.csv**, the source data behind Figure 2E of Murphy et al. (2018), from
Zenodo record 1069104. Released under CC BY-SA 4.0, so it is included here with
attribution, and the files derived from it carry the same licence. 135 bytes,
sha256 beginning `516bfd49deeeae4e`.

**cov_MOTOR**, from the `task_covs` release of Weaver et al. (2026), at
github.com/NicholasJWeaver/BrainCompressibility2025. That repository declares no
licence, so the file is **not redistributed here**. Run

```
python data/external/fetch_external.py
```

to obtain it. The script checks that what arrives is 80,162 bytes with sha256
beginning `e00c39db9c61fada`, which is what the Methods records and what the
paper was computed from.

For convenience `data/derived/cov_MOTOR_matrix.csv` holds the same 100 by 100
array as plain text. It agrees with the original to better than 5e-13 in every
entry. If the upstream licence position changes, the original can be added
alongside it.

## Reproducibility notes

Every random seed is fixed and reported. The simulated ensembles use seeds 42 to
141 at one hundred realisations, and SI Table S8 uses seeds 42 to 101 at sixty.

Where a realisation is disconnected the eigenvector calculation uses the largest
connected component, because the leading eigenvector is not unique otherwise,
while the degree statistics, the heterogeneity and the assortativity use the full
generated graph. The two are read from different objects within one realisation.
This is stated because it affects exact reproduction.

Eigenvector shares are rounded to twelve decimal places before use. The solver
returns values that differ in the last digits between runs, and the figures draw
colours from them, so without the rounding the same script renders two images
that differ by a level of grey. Twelve decimals leaves eight significant figures
on the smallest share in any figure and moves every reported entropy by less than
1e-12 bits.

## Two things a reader should know

**SI Table S8, the scale-free column.** The values printed in that column of the
current SI come from a Barabási–Albert graph grown from a bare seed at three
links per node, which does not hold the link budget the other three columns hold.
`code/s3_matched_nulls.py` prints both that form and the matched form of SI
section 2, side by side, so the difference is on the record rather than buried.

**The cortical threshold sweep.** The compression column of
`cortical_threshold_sweep.csv` is computed against the logarithm of the number of
regions that remain connected at that threshold, not against log2 of 100. The two
agree from 1200 links upward, where every region stays connected. The value the
paper reports, at every one of the 4950 pairs and with no threshold imposed, is
unaffected.

## What is not in here

The code that draws the figures and the movie is not part of this deposit. What
it does is arrange values on a canvas, and every one of those values is computed
by the scripts above and printed by `verify_published_values.py`, so a reader can
check any number in any figure without it. The drawing code is available from the
author on request.

Two SI figures are also waiting on inputs that their owners have not released,
and the analyses those figures belong to are specified in advance in the SI with
their acceptance criteria fixed. They are not in the deposit because there is
nothing yet to compute them from.

## Licence

Code is MIT. See `LICENSE`.

The files derived from `fig2e.csv`, namely `muscle_degree_counts.csv`,
`musculoskeletal_coupling_shares.csv` and
`musculoskeletal_random_null_shares.csv`, inherit CC BY-SA 4.0 from Murphy et al.
(2018). The remaining derived files are CC BY 4.0.

## Citing

See `CITATION.cff`.
