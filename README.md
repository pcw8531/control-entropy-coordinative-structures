Average information across a network measures what a coordinative structure buys in motor control

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21789488.svg)](https://doi.org/10.5281/zenodo.21789488)

Code and data for the paper. Everything printed in the manuscript, the Supplementary Information and the figure captions is recomputed here from the released files.

Chulwook Park · Department of Physical Education, Seoul National University · ORCID 0000-0001-8714-5760

Quick start
bash
pip install -r requirements.txt
python code/verify_published_values.py

That script checks the size and hash of every released file, recomputes every value the paper prints, and exits non-zero if anything fails. It is the gate every other script runs behind.

On a clean checkout it reports 75 of 75 values reproducing.

Reproducing the individual tables
Script	Produces
code/s1_ensembles.py	Table 1, Figure 3, Supplementary Tables 1 to 5
code/s2_measured.py	Supplementary Notes 4 and 5, Supplementary Tables 6 and 7
code/s3_matched_nulls.py	Supplementary Table 8
code/s9_resampling.py	Supplementary Table 9, and the lower row of Figure 4
<!-- CHECK before publishing: the SI now numbers the matched-null table as Supplementary Table 6 (four topology columns, sixty realisations per cell), the measured-structure table as 7 and the outcome-entropy table as 8. Run each script and confirm which table it writes. -->

Every value the paper prints comes out of those scripts. Nothing has to be run in a particular order beyond the verification gate.

Repository layout
Path	What it holds
code/core.py	The equation module, Eq 1 to 13 in the Methods numbering. Every number in the paper traces to one definition here.
code/generators.py	The four matched topology ensembles
code/loaders.py	Readers and the provenance record
code/s*.py	The analyses, named as the Supplementary Information names them
code/verify_published_values.py	The gate. Hashes every input and recomputes every published value.
data/derived/	Nine CSVs, each computed from a released file or transcribed from a published table. Each carries a comment header saying where it came from.
data/external/	Files belonging to other groups. See below.
SOURCES.md	Full provenance for every input file
<!-- CHECK before publishing: the manuscript now labels equations 1 to 19, and cites Eq 14, 16, 18 and 19 by name. Confirm how far core.py actually goes and correct "Eq 1 to 13" if needed. -->
Provenance

Nothing was collected for this paper. Every measured value comes from a published table or a released file. The three team sports are transcribed from printed tables, because none of those studies released a file. The cortical and musculoskeletal structures are read from files their authors released.

The two external files
File	Source	Size	sha256 begins	Licence
fig2e.csv	Source data behind Figure 2E of Murphy et al. (2018), Zenodo 1069104	135 bytes	516bfd49deeeae4e	CC BY-SA 4.0
cov_MOTOR	task_covs release of Weaver et al. (2026), BrainCompressibility2025	80,162 bytes	e00c39db9c61fada	none declared

fig2e.csv is released under CC BY-SA 4.0, so it is included here with attribution, and the files derived from it carry the same licence.

cov_MOTOR comes from a repository that declares no licence, so it is not redistributed here. Fetch it with:

bash
python data/external/fetch_external.py

The script checks that what arrives matches the size and hash above, which is what the Methods record and what the paper was computed from.

For convenience data/derived/cov_MOTOR_matrix.csv holds the same 100 by 100 array as plain text. It agrees with the original to better than 5e-13 in every entry. If the upstream licence position changes, the original can be added alongside it.

Reproducibility notes
Seeds are fixed and reported. The simulated ensembles use seeds 42 to 141 at one hundred realisations, and Supplementary Table 8 uses seeds 42 to 101 at sixty. <!-- CHECK: sixty realisations is the matched-null table, which the current SI numbers as Table 6. -->
Disconnected realisations. The eigenvector calculation uses the largest connected component, because the leading eigenvector is not unique otherwise, while the degree statistics, the heterogeneity and the assortativity use the full generated graph. The two are read from different objects within one realisation. This is stated because it affects exact reproduction.
Rounding. Eigenvector shares are rounded to twelve decimal places before use. The solver returns values that differ in the last digits between runs, and the figures draw colours from them, so without the rounding the same script renders two images that differ by a level of grey. Twelve decimals leaves eight significant figures on the smallest share in any figure and moves every reported entropy by less than 1e-12 bits.
Two things a reader should know

The scale-free column of Supplementary Table 8. <!-- CHECK: table number, see above --> The values printed in that column of the current Supplementary Information come from a Barabasi-Albert graph grown from a bare seed at three links per node, which does not hold the link budget the other three columns hold. code/s3_matched_nulls.py prints both that form and the matched form of Supplementary Note 2 side by side, so the difference is on the record rather than buried.

The cortical threshold sweep. The compression column of cortical_threshold_sweep.csv is computed against the logarithm of the number of regions that remain connected at that threshold, not against log2 of 100. The two agree from 1200 links upward, where every region stays connected. The value the paper reports, at every one of the 4950 pairs and with no threshold imposed, is unaffected.

Not included

The code that draws the figures and the movie is not part of this deposit. What it does is arrange values on a canvas, and every one of those values is computed by the scripts above and printed by verify_published_values.py, so a reader can check any number in any figure without it. The drawing code is available from the author on request.

Two Supplementary figures are waiting on inputs that their owners have not released. The analyses those figures belong to are specified in advance in the Supplementary Information with their acceptance criteria fixed. They are not in the deposit because there is nothing yet to compute them from.

Licence

Code is MIT. See LICENSE.

Files	Licence
muscle_degree_counts.csv, musculoskeletal_coupling_shares.csv, musculoskeletal_random_null_shares.csv	CC BY-SA 4.0, inherited from Murphy et al. (2018)
All remaining derived files	CC BY 4.0
Citing

See CITATION.cff, or cite the deposit directly:

Park, C. (2026). Average information across a network measures what a coordinative structure buys in motor control. Code and data. Zenodo. https://doi.org/10.5281/zenodo.21789488
