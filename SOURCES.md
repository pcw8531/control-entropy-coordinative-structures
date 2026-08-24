# Source manifest

Every input file used in *Average information across a network measures what a coordinative structure buys in motor control*.

**No data were collected for this study.** Each entry below gives the origin, the access route, the file as received, and the licence under which it is used. Reference numbers are those of the manuscript reference list.

---

## Inputs

| Ref | Provides | Source | Accessed as | Licence |
| --- | --- | --- | --- | --- |
| 27 | Football, position centrality, 2014 World Cup knockout stage | Korte & Lames (2018), *Curr. Issues Sport Sci.* **3**, 005 | Transcribed from a printed table. No file was released. | *to be completed* |
| 28 | Basketball, position centrality, 2010 season, published team by team | Fewell, Armbruster, Ingraham, Petersen & Waters (2012), *PLoS ONE* **7**, e47445 | Transcribed from a printed table. No file was released. | *to be completed* |
| 29 | Cortical activity, Human Connectome Project motor task, 100 parcels | Weaver, Faskowitz, Betzel & Lynn (2026), *PNAS* **123**, e2531115123 | `cov_MOTOR` from the `task_covs` release at [BrainCompressibility2025](https://github.com/NicholasJWeaver/BrainCompressibility2025) | none declared at source |
| 30 | Musculoskeletal network, 173 bones and 270 muscles | Murphy, Muldoon, Baker, Lastowka, Bennett, Yang & Bassett (2018), *PLoS Biol.* **16**, e2002811 | `fig2e.csv` from [Zenodo 1069104](https://doi.org/10.5281/zenodo.1069104) | CC BY-SA 4.0 |
| 31 | Expert and novice outcome dispersion, perceptual accuracy task | Park (2026), *Inf. Process. Manag.* **64**, 104951 | Per-participant standard deviations read from the published table. | *to be completed* |
| 44 | Handball, position centrality, positional attack passing analysis | Korte & Lames (2019), *J. Hum. Kinet.* **70**, 209–221 | Transcribed from a printed table. No file was released. | *to be completed* |

Entries marked *to be completed* still need the page or table transcribed, the date of access, and the licence or permission basis.

---

## File-level records for the two released files

| File | Size | sha256 begins | Redistributed here |
| --- | --- | --- | --- |
| `fig2e.csv` | 135 bytes | `516bfd49deeeae4e` | Yes, under CC BY-SA 4.0, with attribution |
| `cov_MOTOR` | 80,162 bytes | `e00c39db9c61fada` | No. The source repository declares no licence. |

Both were taken unmodified.

`cov_MOTOR` is fetched rather than shipped:

```bash
python data/external/fetch_external.py
```

The script verifies the size and hash above before use. For convenience `data/derived/cov_MOTOR_matrix.csv` holds the same 100 by 100 array as plain text under CC BY 4.0, agreeing with the original to better than 5e-13 in every entry.

---

## Derived files

Nine derived CSV files are produced from the inputs above by the scripts in `code/`. Each carries a comment header naming its origin.

```bash
python code/verify_published_values.py
```

That gate checks the size and hash of every released input, recomputes every value the paper prints, and exits non-zero if anything fails. A clean checkout reports **75 of 75 values reproducing**.

---

## Licences

| Files | Licence |
| --- | --- |
| All code | MIT |
| `muscle_degree_counts.csv`, `musculoskeletal_coupling_shares.csv`, `musculoskeletal_random_null_shares.csv` | CC BY-SA 4.0, inherited from Murphy et al. (2018) |
| All remaining derived files | CC BY 4.0 |

---

## Not an input

StatsBomb open data (manuscript ref. 42) is named in the Methods as the release a specified-in-advance analysis will use. No value reported in this paper is computed from it, and it is not part of this deposit.
