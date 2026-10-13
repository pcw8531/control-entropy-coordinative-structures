# Source manifest

Every input file used in *Average information across a network quantifies what a coordinative structure achieves in motor control*.

**No data were collected for this study.** Each entry below gives the origin, the access route, the file as received, and the licence or basis under which it is used.

---

## Inputs

| Provides | Source | Accessed as | Licence or basis of use |
| --- | --- | --- | --- |
| Football, position centrality, 2014 World Cup knockout stage | Korte & Lames (2018), *Curr. Issues Sport Sci.* **3**, 005 | Transcribed from a printed table. No file was released. | Published values, used with citation |
| Basketball, position centrality, 2010 play-offs, published team by team | Fewell, Armbruster, Ingraham, Petersen & Waters (2012), *PLoS ONE* **7**, e47445 | Transcribed from a printed table. No file was released. | Published values, used with citation |
| Cortical activity, Human Connectome Project motor task, 100 parcels | Weaver, Faskowitz, Betzel & Lynn (2026), *PNAS* **123**, e2531115123 | `cov_MOTOR` from the `task_covs` release at [BrainCompressibility2025](https://github.com/NicholasJWeaver/BrainCompressibility2025) | none declared at source |
| Musculoskeletal network, 173 bones and 270 muscles | Murphy, Muldoon, Baker, Lastowka, Bennett, Yang & Bassett (2018), *PLoS Biol.* **16**, e2002811 | `fig2e.csv` from [Zenodo 1069104](https://doi.org/10.5281/zenodo.1069104) | CC BY-SA 4.0 |
| Expert and novice outcome dispersion, perceptual accuracy task | Park (2026a), *Inf. Process. Manag.* **64**, 104951 | Per-participant standard deviations read from Table D.1 of the published article. | Published values, used with citation |
| Handball, position centrality, positional attack passing analysis | Korte & Lames (2019), *J. Hum. Kinet.* **70**, 209–221 | Transcribed from Table 1 of the article and checked against the original on 28 July 2026. No file was released. | Published values, used with citation |

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

Ten derived CSV files are released. Nine are computed from the inputs above or transcribed from them, and each carries a comment header naming its origin. The tenth, `s11_controllability_checks.csv`, is written by the s11 notebook from the generated ensembles of Table 1.

```bash
python code/verify_published_values.py
```

That gate checks the size and hash of every released input, recomputes 75 published values behind Sections 3.1 to 3.3 and the behavioural check in Section 3.4, and exits non-zero if anything fails. A clean checkout reports **75 of 75 values reproducing**.

---

## Licences

| Files | Licence |
| --- | --- |
| All code | MIT |
| `muscle_degree_counts.csv`, `musculoskeletal_coupling_shares.csv`, `musculoskeletal_random_null_shares.csv` | CC BY-SA 4.0, inherited from Murphy et al. (2018) |
| All remaining derived files | CC BY 4.0 |

---

## Not an input

StatsBomb open data (StatsBomb, n.d.) is named in the Methods as the release a specified-in-advance analysis will use. No value reported in this paper is computed from it, and it is not part of this deposit.
