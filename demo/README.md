# Demo — mini reproduction of Figure 6H

Minimal, self-contained example that reproduces a small version of **Figure 6H**
(MILC + |ΔGC3| bubble plot) from **Wu Z. *et al.*** (2026), on a subset of the
data underlying Supplementary Table 18.

The demo runs in **under 30 seconds** on a standard laptop and requires only
three CRAN packages. No external databases, no BLAST, no genome downloads.

---

## Files in this folder

| File | Purpose |
|---|---|
| `demo_input.tsv` | Small (4 targets × 4 references) subset of Supplementary Table 18 in long format: target, reference, reference_type, MILC, \|ΔGC3\| (percentage points). |
| `run_demo.R` | Self-contained R script that reads `demo_input.tsv` and produces the bubble plot. |
| `expected_output.png` | Reference plot the demo should produce (generated once by the authors from a successful run). Compare your run against this file to confirm the code works. |
| `README.md` | This file. |

---

## Software requirements

- **R** ≥ 4.0 (tested on R 4.5.0)
- **R packages** (install once):

```r
install.packages(c("ggplot2", "dplyr", "readr"))
```

No other dependencies. No compiled tools. No internet access needed at run time.

---

## How to run

From the repository root:

```bash
Rscript demo/run_demo.R
```

Or from inside the `demo/` folder:

```bash
cd demo
Rscript run_demo.R
```

Or interactively in RStudio: open `run_demo.R` and hit **Run → Run All**.

**Expected runtime**: < 30 seconds.
**Expected output**: `expected_output.png` — a bubble plot with 4 target rows
(*Msmithii_tetC*, *Mmazei_TEM*, *Mmazei_catA*, *Sislandicus_copR*) and 4 columns
(archaeal host + bacterial comparators). Bubble size encodes MILC; bubble colour
encodes |ΔGC3| (percentage points).

---

## What the demo demonstrates

- The **data model** used throughout the paper: each row is a `(target, reference)` comparison with two composition-divergence metrics (MILC and |ΔGC3|).
- The **plotting logic** used to build Figure 6H (`ggplot2`, bubble size = MILC, fill = |ΔGC3|).
- A quick visual check that the codon-usage divergence framework separates:
  - **host-typical** loci (small bubbles, light fill against their archaeal host — e.g., *copR* vs. *S. islandicus*),
  - **host-atypical / bacterial-affine** loci (small bubbles, light fill against a bacterial reference — e.g., *tetC* vs. *A. salmonicida*).

For the full analysis (all 19 prioritised loci, 21 bacterial comparators, RSCU
PCA, Ma, global rank-based test), see `Genomics/codon analysis.R` and
`GC_MILC.R` in the repository root.

---

## Reference

Wu Z., Godoy-Santos F., Sabino Y.N.V., Huws S.A., Oyama L.B. (2026).
*A domain-wide survey of archaeal antimicrobial resistance genes reveals an
ecologically structured and compositionally heterogeneous resistome.*
DOI: [to be added on publication]
