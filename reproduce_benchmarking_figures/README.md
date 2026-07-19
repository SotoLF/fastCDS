# Reproduce the benchmarking figures (1B, S2, Table S1, S2)


| notebook | figures / tables |
|---|---|
| `notebooks/scaling_and_ram.ipynb` | Fig. 1B, Fig. S2: runtime scaling, throughput, threads x batch-size |
| `notebooks/software_comparison.ipynb` | Table S1: coordinate agreement; Table S2: speed and peak memory |

Each notebook downloads and caches its own inputs, so re-runs reuse what is on disk.
Figures are written to `figures/`.

## Run
```bash
pip install fastCDS pandas numpy scipy matplotlib jupyter
jupyter lab
```

Set the data directory at the top of each notebook (`DATA`, or the `FASTCDS_DATA` env var).

The speed and accuracy comparisons against other tools have their own scripts (and need
R/Bioconductor + TransVar) in `benchmarks/` - see `benchmarks/README.md`.
