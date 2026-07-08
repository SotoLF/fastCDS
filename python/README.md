# fastCDS - Python client + plotter

Python wrapper around the [fastCDS](https://github.com/SotoLF/fastCDS)
C++ mapper. Provides a friendly Python API plus the `fastCDS plot`
subcommand for rendering domain-overlay isoform figures.

```python
import fastCDS as fc

mapper = fc.Mapper("human.idx")
result = mapper.map("ENSP00000269305", aa_start=95, aa_end=288, domain_id="DBD")
result.summary.head()        # one row per query
result.isoform.head()        # plot-ready isoform structure
result.cds_segments.head()   # the exact coding blocks

fc.plot(result, out="TP53_DBD.png")                    # static figure
fc.plot(result, out="TP53_DBD.html")                   # interactive (vanilla JS)
fc.plot(result, out="TP53_DBD.html", engine="plotly")  # interactive (needs plotly)
```

## Requirements

- The `fastCDS` C++ binary must be available on `PATH` as either
  `fastCDS` or `fastCDS-core`, or pointed to via `$FASTCDS_BIN`.
  Build it from the [main repo](https://github.com/SotoLF/fastCDS)
  or install via `mamba install -c bioconda fastCDS`.
- A binary index produced by `fastCDS index --gtf X --out X.idx`.

## Install

```bash
pip install fastCDS
# the vanilla-JS interactive viewer works out of the box;
# add the plotly engine only if you want it:
pip install "fastCDS[plotly]"
```

## See also

- Full documentation, CLI examples, output schemas, and CI/test details:
  <https://github.com/SotoLF/fastCDS#readme>
- Validation against ensembldb / TransVar / Ensembl REST: `PHASES.md` in the repo.

License: MIT.
