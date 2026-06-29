# fastCDS — Python client + plotter

Python wrapper around the [fastCDS](https://github.com/SotoLF/fastCDS)
C++ mapper. Provides a friendly Python API plus the `fastCDS plot`
subcommand for rendering domain-overlay isoform figures.

```python
import fastCDS as p2g

mapper = p2g.Mapper("human.idx")
result = mapper.map("ENSP00000269305", aa_start=95, aa_end=288, domain_id="DBD")
result.coding.head()
result.isoform.head()

p2g.plot(result, out="TP53_DBD.png")          # static figure
p2g.plot(result, html="TP53_DBD.html")        # interactive (needs plotly)
```

## Requirements

- The `fastCDS` C++ binary must be available on `PATH` as either
  `fastCDS` or `fastCDS-core`, or pointed to via `$FASTCDS_BIN`.
  Build it from the [main repo](https://github.com/SotoLF/fastCDS)
  or install via `mamba install -c bioconda fastCDS`.
- A binary index produced by `fastCDS --build-index --gtf X --index X.idx`.

## Install

```bash
pip install fastCDS
# or, with interactive HTML output:
pip install "fastCDS[html]"
```

## See also

- Full documentation, CLI examples, output schemas, and CI/test details:
  <https://github.com/SotoLF/fastCDS#readme>
- Validation against ensembldb / TransVar / Ensembl REST: `PHASES.md` in the repo.

License: MIT.
