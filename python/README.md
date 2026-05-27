# prot2exon — Python client + plotter

Python wrapper around the [prot2exon](https://github.com/SotoLF/Prot2Exon)
C++ mapper. Provides a friendly Python API plus the `prot2exon plot`
subcommand for rendering domain-overlay isoform figures.

```python
import prot2exon as p2g

mapper = p2g.Mapper("human.idx")
result = mapper.map("ENSP00000269305", aa_start=95, aa_end=288, domain_id="DBD")
result.coding.head()
result.isoform.head()

p2g.plot(result, out="TP53_DBD.png")          # static figure
p2g.plot(result, html="TP53_DBD.html")        # interactive (needs plotly)
```

## Requirements

- The `prot2exon` C++ binary must be available on `PATH` as either
  `prot2exon` or `prot2exon-core`, or pointed to via `$PROT2EXON_BIN`.
  Build it from the [main repo](https://github.com/SotoLF/Prot2Exon)
  or install via `mamba install -c bioconda prot2exon`.
- A binary index produced by `prot2exon --build-index --gtf X --index X.idx`.

## Install

```bash
pip install prot2exon
# or, with interactive HTML output:
pip install "prot2exon[html]"
```

## See also

- Full documentation, CLI examples, output schemas, and CI/test details:
  <https://github.com/SotoLF/Prot2Exon#readme>
- Validation against ensembldb / TransVar / Ensembl REST: `PHASES.md` in the repo.

License: MIT.
