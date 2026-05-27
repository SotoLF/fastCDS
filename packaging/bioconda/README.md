# Bioconda recipe — submission checklist

This directory holds the recipe that should be submitted to
[`bioconda/bioconda-recipes`](https://github.com/bioconda/bioconda-recipes).

The recipe is **not yet submitted** — it requires:

1. A tagged release on GitHub (`git tag v2.2.0 && git push --tags`) so the
   `source.url` resolves.
2. Updating `meta.yaml`'s `sha256:` with the actual tarball hash:
   ```bash
   curl -sL https://github.com/SotoLF/Prot2Exon/archive/refs/tags/v2.2.0.tar.gz \
       | shasum -a 256
   ```
3. A fork of `bioconda/bioconda-recipes` and a PR adding this directory at
   `recipes/prot2exon/`.

## Local test

```bash
git clone https://github.com/bioconda/bioconda-recipes
cp -r packaging/bioconda/prot2exon/ bioconda-recipes/recipes/
cd bioconda-recipes
bioconda-utils build --packages prot2exon
```

After the PR merges:

```bash
mamba install -c bioconda prot2exon
prot2exon --help
```
