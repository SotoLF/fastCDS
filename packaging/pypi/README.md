# PyPI upload — dry-run output + checklist

Both distributions have been built and validated locally with `twine check`:

```
prot2exon-2.2.0-py3-none-any.whl  (≈19 KB, pure-Python wheel)
prot2exon-2.2.0.tar.gz            (≈19 KB, sdist)
```

Both `PASSED` `twine check`. A throwaway install from the wheel imports
`prot2exon.Mapper` and `prot2exon.plot` successfully.

## Build the artifacts yourself

```bash
cd python/
python -m pip install --upgrade build twine
rm -rf dist build *.egg-info
python -m build --sdist --wheel
python -m twine check dist/*
```

## Upload (requires your PyPI API token)

1. Create an API token at https://pypi.org/manage/account/token/ (scope it
   to the `prot2exon` project once the first release exists).
2. Test first on TestPyPI:
   ```bash
   python -m twine upload --repository testpypi dist/*
   pip install --index-url https://test.pypi.org/simple/ prot2exon
   ```
3. Real upload:
   ```bash
   python -m twine upload dist/*
   ```

## Why these are tiny

The wheel only contains the **Python wrapper** (`_client.py`, `plot.py`,
`_result.py`, `__init__.py`). The C++ binary is shipped separately via
GitHub releases, `docker pull prot2exon`, or `mamba install -c bioconda
prot2exon`. The wrapper auto-discovers the binary on `PATH` or via
`$PROT2EXON_BIN`.
