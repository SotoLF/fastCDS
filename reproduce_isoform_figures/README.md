# Reproduce the isoform / domain figures (1C, 1D, 1E, 1F)


| step | file | output |
|---|---|---|
| 00 | `00_pull_biomart_pfam.py` | `pfam_human_v115_meta.tsv`, `pfam_human_v115.bed` |
| 01 | `01_pfam_catalytic.py` | `pfam_catalytic.tsv` |
| 02 | `02_build_isoform_universe.py` | `isoforms.tsv`, `genes.tsv` |
| 03 | `03_domain_isoform_master.py` | `domain_isoform_master.tsv`, `domain_genomic_intervals.tsv` |
| - | `Fig1C_intact_skipped_trimmed.ipynb` | Fig. 1C intact / skipped / trimmed bar |
| - | `Fig1D_two_modes_scatter.ipynb` | Fig. 1D skipped fraction vs exon count |
| - | `Fig1E1F_symmetric_exons_and_gsea.ipynb` | Fig. 1E,1F frame-preserving exons + catalytic GSEA |

Needs the `fastCDS` python package. `pfam_catalytic.tsv` is already included (step 01 rebuilds it).

## 1. Download source data
```bash
# Ensembl v115 GTF
wget https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz

# Pfam-A clans (Pfam 38.1)
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.1/Pfam-A.clans.tsv.gz
gunzip Pfam-A.clans.tsv.gz

# Pfam -> GO map (for the catalytic list)
wget https://current.geneontology.org/ontology/external2go/pfam2go
```

## 2. Build the fastCDS index
```bash
fastCDS index --gtf Homo_sapiens.GRCh38.115.gtf.gz --out ensembl_v115_human.idx
```

## 3. Prepare the Pfam inputs
```bash
conda activate base

# Pfam-A domain instances on the human proteome, from Ensembl BioMart.
# The default host serves the current release; to pin release 115 pass its
# archive host, e.g. --host https://<release-115-archive>.ensembl.org/biomart/martservice
python 00_pull_biomart_pfam.py \
    --out-meta pfam_human_v115_meta.tsv --out-bed pfam_human_v115.bed

# catalytic Pfam families (already bundled as pfam_catalytic.tsv; this rebuilds it)
python 01_pfam_catalytic.py --pfam2go pfam2go --out pfam_catalytic.tsv
```

## 4. Build the tables
```bash
python 02_build_isoform_universe.py \
    --gtf Homo_sapiens.GRCh38.115.gtf.gz \
    --out-isoforms isoforms.tsv --out-genes genes.tsv

python 03_domain_isoform_master.py \
    --isoforms isoforms.tsv --genes genes.tsv \
    --pfam-meta pfam_human_v115_meta.tsv --clans Pfam-A.clans.tsv \
    --index ensembl_v115_human.idx --threads 4 \
    --out-master domain_isoform_master.tsv \
    --out-intervals domain_genomic_intervals.tsv
```

## 5. Make the figures
```bash
jupyter lab    # run the three Fig notebooks
```

Set the data path at the top of each notebook to where step 03 wrote the tables.
