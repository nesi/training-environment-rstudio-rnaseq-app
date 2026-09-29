# NeSI training environment RNA-Seq RStudio app

RStudio app for the RNA-Seq workshop running on the NeSI training environment.

## Datasets

The image carries two workshop datasets under `/var/lib/rnaseq/`. The **Workshop dataset** option on the launch form picks one, and `stage-rnaseq-data` copies it into `~/RNA_seq` at session start. Existing files are never overwritten. If `~/RNA_seq` already holds the other dataset, it is moved to `~/RNA_seq_<dataset>_<timestamp>` first.

| Form value | Workshop | Source |
|---|---|---|
| `nz` (default) | [RNA-seq data analysis workflow (NZ spotty wrasse)](https://genomicsaotearoa.github.io/RNA-seq-data-analysis-workflow-NZ/) | `RNA_seq.zip` from the [2026-September release](https://github.com/GenomicsAotearoa/RNA-seq-data-analysis-workflow-NZ/releases/tag/2026-September), downloaded and checksummed at build time |
| `yeast` | [RNA-seq workshop (yeast)](https://genomicsaotearoa.github.io/RNA-seq-workshop/) | bundled in `docker/RNA_seq/` |

The build rebuilds the NZ backup STAR index (`--sjdbOverhang 124 --genomeSAindexNbases 11`), because the one in the zip is 2 GB. It lives at `/var/lib/rnaseq/nz-genomeindex`, and `~/RNA_seq/backup-files/genomeindex` symlinks to it rather than holding a copy in each home directory.
