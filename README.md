# Analysis scripts: FLI embryo assessment

Analysis code for the Advanced Science manuscript. The repository is organised as
one branch per assay; the scripts on each branch sit in its root directory.

| Branch | Content |
|---|---|
| `RNAseq` | upstream processing of the RNA-seq reads: QC, trimming, alignment to GRCm39, gene and transcript counts |
| `Methylation` | CpG data organisation, global methylation, candidate regions, region-level tests, annotation and figures |
| `Proteomics` | quantitative matrix processing, QC, differential abundance, enrichment and figures |

This branch (`main`) holds this README and the MIT licence only.

## Common conventions

Each script sets its data locations in a block at the top of the file, or takes them
as command-line arguments. Run the scripts from the directory that holds the data.
No credential, token or reviewer access code is stored in this repository.

## Reproduction starting points

* RNA-seq starts from raw FASTQ files.
* Methylation starts from existing per-sample CpG counts; alignment and methylation
  calling are not part of this repository.
* Proteomics starts from the processed quantitative matrix. The raw mass-spectrometry
  files have not been provided by the vendor, so identification, quantification and
  acquisition QC cannot be re-checked here.

## Data availability

Sequencing data and the proteomics matrix are shared separately from this code; see
the README on the relevant branch. Data links and accessions are inserted by the
authors once confirmed.

## Licence

MIT, see `LICENSE`. The licence covers this code only and does not extend to the data
or to third-party resources.
