# Example data

`HMP2_taxonomy.tsv` and `HMP2_metadata.tsv` are the example data shipped with
MaAsLin3 (<https://github.com/biobakery/maaslin3>, `inst/extdata`, MIT
license, copyright the MaAsLin3 authors), copied without changes. They come
from the Integrative Human Microbiome Project / Inflammatory Bowel Disease
Multi'omics Database (HMP2 / IBDMDB):

- Lloyd-Price J, Arze C, Ananthakrishnan AN, et al. Multi-omics of the gut
  microbial ecosystem in inflammatory bowel diseases. *Nature* 569, 655–662
  (2019). doi:10.1038/s41586-019-1237-9
- Nickols WA, Kuntz T, Shen J, et al. MaAsLin 3: refining and extending
  generalized multivariable linear models for meta-omic association
  discovery. *Nature Methods* 23, 554–564 (2026).
  doi:10.1038/s41592-025-02923-9

- `HMP2_taxonomy.tsv`: relative abundances of 151 species (columns) in 1527
  stool samples (rows).
- `HMP2_metadata.tsv`: metadata of the same samples, including the
  participant (`participant_id`, 120 participants with repeated samples),
  `diagnosis` (CD, UC, nonIBD), `dysbiosis_state`, `antibiotics`, `age` and
  sequencing depth (`reads`).

`HMP2_bayome_results.tsv` and `HMP2_bayome_lambda.tsv` are the results of
`data-raw/hmp2_example.R` (bayome fitted to these data), read by the vignette.
