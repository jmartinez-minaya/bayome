# bayome: Bayesian Occurrence-Abundance Microbiome Estimation

Bayesian models for microbiome association studies, fitted with
'R-INLA'. The occurrence (presence/absence) and the abundance of each
feature are modelled together: abundance effects are linked to
occurrence effects and information is shared across features. Users
write the same formulas as in `maaslin3()`; bayome translates them into
INLA models. See
[`bayome()`](https://jmartinez-minaya.github.io/bayome/reference/bayome.md),
[`bm_prior()`](https://jmartinez-minaya.github.io/bayome/reference/bm_prior.md),
[`bm_parse_formula()`](https://jmartinez-minaya.github.io/bayome/reference/bm_parse_formula.md)
and
[`bm_inla_formula()`](https://jmartinez-minaya.github.io/bayome/reference/bm_inla_formula.md);
results are plotted with the functions in
[bm_plot](https://jmartinez-minaya.github.io/bayome/reference/bm_plot.md).

## Details

bayome is an independent package, not part of 'MaAsLin3' or 'bioBakery'.
Its interface, pre-processing, two-part structure, median comparison and
output format follow 'MaAsLin3' (MIT license); the Bayesian models, the
joint occurrence-abundance link, the information sharing across features
and the selection rules are new.

## References

Nickols WA, Kuntz T, Shen J, Maharjan S, Mallick H, Franzosa EA,
Thompson KN, Nearing JT, Huttenhower C (2026). MaAsLin 3: refining and
extending generalized multivariable linear models for meta-omic
association discovery. Nature Methods, 23, 554-564.
[doi:10.1038/s41592-025-02923-9](https://doi.org/10.1038/s41592-025-02923-9)

## See also

Useful links:

- <https://jmartinez-minaya.github.io/bayome/>

- <https://github.com/jmartinez-minaya/bayome>

- Report bugs at <https://github.com/jmartinez-minaya/bayome/issues>

## Author

**Maintainer**: Joaquín Martínez-Minaya
<jmartinezminaya@hsph.harvard.edu>

Authors:

- Joaquín Martínez-Minaya <jmartinezminaya@hsph.harvard.edu>
