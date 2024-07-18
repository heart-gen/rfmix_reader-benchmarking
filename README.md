# `RFMix-reader`: Accelerated reading and processing for local ancestry studies

This repository contains the code for the `RFMix-reader` manuscript that
describes the [RFMix-reader](https://github.com/heart-gen/rfmix_reader) software.

## Abstract
**Motivation:** Local ancestry inference is a powerful technique in genetics,
revealing population history and the genetic basis of diseases. It is
particularly valuable for improving eQTL discovery and fine-mapping in admixed
populations. Despite the widespread use of the `RFMix` software for local
ancestry inference, large-scale genomic studies face challenges of high memory
consumption and processing times when handling `RFMix` output files.<br />
**Results:** Here, I present `RFMix-reader`, a new Python-based parsing software,
designed to streamline the analysis of large-scale local ancestry datasets. This
software prioritizes computational efficiency and memory optimization, leveraging
GPUs when available for additional speed boosts. By overcoming these data
processing hurdles, `RFMix-reader` empowers researchers to unlock the full
potential of local ancestry data for understanding human health and health
disparities.<br /> **Availability:** `RFMix-reader` is freely available on PyPI
at
[https://pypi.org/project/rfmix-reader/](https://pypi.org/project/rfmix-reader/),
implemented in Python 3, and supported on Linux, Windows, and Mac OS.<br />
**Contact:** [KynonJade.Benjamin@libd.org](KynonJade.Benjamin@libd.org)<br />
**Supplementary information:** Supplementary data are available at
[https://rfmix-reader.readthedocs.io/en/latest/](https://rfmix-reader.readthedocs.io/en/latest/).


## Data availability
Analysis-ready genotype data will be shared with researchers who obtain database
of Genotypes and Phenotype (dbGaP) access
([phs000979.v3.p2](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000979.v3.p2)).
The 1000 Genomes Project reference data are available at
<http://www.internationalgenome.org/data/>. The simulated data used in the
benchmarking analyses are available on Synapse
[syn61691659](https://www.synapse.org/Synapse:syn61691659). The real data will be
shared with researchers who obtain access to the genotype data (phs000979.v3.p2)
upon request. All scripts generated in this paper are available through GitHub at
<https://github.com/heart-gen/rfmix_reader-benchmarking>.

## Authors
* [Kynon JM Benjamin](https://github.com/Krotosbenjamin)

## Funding

This work was supported by grants from the National Institutes of Health,
National Institute on Minority Health and Health Disparities (NIMHD) K99MD016964.

## Citation

If you the `RFMix-reader` software or anything in this repository please cite the
following pre-print: <https://doi.org/10.1101/2024.07.13.603370>.

Kynon JM Benjamin. "RFMix-reader: Accelerated reading and processing for
local ancestry studies." *bioRxiv*. 2024.
DOI: [10.1101/2024.07.13.603370](https://www.biorxiv.org/content/10.1101/2024.07.13.603370v2).

## License

<img src="https://licensebuttons.net/l/by-nc/3.0/88x31.png" alt width="88"
height="31" scale="0"> Attribution-NonCommercial: CC BY-NC

This license lets others remix, tweak, and build upon our work non-commercially
as long as they acknowledge our work.

[View License Deed](https://creativecommons.org/licenses/by-nc/4.0) | [View Legal
Code](https://creativecommons.org/licenses/by-nc/4.0/legalcode)
