<p align="center">
  <a href="" rel="noopener">
 <img width=200px height=200px src="https://github.com/tf2/CNest/blob/master/src/fluff/cnv_logo.png" alt="Project logo" width="500" height="150"></a>
</p>

<h3 align="center">CNest Main Docs</h3>

<div align="center">

  [![Status](https://img.shields.io/badge/status-active-success.svg)]() 
  [![GitHub Issues](https://img.shields.io/github/issues/kylelobo/The-Documentation-Compendium.svg)](https://github.com/tf2/CNest/issues)
  [![GitHub Pull Requests](https://img.shields.io/github/issues-pr/kylelobo/The-Documentation-Compendium.svg)](https://github.com/tf2/CNest/pulls)
  [![License](https://img.shields.io/badge/license-MIT-blue.svg)](/LICENSE)

</div>

---

<p align="center"> Analysis Methods and Workflows for Large Scale Copy Number Variation Analysis from Next Generation (NGS) Sequencing Data.
    <br> 
</p>

## 📝 Table of Contents
- [About](#about)
- [Usage](#usage)
- [Authors](#authors)
- [Acknowledgments](#acknowledgement)

## 🧐 About <a name = "about"></a>
CNest contains methods that have been speifically developed for large scale analysis of copy number from NGS data. It primarily uses read depth (coverage) information to generate robust copy number estimates for individual samples and is most approprate for use in very large cohorts (minimum of 1000 samples). 

It uses a 'dynamic' reference based approach to generate copy number estimates where the predicted set of ideal internal samples are used to generate a dynamic baseline across the genome for every sample, meaning that each sample has a different reference containing different baseline coverage levels for all target regions across the genome. During this process the reference datasets are selected to minimise several different noise characteristics including the presence and scale of 'genome waves' and an optimisation of the dose response. The number of samples that make up each individual reference can be set at a fixed level (recommended) or derived dynamically for each sample independantly (development).

CNest includes a reliable CNV caller that is based on a Hidden Markov Model (HMM) that calls CNV's for each sample using a techique that could be termed "joint calling" where the models are trained in batches relating to samples that have similar properties. There are also several methods for downstream CNV analysis including CNV merging, population frequency calculations, annotation and PCA. Although CNest will provide highly consistent CNV call sets across very large cohorts the primary reason for its development is the calculation of copy number estimate that are robust enough to allow genome wide association analysis for CNVs from NGS datasets (CNV-GWAS).

The primary 


## 🎈 Usage <a name="usage"></a>
Add notes about how to use the system.

## 🚀 Deployment <a name = "deployment"></a>
Add additional notes about how to deploy this on a live system.

## ⛏️ Built Using <a name = "built_using"></a>
- [Docker](https://docs.docker.com/) - build
- [Nextflow](https://www.nextflow.io/) - workflow
- [WDL](https://github.com/openwdl/wdl) - workflow

## ✍️ Authors <a name = "authors"></a>
- [@tf2](https://github.com/tf2) - Idea, initial work and method development
- [@smshuai](https://github.com/smshuai) - Method development and Nextflow implementation
- [@bshifaw](https://github.com/bshifaw) - WDL implementation

## 🎉 Acknowledgements <a name = "acknowledgement"></a>
- Hat tip to anyone whose code was used
- Inspiration
- References
