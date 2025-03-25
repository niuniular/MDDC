| Category       | Badge                                                                                                                                                                                                                                                                                                                                                                                                                                 |
|--------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Usage**    | [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://github.com/niuniular/MDDC/blob/main/LICENSE.md) [![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/MDDC?color=blue)](https://cran.r-project.org/package=MDDC)                                                                                                                                                     |
| **Release**  | [![R CMD CHECK](https://github.com/niuniular/MDDC/actions/workflows/r.yml/badge.svg)](https://github.com/niuniular/MDDC/actions/workflows/r.yml) [![CRAN](https://www.r-pkg.org/badges/version/MDDC?color=green)](https://cran.r-project.org/package=MDDC) [![CRAN checks](https://badges.cranchecks.info/summary/MDDC.svg)](https://cran.r-project.org/web/checks/check_results_MDDC.html)                                                                                                                                 |
| **Development** | [![codecov](https://codecov.io/gh/niuniular/MDDC/graph/badge.svg?token=YLQNKH2AI6)](https://app.codecov.io/gh/niuniular/MDDC) [![Codacy Badge](https://app.codacy.com/project/badge/Grade/c71ad2e3c52546e7a99ee8fea57d6142)](https://app.codacy.com/gh/niuniular/MDDC/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade) [![CodeFactor](https://www.codefactor.io/repository/github/niuniular/mddc/badge)](https://www.codefactor.io/repository/github/niuniular/mddc) |





# Introduction

* In this package, we present the Modified Detecting Deviating Cells (MDDC) algorithm for adverse event identification.
* For a certain time period, the spontaneous reports can be extracted from the safety database and depicted as an $I \times J$ contingency table, where:
  * $I$ denotes the total number of AEs
  * $J$ denotes the total number of drugs or vaccines
  * With cell counts $n_{ij}$ the total number of reported cases corresponding to the $j$-th drug/vaccine and $i$-th AE
* We are interested in which (AE, drug or vaccine) pairs are signals. The signals refer to potential adverse events that may be caused by a drug/vaccine.
* In the contingency table setting, the signals refer to the cells with $n_{ij}$ abnormally higher than the expected values.
* [**Rousseeuw and Bossche (2018)**](https://wis.kuleuven.be/stat/robust/papers/publications-2018/rousseeuwvandenbossche-ddc-technometrics-2018.pdf) proposed the Detecting Deviating Cells (DDC) algorithm for outlier identification in a multivariate dataset.
* The original DDC algorithm assumes multivariate normality of the data and selects cutoff values based on this assumption. We modify the DDC algorithm to better suit the discrete nature of adverse event data in pharmacovigilance that clearly do not follow a multivariate normal distribution.
* Our Modified Detecting Deviating Cells (MDDC) algorithm has the following characteristics:
  1. It is easy to compute.
  2. It considers AE relationships.
  3. It depends on data-driven cutoffs.
  4. It is independent of the use of ontologies.
* The MDDC algorithm has five steps, with the first two steps identifying univariate outliers via cutoffs, and the next three steps evaluating the signals via the use of AE correlations. The algorithm can be found at **[MDDC algorithm](https://mddc.readthedocs.io/en/latest/user_guide/mddc_algorithm.html)**.

## Authors

* **Anran Liu**
  Email: <anranliu@buffalo.edu>

* **Raktim Mukhopadhyay**
  Email: <raktimmu@buffalo.edu>

* **Marianthi Markatou**
  Email: <markatou@buffalo.edu>

## Maintainers

* **Anran Liu** Email: <anranliu@buffalo.edu>

* **Raktim Mukhopadhyay** Email: <raktimmu@buffalo.edu>

## Documentation

The documentation is hosted at - https://niuniular.github.io/MDDC/index.html

## Citation

If you use this package in your research or work, please cite it as follows:

```
@misc{liu2024mddcrpythonpackage,
      title={MDDC: An R and Python Package for Adverse Event Identification in Pharmacovigilance Data}, 
      author={Anran Liu and Raktim Mukhopadhyay and Marianthi Markatou},
      year={2024},
      eprint={2410.01168},
      archivePrefix={arXiv},
      primaryClass={stat.CO},
      url={https://arxiv.org/abs/2410.01168}, 
}
```

## Funding Information

The work has been supported by Food and Drug Administration, and Kaleida Health Foundation.

## References

Liu, A., Mukhopadhyay, R., and Markatou, M. (2024). MDDC: An R and Python package for adverse event identification in pharmacovigilance data. arXiv preprint. arXiv:2410.01168

Liu, A., Markatou, M., Dang, O., and Ball, R. (2024). Pattern discovery in pharmacovigilance through the Modified Detecting Deviating Cells (MDDC) algorithm. Technical Report, Department of Biostatistics, University at Buffalo.

Rousseeuw, P. J., and Bossche, W. V. D. (2018). Detecting deviating data cells. Technometrics, 60(2), 135-145.
