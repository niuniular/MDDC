# Introduction

- In this package, we present the Modified Detecting Deviating Cells (MDDC) algorithm for adverse event identification.
- For a certain time period, the spontaneous reports can be extracted from the safety database and depicted as an $I \times J$ contingency table, where:
  - $I$ denotes the total number of AEs
  - $J$ denotes the total number of drugs
  - With cell counts $n_{ij}$ the total number of reported cases corresponding to the $j$-th drug and $i$-th AE
- We are interested in which (AE, drug) pairs are signals. The signals refer to potential adverse events that may be caused by a drug.
- In the contingency table setting, the signals refer to the cells with $n_{ij}$ abnormally higher than the expected values.
- [**Rousseeuw and Bossche (2018)**](https://wis.kuleuven.be/stat/robust/papers/publications-2018/rousseeuwvandenbossche-ddc-technometrics-2018.pdf) proposed the Detecting Deviating Cells (DDC) algorithm for outlier identification in a multivariate dataset.
- The original DDC algorithm assumes multivariate normality of the data and selects cutoff values based on this assumption. The foundation of the DDC algorithm lies in detecting deviating data cells within a multivariate dataset. Inspired by this work, we modify the DDC algorithm to better suit the discrete nature of adverse event data in pharmacovigilance.
- Our Modified Detecting Deviating Cells (MDDC) algorithm has the following characteristics:
  1. It is easy to compute.
  2. It considers AE relationships.
  3. It depends on data-driven cutoffs.
- The MDDC algorithm involves five steps, with the first two steps identifying univariate outliers via cutoffs, and the next three steps evaluating the signals via the use of AE correlations. The algorithm can be found at **[MDDC algorithm](https://mddc.readthedocs.io/en/latest/user_guide/mddc_algorithm.html)**.

## Authors

- **Anran Liu** 
  Email: [anranliu@buffalo.edu](mailto:anranliu@buffalo.edu)  

- **Raktim Mukhopadhyay** 
  Email: [raktimmu@buffalo.edu](mailto:raktimmu@buffalo.edu)  

- **Marianthi Markatou** 
  Email: [markatou@buffalo.edu](mailto:markatou@buffalo.edu)  

## Maintainer

- **Anran Liu** 
  Email: [anranliu@buffalo.edu](mailto:anranliu@buffalo.edu) 

- **Raktim Mukhopadhyay**  
  Email: [raktimmu@buffalo.edu](mailto:raktimmu@buffalo.edu)

## Documentation

The documentation is hosted on `Read the Docs` at - https://mddc.readthedocs.io/en/latest/

## Funding Information
The work has been supported by Kaleida Health Foundation, Food and Drug Administration, and Department of Biostatistics, University at Buffalo.

## References

Rousseeuw, P. J., & Bossche, W. V. D. (2018). Detecting deviating data cells. Technometrics, 60(2), 135-145.