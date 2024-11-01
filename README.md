# EpiSmokEr2
EpiSmokEr2 is an R-package for predicting smoking status with whole blood DNA methlyation array data. Based on EpiSmokEr (Bollepalli et al., 2019) which trained on 450K data, EpiSmokEr2 includes an additional model for EPIC/EPICv2 data. The EPIC model was trained using LASSO regression with data from the Young Finns Study (YFS).

# Installation
To install:
```r
devtools::install_github("chloezhu/EpiSmokEr2")
```

# References
Bollepalli S, Korhonen T, Kaprio J, Anders S, Ollikainen M. EpiSmokEr: a robust classifier to determine smoking status from DNA methylation data. Epigenomics. 2019 Oct;11(13):1469-1486. doi: 10.2217/epi-2019-0206. Epub 2019 Aug 30. PMID: 31466478.
