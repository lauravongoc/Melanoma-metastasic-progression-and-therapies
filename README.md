# Characterization of mutational processes affecting metastatic transformation in melanoma and the impact on therapy effectiveness
**__University College London (UCL) -- MSc. Genetics of Human Disease -- 2017-2018__**


## Description
Computational project aiming to investigate the impact of distinct environmental carcinogens and intrinsic mutational processes on disease progression in melanoma. Particular attention is given on the early metastatic transformation of cancer cells. Furthermore, the efficacy of therapies are examined in terms of patient survival, as well as in the context of mutational processes due to risk factors.

The data analyzed was obtained from the following databases:
- [The Cancer Genome Atlas (TCGA)](https://portal.gdc.cancer.gov)
- [International Cancer Genome Consortium (ICGC)](https://dcc.icgc.org)


## Programming language and packages
The analyses in this project were performed in [R](https://www.r-project.org) version 3.4.3 

A number of R packages were used for data access and data analysis:
###### Data access
- [Bioconductor TCGAbiolinks](https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html) version 2.9.2
- [FirebrowseR](https://github.com/mariodeng/FirebrowseR) version 1.1.35

###### Analysis
- Mutational signatures: [deconstructSigs](https://github.com/raerose01/deconstructSigs) version 1.8.0
- Survival: [survival](https://cran.r-project.org/web/packages/survival/index.html) version 2.41.3

###### Plotting and visualization
- [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html) version 3.0.0
- [ggpubr](https://cran.r-project.org/web/packages/ggpubr/index.html) version 0.1.7.999
