## flu-evidence

This repository contains the code used for flu data analysis first published in

Baguelin, Marc, Stefan Flasche, Anton Camacho, Nikolaos Demiris, Elizabeth Miller, and W. John Edmunds. ‘Assessing Optimal Target Populations for Influenza Vaccination Programmes: An Evidence Synthesis and Modelling Study’. PLoS Med 10, no. 10 (8 October 2013): e1001527. doi:10.1371/journal.pmed.1001527.

The code is distributed under a GPLv3 license. Please remember to cite the above manuscript when using this code.

## Install

This code is now mainly developed as an R package and we recommend it to install it as such.

First make sure you have the devtools package install in R by executing

```{r}
install.packages("devtools")
```

from within R. Note that devtools relies on a fairly new R version (>=3.1.0).

Then use devtools to install from fluEvidenceSynthesis by using the following command:

```{r}
library(devtools)
install_github("MJomaba/flu-evidence-synthesis")
```

## Documentation and examples



