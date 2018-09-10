# fluEvidenceSynthesis

This repository contains the fluEvidenceSynthesis package which was described in more detail in:

van Leeuwen, Edwin, Petra Klepac, Dominic Thorrington, Richard Pebody, and Marc Baguelin. “FluEvidenceSynthesis: An R Package for Evidence Synthesis Based Analysis of Epidemiological Outbreaks.” PLOS Computational Biology 13, no. 11 (2017): e1005838. https://doi.org/10.1371/journal.pcbi.1005838.

The code is distributed under a GPLv3 license. Please remember to cite the above manuscript when using this code.

## Install

### Devtools

First make sure you have the devtools package install in R by executing the following in R:

```{r}
install.packages("devtools")
```

Note that devtools relies on a fairly new R version (>=3.1.0). Devtools also relies on a number of packages and one can run into problems installing it. Please refer to their [documentation](https://github.com/hadley/devtools) when encountering problems.

### FluEvidenceSynthesis

Then use devtools to install fluEvidenceSynthesis by using the following command:

```{r}
library(devtools)
install_github("MJomaba/flu-evidence-synthesis", dependencies = TRUE)
```

The above will install the package without the vignettes. To install them as well you first need to install [odin](https://github.com/mrc-ide/odin#installation) following [these instructions](https://github.com/mrc-ide/odin#installation), because one of the vignettes depends on that package. After installing odin you can install the package and its vignettes as follows.

```{r}
library(devtools)
install_github("MJomaba/flu-evidence-synthesis", dependencies = TRUE, build_vignettes = TRUE)
```


## Documentation and examples

All functions provided by the package are documented in R, but we recommend reading the provided vignettes for a good introduction/overview of the provided capabilities. Vignettes on the epidemiological model and parameter inference can be found [here](http://blackedder.github.io/flu-evidence-synthesis/modelling.html), [here](http://blackedder.github.io/flu-evidence-synthesis/adapting-the-transmission-model.html) and [here](http://blackedder.github.io/flu-evidence-synthesis/inference.html). There is also a vignette available on analysing the effects of different [vaccination scenarios](http://blackedder.github.io/flu-evidence-synthesis/vaccination.html).

If you also installed the vignettes locally (see above) you can view them using:

```{r}
# To list the included vignettes:
vignette(package="fluEvidenceSynthesis")

# To open the inference vignette:
vignette("inference",package="fluEvidenceSynthesis")
```

There is also a presentation available with some more details on the package [here](http://blackedder.github.io/flu-evidence-synthesis/RCoursePackage.html).
