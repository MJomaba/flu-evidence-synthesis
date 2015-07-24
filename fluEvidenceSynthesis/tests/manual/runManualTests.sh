#!/bin/bash

# Always run at least once
R -e 'Rcpp::compileAttributes(".",verbose=TRUE)';
# Should actually move needed header in inst/include/fluEvidenceSynthesis.h insted of using echo
R CMD check ../fluEvidenceSynthesis
R CMD INSTALL ../fluEvidenceSynthesis

R -e 'library(testthat);library(fluEvidenceSynthesis);test_dir("tests/manual/")'
