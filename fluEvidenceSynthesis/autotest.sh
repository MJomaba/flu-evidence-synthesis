#!/bin/bash

# Always run at least once
R -e 'Rcpp::compileAttributes(".",verbose=TRUE);library(devtools);devtools::check()'

# This script depends on inotify-hookable
inotify-hookable -t 1000 -i '*.o' -i '\..*' -w data -w src -w tests -f DESCRIPTION -w R -c "R -e 'Rcpp::compileAttributes(\".\",verbose=TRUE);library(devtools);devtools::check()' && sleep 10"
