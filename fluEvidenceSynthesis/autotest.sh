#!/bin/bash

# Always run at least once
R -e 'Rcpp::compileAttributes(".",verbose=TRUE)';
R CMD check ../fluEvidenceSynthesis

# This script depends on inotify-hookable
inotify-hookable -q -w src -w tests -f DESCRIPTION -c "R -e 'Rcpp::compileAttributes(\".\",verbose=TRUE)'; R CMD check ../fluEvidenceSynthesis"
#-c "bin/test_*"
