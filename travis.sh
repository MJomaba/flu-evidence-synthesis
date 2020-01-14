#!/bin/bash

mkdir docs
cd docs
wget https://gitlab.com/BlackEdder/rcourse/raw/master/RIntroduction.html
wget https://gitlab.com/BlackEdder/rcourse/raw/master/RCoursePackage.html
cp ../vignettes/*.html ./

