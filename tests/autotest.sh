#!/bin/bash

# Always run at least once
cmake . -G Ninja && ninja && bin/test_*

# This script depends on inotify-hookable
inotify-hookable -w src -f CMakeLists.txt -w tests -c "ninja && bin/test_*"
#-c "bin/test_*"
