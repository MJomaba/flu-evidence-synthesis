#!/bin/bash

# This script depends on inotify-hookable

inotify-hookable -w src -f CMakeLists.txt -w tests -c "ninja" -c "bin/test_*"
