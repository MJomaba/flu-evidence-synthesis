#!/bin/bash

# This script depends on inotify-hookable

inotify-hookable --watch-files src/*.cc --watch-files src/*.hh --watch-files CMakeLists.txt --watch-files tests/*.cc --watch-files tests/*.hpp --on-modify-command "ninja" --on-modify-command "bin/test_*"
