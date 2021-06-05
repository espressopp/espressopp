#! /bin/bash

git ls-files -- src*.cpp src*.hpp testsuite*.cpp testsuite*.hpp | xargs clang-format -i -style=file
