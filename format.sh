#! /bin/bash

git ls-files -- src*.cpp src*.hpp testsuite*.cpp testsuite*.hpp | xargs clang-format -i -style=file
git ls-files -- src*.py testsuite*.py | xargs autopep8 --in-place --aggressive --aggressive 
