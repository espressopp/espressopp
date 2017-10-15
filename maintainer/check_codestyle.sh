#!/bin/bash

echo "check_codestyle.sh CODESTYLE=$1 TRAVIS_BRANCH=$2"

CODESTYLE=$1
TRAVIS_BRANCH=$2

INCORRECT_FILES=""

function clangformat() {
    for fchanged in $1; do
        if [[ "${fchanged}" =~ .*\.[ch]+pp ]]; then
            clang-format -style=file ${fchanged} | diff -q - ${fchanged} &> /dev/null || INCORRECT_FILES="${INCORRECT_FILES} ${fchanged}"
        fi
    done
}

if [[ $CODESTYLE = ON ]]; then
    FILES="$(git diff --name-only --diff-filter=AM $TRAVIS_BRANCH...HEAD)"
    echo "FILES: ${FILES}"
    clangformat "${FILES}"
    echo "INCORRECT_FILES: ${INCORRECT_FILES}"
    [[ $INCORRECT_FILES ]] && exit 1 || exit 0
fi
exit 0
