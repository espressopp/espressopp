#!/bin/bash

echo "check_codestyle.sh CODESTYLE=$1 TRAVIS_BRANCH=$2 TRAVIS_PULL_REQUEST=$3 TRAVIS_COMMIT=$4 TRAVIS_COMMIT_RANGE=$5"

CODESTYLE=$1
TRAVIS_BRANCH=$2
TRAVIS_PULL_REQUEST=$3
TRAVIS_COMMIT=$4
TRAVIS_COMMIT_RANGE=$5

INCORRECT_FILES=""

function clangformat() {
    for fchanged in $1; do
        if [[ "${fchanged}" =~ .*\.[ch]+pp ]]; then
            clang-format -style=file ${fchanged} | diff -q - ${fchanged} &> /dev/null || INCORRECT_FILES="${INCORRECT_FILES} ${fchanged}"
        fi
    done
}

if [[ $CODESTYLE = ON ]]; then
    if [[ $TRAVIS_PULL_REQUEST != false ]]; then
        FILES=$(git diff --name-only --diff-filter=AM $TRAVIS_BRANCH...HEAD | xargs)
    elif [[ $TRAVIS_COMMIT_RANGE ]]; then
        FILES=$(git diff --name-only --diff-filter=AM $TRAVIS_COMMIT_RANGE | xargs)
    elif [[ $TRAVIS_COMMIT ]]; then
        FILES=$(git show --name-only --no-notes --oneline --diff-filter=AM $TRAVIS_COMMIT | tail -n +2 | xargs)
    else
        FILES=$(git diff --diff-filter=AM --name-only $TRAVIS_BRANCH...HEAD | xargs)
    fi
    echo "FILES: ${FILES}"
    clangformat "${FILES}"
    echo "INCORRECT_FILES: ${INCORRECT_FILES}"
    [[ $INCORRECT_FILES ]] && exit 1 || exit 0
fi
exit 0
