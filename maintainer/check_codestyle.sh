#!/bin/bash
set -uo pipefail

CODESTYLE=$1
TRAVIS_PULL_REQUEST=$2
TRAVIS_COMMIT=$3
TRAVIS_COMMIT_RANGE=$4
TRAVIS_REPO_SLUG=$5

echo "check_codestyle.sh CODESTYLE=${CODESTYLE} TRAVIS_PULL_REQUEST=${TRAVIS_PULL_REQUEST} TRAVIS_COMMIT=${TRAVIS_COMMIT} TRAVIS_COMMIT_RANGE=${TRAVIS_COMMIT_RANGE} TRAVIS_REPO_SLUG=${TRAVIS_REPO_SLUG}"

INCORRECT_FILES=""

function clangformat() {
    for fchanged in $1; do
        if [[ "${fchanged}" =~ .*\.[ch]+pp ]]; then
            clang-format -style=file ${fchanged} | diff -q - ${fchanged} &> /dev/null || INCORRECT_FILES="${INCORRECT_FILES} ${fchanged}"
        fi
    done
}

#function pep8format() {
#}

if [ "X$CODESTYLE" = "XON" ]; then
    if [ "X$TRAVIS_PULL_REQUEST" != "Xfalse" ]; then
        # Gets list of files from pull request.
        FILES="`curl -s "https://api.github.com/repos/${TRAVIS_REPO_SLUG}/pulls/${TRAVIS_PULL_REQUEST}/files" | grep '"filename":' | grep '"filename":' | sed -e 's/.*"filename": "\(.*\)",/\1/g' | xargs`"
    elif [ "X$TRAVIS_COMMIT_RANGE" != "X" ]; then
        FILES="`git diff --name-only ${TRAVIS_COMMIT_RANGE} | xargs`"
    else
        # Gets list from single commit.
        FILES="`git show --no-notes --oneline --name-only ${TRAVIS_COMMIT} | tail -n +2 | xargs`"
    fi
    echo "FILES: ${FILES}"
    clangformat "${FILES}"
    echo "INCORRECT_FILES: ${INCORRECT_FILES}"
    [ "$INCORRECT_FILES" ] && exit 1 || exit 0
fi
exit 0
