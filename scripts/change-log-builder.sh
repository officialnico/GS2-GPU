#!/usr/bin/env bash
#
# Script for building a change log from merged pull requests.

set -e

readonly PROGRAM=$(basename $0)

readonly ERR_WRONG_NUMBER_ARGS=2
readonly ERR_BAD_ARGUMENT=3

function usage() {
    echo "usage: $PROGRAM [START END]"
    echo "       $PROGRAM --help|-h"
    echo
    echo "    Script for building a change log from merged pull requests"
    echo
    echo "    If START and END not given, uses all tags instead"
}

# This does the heavy lifting, transforming the git log into nice markdown
function format_commit_range {
    old=$1
    new=$2
    tag_date=$(git log -1 --pretty=format:'%ad' --date=short ${previous_tag})
    printf "## ${new} (${tag_date})\n\n"
    git log ${old}..${new} --grep "pull request" --pretty=format:'%aN %s %b' |\
        grep "Merged" |\
        sed 's;^\(.*\) Merged.*pull request #\([0-9]*\)) \(.*\);* \3 [PR #\2 | \1](https://bitbucket.org/gyrokinetics/gs2/pull-requests/\2);'
    printf "\n\n"
}

# Check that the argument is something git will understand
function check_is_treeish {
    if ! git rev-parse $1 2>/dev/null >/dev/null; then
        echo "$1 is not a commit, tag or branch"
        exit $ERR_BAD_ARGUMENT
    fi
}

# Check the inputs
if [[ $1 == "--help" || $1 == "-h" ]]; then
    usage
    exit 0
fi

if [[ $# != 0 && $# != 2 ]]; then
    usage
    exit $ERR_WRONG_NUMBER_ARGS
fi

# We either got two arguments, in which case give the changelog
# between those two, or we got no arguments, in which case, give the
# changelog for each tag
if [[ $# == 2 ]]; then
    check_is_treeish $1
    check_is_treeish $2
    format_commit_range $1 $2
else
    previous_tag=0
    for current_tag in $(git tag --sort=-creatordate)
    do
        if [ "$previous_tag" != 0 ];then
            format_commit_range $current_tag $previous_tag
        fi
        previous_tag=${current_tag}
    done
fi
