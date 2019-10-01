#!/usr/bin/env bash
: '
 Script for building a change log from merged pull requests. 
 Usage:
 $ sh change-log-builder.sh > CHANGELOG.MD 
'

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

previous_tag=0
for current_tag in $(git tag --sort=-creatordate)
do
    if [ "$previous_tag" != 0 ];then
        format_commit_range $current_tag $previous_tag
    fi
    previous_tag=${current_tag}
done

