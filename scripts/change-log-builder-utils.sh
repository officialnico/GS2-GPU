#!/usr/bin/env bash
: '
 Script for building a change log for utils from merged pull requests. 
 Usage:
 From within the utils folder:
 $ sh ../scripts/change-log-builder-utils.sh > CHANGELOG.MD 
'
previous_tag=0
for current_tag in $(git tag --sort=-creatordate)
do

if [ "$previous_tag" != 0 ];then
    tag_date=$(git log -1 --pretty=format:'%ad' --date=short ${previous_tag})
    printf "## ${previous_tag} (${tag_date})\n\n"
    git log ${current_tag}..${previous_tag} --grep "pull request" --pretty=format:'%aN %s %b' |\
    grep "Merged" |\
    sed 's;^\(.*\) Merged.*pull request #\([0-9]*\)) \(.*\);* \3 [PR #\2 | \1](https://bitbucket.org/gyrokinetics/utils/pull-requests/\2);'
    printf "\n\n"
fi
previous_tag=${current_tag}
done

