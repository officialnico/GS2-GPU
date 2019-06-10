#!/bin/sh

FILE=$6/$1".xml"

/bin/cat <<EOM >$FILE
<?xml version="1.0" encoding="UTF-8"?>
<testsuite name="$1" tests="$2" errors="$3" failures="$4" skip="$5">
    <testcase classname="$1"
              name="$1" time="0">
EOM

if [ $3 != 0 ]
then
/bin/cat <<EOM >>$FILE
        <error type="test failure" message="$1 failed">
        </error>
EOM

fi

/bin/cat <<EOM >>$FILE
    </testcase>
</testsuite>
EOM


