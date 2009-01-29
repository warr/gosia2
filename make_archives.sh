#!/bin/bash

list=$(git tag -l |grep ^200)

for t in $list ; do
    echo "Creating archive for gosia2 $t"
    git archive --format=tar --prefix=gosia2-$t/ $t | gzip -9 > /tmp/gosia2_$t.tgz
done
