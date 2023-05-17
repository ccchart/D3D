#!/bin/bash

curdir=`pwd`
scriptdirname=`readlink \-f \$0`
scriptdir=`dirname $scriptdirname`
srcroot=$scriptdir/../..

echo Script      : $scriptdirname
echo In directory: $curdir
echo Executing   : $scriptdir/waqbininputdump $1 $2 >$3
                   $scriptdir/waqbininputdump $1 $2 >$3
