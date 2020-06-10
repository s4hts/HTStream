#!/usr/bin/env bash
set -euo pipefail
BUILDDIR=$1
echo $0

if [[ "$OSTYPE" == "darwin"* ]]; then
    findarg="-perm +111"
else
    findarg=-executable
fi

## note you must run this from the regression dir.
## first arg is build dir second is either test or regen

fastqr1=./test_r1.fastq.gz
fastqr2=./test_r2.fastq.gz
fastq_rawr1=./test_r1.fastq
fastq_rawr2=./test_r2.fastq

gunzip -kf $fastqr1
gunzip -kf $fastqr2

regenerate() {
    for prog in `find $BUILDDIR -maxdepth 2 -name "hts*" -type f $findarg -not -name "*_test"`; do
        out=${prog##*/}
        echo running $prog output to $out
        $prog -1 $fastqr1 -2 $fastqr2 -t $out -F -L $out.json
        $prog -1 $fastq_rawr1 -2 $fastq_rawr2 -t $out -F -u

        if [ ${prog##*/} == 'hts_SuperDeduper' ]
        then
           mv $out.tab6.gz $out.tmp && gzip -dc $out.tmp | sort | gzip > $out.tab6.gz
           mv $out.tab6 $out.tmp && cat $out.tmp | sort > $out.tab6
           rm $out.tmp
        fi
    done

    ## run piped commands to make sure log output is good
    superd=`find $BUILDDIR -maxdepth 2 -name "hts_SuperDeduper" -type f $findarg -not -name "*_test"`
    stats=`find $BUILDDIR -maxdepth 2 -name "hts_Stats" -type f $findarg -not -name "*_test"`

    $superd -1 $fastqr1 -2 $fastqr2 -L chain.json | $stats -A chain.json > /dev/null
}

testrun() {
    for prog in `find $BUILDDIR  -maxdepth 2 -name "hts*" -type f $findarg -not -name "*_test"`; do
        mkdir -p test
        out=test/${prog##*/}
        echo running $prog output to $out
        $prog -1 $fastqr1 -2 $fastqr2 -t $out -F -L $out.json
        $prog -1 $fastq_rawr1 -2 $fastq_rawr2 -t $out -F -u

        orig=${prog##*/}
        if [ ${prog##*/} == 'hts_SuperDeduper' ]
        then
            echo sorting superDeduper because its output is non-deterministic
            cp $out.tab6.gz $out.tmp && gzip -dc  $out.tmp | sort | gzip > $out.tab6.gz
            echo "cp $out.tab6 $out.tmp && cat $out.tmp | sort > $out.tab6"
            cp $out.tab6 $out.tmp && cat $out.tmp | sort > $out.tab6
            rm $out.tmp

            echo diff $out.tab6 $orig.tab6
            diff $out.tab6 $orig.tab6
            echo zdiff $out.tab6.gz $orig.tab6.gz
            diff <(gzip -dc $out.tab6.gz) <(gzip -dc $orig.tab6.gz)
            rm $out.tab6 $out.tab6.gz
        else
            echo zdiff $out.tab6.gz $orig.tab6.gz
            diff <(gzip -dc $out.tab6.gz) <(gzip -dc $orig.tab6.gz)
            echo diff $out.tab6 $orig.tab6
            diff $out.tab6 $orig.tab6
            rm $out.tab6.gz $out.tab6
        fi
        echo diff $out.json $orig.json
        diff <(cat $out.json | jq "del(.[0].Program_details.version) | del(.[0].Program_details.options[\"stats-file\"]) | del(.[0].Program_details.options[\"tab-output\"])") \
            <(cat $orig.json | jq "del(.[0].Program_details.version) | del(.[0].Program_details.options[\"stats-file\"]) | del(.[0].Program_details.options[\"tab-output\"])")

        ## run piped commands to make sure log output is good
        superd=`find $BUILDDIR -maxdepth 2 -name "hts_SuperDeduper" -type f $findarg -not -name "*_test"`
        stats=`find $BUILDDIR -maxdepth 2 -name "hts_Stats" -type f $findarg -not -name "*_test"`

        $superd -1 $fastqr1 -2 $fastqr2 -L test/chain.json | $stats -A test/chain.json > /dev/null

        echo diff test/chain.json chain.json
        diff <(cat test/chain.json | jq "del(.[0].Program_details.version) | del(.[0].Program_details.options[\"stats-file\"]) | del(.[0].Program_details.options[\"append-stats-file\"]) | del(.[1].Program_details.version) | del(.[1].Program_details.options[\"stats-file\"]) | del(.[1].Program_details.options[\"append-stats-file\"])") \
            <(cat chain.json | jq "del(.[0].Program_details.version) | del(.[0].Program_details.options[\"stats-file\"]) | del(.[0].Program_details.options[\"append-stats-file\"]) | del(.[1].Program_details.version) | del(.[1].Program_details.options[\"stats-file\"]) | del(.[1].Program_details.options[\"append-stats-file\"])")

        rm -rf test
    done
}

if [ $2 == 'test' ]
then
    testrun
elif [ $2 == 'regen' ]
then
    regenerate
fi

rm $fastq_rawr1 $fastq_rawr2
rm stats.log
