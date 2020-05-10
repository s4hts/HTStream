#!/usr/bin/env bash
set -euo pipefail
export BUILDDIR=$1
echo $0

## note you must run this from the regression dir.
## first arg is build dir second is either test or regen

fastqr1=./test_r1.fastq.gz
fastqr2=./test_r2.fastq.gz

regenerate() {
    for prog in `find $BUILDDIR -maxdepth 2 -name "hts*" -type f -executable -not -name "*_test"`; do
        out=${prog##*/}
        echo running $prog output to $out
        $prog -1 $fastqr1 -2 $fastqr2 -t $out -F
    done
}

testrun() {
    for prog in `find $BUILDDIR  -maxdepth 2 -name "hts*" -type f -executable -not -name "*_test"`; do
        out=${prog##*/}.test
        echo running $prog output to $out
        $prog -1 $fastqr1 -2 $fastqr2 -t $out -F

        if [ ${prog##*/} == 'hts_SuperDeduper' ]
        then
            echo sorting superDeduper because its output is non-deterministic
            mv $out.tab6.gz $out.tmp && zcat $out.tmp | sort | gzip -c > $out.tab6.gz
            rm $out.tmp

            orig=${out%%.*}
            zcat $orig.tab6.gz | sort | gzip -c > $orig.sorted.tab6.gz
            echo zdiff $out.tab6.gz $orig.sorted.tab6.gz
            zdiff $out.tab6.gz $orig.sorted.tab6.gz
            rm $out.tab6.gz
            rm $orig.sorted.tab6.gz
        else
            echo zdiff $out.tab6.gz ${out%%.*}.tab6.gz
            zdiff $out.tab6.gz ${out%%.*}.tab6.gz
            rm $out.tab6.gz
        fi
    done
}

if [ $2 == 'test' ]
then
    testrun
elif [ $2 == 'regen' ]
then
    regenerate
fi
