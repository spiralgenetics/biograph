#!/bin/sh

BIOGRAPH=modules/biograph/biograph
TINYHUMAN=/share/datasets/tinyhuman-rand
TMPDIR=`mktemp -d`

trap "echo rm -rf ${TMPDIR}/tinyhuman.bg" EXIT

set -e

if [ -f /scratch/tinyhuman-loop.bg ]
then
    cp -sr /scratch/tinyhuman-loop.bg ${TMPDIR}/tinyhuman.bg
else
    cp -sr ${TINYHUMAN}/tinyhuman-loop.bg ${TMPDIR}/tinyhuman.bg
fi
chmod -R u+w ${TMPDIR}/tinyhuman.bg

REFMAP=${TMPDIR}/tinyhuman.bg.refmap
if [ -d /home/${USER}/datasets/tinyhuman-loop.bg ]
then
    REFMAP=/home/${USER}/datasets/tinyhuman-loop.bg.refmap
elif [ -d /scratch/tinyhuman-loop.bg ] || [ -e /scratch/tinyhuman-loop.bg.refmap ]
then
    REFMAP=/scratch/tinyhuman-loop.bg.refmap
fi

echo 'Calling variants on '${TMPDIR}/tinyhuman.bg
rm -f /tmp/compare-tinyhuman-rand.prof
time env CPUPROFILE=/tmp/compare-tinyhuman-rand.prof \
${BIOGRAPH} variants \
	    --in ${TMPDIR}/tinyhuman.bg \
	    --ref /reference/hs37d5 \
	    --ref-map ${REFMAP} \
	    --out ${TMPDIR}/tinyhuman-biograph-unfiltered.vcf \
	    --bed ${TINYHUMAN}/tinyhuman-reads.bed \
	    --assemblies-out ${TMPDIR}/assemblies-out.csv \
	    --verbose-trace-work=true
echo
echo 'Filtering variants present in tinyhuman-calls.bed'
vcftools --vcf ${TMPDIR}/tinyhuman-biograph-unfiltered.vcf \
	 --bed ${TINYHUMAN}/tinyhuman-calls.bed \
	 --recode --recode-INFO-all \
	 --out ${TMPDIR}/tinyhuman-biograph
# Join alleles, and fix the genotyping.  "bcftools norm" joins two
# "0/1" genotypes into a "0/2" genotype when we'd really them to be
# joined as a "1/2" genotype instead.
bcftools norm -f /reference/hs37d5/source.fasta -m +any ${TMPDIR}/tinyhuman-biograph.recode.vcf |
    sed 's@\t0/2:@\t1/2:@' > ${TMPDIR}/tinyhuman-biograph-norm.vcf
vcf-sort < ${TMPDIR}/tinyhuman-biograph-norm.vcf > ${TMPDIR}/tinyhuman-biograph.vcf
bgzip ${TMPDIR}/tinyhuman-biograph.vcf
tabix -p vcf ${TMPDIR}/tinyhuman-biograph.vcf.gz
echo
echo 'Evaluating VCF quality versus ground truth'
ln -s ${TINYHUMAN}/tinyhuman-ground-truth.vcf.gz ${TMPDIR}/tinyhuman-ground-truth.vcf.gz
rtg vcfeval -b ${TINYHUMAN}/tinyhuman-ground-truth.vcf.gz \
    -c ${TMPDIR}/tinyhuman-biograph.vcf.gz \
    -t /share/datasets/HG001/hg19.sdf --all-records \
    -o ${TMPDIR}/tinyhuman-rand-comparison --squash-ploidy
echo
echo 'Random selection of false negatives:'
zcat ${TMPDIR}/tinyhuman-rand-comparison/fn.vcf.gz |shuf -n 10
echo 'Random selection of false positives:'
zcat ${TMPDIR}/tinyhuman-rand-comparison/fp.vcf.gz |shuf -n 10
exit 0
