./../src/strie freq \
  --out "isize.txt" \
  --bam "/media/Dagge/UCC other/1000GP/Data/NA12878/NA12878.chrom22.ILLUMINA.bwa.CEU.high_coverage.20100311.bam" \
  --ival 36-450 \
  --mapq 30

./../src/strie est \
  --out "est.txt" \
  --bam "/media/Dagge/UCC other/1000GP/Data/NA12878/NA12878.chrom22.ILLUMINA.bwa.CEU.high_coverage.20100311.bam" \
  --loc "str_locus_coordinates_chr22_15-mers.txt" \
  --prior "indel_priors_15mers.txt" \
  --isize "isize.txt" \
  --indelsizemax 150 \
  --indelstepsize 15 \
  --refsize 10-250 \
  --ival 36-450 \
  --mapq 30

