###############
### Compilation
###############

# Download STRIE
git clone git://github.com/sunnybutterfly/strie
cd strie

# If you have the BAMTOOLS C++ library installed, enter the paths here if 
# necessary (e.g. BAMTOOLS_SRC_DIR="/usr/include/", 
# BAMTOOLS_LIB_DIR="/usr/include/").
BAMTOOLS_SRC_DIR=""
BAMTOOLS_LIB_DIR=""

# Otherwise, install as follows.
BAMTOOLS_SRC_DIR="../bamtools"
BAMTOOLS_LIB_DIR="../bamtools"

git clone https://github.com/pezmaster31/bamtools
cd bamtools; mkdir build; cd build; cmake ..; make
cd ..; cd ..

# If you have the BOOST C++ libraries installed, enter the path here if 
# necessary (e.g. "/usr/include/").
BOOST_SRC_DIR=""
BOOST_LIB_DIR=""

# Otherwise, install as follows.
BOOST_SRC_DIR="../boost-trunk"
BOOST_LIB_DIR="../boost-trunk"

svn co http://svn.boost.org/svn/boost/trunk boost-trunk
cd boost-trunk; ./bootstrap.sh --prefix=./; ./b2 install
cd ..

# Compile STRIE
cd src
g++-4.7 -o strie main.cpp -std=c++11 -std=gnu++11 -g -O2 -fPIC -Wall ${BAMTOOLS_LIB_DIR}/lib/libbamtools.a -I${BAMTOOLS_SRC_DIR}/include -I{BOOST_SRC_DIR} -L{BOOST_LIB_DIR} -l boost_filesystem -l boost_system -l boost_program_options-mt -lz


################
### Introduction
################

STRIE (Short Tandem Repeat Indel Estimator) is a software package for estimating
genotypes of STR loci. It estimates genotypes using a statistical algorithm of
based on estimated insert sizes of paired-end (PE) reads mapped to a reference 
genome. PE reads that have one read mapped to the left side of an STR locus and
the other read mapped to the right side have a different insert size 
distribution than other reads if there exist indels in the STR locus, nearby 
which they are mapped. Such PE reads are called SRPs (Spanning Read Pairs)

An insert size distribution of a PE read library is unknown. It can be 
estimated using the distance computed after mapping to a reference genome. This
is called an MPERS (Mapped Paired-End Read Separations; my terminology) 
distribution.

STRIE was originally developed for repeat locus lengths 1-6 nt for use with the
Illumina read libraries in the 1,000 Genomes Project. Two trios were sequenced, 
a CEU trio (NA12878, NA12877, NA12882) and a YRI trio (NA19238, NA19239, 
NA19240). It turned out in simulations that the sequencing coverage was far too 
low for the intended use. In other simulations with 15-mers, it appeared that 
NA12878 has got the right sequencing coverage to be used. The repeat loci 
should have a reference length of 30-250 nt and contain indels no longer than 
(+/-) 30 nt.

STRIE contains two algorithms. First, an MPERS distribution of each read 
library must be created. This is STRIE FREQ. The other algorithm, STRIE EST, 
estimates genotypes of repeat loci. The statistical algorithm in STRIE EST 
requires a set of prior probabilities located in "test/indel_priors_15mers.txt".
STR locus coordinates of 15-mers are located in 
"test/str_locus_coordinates_chr22_15-mers.txt". An example of how to run STRIE
can be found in the file "test/test.sh" and it requires the BAM file of chr. 22
of NA12878.

The BAM files of NA12878 from the 1,000 Genomes Project are 
NA12878.chromX.ILLUMINA.bwa.CEU.high_coverage.20100311.bam, where X is a 
chromosome, located at 
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/pilot2_high_cov_GRCh37_bams/data/NA12878/alignment/.
To estimate the genotypes, run first STRIE FREQ and then STRIE EST.


################
### Availability
################

STRIE is relased under [GPLv3][1] and can be downloaded from [2] at Github.


################
### Seeking help
################

Please write questions to [3].


##############
### Citing BWA
##############

STRIE was developed in part by PhD student Weldon Whitener in Richard Durbin's 
group at the Wellcome Trust Sanger Institute, Cambridge, UK. I later rewrote it
in C++ as a PhD student at UCC, Cork, Ireland under Avril Coghlan, who works at
Wellcome Trust Sanger Institute. STRIE remains as of yet unpublished.

Dag Lyberg
November 20/11/2013.


[1]: http://en.wikipedia.org/wiki/GNU_General_Public_License
[2]: https://github.com/sunnybutterfly/strie
[3]: https://github.com/sunnybutterfly/strie/issues
