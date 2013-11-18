BAMTOOLS_DIR="/home/daggy/lib/bamtools-2.2.2"
#BAMTOOLS_DIR="/home/daggy/usr/lib/bamtools-2.6.4"
g++-4.7 -o strie main.cpp -std=c++11 -std=gnu++11 -g -O2 -fPIC -Wall ${BAMTOOLS_DIR}/lib/libbamtools.a -I${BAMTOOLS_DIR}/include -lz -I"/usr/include/" -l boost_filesystem -l boost_system -l boost_program_options-mt
