#SAMTOOLS="/media/DAGGE/Universitet/UCC/Projects/SAMtools/samtools-0.1.10"
#SAMDIR="~/lib/samtools-0.1.12a"


# Old boole.
#g++ -o strie main.cpp -g -O2 -fPIC -Wall ${BAMTOOLS_DIR}/lib/libbamtools.a -I${BAMTOOLS_DIR}/include -lz
#g++ -o strie main.cpp -g -O2 -fPIC -Wall ${BAMTOOLS_DIR}/lib/libbamtools.a -I${BAMTOOLS_DIR}/include -lz -I${HOME}/usr/include -I${HOME}/usr/local/include -L/home/bcri/dlyberg/usr/local/lib -lPocoDataSQLite -Wl,-rpath,/home/bcri/dlyberg/usr/local/lib

# New boole (test).
g++ -o strie_weldon main.cpp -g -O2 -fPIC -Wall ${BAMTOOLS_DIR}/lib/libbamtools.a -I${BAMTOOLS_DIR}/include -lz -I${HOME}/usr/include -I${HOME}/usr/local/include -L/home/bcri/dlyberg/usr/local/lib -lPocoFoundation -lPocoDataSQLite -Wl,-rpath,/home/bcri/dlyberg/usr/local/lib

# On office computer.
#g++ -o strie main.cpp -g -O2 -fPIC -Wall ${BAMTOOLS_DIR}/lib/libbamtools.a -I${BAMTOOLS_DIR}/include -lz -I"/usr/include/" -lPocoData -lPocoFoundation -lPocoSQLite  #-I"/usr/include/qt4/QtSql/" -I"/usr/include/qt4/" -lQtSql


#cp strie ${HOME}/usr/bin/strie
ln -s -f ${HOME}"/Project/STRIE for paper/"strie_weldon ${HOME}/usr/bin/strie_weldon

