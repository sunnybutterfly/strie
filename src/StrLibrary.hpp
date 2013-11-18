#ifndef STRLIBRARY_HPP_
#define STRLIBRARY_HPP_

// STL libraries
#include <cmath>        // log(), exp()
#include <cstdlib>      // abs()
#include <fstream>      // open(), close()
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>


// External libraries
#include <boost/multi_array.hpp>


// Namespace
using namespace std;
using namespace boost;


// Class
class StrLibrary
{
// Functions
public:
    StrLibrary();
    ~StrLibrary();
    StrLibrary(const StrLibrary & rhs);

    bool IsEmptyMpersPdf();

// Variables
    map<int,int>        insertSizeFreqMap;
    multi_array<T,2>    mpersPdfArr;
    map<pair<int,int>,int> overlapFfMap;
    multi_array<int,2>  rLenFfArr;
    multi_array<int,2>  rLenPdfArr;
    map<int, int>       rLenNRead;

    string libTag;

    bool isPicked;
    int inUse;
    int nRead;
    int libInsertSizeMax;

    T pct0;
    T pct2p5;
    T pct25;
    T pct50;
    T pct75;
    T pct97p5;
    T pct100;

    T mean;
    T stDev;

private:

};



StrLibrary::StrLibrary()
{
    mpersPdfArr.resize(extents[0][0]);

    isPicked = true;
    libTag = "";
    nRead = 0;
    libInsertSizeMax = 450;

    pct0 = 0;
    pct2p5 = 0;
    pct25 = 0;
    pct50 = 0;
    pct75 = 0;
    pct97p5 = 0;
    pct100 = 0;

    mean = 0;
    stDev = 0;
}



StrLibrary::~StrLibrary()
{

}



StrLibrary::StrLibrary(const StrLibrary & rhs)
{
    this->mpersPdfArr.resize(extents[rhs.mpersPdfArr.shape()[0]] [rhs.mpersPdfArr.shape()[1]]);

    for (int i0 = 0; i0 < int(rhs.mpersPdfArr.shape()[0]); i0++)
        for (int i1 = 0; i1 < int(rhs.mpersPdfArr.shape()[1]); i1++)
            this->mpersPdfArr[i0][i1] = rhs.mpersPdfArr[i0][i1];

    this->insertSizeFreqMap = rhs.insertSizeFreqMap;
    this->libTag = rhs.libTag;
    this->isPicked = rhs.isPicked;
    this->nRead = rhs.nRead;
    this->libInsertSizeMax = rhs.libInsertSizeMax;

    this->pct2p5 = rhs.pct2p5;
    this->pct25 = rhs.pct25;
    this->pct50 = rhs.pct50;
    this->pct75 = rhs.pct75;
    this->pct97p5 = rhs.pct97p5;
    this->mean = rhs.mean;
    this->stDev = rhs.stDev;
}



bool StrLibrary::IsEmptyMpersPdf()
{
    if (mpersPdfArr.shape()[0] == 0 && mpersPdfArr.shape()[1] == 0)
        return true;

    return false;
}


#endif /* STRLIBRARY_HPP_ */
