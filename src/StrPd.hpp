#ifndef STRPD_HPP_
#define STRPD_HPP_

// STL libraries
#include <algorithm>    // swap(), max_element
#include <cmath>        // log(), exp()
#include <ctime>
#include <cstdlib>      // abs(), itoa()
#include <fstream>      // open(), close()
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

// External libraries

// Own files
#include "type.hpp"

// Namespace
using namespace std;
using namespace boost;

// Class
class StrPd
{
public:
  StrPd();
  ~StrPd();
  StrPd(const StrPd &rhs);

  int CalcPriorPdf();
  //int CalcPriorPdfMax();

  T GetPriorPdfMax();

  int GenPriorUnif(int indelStepSizeNew, int nIndelStepNew);

  std::vector<int> GetIndelSize();
  int GetIndelSize(
      int & indelStepSize,
      int & nIndelStep,
      int & indelSizeMax,
      int & nIndelElem);
  int GetPrior(
      T * priorCur,
      int & indelStepSize,
      int & nIndelStep,
      int & indelSizeMax,
      int & nIndelElem);
  T GetPwrPar();

  int SetPrior(const T * priorNew, int indelStepSizeNew, int nIndelStepNew);
  int SetPrior(
      const list<T> & priorListNew,
      int indelSizeMaxNew,
      int indelStepSizeNew);
  int SetPwrPar(const T pwrParNew);

protected:

  int CalcPriorPdf(
      list<T> priorList,
      int indelSizeMaxNew,
      int indelStepSizeNew);
  int CalcPriorPdfMax();

  int NormalizePrior();
  int PowerTransform();

  multi_array<T, 1> prior;      // Array of indel probabilities of alleles.
  multi_array<T, 2> priorPdf;   // Array of indel probabilities of genotypes.
  std::vector<int> indelShift; // Container of indel sizes.

  int indelSizeMax;   // Maximum indel size.
  int indelStepSize;  // Size of indel steps.
  int nIndelElem;     // No. of indel elements.
  int nIndelStep;     // No. of indel steps.
  T priorPdfMax;    // Maximum value of PDF of priors.

private:

  T StrToT(const string & sth);
  T pwrPar;         // Power transform parameter.
};

StrPd::StrPd()
{
  indelStepSize = 3;
  nIndelStep = 10;
  indelSizeMax = indelStepSize * nIndelStep;
  nIndelElem = 2 * nIndelStep + 1;
  pwrPar = 1.0;

  // Generate uniform prior probabilities.
  GenPriorUnif(indelStepSize, nIndelStep);
  CalcPriorPdfMax();
}

StrPd::~StrPd()
{

}

int StrPd::CalcPriorPdf()
{
  T sum = 0;

  // Change indelShift.
  indelShift.resize(nIndelElem);

  for (int i0 = 0; i0 < nIndelElem; i0++)
    indelShift[i0] = -indelSizeMax + indelStepSize * i0;

  // Normalize values.
  for (int i0 = 0; i0 < nIndelElem; i0++)
    sum += prior[i0];

  for (int i0 = 0; i0 < nIndelElem; i0++)
    prior[i0] /= sum;

  // Compute prior probabilities
  priorPdf.resize(extents[nIndelElem][nIndelElem]);

  for (int i0 = 0; i0 < nIndelElem; i0++)
  {
    int indelNo1 = indelShift[i0];

    for (int i1 = 0; i1 < i0 + 1; i1++)
    {
      int indelNo2 = indelShift[i1];

      // Homozygous reference.
      if (indelNo1 == 0)
      {
        if (indelNo2 == 0)
          priorPdf[i0][i1] = 0;
        else
          priorPdf[i0][i1] = std::log(prior[i1]);
      }
      // Heterozygous one reference.
      else if (indelNo2 == 0)
      {
        if (indelNo1 == 0)
          priorPdf[i0][i1] = 0;
        else
          priorPdf[i0][i1] = std::log(prior[i0]);
      }
      else
      {
        // Homozygous non-reference.
        if (indelNo1 == indelNo2)
        {
          priorPdf[i0][i1] = std::log(prior[i0]) - std::log(2);
        }
        // Heterozygous non-reference.
        else if (abs(indelNo1) < abs(indelNo2))
        {
          T val1 = std::log(prior[i0]);
          T val2 = std::log(prior[i1]);

          priorPdf[i0][i1] = val1 / 2 + val2;
        }
        else if (abs(indelNo1) > abs(indelNo2))
        {
          T val1 = std::log(prior[i0]);
          T val2 = std::log(prior[i1]);

          priorPdf[i0][i1] = val1 + val2 / 2;
        }
        else
        {
          T val1 = std::log(prior[i0]);
          T val2 = std::log(prior[i1]);
          priorPdf[i0][i1] = (val1 + val2) * 3 / T(4);
        }
      }
    }
  }

  NormalizePrior();
  PowerTransform();
  NormalizePrior();
  CalcPriorPdfMax();

  return 0;
}

int StrPd::CalcPriorPdfMax()
{
  T valMax = priorPdf[0][0];

  for (int i0 = 0; i0 < nIndelElem; i0++)
    for (int i1 = 0; i1 < i0 + 1; i1++)
      valMax = max(priorPdf[i0][i1], valMax);

  priorPdfMax = valMax;

  return 0;
}

int StrPd::GenPriorUnif(int indelStepSizeNew, int nIndelStepNew)
{
  if (indelStepSizeNew < 1)
  {
    cout << "The minimum of indel step size is 1." << endl;
    return 1;
  }

  if (nIndelStepNew < 1)
  {
    cout << "The minimum of no. of indel steps is 1." << endl;
    return 2;
  }

  indelStepSize = indelStepSizeNew;
  nIndelStep = nIndelStepNew;
  indelSizeMax = indelStepSize * nIndelStep;
  nIndelElem = nIndelStep * 2 + 1;

  T prob = 1 / double(nIndelElem);
  T probLog = std::log(prob);

  prior.resize(extents[nIndelElem]);
  priorPdf.resize(extents[nIndelElem][nIndelElem]);

  for (int i0 = 0; i0 < nIndelElem; i0++)
    prior[i0] = prob;

  for (int i0 = 0; i0 < nIndelElem; i0++)
    for (int i1 = 0; i1 < i0 + 1; i1++)
      priorPdf[i0][i1] = probLog;

  indelShift.resize(nIndelElem);

  for (int i0 = 0; i0 < nIndelElem; i0++)
    indelShift[i0] = (i0 - nIndelStep) * indelStepSize;

  return 0;
}

std::vector<int> StrPd::GetIndelSize()
{
  std::vector<int> indelSize;
  int nIndelShift = indelShift.size();

  if (nIndelShift == 0)
    return indelSize;

  indelSize.push_back(indelShift[0]);

  if (nIndelShift == 1)
    indelSize.push_back(indelShift[0]);
  else
    indelSize.push_back(indelShift[nIndelShift - 1]);

  return indelSize;
}

int StrPd::GetIndelSize(
    int & indelStepSize,
    int & nIndelStep,
    int & indelSizeMax,
    int & nIndelElem)
{
  indelStepSize = this->indelStepSize;
  nIndelStep = this->nIndelStep;
  indelSizeMax = this->indelSizeMax;
  nIndelElem = this->nIndelElem;

  return 0;
}

int StrPd::GetPrior(
    T * priorCur,
    int & indelStepSize,
    int & nIndelStep,
    int & indelSizeMax,
    int & nIndelElem)
{
  if (prior.shape()[0] == 0)
  {
    cout << "The size of the vector of prior probabilities is zero." << endl;
    return 1;
  }

  //priorCur = new T[nIndelElem];

  for (int i0 = 0; i0 < nIndelElem; i0++)
    priorCur[i0] = prior[i0];

  indelStepSize = this->indelStepSize;
  nIndelStep = this->nIndelStep;
  indelSizeMax = this->indelSizeMax;
  nIndelElem = this->nIndelElem;

  return 0;
}

T StrPd::GetPriorPdfMax()
{
  return priorPdfMax;
}

T StrPd::GetPwrPar()
{
  return pwrPar;
}

int StrPd::NormalizePrior()
{
  T sum = 0.0;
  T sumLog;

  // Calculate sum of existing values.
  for (int i0 = 0; i0 < nIndelElem; i0++)
    for (int i1 = 0; i1 < i0 + 1; i1++)
      sum += std::exp(priorPdf[i0][i1]);

  sumLog = std::log(sum);

  // Normalize such that sum is 1.
  for (int i0 = 0; i0 < nIndelElem; i0++)
    for (int i1 = 0; i1 < i0 + 1; i1++)
      priorPdf[i0][i1] -= sumLog;

  return 0;
}

int StrPd::PowerTransform()
{
  // Calculate sum of existing values.
  for (int i0 = 0; i0 < nIndelElem; i0++)
    for (int i1 = 0; i1 < i0 + 1; i1++)
      priorPdf[i0][i1] *= pwrPar;

  return 0;
}

int StrPd::SetPrior(const T * priorNew, int indelStepSizeNew, int nIndelStepNew)
{
  if (indelStepSizeNew < 1)
  {
    cout
        << "The step size of the vector of new prior probabilities must be greater than zero."
        << endl;
    return 1;
  }

  if (nIndelStepNew < 1)
  {
    cout
        << "The no. of steps of the vector of new prior probabilities must be greater than zero."
        << endl;
    return 2;
  }

  //for (int i0 = 0; i0 < nIndelStepNew * 2 + 1; i0++)
  //    cout << priorNew[i0] << endl;

  for (int i0 = 0; i0 < nIndelStepNew * 2 + 1; i0++)
  {
    if (priorNew[i0] < 0)
    {
      cout << "All prior probabilities must be positive." << endl;
      return 3;
    }
  }

  //T sum;

  indelStepSize = indelStepSizeNew;
  nIndelStep = nIndelStepNew;
  indelSizeMax = indelStepSize * nIndelStep;
  nIndelElem = nIndelStep * 2 + 1;

  // Copy values.
  prior.resize(extents[nIndelElem]);

  for (int i0 = 0; i0 < nIndelElem; i0++)
    prior[i0] = priorNew[i0];

  // Update Prior PDF.
  CalcPriorPdf();

  return 0;
}

int StrPd::SetPrior(
    const list<T> & priorListNew,
    int indelSizeMaxNew,
    int indelStepSizeNew)
{
  list<T>::const_iterator cit0;
  int count = 0;

  indelSizeMax = indelSizeMaxNew;
  indelStepSize = indelStepSizeNew;
  nIndelStep = indelSizeMax / indelStepSize;
  nIndelElem = 2 * nIndelStep + 1;

  prior.resize(extents[nIndelElem]);
  priorPdf.resize(extents[nIndelElem][nIndelElem]);

  for (cit0 = priorListNew.begin(); cit0 != priorListNew.end(); cit0++)
  {
    prior[count] = *cit0;
    count++;
  }

  CalcPriorPdf();

  return 0;
}

int StrPd::SetPwrPar(const T pwrParNew)
{
  if (pwrParNew < 0.1 || pwrParNew > 2)
  {
    cout << "The power transform parameter must lie in the range [0.1,2]."
        << endl;
    return 1;
  }

  pwrPar = pwrParNew;

  // Normalize, compute power transform and find largest value in PDF of prior.
  NormalizePrior();
  PowerTransform();
  NormalizePrior();
  CalcPriorPdfMax();

  return 0;
}

T StrPd::StrToT(const string & sth)
{
  istringstream instr(sth);
  T val;
  instr >> val;

  return val;
}

#endif /* STRPD_HPP_ */
