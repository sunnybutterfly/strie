#include <cstring>
#include <cctype>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

//#include "boost/lexical_cast.hpp"
#include "boost_libraries.hpp"

#include <sys/stat.h>
#include <sys/types.h>

#include "MenuAuxiliary.cpp"
#include "MenuEst.cpp"
#include "MenuFreq.cpp"
#include "StrIndelSize.hpp"

using namespace std;
using namespace boost;

int OptionEst(const std::vector<string> &argVec, int nArg);
int OptionFreq(const std::vector<string> &argVec, int nArg);
void ShowMenu();

const static string VERSION = "1.0.0";

int main(int argc, char *argv[])
{

  po::options_description desc("Allowed options");
  //int numThreads, tid;
  string command_choice;

  desc.add_options()
    ("help", "produce help message");

  if (argc < 2)
  {
    ShowMenu();
    return 0;
  }

  command_choice.assign(argv[1]);

  if (command_choice == "est")
    MenuEstSrp(argc-1, argv+1);
  else if (command_choice == "freq")
    MenuFreqSrp(argc-1, argv+1);
  else
    ShowMenu();

  return 0;
}

void ShowMenu()
{
  cout << "\nProgram: STRIE (Short Tandem Repeat Indel Estimator)";
  cout << "\nVersion: " << VERSION;

  cout << endl << endl;
  cout << "Usage: flankid <command> [options]";

  cout << endl << endl;
  cout << "Command:";
  cout << "\n\test" << "\tGenotype estimation";
  cout << "\n\tfreq" << "\tFrequency sampling of MPERS";

  cout << endl << endl;
}
