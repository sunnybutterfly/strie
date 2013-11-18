// STL libraries
//#include <cstdlib>      // abs(), rand()
#include <fstream>      // open(), close()
#include <iostream>
#include <string>
#include <vector>

// External libraries
#include <sys/time.h>   // gettimeofday
// Own files

//#include "boost_libraries.hpp"

#include "StrIndelSize.hpp"
#include "MenuAuxiliary.hpp"
#include "MenuFreq.hpp"

void SetParametersFreqSrp(po::variables_map &vm)
{
  const string out_path = vm["out"].as<string>();
  filesystem::path out_fp = filesystem::path(out_path);
  const string bam_path = vm["bam"].as<string>();
  filesystem::path bam_fp = filesystem::path(bam_path);

  StrIndelSize str = StrIndelSize();
  string bamPath;
  string outPath;

  const int mapq_mpers_min = vm["mapq"].as<int>();
  const string insert_size_str = vm["ival"].as<string>();
  std::vector<string> insert_sizes;
  int insert_size_min;
  int insert_size_max;
  const Region region = { 0, "", 0, 0 };
  const int n_read_max = numeric_limits<int>::max();

  // Check that file paths exist.
  if (!(boost::filesystem::exists(bam_fp) && boost::filesystem::is_regular_file(bam_fp)))
  {
    throw "BAM file path not correct";
  }

  // Check that BAM file can be opened (exit otherwise) .
  str.OpenBam(bam_fp);

  // Read BAM header.
  str.ReadBamHeader();

  // Set remaining options (exit if wrong).
  boost::algorithm::split(insert_sizes, insert_size_str, is_any_of("-"));
  try
  {
    if (insert_sizes.size() != 2)
      throw "";
    insert_size_min = boost::lexical_cast<int>(insert_sizes[0]);
    insert_size_max = boost::lexical_cast<int>(insert_sizes[1]);
  }
  catch(std::exception &e)
  {
    throw "Insert size argument is incorrect.";
  }

  str.SetMapQMpersMin(mapq_mpers_min);
  str.SetInsertSize(insert_size_min, insert_size_max);
  str.SetReadCountMax(n_read_max);
  // Insert size interval
  boost::algorithm::split(insert_sizes, insert_size_str, is_any_of("-"));
  try
  {
    if (insert_sizes.size() != 2)
      throw "";
    insert_size_min = boost::lexical_cast<int>(insert_sizes[0]);
    insert_size_max = boost::lexical_cast<int>(insert_sizes[1]);
  }
  catch(std::exception &e)
  {
    throw "Insert size argument is incorrect.";
  }
  // Read insert size frequencies.
  str.ReadInsertSize(region.refIdStr, region.refBeg, region.refEnd);

  // Compute MPERS distributions.
  str.CalcMpersPdf();

  // Print some information.
  str.PrintLibInfo(0);

  // Write insert size frequencies to file.
  str.WriteInsertSize(out_fp);

  return;
}

void MenuFreqSrp(int argc, char *argv[])
{
  po::options_description desc("Options");
  po::variables_map vm;

  // add options to variables
  desc.add_options()("out", po::value<string>()->required(),
      "text output file");
  desc.add_options()("bam", po::value<string>()->required(), "BAM input file");

  desc.add_options()("ival", po::value<string>()->default_value("36-450"),
      "insert size interval of all libraries, min-max [36-450]");
  desc.add_options()("mapq", po::value<int>()->default_value(30),
      "MAPQ score minimum, INT [30]");

  if (argc <= 1)
  {
    ShowMenuFreqSrp(desc);
    return;
  }

  try
  {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // start finding random reads
    SetParametersFreqSrp(vm);
  } catch (std::exception &e)
  {
    cout << "hoopla" << endl;
    cout << e.what() << "\n";
    return;
  }

}

void ShowMenuFreqSrp(const po::options_description &desc)
{
  cout << endl;
  cout << "Usage: strie freq [options]" << endl << endl;
  cout << desc;
  cout << endl;
}

