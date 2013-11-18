#ifndef STR_MENU_EST_C
#define STR_MENU_EST_C

// STL libraries
//#include <cstdlib>      // abs(), rand()
#include <fstream>      // open(), close()
#include <iostream>
#include <string>
#include <vector>

// External libraries
#include <sys/time.h>   // gettimeofday
// Own files
#include "StrIndelSize.hpp"
#include "MenuAuxiliary.hpp"
#include "MenuEst.hpp"

// Namespace
using namespace std;

void SetParametersEstSrp(po::variables_map &vm)
{
  StrIndelSize str = StrIndelSize();
  const string out_path = vm["out"].as<string>();
  const filesystem::path out_fp = filesystem::path(out_path);
  const string bam_path = vm["bam"].as<string>();
  const filesystem::path bam_fp = filesystem::path(bam_path);
  const string loc_path = vm["loc"].as<string>();
  const filesystem::path loc_fp = filesystem::path(loc_path);
  const string prior_path = vm["prior"].as<string>();
  const filesystem::path prior_fp = filesystem::path(prior_path);
  const string isize_path = vm["isize"].as<string>();
  const filesystem::path isize_fp = filesystem::path(isize_path);

  const string insert_size_str = vm["ival"].as<string>();
  std::vector<string> insert_sizes;
  int insert_size_min, insert_size_max;
  const int indel_size_max = vm["indelsizemax"].as<int>();
  const int indel_step_size = vm["indelstepsize"].as<int>();
  const string ref_size_str = vm["refsize"].as<string>();
  std::vector<string> ref_sizes;
  int ref_size_min, ref_size_max;
  const bool check_neighbor = vm["neighbor"].as<string>() == "" ? true : false;

  const int anchor_len = vm["alen"].as<int>();
  const int mapq_min = vm["mapq"].as<int>();

  // Check that file paths exist.
  if (!(boost::filesystem::exists(bam_fp) && boost::filesystem::is_regular_file(bam_fp)))
  {
    throw "BAM file path not correct";
  }

  // Check that file paths exist.
  if (!(boost::filesystem::exists(loc_fp) && boost::filesystem::is_regular_file(loc_fp)))
  {
    throw "Path not correct of file containing STR locus coordinates";
  }

  // Check that file paths exist.
  if (!(boost::filesystem::exists(prior_fp) && boost::filesystem::is_regular_file(prior_fp)))
  {
    throw "Path not correct of file containing prior probabilities";
  }

  // Check that file paths exist.
  if (!(boost::filesystem::exists(isize_fp) && boost::filesystem::is_regular_file(isize_fp)))
  {
    throw "Path not correct of file containing insert size frequencies";
  }

  // Check that BAM file can be opened (exit otherwise) .
  str.OpenBam(bam_fp);

  // Read BAM header.
  cout << "Reading BAM header." << endl;
  str.ReadBamHeader();

  // Check that files can be opened (exit otherwise) .
  cout << "Reading locus coordinate input.";
  str.ReadLocus(loc_fp);

  cout << "\nReading insert size input.";
  str.ReadInsertSize(isize_fp);

  cout << "\nReading prior probability input.";
  str.ReadPrior(prior_fp, 0, indel_size_max, indel_step_size);

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

  str.SetInsertSize(insert_size_min, insert_size_max);

  // Reference size interval
  boost::algorithm::split(ref_sizes, ref_size_str, is_any_of("-"));
  try
  {
    if (ref_sizes.size() != 2)
      throw "";
    ref_size_min = boost::lexical_cast<int>(ref_sizes[0]);
    ref_size_max = boost::lexical_cast<int>(ref_sizes[1]);
  }
  catch(std::exception &e)
  {
    throw "Reference interval size argument is incorrect.";
  }

  str.SetRefSize(ref_size_min, ref_size_max);

  // Set minimum MAPQ score.
  str.SetMapQMin(mapq_min);
  str.SetAnchorLen(anchor_len);
  str.SetCheckNeighbor(check_neighbor);

  // Compute MPERS distributions.
  str.CalcMpersPdf();

  // Open locus output file.
  cout << "\nOpening output file.";
  str.OpenLocOutput(out_fp);

  // Run files.
  str.LaunchStrEstimation();

  return;
}

void MenuEstSrp(int argc, char *argv[])
{
  po::options_description desc("Options");
  po::variables_map vm;

  // add options to variables
  desc.add_options()("out", po::value<string>()->required(), "text output file");
  desc.add_options()("bam", po::value<string>()->required(), "BAM input file");
  desc.add_options()("loc", po::value<string>()->required(), "STR locus coordinates");
  desc.add_options()("prior", po::value<string>()->required(), "prior probabilities");
  desc.add_options()("isize", po::value<string>()->required(), "insert size probabilities (from \"strie freq\")");

  desc.add_options()(
    "indelsizemax",
    po::value<int>()->default_value(45),
    "indel size maximum, INT [45]");
  desc.add_options()(
    "indelstepsize",
    po::value<int>()->default_value(15),
    "indel step size, INT [15]");
  desc.add_options()(
    "refsize",
    po::value<string>()->default_value("10-250"),
    "reference size interval, min-max [10-250]");
  desc.add_options()(
    "skiplib",
    po::value<string>()->default_value(""),
    "skip library, STR");
  desc.add_options()(
    "ival",
    po::value<string>()->default_value("36-450"),
    "insert size interval of all libraries, min-max [36-450]");
  desc.add_options()(
    "libival",
    po::value<string>()->default_value(""),
    "insert size interval of a single library (overrides -ival), STR min-max");
  desc.add_options()(
    "alen",
    po::value<int>()->default_value(3),
    "anchor length, INT [3]");
  desc.add_options()(
    "mapq",
    po::value<int>()->default_value(30),
    "MAPQ score minimum, INT [30]");
  desc.add_options()(
    "neighbor",
    po::value<string>()->default_value("1"),
    "if read overlaps two STR loci, use it");
  //desc.add_options()(
  //  "locdist",
  //  po::value<int>()->default_value(50),
  //  "minimum distance between STR loci");

  if (argc <= 1)
  {
    ShowMenuEstSrp(desc);
    return;
  }

  try
  {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // start finding random reads
    SetParametersEstSrp(vm);
  }
  catch (std::exception &e)
  {
    cout << e.what() << "\n";
    return;
  }

}

void ShowMenuEstSrp(const po::options_description &desc)
{
  cout << endl;
  cout << "Usage: strie est [options]" << endl << endl;
  cout << desc;
  cout << endl;
  cout << "Example:" << endl;
  cout << "  strie est";
  cout << " --out out.txt";
  cout << " --bam file.bam";
  cout << " --loc locus.txt";
  cout << " --prior prior.txt";
  cout << " --isize isize.txt";
  cout << " --indelsizemax 45";
  cout << " --indelstepsize 15";
  cout << " --refsize 10-250";
  cout << " --skiplib lib1";
  cout << " --ival 36-450";
  cout << " --mapq 30";
  cout << endl;
  cout << endl;


}

#endif // STR_MENU_EST_C
