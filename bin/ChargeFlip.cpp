#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"
#include "CombineHarvester/CombineTools/interface/Systematics.h"
#include "CombineHarvester/CombineTools/interface/BinByBin.h"

#include <boost/filesystem.hpp> // boost::filesystem::
#include <boost/program_options.hpp> // boost::program_options::
#include <boost/range/adaptor/map.hpp> // boost::adaptors::map_keys
#include <boost/range/algorithm/copy.hpp> // boost::copy()
#include <boost/algorithm/string/join.hpp> // boost::algorithm::join()
#include <boost/algorithm/string/case_conv.hpp> // boost::algorithm::to_lower_copy()

#include <algorithm> // std::transform()

/*! @file ChargeFlip.cpp
    @brief Generate datacards for charge flip fit

    @author  Andres Tiko <andres.tiko@cern.ch>    
             Karl Ehatäht <karl.ehataht@cern.ch>
*/

/** TODO: support for muons */
/** TODO: option to omit certain combinations of charge and category for signal/background */
/** TODO: fix boost program_option notifiers */

std::vector<std::string> SHAPESYST_NAMES = {
  "CMS_ttHl_electronER",
  "CMS_ttHl_electronESBarrel",
  "CMS_ttHl_electronESEndcap",
};

struct DatacardParams
{
  DatacardParams()
    : year(-1)
  {}
  std::string input_file;
  std::string output_dir;
  std::vector<std::string> charges;
  std::vector<std::string> signal_processes;
  std::vector<std::string> background_processes;
  std::string prefix;
  std::vector<std::string> shape_systs_signal;
  std::vector<std::string> shape_systs_background;
  int year;
};

std::ostream &
operator<<(std::ostream & os,
           const DatacardParams & params)
{
  os
    << "Input:            " << params.input_file                                           << "\n"
       "Output:           " << params.output_dir                                           << "\n"
       "Charges:          " << boost::algorithm::join(params.charges, ", ")                << "\n"
       "Signal:           " << boost::algorithm::join(params.signal_processes, ", ")       << "\n"
       "Background:       " << boost::algorithm::join(params.background_processes, ", ")   << "\n"
       "Signal shape:     " << boost::algorithm::join(params.shape_systs_signal, ", ")     << "\n"
       "Background shape: " << boost::algorithm::join(params.shape_systs_background, ", ") << "\n"
       "Year:             " << params.year                                                 << "\n"
       "Prefix:           " << params.prefix                                               << '\n'
  ;
  return os;
}

std::map<int, double> LUMINOSITY_UNC = { // in %
  { 2016, 2.5 },
  { 2017, 2.3 },
  { 2018, 2.5 },
};

std::vector<std::string> ERAS = {
  "13TeV",
};

std::vector<std::pair<int, std::string>> BINNING = {
  { 0,  "BB_LL" },
  { 1,  "BB_ML" },
  { 2,  "BB_MM" },
  { 3,  "BB_HL" },
  { 4,  "BB_HM" },
  { 5,  "BB_HH" },
  { 6,  "EE_LL" },
  { 7,  "EE_ML" },
  { 8,  "EE_MM" },
  { 9,  "EE_HL" },
  { 10, "EE_HM" },
  { 11, "EE_HH" },
  { 12, "BE_LL" },
  { 13, "BE_ML" },
  { 14, "EB_ML" },
  { 15, "BE_MM" },
  { 16, "BE_HL" },
  { 17, "EB_HL" },
  { 18, "BE_HM" },
  { 19, "EB_HM" },
  { 20, "BE_HH" },
};

std::vector<std::string> CHARGES = {
  "OS",
  "SS",
};

std::map<std::string, double> SIGNALS_UNC = { // in %
  { "DY", 50.},
};

std::map<std::string, double> BACKGROUNDS_UNC = { // in %
  { "DY_fake",    50. },
  { "WJets",     400. }, // Exclude by default as not enough statistics for a reasonable shape?
  { "Singletop",  50. },
  { "Diboson",    50. },
  { "TTbar",      50. },
};

/*!
    indir - directory with histograms for data and pseudodata created by create_pseudodata_datacard.py
    outdir - where to put output, relative name
    type: data or pseudodata
*/
bool
create_datacard(const DatacardParams & params)
{
  // Create an empty CombineHarvester instance that will hold all of the
  // datacard configuration and histograms etc.
  ch::CombineHarvester cb;

  std::cout << ">> Creating processes and observations...\n";
  for(const std::string & era: ERAS)
  {
    for(const std::string & charge: params.charges)
    {
      const std::string key = Form("%s_%s", charge.data(), era.data());
      cb.AddObservations({ "*" }, { "htt" }, { era }, { charge }, BINNING);
      cb.AddProcesses   ({ "*" }, { "htt" }, { era }, { charge }, params.background_processes, BINNING, false);
      cb.AddProcesses   ({ "*" }, { "htt" }, { era }, { charge }, params.signal_processes,     BINNING, true);
    }
  }

  assert(LUMINOSITY_UNC.count(params.year));
  const double lumi_unc = 1. + LUMINOSITY_UNC.at(params.year) / 100.;
  for(const std::string & charge: params.charges)
  {
    //syst on luminosity
    cb.cp().channel({ charge }).signals()
      .AddSyst(cb, "lumi", "lnN", ch::syst::SystMap<>::init(lumi_unc))
    ;
    cb.cp().channel({ charge }).backgrounds()
      .AddSyst(cb, "lumi", "lnN", ch::syst::SystMap<>::init(lumi_unc))
    ;

    //systs on normalization
    for(const std::string & signal_process: params.signal_processes)
    {
      assert(SIGNALS_UNC.count(signal_process));
      cb.cp().channel({ charge }).signals()
        .AddSyst(
          cb, Form("%s_norm", signal_process.data()), "lnN",
          ch::syst::SystMap<>::init(1. + SIGNALS_UNC.at(signal_process) / 100.)
        )
      ;
    }
    for(const std::string & background_process: params.background_processes)
    {
      assert(BACKGROUNDS_UNC.count(background_process));
      const std::string syst_label = background_process == "DY_fake" ?
        "fake"                                                       :
        boost::algorithm::to_lower_copy(background_process)
      ;
      cb.cp().channel({ charge }).process({ background_process })
        .AddSyst(
          cb, Form("%s_norm", syst_label.data()), "lnN",
          ch::syst::SystMap<>::init(1. + BACKGROUNDS_UNC.at(background_process) / 100.)
        )
      ;
    }
  }

  for(const std::string & charge: params.charges)
  {
    const std::set<std::string> bins = cb.cp().channel({ charge }).bin_set();
    for(const std::string & bin: bins)
    {

      for(const std::string & shape_syst_sig: params.shape_systs_signal)
      {
        // skip if(bin == "EE_LL" && charge == "SS") ?
        cb.cp().channel({ charge }).bin({ bin }).signals().AddSyst(
          cb, shape_syst_sig, "shape", ch::syst::SystMap<>::init(1.00)
        );
      }
      for(const std::string & shape_syst: params.shape_systs_background)
      {
        cb.cp().channel({ charge }).bin({ bin }).backgrounds().AddSyst(
          cb, shape_syst, "shape", ch::syst::SystMap<>::init(1.00)
        );
      }
    }
  }

  std::cout << ">> Extracting histograms from input root files...\n";
  for(const std::string & era : ERAS)
  {
    for(const std::string & charge : params.charges)
    {
      std::cout
        << "Input file: " << params.input_file << "\n, "
           "era: "        << era               << ", "
           "charge: "     << charge            << '\n'
      ;
      const std::string shape_central = Form(
        "%s_%s_$BIN/rebinned/$PROCESS_rebinned",             params.prefix.data(), charge.data()
      );
      const std::string shape_sys     = Form(
        "%s_%s_$BIN/rebinned/$PROCESS_$SYSTEMATIC_rebinned", params.prefix.data(), charge.data()
      );

      cb.cp().channel({ charge }).backgrounds().ExtractShapes(params.input_file, shape_central, shape_sys);
      cb.cp().channel({ charge }).signals().    ExtractShapes(params.input_file, shape_central, shape_sys);
    }
  }

  ch::BinByBinFactory bbb = ch::BinByBinFactory()
    .SetAddThreshold(0.05)
    .SetMergeThreshold(0.5)
    .SetFixNorm(true)
  ;

  bbb.AddBinByBin(cb.cp().backgrounds(), cb);
  bbb.AddBinByBin(cb.cp().signals(), cb);

  // This function modifies every entry to have a standardised bin name of
  // the form: {analysis}_{channel}_{bin_id}_{era}
  // which is commonly used in the htt analyses
  ch::SetStandardBinNames(cb);

  // First we generate a set of bin names:
  const std::set<std::string> bins = cb.bin_set();
  // This method will produce a set of unique bin names by considering all
  // Observation, Process and Systematic entries in the CombineHarvester
  // instance.


  for(const std::string & charge: params.charges)
  {
    // ! Where to write output datacards, update accordingly
    const std::string folder = Form("%s/cards/%scards", params.output_dir.data(), charge.data());
    const std::string folder_common = Form("%s/common", folder.data());
    const std::string output_filename = Form("%s/htt_%s.input.root", folder_common.data(), charge.data());

    if(! boost::filesystem::create_directories(folder_common))
    {
      std::cerr << "Unable to create directory: " << folder_common << '\n';
      return false;
    }

    TFile output(output_filename.data(), "recreate");
    const std::set<std::string> bins = cb.cp().channel({ charge }).bin_set();
    for(const std::string & bin : bins)
    {
      const std::string datacard = Form("%s/%s.txt", folder.data(), bin.data());
      std::cout << ">> Writing datacard for bin " << bin << " to " << datacard << '\r' << std::flush;
      cb.cp().channel({ charge }).bin({ bin }).WriteDatacard(datacard, output);
    }
    output.Close();
    std::cout << ">> Created output file: " << output_filename << '\n';
  }
  return true;
}

template<typename T>
void
notifier(const std::vector<T> & choices,
         const T & selection,
         const std::string & option_name)
{
  if(std::find(choices.cbegin(), choices.cend(), selection) == choices.cend())
  {
    throw boost::program_options::validation_error(
      boost::program_options::validation_error::invalid_option_value, option_name
    );
  }
}

template<typename T>
void
notifierv(const std::vector<T> & choices,
          const std::vector<T> & selections,
          const std::string & option_name)
{
  std::vector<T> invalid_selections;
  for(const T & selection: selections)
  {
    if(std::find(choices.cbegin(), choices.cend(), selection) == choices.cend())
    {
      invalid_selections.push_back(selection);
    }
  }
  if(! invalid_selections.empty())
  {
    throw boost::program_options::validation_error(
      boost::program_options::validation_error::invalid_option_value, option_name
    );
  }
}

template<typename T,
         typename U>
std::vector<T>
get_map_keys(const std::map<T, U> & map)
{
  std::vector<T> keys;
  boost::copy(map | boost::adaptors::map_keys, std::back_inserter(keys));
  return keys;
}

int
main(int argc,
     const char ** argv)
{
  const std::vector<std::string> signal_keys = get_map_keys(SIGNALS_UNC);
  const std::vector<std::string> background_keys = get_map_keys(BACKGROUNDS_UNC);
  const std::vector<int> years = get_map_keys(LUMINOSITY_UNC);

  const std::string charges_joined = boost::algorithm::join(CHARGES, ", ");
  const std::string signals_joined = boost::algorithm::join(signal_keys, ", ");
  const std::string backgrounds_joined = boost::algorithm::join(background_keys, ", ");
  const std::string shapes_joined = boost::algorithm::join(SHAPESYST_NAMES, ", ");

  std::vector<std::string> years_str;
  std::transform(
    years.cbegin(), years.cend(), std::back_inserter(years_str),
    [](int year) -> std::string { return std::to_string(year); }
  );

  const std::vector<std::string> shape_syst_signal_default = SHAPESYST_NAMES;
  const std::vector<std::string> shape_syst_background_default = {
    "CMS_ttHl_electronESBarrel",
    "CMS_ttHl_electronESEndcap",
  };
  const std::string shape_syst_signal_default_joined = boost::algorithm::join(shape_syst_signal_default, ", ");
  const std::string shape_syst_background_default_joined = boost::algorithm::join(shape_syst_background_default, ", ");

  DatacardParams params;
  bool force = false;
  boost::program_options::options_description desc("Options");
  desc.add_options()
    (
      "help,h",
      "Help"
    )
    (
      "input,i",
      boost::program_options::value<std::string>(&params.input_file),
      "Input pseudo/data datacard"
    )
    (
      "output,o",
      boost::program_options::value<std::string>(&params.output_dir),
      "Output directory"
    )
    (
      "charge,c",
      boost::program_options::value<std::vector<std::string>>(&params.charges)
        ->multitoken()
        ->default_value(CHARGES, charges_joined)
        ->notifier([](const std::vector<std::string> & charge) -> void {
            return notifierv(CHARGES, charge, "charge");
          })
      ,
      Form("Charge sums (choices: %s)", charges_joined.data())
    )
    (
      "sig,s",
      boost::program_options::value<std::vector<std::string>>(&params.signal_processes)
        ->multitoken()
        ->default_value(signal_keys, signals_joined)
        ->notifier([&signal_keys](const std::vector<std::string> & sig_procs) -> void {
            return notifierv(signal_keys, sig_procs, "sig");
          })
      ,
      Form("Signal processes (choices: %s)", signals_joined.data())
    )
    (
      "bkg,b",
      boost::program_options::value<std::vector<std::string>>(&params.background_processes)
        ->multitoken()
        ->default_value(background_keys, backgrounds_joined)
        ->notifier([&background_keys](const std::vector<std::string> & bkg_procs) -> void {
            return notifierv(background_keys, bkg_procs, "bkg");
          })
      ,
      Form("Background processes (choices: %s)", backgrounds_joined.data())
    )
    (
      "year,y",
      boost::program_options::value<int>(&params.year)
        ->notifier([&years](int year) -> void {
            return notifier(years, year, "year");
          })
      ,
      Form("Year (choices: %s)", boost::algorithm::join(years_str, ", ").data())
    )
    (
      "shape-signal,S",
      boost::program_options::value<std::vector<std::string>>(&params.shape_systs_signal)
        ->multitoken()
        ->default_value(shape_syst_signal_default, shape_syst_signal_default_joined)
        ->notifier([&shape_syst_signal_default](const std::vector<std::string> & syst) -> void {
            return notifierv(shape_syst_signal_default, syst, "shape-signal");
          })
      ,
      Form("Shape systematics for signal (choices: %s)", shape_syst_signal_default_joined.data())
    )
    (
      "shape-background,B",
      boost::program_options::value<std::vector<std::string>>(&params.shape_systs_background)
        ->multitoken()
        ->default_value(shape_syst_background_default, shape_syst_background_default_joined)
        ->notifier([&shape_syst_background_default](const std::vector<std::string> & syst) -> void {
            return notifierv(shape_syst_background_default, syst, "shape-background");
          })
      ,
      Form("Shape systematics for background (choices: %s)", shape_syst_background_default_joined.data())
    )
    (
      "prefix,p",
      boost::program_options::value<std::string>(&params.prefix)
        ->default_value("ttH_charge_flip")
      ,
      Form("Prefix in histogram names")
    )
    (
      "force,f",
      boost::program_options::bool_switch(&force),
      "Force the creation of output directory"
    )
  ;

  try
  {
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
    boost::program_options::notify(vm);
    if(vm.count("help"))
    {
      std::cout << desc << '\n';
      return EXIT_SUCCESS;
    }
    for(const std::string & vm_key: { "input", "output", "year" })
    {
      if(! vm.count(vm_key))
      {
        std::cerr << "Option '" << vm_key << "' was not specified\n";
        return EXIT_FAILURE;
      }
    }
  }

  catch(const boost::program_options::error & err)
  {
    std::cerr << err.what() << '\n';
  }


  if(! boost::filesystem::exists(params.input_file))
  {
    std::cerr << "No such input path: " << params.input_file << '\n';
    return EXIT_FAILURE;
  }
  else if(! boost::filesystem::is_regular_file(params.input_file))
  {
    std::cerr << "No such input file: " << params.input_file << '\n';
    return EXIT_FAILURE;
  }

  if(! boost::filesystem::exists(params.output_dir))
  {
    if(! force)
    {
      std::cerr << "No such output directory: " << params.output_dir << " (use -f to create it)\n";
      return EXIT_FAILURE;
    }
    else
    {
      if(boost::filesystem::create_directories(params.output_dir))
      {
        std::cerr << "Created output directory: " << params.output_dir << '\n';
      }
      else
      {
        std::cerr << "Unable to create output directory: " << params.output_dir << '\n';
        return EXIT_FAILURE;
      }
    }
  }
  else if(! boost::filesystem::is_directory(params.output_dir))
  {
    std::cerr << "Output path is not a directory: " << params.output_dir << '\n';
    return EXIT_FAILURE;
  }

  std::cout << params << '\n';

//  if(! create_datacard(params))
//  {
//    return EXIT_FAILURE;
//  }
  return EXIT_SUCCESS;
}
