#!/usr/bin/env python

from tthAnalysis.ChargeFlipEstimation.utils import read_category_ratios, readMisIDRatios, BIN_NAMES_COMPOSITE, \
                                                   BIN_NAMES_COMPOSITE_NICE, read_exclude_bins, bin_names_single, \
                                                   get_bin_nr, SmartFormatter
from tthAnalysis.ChargeFlipEstimation.matrix_solver import calculate_solution, print_ratios_latex, calculate
from tthAnalysis.ChargeFlipEstimation.plot_pulls import readMisIDRatiosGen, readCategoryRatiosGen, make_pull_plot_21, \
                                                        makeCatRatiosFrom6, compare_misIdRatios

import ROOT
import math
import argparse
import os

ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2()

"""@file docstring
Script for plotting pulls comparing result from 21 and 6 categories
Selects which categories of 21 to drop because of correlations and solves the equations for different cases

@author Andres Tiko <andres.tiko@cern.ch>
"""

# Number of sigmas difference to consider fit results not compatible
NSIGMAS = 1.2815515  # Corresponds to p-value of 0.1
# File to store list of excluded categories
EXCLUDED_FILE = "/home/ssawant/VHbbNtuples_10_x/CMSSW_10_2_13/src/tthAnalysis/ChargeFlipEstimation/bin/../data/excluded_categories.txt"


# Writes to the specified file list of categories to exclude
def select_categories(chi2s, catRatios):
  f = open(EXCLUDED_FILE, "w")
  nans = []
  print("Checking which categories to drop:")
  for (k, v) in chi2s.items():
    if v > NSIGMAS or math.isnan(catRatios[k][0]):
      f.write("%s\n" % k)
    print("Category %s, with a difference of %.2f sigmas, %s." % (k, v, "dropped" if (v > NSIGMAS or math.isnan(catRatios[k][0])) else "not dropped"))
    if math.isnan(catRatios[k][0]):
      print(" (NaN in fit results)")
      nans.append(k)
    else:
      print("")
  f.close()
  return nans

RATE_BINS = {
  0 : 'BL',
  1 : 'BM',
  2 : 'BH',
  3 : 'EL',
  4 : 'EM',
  5 : 'EH',
}

if __name__ == "__main__":
  parser = argparse.ArgumentParser(
    formatter_class = lambda prog: SmartFormatter(prog, max_help_position = 35)
  )
  parser.add_argument('-f', '--fits-data',
    type = str, dest = 'fits_data', metavar = 'directory', required = True,
    help = 'R|Input directory',
  )
  parser.add_argument('-F', '--fits-pseudodata',
    type = str, dest = 'fits_pseudodata', metavar = 'directory', required = True,
    help = 'R|Input directory',
  )
  parser.add_argument('-i', '--input-hadd',
    type = str, dest = 'input_hadd', metavar = 'file', required = True,
    help = 'R|Input hadd stage2 file',
  )
  parser.add_argument('-o', '--output',
    type = str, dest = 'output', metavar = 'directory', required = True,
    help = 'R|Output directory',
  )
  parser.add_argument('-e', '--exclude',
    type = int, dest = 'exclude', metavar = 'bin', required = False, nargs = '+', default = [],
    help = 'R|Exclude bins',
  )
  args = parser.parse_args()
  #TODO remove SS BB_LL

  input_hadd_stage2 = args.input_hadd
  output_dir = os.path.abspath(args.output)
  exclude_bins = args.exclude

  # Check 1] Check MC closure to solve 21 linear equations using dummy 6 rates::
  #   6 dummy rates ..>
  #   Add them to get 21 rations -->
  #   Solve 21 equation for 6 rates ->
  #   6 rates (caculated)
  print("Checking MC closure to solve 21 linear equations using dummy 6 rates")
  rate_dummy = 2.e-5
  eRate_dummy = 1.e-6
  misIDRatiosNum = [ (rate_dummy, eRate_dummy) for _ in RATE_BINS  ]
  misIDRatios = { bin_names_single[i] : (rate_dummy, eRate_dummy, eRate_dummy) for i in RATE_BINS }
  print("misIDRatios: (6 dummy rates)")
  for bin_idx, bin in RATE_BINS.items():
    print("  {} {}:  {}".format(bin_idx, bin, misIDRatios[bin]))

  # Calculate 21 rates from misIdRatios by adding the corresponding of the 6 rates
  catRatiosNum_genSum, catRatios_genSum = makeCatRatiosFrom6(misIDRatios, exclude_bins)
  print("catRatios_genSum: (6 dummy rates ..> Add them to get 21 rations) ")
  for key, value in catRatios_genSum.items():
    print("  {}:  {}".format(key, value))

  print("21Cat -> 6Rates:")
  (rates, uncs) = calculate(catRatiosNum_genSum, exclude_bins)
  print("CalculatedRates:")
  for i in range(0, len(rates)):
    print("\t {}:  {} +- {}".format(i, rates[i], uncs[i]))
  compare_misIdRatios(misIDRatiosNum, rates, uncs, name = "gen_dummy_closure", mydir = output_dir, outFileName = "")

  sys.exit(0)

  # Check 2.1] Check MC clousre using the gen rates::
  #   6 gen rates ..>
  #   Add them to get 21 rations -->
  #   Solve 21 equation for 6 rates -->
  #   6 gen rates (caculated)
  print("Checking MC clousre using the gen rates")
  (misIDRatiosNum, misIDRatios) = readMisIDRatiosGen(input_hadd_stage2)  # read 6 eMisId w.r.t. gen-pT-eta from stage2.root
  print("misIDRatios:: (6 gen rates) ")
  for i in range(6):
    print("\t {} {}:  {}".format(i, RATE_BINS[i], misIDRatios[RATE_BINS[i]]))

  # Calculate 21 rates from misIdRatios by adding the corresponding of the 6 rates
  catRatiosNum_genSum, catRatios_genSum = makeCatRatiosFrom6(misIDRatios)
  print("catRatios_genSum:: (6 gen rates ..> Add them to get 21 rations) ")
  for key, value in catRatios_genSum.items():
    print("\t {}:  {}".format(key, value))

  print("21Cat --> 6Rates:: ( gen rates ..> Add them to get 21 rations --> Solve 21 equation for 6 rates) ");
  (rates, uncs) = calculate(catRatiosNum_genSum, exclude_bins)
  print("CalculatedRates: misIDRatios(gen) --> catRatios_genSum  --> Rates::")
  for i in range(0, len(rates)):
    print("\t {}:  {} +- {}".format(i, rates[i], uncs[i]))

  compare_misIdRatios(misIDRatiosNum, rates, uncs, name = "gen_closure", mydir = output_dir, outFileName = "")

  print("Check 2.2] Check MC clousre using the gen rates:: 6 gen_rec rates ..> Add them to get 21 rations --> Solve 21 equation for 6 rates --> 6 gen_rec rates (caculated) :: ")
  (misIDRatiosNum, misIDRatios) = readMisIDRatiosGen(input_hadd_stage2,
                                                     rec = "_rec")  # read 6 eMisId w.r.t. gen-pT-eta from stage2.root
  print("misIDRatios:: (6 gen_rec rates) ")
  for i in range(6):
    print("\t {} {}:  {}".format(i, RATE_BINS[i], misIDRatios[RATE_BINS[i]]))

  # Calculate 21 rates from misIdRatios by adding the corresponding of the 6 rates
  catRatiosNum_genSum, catRatios_genSum = makeCatRatiosFrom6(misIDRatios)
  print("catRatios_genSum:: (6 gen_rec rates ..> Add them to get 21 rations) ")
  for key, value in catRatios_genSum.items():
    print("\t {}:  {}".format(key, value))

  print("21Cat --> 6Rates:: ( gen_rec rates ..> Add them to get 21 rations --> Solve 21 equation for 6 rates) ");
  (rates, uncs) = calculate(catRatiosNum_genSum, exclude_bins)
  print("CalculatedRates: misIDRatios(gen_rec) --> catRatios_genSum  --> Rates::")
  for i in range(0, len(rates)):
    print("\t {}:  {} +- {}".format(i, rates[i], uncs[i]))

  compare_misIdRatios(misIDRatiosNum, rates, uncs, name = "gen_rec_closure", mydir = output_dir, outFileName = "")

  print("Check 3] Compare 21 ratios (gen) read from stage2_mass_ll OS/SS histograms and 21 ratios calculated from 6 (gen) rates :: ")
  (misIDRatiosNum, misIDRatios) = readMisIDRatiosGen(input_hadd_stage2)  # read 6 eMisId w.r.t. gen-pT-eta from stage2.root
  print_ratios_latex(misIDRatiosNum, "gen");
  catRatiosNum, catRatios = readCategoryRatiosGen(input_hadd_stage2)  # calculate 21 r=SS/(OS+SS) from gen_massll from stage2.root
  chi2s = make_pull_plot_21(misIDRatios, catRatios, mydir = output_dir, y_range = (-0.001, 0.011),
                            outFileName = "fit_output_pseudodata_%s/fit_res.root" % FITNAME)

  print("Check 4] Compare 21 ratios (gen_rec w.r.t reconstructed pT, eta) read from stage2_mass_ll OS/SS histograms and 21 ratios calculated from 6 (gen) rates :: ")
  (misIDRatiosNum, misIDRatios) = readMisIDRatiosGen(input_hadd_stage2, rec = "_rec")
  print_ratios_latex(misIDRatiosNum, "gen_rec");
  catRatiosNum, catRatios = readCategoryRatiosGen(input_hadd_stage2, gen = "gen_rec")
  chi2s = make_pull_plot_21(misIDRatios, catRatios, mydir = output_dir, y_range = (-0.001, 0.011), name = "gen_rec",
                            outFileName = "fit_output_pseudodata_%s/fit_res.root" % FITNAME)

  print("Check 5] Closure test to see if we get 6 number back if we construct the 21 from the 6 (gen_rec) rates and then fit :: ")
  print "Turns out this underestimates uncertainty (due to correlations)"
  (exclude_bins, exclude_bins_num) = read_exclude_bins(EXCLUDED_FILE)
  catRatiosNum, catRatios = makeCatRatiosFrom6(misIDRatios, exclude_bins)
  calculate_solution(catRatiosNum, exclude_bins_num, FITNAME, "closure", "pseudodata")

  print("Step 6] Actual calculation starts from here. Read 21 category ratios. Drop categories which we don't want to consider or those have large chi2 in (gen_rec) comparison of 21 rates (from mass_ll histogram) and (from sum of 6 gen_rec rates) (Check 4) ")
  catRatiosNum, catRatios = read_category_ratios("fit_output_pseudodata_%s/results_cat.txt" % (FITNAME))
  # Selects categories to drop and writes them to file
  # Comment out the following line and edit the file manually to specify the categories to be dropped yourself
  nans = []
  nans = select_categories(chi2s, catRatios)

  print("Step 6.1] Make pull plots for both cases of not dropping and dropping categories")
  for exclude in [False, True]:
    fittypestring = "_gen_rec"
    file_misId = "fit_output_pseudodata_%s/fit_res%s.root" % (FITNAME, fittypestring)
    name = "gen_rec_fit"
    exclude_bins, exclude_bins_num = [], []
    if exclude:
      (exclude_bins, exclude_bins_num) = read_exclude_bins(EXCLUDED_FILE)
      name += "_exclusions"
      fittypestring += "_exclusions"
    catRatiosNum, catRatios = readCategoryRatiosGen(input_hadd_stage2, exclude_bins, gen = "gen_rec")
    calculate_solution(catRatiosNum, exclude_bins_num, FITNAME, fittypestring, "pseudodata")
    misIDRatios = readMisIDRatios(file_misId)
    make_pull_plot_21(misIDRatios, catRatios, mydir = output_dir, name = name, y_range = (-0.001, 0.011),
                      excluded = exclude_bins)

  print("Step 6.2] Fit results for pseudodata and data first without and then with excluding some categories")
  FITTYPE = ""  # can use also "shapes" or "hybrid" here if the fit results exist
  # fCompare = TFile("CompareResults.root","recreate")
  for datastring in ["pseudodata", "data"]:
    fittypestring = FITTYPE
    for exclude in [False, True]:
      if len(FITTYPE) > 0: fittypestring = "_" + FITTYPE
      file_cats = "fit_output_%s_%s/results_cat%s.txt" % (datastring, FITNAME, fittypestring)
      print "Data: ", datastring, ", exclude: ", exclude, ", i/p result file: ", file_cats;
      name = datastring
      # Still exclude NaN fit results from non_exclusion results
      exclude_bins, exclude_bins_num = nans, map(get_bin_nr, nans)
      if exclude:
        (exclude_bins, exclude_bins_num) = read_exclude_bins(EXCLUDED_FILE)
        name += "_exclusions"
        fittypestring += "_exclusions"
      file_misId = "fit_output_%s_%s/fit_res%s.root" % (datastring, FITNAME, fittypestring)
      catRatiosNum, catRatios = read_category_ratios(file_cats, exclude_bins)
      calculate_solution(catRatiosNum, exclude_bins_num, FITNAME, fittypestring, datastring)
      misIDRatios = readMisIDRatios(file_misId)

      make_pull_plot_21(misIDRatios, catRatios, mydir = output_dir, name = name, y_range = (-0.001, 0.011),
                        excluded = exclude_bins, outFileName = file_misId)

  print("End..")

  # Read generator-level ratios (misIDRatios: 6 for single electrons, catRatios for 21 categories of double electrons)
  # Numeric values are used for matrix solver
  (misIDRatiosNum, misIDRatios) = readMisIDRatiosGen(input_hadd_stage2)  # read 6 eMisId w.r.t. gen-pT-eta from stage2.root
  catRatiosNum, catRatios = readCategoryRatiosGen(input_hadd_stage2)  # calculate 21 r=SS/(OS+SS) from gen_massll from stage2.root
  # Print latex results for gen-level ratios
  print_ratios_latex(misIDRatiosNum, "gen")

  # Makes a pull plot comparing the 21 numbers to sums of respective ones from 6
  print("The ratios for gen-level electrons with gen-level pT and eta:")
  # chi2s = make_pull_plot_21(misIDRatios, catRatios, mydir = mydir1, y_range = (-0.001, 0.011))
  chi2s = make_pull_plot_21(misIDRatios, catRatios, mydir = output_dir, y_range = (-0.001, 0.011),
                            outFileName = "fit_output_pseudodata_%s/fit_res.root" % FITNAME)

  print("The ratios for gen-level electrons with reconstructed pT and eta:")
  (misIDRatiosNum, misIDRatios) = readMisIDRatiosGen(input_hadd_stage2, rec = "_rec")
  catRatiosNum, catRatios = readCategoryRatiosGen(input_hadd_stage2, gen = "gen_rec")
  print_ratios_latex(misIDRatiosNum, "gen_rec")
  # chi2s = make_pull_plot_21(misIDRatios, catRatios, mydir = mydir1, y_range = (-0.001, 0.011), name = "gen_rec")
  chi2s = make_pull_plot_21(misIDRatios, catRatios, mydir = output_dir, y_range = (-0.001, 0.011), name = "gen_rec",
                            outFileName = "fit_output_pseudodata_%s/fit_res.root" % FITNAME)

  print("Closure test to see if we get 6 number back if we construct the 21 from the 6 and then fit")
  print("Turns out this underestimates uncertainty (due to correlations)")
  (exclude_bins, exclude_bins_num) = read_exclude_bins(EXCLUDED_FILE)
  catRatiosNum, catRatios = makeCatRatiosFrom6(misIDRatios, exclude_bins)
  calculate_solution(catRatiosNum, exclude_bins_num, FITNAME, "closure", "pseudodata")

  catRatiosNum, catRatios = read_category_ratios("fit_output_pseudodata_%s/results_cat.txt" % (FITNAME))
  # Selects categories to drop and writes them to file
  # Comment out the following line and edit the file manually to specify the categories to be dropped yourself
  nans = []
  nans = select_categories(chi2s, catRatios)

  print("Make pull plots for both cases of not dropping and dropping categories")
  for exclude in [False, True]:
    fittypestring = "_gen_rec"
    file_misId = "fit_output_pseudodata_%s/fit_res%s.root" % (FITNAME, fittypestring)
    name = "gen_rec_fit"
    exclude_bins, exclude_bins_num = [], []
    if exclude:
      (exclude_bins, exclude_bins_num) = read_exclude_bins(EXCLUDED_FILE)
      name += "_exclusions"
      fittypestring += "_exclusions"
    catRatiosNum, catRatios = readCategoryRatiosGen(input_hadd_stage2, exclude_bins, gen = "gen_rec")
    calculate_solution(catRatiosNum, exclude_bins_num, FITNAME, fittypestring, "pseudodata")
    misIDRatios = readMisIDRatios(file_misId)
    make_pull_plot_21(misIDRatios, catRatios, mydir = output_dir, name = name, y_range = (-0.001, 0.011),
                      excluded = exclude_bins)

  print("Fit results for pseudodata and data first without and then with excluding some categories")
  FITTYPE = ""  # can use also "shapes" or "hybrid" here if the fit results exist
  # fCompare = TFile("CompareResults.root","recreate")
  for datastring in ["pseudodata", "data"]:
    fittypestring = FITTYPE
    for exclude in [False, True]:
      if len(FITTYPE) > 0: fittypestring = "_" + FITTYPE
      file_cats = "fit_output_%s_%s/results_cat%s.txt" % (datastring, FITNAME, fittypestring)
      print("Data: ", datastring, ", exclude: ", exclude, ", i/p result file: ", file_cats)
      name = datastring
      # Still exclude NaN fit results from non_exclusion results
      exclude_bins, exclude_bins_num = nans, map(get_bin_nr, nans)
      if exclude:
        (exclude_bins, exclude_bins_num) = read_exclude_bins(EXCLUDED_FILE)
        name += "_exclusions"
        fittypestring += "_exclusions"
      file_misId = "fit_output_%s_%s/fit_res%s.root" % (datastring, FITNAME, fittypestring)
      catRatiosNum, catRatios = read_category_ratios(file_cats, exclude_bins)
      calculate_solution(catRatiosNum, exclude_bins_num, FITNAME, fittypestring, datastring)
      misIDRatios = readMisIDRatios(file_misId)

      make_pull_plot_21(misIDRatios, catRatios, mydir = output_dir, name = name, y_range = (-0.001, 0.011),
                        excluded = exclude_bins, outFileName = file_misId)

  print("End..")