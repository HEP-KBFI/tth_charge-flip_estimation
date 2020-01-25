#!/usr/bin/env python

from tthAnalysis.ChargeFlipEstimation.utils import read_category_ratios, readMisIDRatios, read_exclude_bins, \
                                                   BIN_NAMES_SINGLE, get_bin_nr, SmartFormatter
from tthAnalysis.ChargeFlipEstimation.matrix_solver import calculate_solution, print_ratios_latex, calculate
from tthAnalysis.ChargeFlipEstimation.plot_pulls import readMisIDRatiosGen, readCategoryRatiosGen, make_pull_plot_21, \
                                                        makeCatRatiosFrom6, compare_misIdRatios
import scipy.stats
import ROOT
import math
import argparse
import os
import sys

ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2()

"""@file docstring
Script for plotting pulls comparing result from 21 and 6 categories
Selects which categories of 21 to drop because of correlations and solves the equations for different cases

@author Andres Tiko <andres.tiko@cern.ch>
"""

# Number of sigmas difference to consider fit results not compatible
PVALUE = 0.1
NSIGMAS = scipy.stats.norm.ppf(1. - PVALUE)

# File to store list of excluded categories
EXCLUDED_FILE = "/home/ssawant/VHbbNtuples_10_x/CMSSW_10_2_13/src/tthAnalysis/ChargeFlipEstimation/bin/../data/excluded_categories.txt"

EXCLUDED_CATEGORIES_2016 = [ "BM_BL", "EM_EM", "EH_EH", "EM_BL", "BH_EH", "EH_EL", "BL_EL", "EL_EL", "BH_BL", "BH_BM" ]
EXCLUDED_CATEGORIES_2017 = [ "EM_EM", "EM_EL", "EH_EH", "BL_BL", "EH_BM", "EH_EL", "EL_EL", "EH_BL", "BH_BM" ]
EXCLUDED_CATEGORIES_2018 = [ "BM_BL", "EM_EM", "EM_EL", "BH_EM", "BL_BL", "EH_BM", "EH_EL", "BH_BH", "BL_EL", "EL_EL", "BH_BL", "BH_BM" ]

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

RATE_BINS = { cat_idx : cat for cat_idx, cat in enumerate(BIN_NAMES_SINGLE) }

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

  input_hadd_stage2 = args.input_hadd
  output_dir = os.path.abspath(args.output)
  exclude_bins = args.exclude

  # Check 1: MC closure to solve 21 linear equations using dummy 6 rates:
  #   6 dummy rates ->
  #   Add them to get 21 rations ->
  #   Solve 21 equation for 6 rates ->
  #   6 rates (caculated)
  print('=' * 120)
  print("Checking MC closure to solve 21 linear equations using dummy 6 rates")
  rate_dummy = 2.e-5
  eRate_dummy = 1.e-6
  misIDRatiosNum_dummy = [ (rate_dummy, eRate_dummy) for _ in RATE_BINS  ]
  misIDRatios_dummy = { BIN_NAMES_SINGLE[bin_idx] : (rate_dummy, eRate_dummy, eRate_dummy) for bin_idx in RATE_BINS }
  print("Dummy mis-identification ratios:")
  for bin_idx, bin in RATE_BINS.items():
    print("  {} {}:  {} +/- {}".format(bin_idx, bin, misIDRatios_dummy[bin][0], misIDRatios_dummy[bin][1]))

  # Calculate 21 rates from misIDRatios_dummy by adding the corresponding of the 6 rates
  catRatiosNum_dummySum, catRatios_dummySum = makeCatRatiosFrom6(misIDRatios_dummy, exclude_bins)
  print("Ratios by category for 6 dummy rates:")
  for bin_idx, value in catRatios_dummySum.items():
    print("  {}: {} +/- {}".format(bin_idx, value[0], value[1]))

  rates_dummy, uncs_dummy = calculate(catRatiosNum_dummySum, exclude_bins)
  print("Calculated dummy rates:")
  for bin_idx, rate in enumerate(rates_dummy):
    print("  {}: {} +- {}".format(bin_idx, rate, uncs_dummy[bin_idx]))
  compare_misIdRatios(misIDRatiosNum_dummy, rates_dummy, uncs_dummy, name = "dummy_closure", output_dir = output_dir)
  print('=' * 120 + '\n')

  # Check 2: MC clousre using the gen rates:
  #   6 gen rates ->
  #   Add them to get 21 rations ->
  #   Solve 21 equation for 6 rates ->
  #   6 gen rates (caculated)
  print('=' * 120)
  print("Checking MC closure to solve 21 linear equations using generator-level rates")
  misIDRatiosNum_gen, misIDRatios_gen = readMisIDRatiosGen(input_hadd_stage2, rec = False)
  print("Generator-level mis-identification ratios:")
  for bin_idx, bin in RATE_BINS.items():
    print("  {} {}:  {} +/- {}".format(bin_idx, bin, misIDRatios_gen[bin][0], misIDRatios_gen[bin][1]))

  # Calculate 21 rates from misIDRatios_gen by adding the corresponding of the 6 rates
  catRatiosNum_genSum, catRatios_genSum = makeCatRatiosFrom6(misIDRatios_gen, exclude_bins)
  print("Ratios by category for 6 generator-level rates:")
  for bin_idx, value in catRatios_genSum.items():
    print("  {}: {} +/- {}".format(bin_idx, value[0], value[1]))

  rates_gen, uncs_gen = calculate(catRatiosNum_genSum, exclude_bins)
  print("Calculated generator-level rates:")
  for bin_idx, rate in enumerate(rates_gen):
    print("  {}: {} +- {}".format(bin_idx, rate, uncs_gen[bin_idx]))
  compare_misIdRatios(misIDRatiosNum_gen, rates_gen, uncs_gen, name = "gen_closure", output_dir = output_dir)
  print('=' * 120 + '\n')

  # Check 3: MC clousre using gen vs rec rates:
  #   6 gen vs rec rates ->
  #   Add them to get 21 rations ->
  #   Solve 21 equation for 6 rates ->
  #   6 gen vs rec rates (caculated)
  print('=' * 120)
  print("Checking MC closure to solve 21 linear equations using generator vs reconstruction level rates")
  misIDRatiosNum_genRec, misIDRatios_genRec = readMisIDRatiosGen(input_hadd_stage2, rec = True)
  print("Mis-identification ratios for generator vs reconstruction level:")
  for bin_idx, bin in RATE_BINS.items():
    print("  {} {}:  {} +/- {}".format(bin_idx, bin, misIDRatios_genRec[bin][0], misIDRatios_genRec[bin][1]))

  # Calculate 21 rates from misIDRatios_genRec by adding the corresponding of the 6 rates
  catRatiosNum_genRecSum, catRatios_genRecSum = makeCatRatiosFrom6(misIDRatios_genRec)
  print("Ratios by category for 6 generator vs reconstruction level rates:")
  for bin_idx, value in catRatios_genRecSum.items():
    print("  {}: {} +/- {}".format(bin_idx, value[0], value[1]))

  rates_genRec, uncs_genRec = calculate(catRatiosNum_genRecSum, exclude_bins)
  print("Calculated generator vs reconstruction level rates:")
  for bin_idx, rate in enumerate(rates_genRec):
    print("  {}: {} +- {}".format(bin_idx, rate, uncs_gen[bin_idx]))
  compare_misIdRatios(misIDRatiosNum_genRec, rates_genRec, uncs_genRec, name = "genRec_closure", output_dir = output_dir)
  print('=' * 120 + '\n')

  sys.exit(0)

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