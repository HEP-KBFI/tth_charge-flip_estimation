#!/usr/bin/env python

from tthAnalysis.ChargeFlipEstimation.utils import read_category_ratios, readMisIDRatios, fit_results_to_file, \
                                                   BIN_NAMES_SINGLE, BIN_NAMES_COMPOSITE_NICE, get_bin_nr, SmartFormatter
from tthAnalysis.ChargeFlipEstimation.matrix_solver import calculate_solution, calculate
from tthAnalysis.ChargeFlipEstimation.plot_pulls import readMisIDRatiosGen, readCategoryRatiosGen, make_pull_plot_21, \
                                                        makeCatRatiosFrom6, compare_misIdRatios
import scipy.stats
import ROOT
import math
import argparse
import os
import copy
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
RATE_BINS = { cat_idx : cat for cat_idx, cat in enumerate(BIN_NAMES_SINGLE) }

def read_fits(input_dir):
  fit_results = os.path.join(input_dir, "results_cat.txt")
  catRatiosNum, catRatios = read_category_ratios(fit_results)
  nans = []
  for bin_name in catRatios:
    if math.isnan(catRatios[bin_name][0]):
      nans.append(bin_name)
  nans_num = [ get_bin_nr(bin_name) for bin_name in nans ]
  print("Read the following rates from file {}".format(fit_results))
  for bin_name in BIN_NAMES_COMPOSITE_NICE:
    catRatio = catRatios[bin_name]
    print("  {}: {} + {} - {}".format(bin_name, catRatio[0], catRatio[1], catRatio[2]))
  return catRatios, catRatiosNum, nans, nans_num

def merge_excludable_bins(exclude_by_singificance = None, exclude_by_nan = None, exclude_by_request = None):
  exclude_result = []
  if exclude_by_singificance:
    assert(all(bin_name in BIN_NAMES_COMPOSITE_NICE for bin_name in exclude_by_singificance))
    for bin_name in exclude_by_singificance:
      if bin_name not in exclude_result:
        exclude_result.append(bin_name)
    print("Excluding because not passing the signficance test: {}".format(", ".join(exclude_by_singificance)))
  if exclude_by_nan:
    assert (all(bin_name in BIN_NAMES_COMPOSITE_NICE for bin_name in exclude_by_nan))
    for bin_name in exclude_by_nan:
      if bin_name not in exclude_result:
        exclude_result.append(bin_name)
    print("Excluding because fit resulted in NaN: {}".format(", ".join(exclude_by_nan)))
  if exclude_by_request:
    assert (all(bin_name in BIN_NAMES_COMPOSITE_NICE for bin_name in exclude_by_request))
    for bin_name in exclude_by_nan:
      if bin_name not in exclude_result:
        exclude_result.append(bin_name)
    print("Excluding because explicitly requested to skip: {}".format(", ".join(exclude_by_request)))
  print("Verdit -> excluding the following categories: {}".format(", ".join(exclude_result)))
  exclude_result_num = [ get_bin_nr(bin_name) for bin_name in exclude_result ]
  return exclude_result, exclude_result_num

def exclude_cats(catRatiosNum, catRatios, exclude_bins, exclude_bins_num):
  assert(len(exclude_bins) == len(exclude_bins_num))
  catRatios_excl = copy.deepcopy(catRatios)
  for bin_name in exclude_bins:
    assert(bin_name in catRatios_excl)
    del catRatios_excl[bin_name]
  catRatiosNum_excl = []
  for bin_idx, ratios in enumerate(catRatiosNum):
    if bin_idx not in exclude_bins_num:
      catRatiosNum_excl.append(ratios)
  assert(len(catRatios_excl) == len(catRatiosNum_excl))
  return catRatios_excl, catRatiosNum_excl

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
    type = str, dest = 'exclude', metavar = 'bin', required = False, nargs = '+', default = [],
    help = 'R|Additional bins to exclude (in nice nomenclature, eg BL_BL not BB_LL)',
  )
  args = parser.parse_args()

  input_hadd_stage2 = args.input_hadd
  input_data_dir = args.fits_data
  input_pseudodata_dir = args.fits_pseudodata
  output_dir = os.path.abspath(args.output)
  exclude_bins_additional = args.exclude
  assert(os.path.isfile(input_hadd_stage2))
  assert(os.path.isdir(input_data_dir))
  assert (os.path.isdir(input_pseudodata_dir))
  assert(all(bin in BIN_NAMES_COMPOSITE_NICE for bin in exclude_bins_additional))

  print("Input hadd stage2 file:     {}".format(input_hadd_stage2))
  print("Input data directory:       {}".format(input_data_dir))
  print("Input pseudodata directory: {}".format(input_pseudodata_dir))
  print("Output directory:           {}".format(output_dir))
  print("Bins to explicitly exclude: {}".format(exclude_bins_additional))

  # Check 1: MC closure to solve 21 linear equations using dummy 6 rates:
  #   6 dummy rates ->
  #   Add them to get 21 rations ->
  #   Solve 21 equation for 6 rates ->
  #   6 rates (caculated)
  print('=' * 120)
  print("Checking MC closure to solve 21 linear equations using dummy 6 rates")
  exclude_bins_dummy = []
  rate_dummy = 2.e-5
  eRate_dummy = 1.e-6
  misIDRatiosNum_dummy = [ (rate_dummy, eRate_dummy) for _ in RATE_BINS  ]
  misIDRatios_dummy = { BIN_NAMES_SINGLE[bin_idx] : (rate_dummy, eRate_dummy, eRate_dummy) for bin_idx in RATE_BINS }
  print("Dummy mis-identification ratios:")
  for bin_idx, bin in RATE_BINS.items():
    print("  {} {}:  {} +/- {}".format(bin_idx, bin, misIDRatios_dummy[bin][0], misIDRatios_dummy[bin][1]))

  catRatiosNum_dummySum, catRatios_dummySum = makeCatRatiosFrom6(misIDRatios_dummy, exclude_bins_dummy)
  print("Ratios by category for 6 dummy rates:")
  for bin_idx, value in catRatios_dummySum.items():
    print("  {}: {} +/- {}".format(bin_idx, value[0], value[1]))

  rates_dummy, uncs_dummy = calculate(catRatiosNum_dummySum, exclude_bins_dummy)
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
  exclude_bins_gen = []
  misIDRatiosNum_gen, misIDRatios_gen = readMisIDRatiosGen(input_hadd_stage2, rec = False)
  print("Generator-level mis-identification ratios:")
  for bin_idx, bin in RATE_BINS.items():
    print("  {} {}:  {} +/- {}".format(bin_idx, bin, misIDRatios_gen[bin][0], misIDRatios_gen[bin][1]))

  catRatiosNum_genSum, catRatios_genSum = makeCatRatiosFrom6(misIDRatios_gen, exclude_bins_gen)
  print("Ratios by category for 6 generator-level rates:")
  for bin_idx, value in catRatios_genSum.items():
    print("  {}: {} +/- {}".format(bin_idx, value[0], value[1]))

  rates_gen, uncs_gen = calculate(catRatiosNum_genSum, exclude_bins_gen)
  print("Calculated generator-level rates:")
  for bin_idx, rate in enumerate(rates_gen):
    print("  {}: {} +- {}".format(bin_idx, rate, uncs_gen[bin_idx]))
  compare_misIdRatios(misIDRatiosNum_gen, rates_gen, uncs_gen, name = "gen_closure", output_dir = output_dir)

  print("Comparing 21 ratios from di-lepton mass histograms to the calculated gen rates")
  catRatiosNum_gen, catRatios_gen = readCategoryRatiosGen(input_hadd_stage2, is_gen = True)
  chi2s_gen = make_pull_plot_21(
    misIDRatios_gen, catRatios_gen, name = "gen", output_dir = output_dir,
    y_range = (-0.001, 0.011), excluded = exclude_bins_gen,
  )
  print("chi2:")
  for bin, chi2 in chi2s_gen.items():
    print("  {} {:.3f}{}".format(bin, chi2, " > {:.3f} => exclude".format(NSIGMAS) if chi2 > NSIGMAS else ""))
  exclude_gen = [ bin for bin, chi2 in chi2s_gen.items() if chi2 > NSIGMAS ]
  nof_exclude_gen = len(exclude_gen)
  print("To keep = {} / exclude = {}".format(len(chi2s_gen) - nof_exclude_gen, nof_exclude_gen))
  print('=' * 120 + '\n')

  # Check 3: MC clousre using gen vs rec rates:
  #   6 gen vs rec rates ->
  #   Add them to get 21 rations ->
  #   Solve 21 equation for 6 rates ->
  #   6 gen vs rec rates (caculated)
  print('=' * 120)
  print("Checking MC closure to solve 21 linear equations using generator vs reconstruction level rates")
  exclude_bins_genRec = []
  misIDRatiosNum_genRec, misIDRatios_genRec = readMisIDRatiosGen(input_hadd_stage2, rec = True)
  print("Mis-identification ratios for generator vs reconstruction level:")
  for bin_idx, bin in RATE_BINS.items():
    print("  {} {}:  {} +/- {}".format(bin_idx, bin, misIDRatios_genRec[bin][0], misIDRatios_genRec[bin][1]))

  catRatiosNum_genRecSum, catRatios_genRecSum = makeCatRatiosFrom6(misIDRatios_genRec)
  print("Ratios by category for 6 generator vs reconstruction level rates:")
  for bin_idx, value in catRatios_genRecSum.items():
    print("  {}: {} +/- {}".format(bin_idx, value[0], value[1]))

  rates_genRec, uncs_genRec = calculate(catRatiosNum_genRecSum, exclude_bins_genRec)
  print("Calculated generator vs reconstruction level rates:")
  for bin_idx, rate in enumerate(rates_genRec):
    print("  {}: {} +- {}".format(bin_idx, rate, uncs_genRec[bin_idx]))
  compare_misIdRatios(misIDRatiosNum_genRec, rates_genRec, uncs_genRec, name = "genRec_closure", output_dir = output_dir)

  print("Comparing 21 ratios from di-lepton mass histograms to the calculated generator vs reconstruction level rates")
  catRatiosNum_genRec, catRatios_genRec = readCategoryRatiosGen(input_hadd_stage2, is_gen = False)
  chi2s_genRec = make_pull_plot_21(
    misIDRatios_genRec, catRatios_genRec, name = "genRec", output_dir = output_dir,
    y_range = (-0.001, 0.011), excluded = exclude_bins_genRec,
  )
  print("chi2:")
  for bin, chi2 in chi2s_genRec.items():
    print("  {} {:.3f}{}".format(bin, chi2, " > {:.3f} => exclude".format(NSIGMAS) if chi2 > NSIGMAS else ""))
  exclude_genRec = [ bin for bin, chi2 in chi2s_genRec.items() if chi2 > NSIGMAS ]
  nof_exclude_genRec = len(exclude_genRec)
  print("To keep = {} / exclude = {}".format(len(chi2s_genRec) - nof_exclude_genRec, nof_exclude_genRec))

  # Check 4: test to see if we get 6 number back if we construct the 21 rates from
  #          the 6 generator vs reconstruction level rates and then fit
  exclude_testGenRec, exclude_testGenRec_num = merge_excludable_bins(
    exclude_by_singificance = exclude_genRec, exclude_by_request = exclude_bins_additional
  )
  catRatiosNum_testGenRec, catRatios_testGenRec = makeCatRatiosFrom6(misIDRatios_genRec, exclude_testGenRec_num)
  rates_testGenRec, uncs_testGenRec = calculate(catRatiosNum_testGenRec, exclude_testGenRec_num)
  print("Re-calculated generator vs reconstruction level rates:")
  for bin_idx, rate in enumerate(rates_testGenRec):
    print("  {}: {} +- {}".format(bin_idx, rate, uncs_testGenRec[bin_idx]))
  fitResult_genRec = os.path.join(output_dir, "fit_result_genRec.root")
  fit_results_to_file(rates_testGenRec, uncs_testGenRec, fitResult_genRec)
  print('=' * 120 + '\n')

  # Actual calculation starts from here:
  #   - read 21 category ratios
  #   - drop categories which we don't want to consider or those that have large chi2 in Check 3/4
  print('=' * 120)
  catRatios_pseudo, catRatiosNum_pseudo, nans_pseudo, nans_pseudo_num = read_fits(input_pseudodata_dir)

  # Creating pull plots for both cases of not dropping (categories except the ones with NaN ratio) and dropping categories
  for exclude in [ False, True ]:
    exclude_suffix = "_exclusions" if exclude else ""
    output_fit_pseudo_excl = os.path.join(output_dir, "fit_result_genRec_pseudo{}.root".format(exclude_suffix))

    exclude_bins_genRec_excl, exclude_bins_num_genRec_excl = merge_excludable_bins(
      exclude_by_singificance = exclude_genRec if exclude else None,
      exclude_by_nan          = nans_pseudo,
      exclude_by_request      = exclude_bins_additional
    )
    catRatiosNum_genRec_excl, catRatios_genRec_excl = readCategoryRatiosGen(
      input_hadd_stage2, is_gen = False, exclude_bins = exclude_bins_genRec_excl
    )
    calculate_solution(catRatiosNum_genRec_excl, exclude_bins_num_genRec_excl, output_fit_pseudo_excl)
    misIDRatios_genRec_excl = readMisIDRatios(output_fit_pseudo_excl)
    chi2s_genRec_excl = make_pull_plot_21(
      misIDRatios_genRec_excl, catRatios_genRec_excl, name = "genRec_fit{}".format(exclude_suffix),
      output_dir = output_dir, y_range = (-0.001, 0.011), excluded = exclude_bins_num_genRec_excl
    )
    print("chi2:")
    for bin, chi2 in chi2s_genRec_excl.items():
      print("  {} {:.3f}{}".format(bin, chi2, " > {:.3f} => exclude".format(NSIGMAS) if chi2 > NSIGMAS else ""))
  print('=' * 120 + '\n')

  # Fit results for pseudodata and data first without and then with excluding some categories
  print('=' * 120)
  catRatios_data, catRatiosNum_data, nans_data, nans_data_num = read_fits(input_data_dir)
  for is_data in [ False, True ]:
    for exclude in [ False, True ]:
      name = "{}data{}".format("" if is_data else "pseudo", "_exclusions" if exclude else "")
      # read the fit ratios first!
      exclude_bins_excl, exclude_bins_num_excl = merge_excludable_bins(
        exclude_by_singificance = exclude_genRec if exclude else None,
        exclude_by_nan          = nans_data if is_data else nans_pseudo,
        exclude_by_request      = exclude_bins_additional
      )
      output_fit_excl = os.path.join(output_dir, "fit_result_{}.root".format(name))
      catRatios_excl, catRatiosNum_excl = exclude_cats(
        catRatiosNum_data if is_data else catRatiosNum_pseudo,
        catRatios_data    if is_data else catRatios_pseudo,
        exclude_bins_excl, exclude_bins_num_excl
      )
      calculate_solution(catRatiosNum_excl, exclude_bins_num_excl, output_fit_excl)
      misIDRatios_excl = readMisIDRatios(output_fit_excl)
      chi2s_excl = make_pull_plot_21(
        misIDRatios_excl, catRatios_excl, name = name, output_dir = output_dir,
        y_range = (-0.001, 0.011), excluded = exclude_bins_num_excl
      )
      print("chi2:")
      for bin, chi2 in chi2s_excl.items():
        print("  {} {:.3f}{}".format(bin, chi2, " > {:.3f} => exclude".format(NSIGMAS) if chi2 > NSIGMAS else ""))
  print('=' * 120 + '\n')
  sys.exit(0)

  # Read generator-level ratios (misIDRatios: 6 for single electrons, catRatios for 21 categories of double electrons)
  # Numeric values are used for matrix solver
  (misIDRatiosNum, misIDRatios) = readMisIDRatiosGen(input_hadd_stage2)  # read 6 eMisId w.r.t. gen-pT-eta from stage2.root
  catRatiosNum, catRatios = readCategoryRatiosGen(input_hadd_stage2)  # calculate 21 r=SS/(OS+SS) from gen_massll from stage2.root
  # Print latex results for gen-level ratios

  # Makes a pull plot comparing the 21 numbers to sums of respective ones from 6
  print("The ratios for gen-level electrons with gen-level pT and eta:")
  # chi2s = make_pull_plot_21(misIDRatios, catRatios, mydir = mydir1, y_range = (-0.001, 0.011))
  chi2s = make_pull_plot_21(misIDRatios, catRatios, mydir = output_dir, y_range = (-0.001, 0.011),
                            outFileName = "fit_output_pseudodata_%s/fit_res.root" % FITNAME)

  print("The ratios for gen-level electrons with reconstructed pT and eta:")
  (misIDRatiosNum, misIDRatios) = readMisIDRatiosGen(input_hadd_stage2, rec = "_rec")
  catRatiosNum, catRatios = readCategoryRatiosGen(input_hadd_stage2, gen = "gen_rec")
  # chi2s = make_pull_plot_21(misIDRatios, catRatios, mydir = mydir1, y_range = (-0.001, 0.011), name = "gen_rec")
  chi2s = make_pull_plot_21(misIDRatios, catRatios, mydir = output_dir, y_range = (-0.001, 0.011), name = "gen_rec",
                            outFileName = "fit_output_pseudodata_%s/fit_res.root" % FITNAME)

  print("Closure test to see if we get 6 number back if we construct the 21 from the 6 and then fit")
  print("Turns out this underestimates uncertainty (due to correlations)")
  #(exclude_bins, exclude_bins_num) = read_exclude_bins(EXCLUDED_FILE)
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
      #(exclude_bins, exclude_bins_num) = read_exclude_bins(EXCLUDED_FILE)
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
        #(exclude_bins, exclude_bins_num) = read_exclude_bins(EXCLUDED_FILE)
        name += "_exclusions"
        fittypestring += "_exclusions"
      file_misId = "fit_output_%s_%s/fit_res%s.root" % (datastring, FITNAME, fittypestring)
      catRatiosNum, catRatios = read_category_ratios(file_cats, exclude_bins)
      calculate_solution(catRatiosNum, exclude_bins_num, FITNAME, fittypestring, datastring)
      misIDRatios = readMisIDRatios(file_misId)

      make_pull_plot_21(misIDRatios, catRatios, mydir = output_dir, name = name, y_range = (-0.001, 0.011),
                        excluded = exclude_bins, outFileName = file_misId)

  print("End..")