#!/usr/bin/env python
# -*- coding: utf-8 -*-

from tthAnalysis.ChargeFlipEstimation.utils import read_category_ratios, readMisIDRatios, fit_results_to_file, \
                                                   BIN_NAMES_SINGLE, BIN_NAMES_COMPOSITE_NICE, get_bin_nr, SmartFormatter
from tthAnalysis.ChargeFlipEstimation.matrix_solver import calculate_solution, calculate, get_solution_latex, \
                                                           calculate_solution_wToys, calculate_wToys
from tthAnalysis.ChargeFlipEstimation.plot_pulls import readMisIDRatiosGen, readCategoryRatiosGen, make_pull_plot_21, \
                                                        makeCatRatiosFrom6, compare_misIdRatios, plot_ratios, plot_rates
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
@author Karl Ehatäht <karl.ehataht@cern.ch>
"""

# Number of sigmas difference to consider fit results not compatible
PVALUE = 0.001 
NSIGMAS = scipy.stats.norm.ppf(1. - PVALUE)
RATE_BINS = { cat_idx : cat for cat_idx, cat in enumerate(BIN_NAMES_SINGLE) }
USE_GENREC_TO_EXCLUDE_CATEGORIES = True

def read_fits(fit_results):
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
  parser.add_argument('-l', '--latex',
    dest = 'latex', action = 'store_true', default = False,
    help = 'R|Print results in LaTeX-formatted table',
  )
  parser.add_argument('-s', '--seed',
    type = int, dest = 'seed', metavar = 'seed', required = False, default = 12345,
    help = 'R|Seed',
  )
  parser.add_argument('-t', '--toys',
    type = int, dest = 'toys', metavar = 'toys', required = False, default = -1,
    help = 'R|Use toys',
  )
  parser.add_argument('-p', '--placeholder',
    type = float, dest = 'placeholder', default = 3.5e-5, required = False,
    help = 'R|Value to use for negative flip rates',
  )
  args = parser.parse_args()

  input_hadd_stage2 = args.input_hadd
  input_data_file = args.fits_data
  input_pseudodata_file = args.fits_pseudodata
  output_dir = os.path.abspath(args.output)
  exclude_bins_additional = args.exclude
  print_in_latex = args.latex
  fallback_value = args.placeholder
  nof_toys = args.toys
  use_toys = nof_toys > 0
  seed = args.seed
  assert(os.path.isfile(input_hadd_stage2))
  assert(os.path.isfile(input_data_file))
  assert(os.path.isfile(input_pseudodata_file))
  assert(all(bin in BIN_NAMES_COMPOSITE_NICE for bin in exclude_bins_additional))

  print("Input hadd stage2 file:     {}".format(input_hadd_stage2))
  print("Input data file:            {}".format(input_data_file))
  print("Input pseudodata file:      {}".format(input_pseudodata_file))
  print("Output directory:           {}".format(output_dir))
  print("Bins to explicitly exclude: {}".format(exclude_bins_additional))
  print("Use toys:                   {}".format(use_toys))
  print("Number of toys:             {}".format(nof_toys))
  print("PVALUE:                     {}".format(PVALUE))
  print("  --> NSIGMAS:              {}".format(NSIGMAS))
  print("USE_GENREC_TO_EXCLUDE_CATEGORIES:  {}".format(USE_GENREC_TO_EXCLUDE_CATEGORIES))

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

  if use_toys:
    rates_dummy, uncs_dummy = calculate_wToys(catRatios_dummySum, exclude_bins_dummy, seed, nof_toys)
  else:
    rates_dummy, uncs_dummy = calculate(catRatiosNum_dummySum, exclude_bins_dummy)
  print("Calculated dummy rates:")
  for bin_idx, rate in enumerate(rates_dummy):
    print("  {}: {} +- {}".format(bin_idx, rate, uncs_dummy[bin_idx]))
  sys.stdout.flush()  
  compare_misIdRatios(misIDRatiosNum_dummy, rates_dummy, uncs_dummy, name = "dummy_closure", output_dir = output_dir)
  sys.stdout.flush()
  if print_in_latex:
    print(get_solution_latex(rates_dummy, uncs_dummy, "dummy"))
  print('=' * 120 + '\n')
  sys.stdout.flush()
  
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
  sys.stdout.flush()
  
  if use_toys:
    rates_gen, uncs_gen = calculate_wToys(catRatios_genSum, exclude_bins_gen, seed, nof_toys)
  else:
    rates_gen, uncs_gen = calculate(catRatiosNum_genSum, exclude_bins_gen)
  print("Calculated generator-level rates:")
  for bin_idx, rate in enumerate(rates_gen):
    print("  {}: {} +- {}".format(bin_idx, rate, uncs_gen[bin_idx]))
  if print_in_latex:
    print(get_solution_latex(rates_gen, uncs_gen, "generator-level"))
  sys.stdout.flush()
  compare_misIdRatios(misIDRatiosNum_gen, rates_gen, uncs_gen, name = "gen_closure", output_dir = output_dir)
  sys.stdout.flush()
  
  print("Comparing 21 ratios from di-lepton mass histograms to the calculated gen rates")
  catRatiosNum_gen, catRatios_gen = readCategoryRatiosGen(input_hadd_stage2, is_gen = True)
  sys.stdout.flush()
  chi2s_gen = make_pull_plot_21(
    misIDRatios_gen, catRatios_gen, name = "gen", output_dir = output_dir,
    y_range = (-0.001, 0.011), excluded = exclude_bins_gen,
  )
  sys.stdout.flush()
  print("chi2:")
  for bin, chi2 in chi2s_gen.items():
    print("  {} {:.3f}{}".format(bin, chi2, " > {:.3f} => exclude".format(NSIGMAS) if chi2 > NSIGMAS else ""))
  exclude_gen = [ bin for bin, chi2 in chi2s_gen.items() if chi2 > NSIGMAS ]
  nof_exclude_gen = len(exclude_gen)
  print("To keep = {} / exclude = {}".format(len(chi2s_gen) - nof_exclude_gen, nof_exclude_gen))
  print('=' * 120 + '\n')
  sys.stdout.flush()
  
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

  if use_toys:
    rates_genRec, uncs_genRec = calculate_wToys(catRatios_genRecSum, exclude_bins_genRec, seed, nof_toys)
  else:
    rates_genRec, uncs_genRec = calculate(catRatiosNum_genRecSum, exclude_bins_genRec)
  print("Calculated generator vs reconstruction level rates:")
  for bin_idx, rate in enumerate(rates_genRec):
    print("  {}: {} +- {}".format(bin_idx, rate, uncs_genRec[bin_idx]))
  if print_in_latex:
    print(get_solution_latex(rates_genRec, uncs_genRec, "generator vs reconstruction level"))
  sys.stdout.flush()
  compare_misIdRatios(misIDRatiosNum_genRec, rates_genRec, uncs_genRec, name = "genRec_closure", output_dir = output_dir)
  sys.stdout.flush()
  
  print("Comparing 21 ratios from di-lepton mass histograms to the calculated generator vs reconstruction level rates")
  catRatiosNum_genRec, catRatios_genRec = readCategoryRatiosGen(input_hadd_stage2, is_gen = False)
  sys.stdout.flush()
  chi2s_genRec = make_pull_plot_21(
    misIDRatios_genRec, catRatios_genRec, name = "genRec", output_dir = output_dir,
    y_range = (-0.001, 0.011), excluded = exclude_bins_genRec,
  )
  sys.stdout.flush()
  print("catRatios_genRec vs catRatios_genRecSum: ")
  for bin, r in catRatios_genRec.items():
    print(" {}: {:.5f} +- {:.5f}   {:.5f} +- {:.5f}".format(bin, catRatios_genRec[bin][0],catRatios_genRec[bin][1],  catRatios_genRecSum[bin][0],catRatios_genRecSum[bin][1]))
  print("catRatios_genRec:    {}".format(catRatios_genRec))
  print("catRatios_genRecSum: {}".format(catRatios_genRecSum))
  sys.stdout.flush()
  print("chi2:")
  for bin, chi2 in chi2s_genRec.items():
    print("  {} {:.3f}{}".format(bin, chi2, " > {:.3f} => exclude".format(NSIGMAS) if chi2 > NSIGMAS else ""))
  exclude_genRec = [ bin for bin, chi2 in chi2s_genRec.items() if chi2 > NSIGMAS ]
  nof_exclude_genRec = len(exclude_genRec)
  print("To keep = {} / exclude = {}".format(len(chi2s_genRec) - nof_exclude_genRec, nof_exclude_genRec))
  sys.stdout.flush()
  print("exclude_genRec ({}): {}".format(len(exclude_genRec), exclude_genRec))
  sys.stdout.flush()
  print("Siddh")
  sys.stdout.flush()
  
  # Check 4: test to see if we get 6 number back if we construct the 21 rates from
  #          the 6 generator vs reconstruction level rates and then fit
  if USE_GENREC_TO_EXCLUDE_CATEGORIES:
    exclude_testGenRec, exclude_testGenRec_num = merge_excludable_bins(
      exclude_by_singificance = exclude_genRec, exclude_by_request = exclude_bins_additional
    )
    catRatiosNum_testGenRec, catRatios_testGenRec = makeCatRatiosFrom6(misIDRatios_genRec, exclude_testGenRec_num)
    if use_toys:
      rates_testGenRec, uncs_testGenRec = calculate_wToys(catRatios_testGenRec, exclude_testGenRec_num, seed, nof_toys)
    else:
      rates_testGenRec, uncs_testGenRec = calculate(catRatiosNum_testGenRec, exclude_testGenRec_num)
    print("Re-calculated generator vs reconstruction level rates:")
    for bin_idx, rate in enumerate(rates_testGenRec):
      print("  {}: {} +- {}".format(bin_idx, rate, uncs_testGenRec[bin_idx]))
    if print_in_latex:
      print(get_solution_latex(rates_testGenRec, uncs_testGenRec, "generator vs reconstruction level (recomputed)"))
    fitResult_genRec = os.path.join(output_dir, "fit_result_genRec.root")
    fit_results_to_file(rates_testGenRec, uncs_testGenRec, fitResult_genRec, fallback_value)
    print('=' * 120 + '\n')
  sys.stdout.flush()
  
  # Actual calculation starts from here:
  #   - read 21 category ratios
  #   - drop categories which we don't want to consider or those that have large chi2 in Check 3/4
  print('=' * 120)
  catRatios_pseudo, catRatiosNum_pseudo, nans_pseudo, nans_pseudo_num = read_fits(input_pseudodata_file)
  sys.stdout.flush()
  
  # Creating pull plots for both cases of not dropping (categories except the ones with NaN ratio) and dropping categories
  for exclude in [ False, True ]:
    if not USE_GENREC_TO_EXCLUDE_CATEGORIES:
      break
    
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
    datastring = "gen vs reco ({} bins)".format("excluded" if exclude else "all")
    if use_toys:      
      calculate_solution_wToys(
        catRatios_genRec_excl, exclude_bins_num_genRec_excl, seed, nof_toys, output_fit_pseudo_excl, fallback_value,
        datastring if print_in_latex else ""
      )
    else:
      calculate_solution(
        catRatiosNum_genRec_excl, exclude_bins_num_genRec_excl, output_fit_pseudo_excl, fallback_value,
        datastring if print_in_latex else ""
      )
    misIDRatios_genRec_excl = readMisIDRatios(output_fit_pseudo_excl)
    sys.stdout.flush()
    chi2s_genRec_excl = make_pull_plot_21(
      misIDRatios_genRec_excl, catRatios_genRec_excl, name = "genRec_fit{}".format(exclude_suffix),
      output_dir = output_dir, y_range = (-0.001, 0.011), excluded = exclude_bins_num_genRec_excl
    )
    sys.stdout.flush()
    print("chi2:")
    for bin, chi2 in chi2s_genRec_excl.items():
      print("  {} {:.3f}{}".format(bin, chi2, " > {:.3f} => exclude".format(NSIGMAS) if chi2 > NSIGMAS else ""))      
    sys.stdout.flush()
  print('=' * 120 + '\n')
  sys.stdout.flush()
  
  # Fit results for pseudodata and data first without and then with excluding some categories
  print('=' * 120)
  catRatios_data, catRatiosNum_data, nans_data, nans_data_num = read_fits(input_data_file)
  exclude_dict = {}
  for exclude in [ False, True ]:
    catRatios_excl_dict, misIDRatios_excl_dict = {}, {}
    for is_data in [ False, True ]:
      name = "{}data{}".format("" if is_data else "pseudo", "_exclusions" if exclude else "")
      print("\n\nRunning for {}:".format(name))
      if exclude:
        assert(is_data in exclude_dict)
      exclude_bins_excl, exclude_bins_num_excl = merge_excludable_bins(
        #exclude_by_singificance = exclude_dict[is_data] if exclude else None,
        #exclude_by_singificance = exclude_genRec if exclude else None,
        exclude_by_singificance = (exclude_genRec if USE_GENREC_TO_EXCLUDE_CATEGORIES else exclude_dict[is_data]) if exclude else None,
        exclude_by_nan          = nans_data if is_data else nans_pseudo,
        exclude_by_request      = exclude_bins_additional
      )
      output_fit_excl = os.path.join(output_dir, "fit_result_{}.root".format(name))
      catRatios_excl, catRatiosNum_excl = exclude_cats(
        catRatiosNum_data if is_data else catRatiosNum_pseudo,
        catRatios_data    if is_data else catRatios_pseudo,
        exclude_bins_excl, exclude_bins_num_excl
      )
      if use_toys:
        calculate_solution_wToys(
          catRatios_excl, exclude_bins_num_excl, seed, nof_toys, output_fit_pseudo_excl, fallback_value,
          name if print_in_latex else ""
        )
      else:
        calculate_solution(
          catRatiosNum_excl, exclude_bins_num_excl, output_fit_excl, fallback_value, name if print_in_latex else ""
        )
      misIDRatios_excl = readMisIDRatios(output_fit_excl)
      sys.stdout.flush()
      chi2s_excl = make_pull_plot_21(
        misIDRatios_excl, catRatios_excl, name = name, output_dir = output_dir,
        y_range = (-0.001, 0.011), excluded = exclude_bins_num_excl
      )
      sys.stdout.flush()
      to_exclude = []
      print("chi2:")
      for bin, chi2 in chi2s_excl.items():
        exclude_by_significance = chi2 > NSIGMAS
        print("  {} {:.3f}{}".format(bin, chi2, " > {:.3f} => exclude".format(NSIGMAS) if exclude_by_significance else ""))
        if exclude_by_significance:
          to_exclude.append(bin)
      if not exclude:
        assert(is_data not in exclude_dict)
        exclude_dict[is_data] = to_exclude
      catRatios_excl_dict[is_data] = catRatios_excl
      misIDRatios_excl_dict[is_data] = misIDRatios_excl
    sys.stdout.flush()
    print("exclude_genRec: {}".format(exclude_genRec)) 
    print("exclude_dict: {}".format(exclude_dict))
    sys.stdout.flush()
    title = "With{} excluding the categories".format("" if exclude else "out")
    plot_ratio_filenames = [
      os.path.join(output_dir, "CompareRatios_w{}ExcludedBins.{}".format("" if exclude else "o", ext)) \
      for ext in [ "png", "pdf" ]
    ]
    plot_ratios(catRatios_excl_dict[True], catRatios_excl_dict[False], plot_ratio_filenames, title)
    plot_rate_filenames = [
      os.path.join(output_dir, "CompareRates_w{}ExcludedBins.{}".format("" if exclude else "o", ext)) \
      for ext in [ "png", "pdf" ]
    ]
    plot_rates(misIDRatios_excl_dict[True], misIDRatios_excl_dict[False], misIDRatios_gen, plot_rate_filenames, title)
  print('=' * 120 + '\n')
