#!/usr/bin/env python
# -*- coding: utf-8 -*-

from tthAnalysis.ChargeFlipEstimation.utils import mkdir_p, SmartFormatter

import ROOT

import subprocess
import os
import math
import argparse
import shutil

ROOT.gROOT.SetBatch(True)

"""@file docstring
Script for running electron charge flip estimation fit

@author Andres Tiko <andres.tiko@cern.ch>
@author Karl Ehatäht <karl.ehataht@cern.ch>
"""

#TODO parallelize fits

OBSERVATION_STR = 'observation'

COMBINE_SETTINGS = {
  'data' : {
    'electron' : {
      '2016' : {
        13 : "--cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerType Minuit",
        17 : "--cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerType Minuit",
      },
      '2017' : {
        3  : "--cminDefaultMinimizerStrategy 0 --rMin -4.0 --rMax 20.0 --cminDefaultMinimizerTolerance 2000",
        4  : "--cminDefaultMinimizerStrategy 0 --rMin -4.0 --rMax 20.0 --cminDefaultMinimizerTolerance 1000",
        5  : "--cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerType Minuit --rMin -4.0 --rMax 20.0 --cminDefaultMinimizerTolerance 2000",
        8  : "--cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerType Minuit --rMin -4.0 --rMax 20.0 --cminDefaultMinimizerTolerance 1000",
        15 : "--cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerType Minuit",
      },
      '2018' : {
        2  : "--cminDefaultMinimizerType Minuit --cminDefaultMinimizerStrategy 0  --cminDefaultMinimizerTolerance 25 --rMin -4.0 --rMax 20.0",
        5  : "--cminDefaultMinimizerStrategy 0 --rMin -4.0 --rMax 20.0 --cminDefaultMinimizerTolerance 25",
        11 : "--cminDefaultMinimizerStrategy 0 --rMin -4.0 --rMax 20.0 --cminDefaultMinimizerTolerance 1000",
        16 : "--cminDefaultMinimizerType Minuit --cminDefaultMinimizerStrategy 0 --rMin -4.0 --rMax 20.0",
      },
    },
  },
  'pseudodata' : {
    'electron' : {
      '2016' : {
        2  : "--cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerType Minuit",
      },
    },
  },
}

def find_zero_bins(input_datacard, skip_bins):
  for bin in range(21):
    input_datacard_filename = os.path.join(input_datacard, "cards", "SScards", "htt_SS_{}_13TeV.txt".format(bin))
    if not os.path.isfile(input_datacard_filename):
      raise ValueError("No such file: %s" % input_datacard_filename)

    found_observed = -1
    with open(input_datacard_filename, 'r') as input_datacard_fptr:
      for line in input_datacard_fptr:
        line_stripped = line.rstrip()
        if line_stripped.startswith(OBSERVATION_STR):
          line_split = line_stripped.split()
          assert(len(line_split) == 2)
          found_observed = int(line_split[1])
          break
    if found_observed < 0:
      raise RuntimeError(
        "Unable to line that starts with '%s' from file %s" % (OBSERVATION_STR, input_datacard_filename)
      )
    elif found_observed == 0 and bin not in skip_bins:
      print("Skipping bin {} because zero observations found")
      skip_bins.append(skip_bins)

  return skip_bins

def read_fit_result(fit_file, postfit_file, bin):
  postfit_file_ptr = ROOT.TFile(postfit_file, 'read')
  fail_h = postfit_file_ptr.Get("fail_prefit/DY")
  pass_h = postfit_file_ptr.Get("pass_prefit/DY")

  fit_file_ptr = ROOT.TFile(fit_file, 'read')
  tree = fit_file_ptr.Get("tree_fit_sb")

  tree.Draw("fit_status>>histStatus")
  fit_status = ROOT.gDirectory.Get("histStatus").GetMean()
  if fit_status != 0:
    # raise AssertionError("Input fit %d did not converge! Fit status %d." % (bin, fit_status))
    print("Input fit {} did not converge! Fit status {}".format(bin, fit_status))
    return ()

  tree.Draw("SF>>hist")
  mu = (ROOT.gDirectory.Get("hist")).GetMean()

  tree.Draw("SFErr>>histErr")
  muErr = (ROOT.gDirectory.Get("histErr")).GetMean()

  tree.Draw("SFLoErr>>histLo")
  muLoErr = (ROOT.gDirectory.Get("histLo")).GetMean()

  tree.Draw("SFHiErr>>histHi")
  muHiErr = (ROOT.gDirectory.Get("histHi")).GetMean()

  print("pass: {},  fail: {},  mu: {}, muErr: {}, muLoErr: {}, muHiErr: {} ".format(
    pass_h.Integral(), fail_h.Integral(), mu, muErr, muLoErr, muHiErr, muHiErr / muLoErr)
  )
  if muLoErr > 0.:
    print("\t\t muHiErr/muLoErr: {}".format(muHiErr / muLoErr))

  try:
    # postFit distributions are scaled to scale factor 1, need to multiply by fitted number
    bestFit = mu * pass_h.Integral() / (fail_h.Integral() + pass_h.Integral())
  except AttributeError:
    if fail_h.Integral() > 0:
      bestFit = 0.
    else:
      # raise
      print("fail_histo integral = 0")
      return ()

  if abs(muHiErr / muLoErr) > 2 or abs(muHiErr / muLoErr) < 1. / 2:
    # In this case probably failed to find crossing, use symmetric error
    print("Strange asymmetric errors! Probably failed to find crossing.")
    return ()

  fitHiErr = muHiErr / mu * bestFit
  fitLoErr = muLoErr / mu * bestFit
  print("r: {},  fitHiErr: {}, fitLoErr: {}".format(bestFit, fitHiErr, fitLoErr))

  # Use Poisson mean of getting 0 events for uncertainty if no observed events
  if fitHiErr == 0.:
    lambda_poisson0 = -math.log(0.32)
    fitHiErr = lambda_poisson0 / (lambda_poisson0 + fail_h.Integral())
    fitLoErr = fitHiErr
    print("\t\t changed >>> r: {},  fitHiErr: {}, fitLoErr: {}".format(bestFit, fitHiErr, fitLoErr))

  return (int(bin), bestFit, fitHiErr, fitLoErr)

def failed_result(bin):
  return (int(bin), float('nan'), float('nan'), float('nan'))

def make_fits(input_dir, data_type, lepton_type, era, skip_bins = None):
  assert(os.path.isdir(input_dir))

  datacard_dir = os.path.join(input_dir, "cards")
  workspace_dir = os.path.join(input_dir, "workspaces")
  mkdir_p(datacard_dir)
  mkdir_p(workspace_dir)

  fit_results = []
  fit_results_failed = []

  fit_dir = os.path.join(input_dir, "fit")
  mkdir_p(fit_dir)

  for bin in range(21):
    print("\n\n{} bin {}: ------------------------- \n".format(input_dir, bin))

    current_card_base = "card_{}.txt".format(bin)
    current_card = os.path.join(datacard_dir, current_card_base)
    current_workspace = os.path.join(datacard_dir, "workspace_{}.root".format(bin))

    command_combineCards = "combineCards.py " \
      "pass={datacard_dir}/SScards/htt_SS_{bin}_13TeV.txt " \
      "fail={datacard_dir}/OScards/htt_OS_{bin}_13TeV.txt > {current_card}".format(
        datacard_dir = datacard_dir,
        bin          = bin,
        current_card = current_card,
    )
    print("Running: {}".format(command_combineCards))

    # 1. step: combine SS and OS datacatds
    if skip_bins and bin in skip_bins:
      fit_results.append(failed_result(bin))
      continue
    else:
      subprocess.call(command_combineCards, shell = True)

    # Hack to prevent PostFitShapesFromWorkspace messing up directory structure - copy datacard to current directory
    current_card_backup = '{}.{}'.format(current_card, 'bak')
    shutil.copyfile(current_card, current_card_backup)

    # 2. Make Roofit workspace from datacard
    command_text2ws = "text2workspace.py {} -o {} -P HiggsAnalysis.CombinedLimit.TagAndProbeModel:tagAndProbe".format(
      current_card, current_workspace
    )
    print("Running: {}".format(command_text2ws))
    subprocess.call(command_text2ws, shell = True)

    # Specify output directory for fit results
    fit_bin_dir = os.path.join(input_dir, "bin{}".format(bin))
    mkdir_p(fit_bin_dir)

    # 3. Perform fit with combine
    # Default fit settings specified in else clause
    # But always do not give convergence - settings for specific fits defined here, adjust as necessary:
    specific_settings = "--robustFit 1 "
    if data_type    in COMBINE_SETTINGS and \
        lepton_type in COMBINE_SETTINGS[data_type] and \
        era         in COMBINE_SETTINGS[data_type][lepton_type] and \
        bin         in COMBINE_SETTINGS[data_type][lepton_type][era]:
      specific_settings += COMBINE_SETTINGS[data_type][lepton_type][era][bin]

    commandCombine = "combine -v 0 -M FitDiagnostics {} --out {} --plots --saveNormalizations --skipBOnlyFit --saveShapes " \
      "--saveWithUncertainties --maxFailedSteps 20 {}".format(current_workspace, fit_bin_dir, specific_settings)
    print("Running: {}".format(commandCombine))
    subprocess.call(commandCombine, shell = True)

    # 4. Create postfit plots
    fit_file = os.path.join(fit_bin_dir, "fitDiagnostics.root")
    postfit_file = os.path.join(fit_bin_dir, "output_postfit.root")
    command_postFit = "PostFitShapesFromWorkspace -d {} -w {} -o {} -f {}:fit_s " \
                      "--postfit --sampling".format(current_card_backup, current_workspace, postfit_file, fit_file)
    print("Running: {}".format(command_postFit))
    subprocess.call(command_postFit, shell = True)

    # 5. Add to list of fit results
    fit_status = ()
    if os.path.exists(fit_file) and os.path.exists(postfit_file):
      fit_status = read_fit_result(fit_file, postfit_file, bin = bin)

    print("fit_status: {}".format(fit_status))
    if fit_status:
      fit_results.append(fit_status)
      print(">>> fitResults appended")
    else:
      print("{} bin {} fit did not converge: {}".format(input_dir, bin, fit_status))
      fit_results_failed.append(bin)

  if fit_results_failed:
    with open(os.path.join(fit_dir, "failedFits.txt"), "w") as fFailedFits:
      for bin in fit_results_failed:
        fFailedFits.write("{} \t {} \n".format(input_dir, bin))

  # Output fit results
  with open(os.path.join(fit_dir, "results_cat.txt"), "w") as f:
    for i, fr in enumerate(fit_results):
      print("RES: %d   %d, %.8f + %.8f - %.8f" % (i, fr[0], fr[1], fr[2], fr[3]))
      f.write("%d, %.8f, %.8f, %.8f\n" % (fr[0], fr[1], fr[2], fr[3]))

if __name__ == "__main__":
  parser = argparse.ArgumentParser(
    formatter_class = lambda prog: SmartFormatter(prog, max_help_position = 35)
  )
  parser.add_argument('-i', '--input',
    type = str, dest = 'input_data', metavar = 'directory', required = True,
    help = 'R|Input directory',
  )
  parser.add_argument('-l', '--lepton-type',
    type = str, dest = 'lepton_type', metavar = 'type', required = True, choices = [ 'electron', 'muon' ],
    help = 'R|Lepton type',
  )
  parser.add_argument('-t', '--data-type',
    type = str, dest = 'data_type', metavar = 'type', required = True, choices = [ 'data', 'pseudodata' ],
    help = 'R|Datacard type',
  )
  parser.add_argument('-e', '--era',
    type = int, dest = 'era', metavar = 'year', required = True, choices = [ 2016, 2017, 2018 ],
    help = 'R|Era',
  )
  parser.add_argument('-s', '--skip',
    type = int, dest = 'skip', metavar = 'bin', required = False, nargs = '+', default = [],
    help = 'R|Skip bins',
  )
  parser.add_argument('-S', '--skip-automatically',
    type = bool, dest = 'skip_automatically', required = False, default = True,
    help = 'R|Skip bins that have no data and pseudodata',
  )
  args = parser.parse_args()

  skip_bins = args.skip
  if args.skip_automatically:
    skip_bins = find_zero_bins(args.input_data, skip_bins)

  print("Input:              {}".format(args.input_data))
  print("Data type:          {}".format(args.data_type))
  print("Lepton type:        {}".format(args.lepton_type))
  print("Era:                {}".format(args.era))
  print("Output:             {}".format(args.output))
  print("Skipping bins:      {}".format(args.skip))
  print("Skipping zero bins: {}".format(args.skip_automatically))

  make_fits(
    input_dir   = args.input_data,
    data_type   = args.data_type,
    lepton_type = args.lepton_type,
    era         = args.era,
    skip_bins   = skip_bins,
  )
