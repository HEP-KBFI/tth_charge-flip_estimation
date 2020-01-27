#!/usr/bin/env python
# -*- coding: utf-8 -*-

from tthAnalysis.ChargeFlipEstimation.utils import mkdir_p, SmartFormatter

import ROOT

import subprocess
import os
import math
import argparse
import multiprocessing
import sys
import signal

ROOT.gROOT.SetBatch(True)

"""@file docstring
Script for running electron charge flip estimation fit

@author Andres Tiko <andres.tiko@cern.ch>
@author Karl Ehatäht <karl.ehataht@cern.ch>
"""

fit_results = {}
fit_results_failed = []

OBSERVATION_STR = 'observation'

COMBINE_SETTINGS = [
  "",
  "--cminDefaultMinimizerStrategy 0",
  "--cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerType Minuit",
  "--cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerType Minuit --rMin -4.0 --rMax 20.0",
  "--cminDefaultMinimizerStrategy 0 --rMin -4.0 --rMax 20.0 --cminDefaultMinimizerTolerance 100",
  "--cminDefaultMinimizerStrategy 0 --rMin -4.0 --rMax 20.0 --cminDefaultMinimizerTolerance 1000",
  "--cminDefaultMinimizerStrategy 0 --rMin -4.0 --rMax 20.0 --cminDefaultMinimizerTolerance 2000",
  "--cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerType Minuit --rMin -4.0 --rMax 20.0 --cminDefaultMinimizerTolerance 100",
  "--cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerType Minuit --rMin -4.0 --rMax 20.0 --cminDefaultMinimizerTolerance 1000",
  "--cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerType Minuit --rMin -4.0 --rMax 20.0 --cminDefaultMinimizerTolerance 2000",
]

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
          found_observed = int(float(line_split[1]))
          break
    if found_observed < 0:
      raise RuntimeError(
        "Unable to line that starts with '%s' from file %s" % (OBSERVATION_STR, input_datacard_filename)
      )
    elif found_observed == 0 and bin not in skip_bins:
      print("Skipping bin {} because zero observations found".format(bin))
      skip_bins.append(bin)

  return skip_bins

def failed_result(bin, status):
  return (int(bin), float('nan'), float('nan'), float('nan'), status)

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
    return None

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
      return None

  if abs(muHiErr / muLoErr) > 2 or abs(muHiErr / muLoErr) < 1. / 2:
    # In this case probably failed to find crossing, use symmetric error
    print("Strange asymmetric errors! Probably failed to find crossing.")
    return None

  fitHiErr = muHiErr / mu * bestFit
  fitLoErr = muLoErr / mu * bestFit
  print("r: {},  fitHiErr: {}, fitLoErr: {}".format(bestFit, fitHiErr, fitLoErr))

  # Use Poisson mean of getting 0 events for uncertainty if no observed events
  if fitHiErr == 0.:
    lambda_poisson0 = -math.log(0.32)
    fitHiErr = lambda_poisson0 / (lambda_poisson0 + fail_h.Integral())
    fitLoErr = fitHiErr
    print("Changed >>> r: {},  fitHiErr: {}, fitLoErr: {}".format(bestFit, fitHiErr, fitLoErr))

  if bestFit < 0.:
    print("WARNING: got negative fit result => considering it as invalid")
    return None

  return (int(bin), bestFit, fitHiErr, fitLoErr, 0)

def update_fit_results(results):
  bin, bestFit, fitHiErr, fitLoErr, status = results
  if status > 0:
    fit_results_failed.append(bin)
  fit_results[bin] = (bestFit, fitHiErr, fitLoErr)

def fit_bin(input_dir, bin, skip_bins):
  print("{} bin {}: ------------------------- ".format(input_dir, bin))

  datacard_dir = os.path.join(input_dir, "cards")
  current_card_base = "card_{}.txt".format(bin)
  current_card = os.path.join(datacard_dir, current_card_base)
  current_workspace = os.path.join(datacard_dir, "workspace_{}.root".format(bin))

  command_combineCards = "combineCards.py " \
    "pass=SScards/htt_SS_{bin}_13TeV.txt " \
    "fail=OScards/htt_OS_{bin}_13TeV.txt > {current_card}".format(
      datacard_dir = datacard_dir,
      bin          = bin,
      current_card = current_card_base,
  )

  # 1. step: combine SS and OS datacatds
  if skip_bins:
    return failed_result(bin, 1)
  else:
    print("Running: {}".format(command_combineCards))
    subprocess.call(command_combineCards, shell = True, cwd = datacard_dir)

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
  combine_out = os.path.join(datacard_dir, "combine_out_{}.log".format(bin))
  fit_file = os.path.join(fit_bin_dir, "fitDiagnostics.root")
  postfit_file = os.path.join(fit_bin_dir, "output_postfit.root")
  fit_out = os.path.join(fit_bin_dir, "out.log")

  fit_status = ()
  common_settings = "--robustFit 1 "
  for settings in COMBINE_SETTINGS:
    specific_settings = common_settings + settings
    commandCombine = "combine -v 0 -M FitDiagnostics {} --out {} --plots --saveNormalizations --skipBOnlyFit " \
      "--saveShapes --saveWithUncertainties --maxFailedSteps 20 {} &> {}".format(
        current_workspace, fit_bin_dir, specific_settings, combine_out
    )
    print("Running: {}".format(commandCombine))
    subprocess.call(commandCombine, shell = True, cwd = fit_bin_dir)

    # 4. Create postfit plots
    command_postFit = "PostFitShapesFromWorkspace -d {} -w {} -o {} -f {}:fit_s --postfit --sampling &> {}".format(
      current_card, current_workspace, postfit_file, fit_file, fit_out
    )
    print("Running: {}".format(command_postFit))
    subprocess.call(command_postFit, shell = True)

    # 5. Add to list of fit results
    if os.path.exists(fit_file) and os.path.exists(postfit_file):
      fit_status = read_fit_result(fit_file, postfit_file, bin = bin)

    print("fit_status: {}".format(fit_status))
    if fit_status:
      break

  if fit_status:
    return fit_status
  else:
    print("{} bin {} fit did not converge: {}".format(input_dir, bin, fit_status))
    return failed_result(bin, 2)

def make_fits(input_dir, output_file, whitelist = None, skip_bins = None, jobs = -1, force = False):
  assert(os.path.isdir(input_dir))
  output_dir = os.path.dirname(output_file)
  if not os.path.isdir(output_dir):
    if os.path.isfile(output_dir):
      raise ValueError("Cannot create directory %d because it's a file" % output_dir)
    if not force:
      raise ValueError("Use -f/--force to create output directory %d" % output_dir)
    os.makedirs(output_dir)

  bins = list(range(21)) if not whitelist else whitelist
  original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
  fit_pool = multiprocessing.Pool(jobs if jobs > 0 else (multiprocessing.cpu_count() - 1))
  signal.signal(signal.SIGINT, original_sigint_handler)
  for bin in bins:
    try:
      fit_pool.apply_async(
        fit_bin,
        args = (input_dir, bin, skip_bins and bin in skip_bins),
        callback = update_fit_results,
      )
    except KeyboardInterrupt:
      fit_pool.terminate()
      sys.exit(1)

  fit_pool.close()
  fit_pool.join()

  # Output fit results
  with open(output_file, "w") as f:
    for bin in sorted(fit_results.keys()):
      fr = fit_results[bin]
      print("Bin = %d, r = %.8f + %.8f - %.8f" % (bin, fr[0], fr[1], fr[2]))
      f.write("%d, %.8f, %.8f, %.8f\n" % (bin, fr[0], fr[1], fr[2]))

  if fit_results_failed:
    output_file_failed = '{}.failed'.format(output_file)
    with open(output_file_failed, "w") as fFailedFits:
      for bin in fit_results_failed:
        fFailedFits.write("{}\n".format(bin))

if __name__ == "__main__":
  parser = argparse.ArgumentParser(
    formatter_class = lambda prog: SmartFormatter(prog, max_help_position = 35)
  )
  parser.add_argument('-i', '--input',
    type = str, dest = 'input_data', metavar = 'directory', required = True,
    help = 'R|Input directory',
  )
  parser.add_argument('-o', '--output',
    type = str, dest = 'output', metavar = 'file', required = True,
    help = 'R|Output file',
  )
  parser.add_argument('-s', '--skip',
    type = int, dest = 'skip', metavar = 'bin', required = False, nargs = '+', default = [],
    help = 'R|Skip bins',
  )
  parser.add_argument('-j', '--jobs',
    type = int, dest = 'jobs', metavar = 'bin', required = False, default = 16,
    help = 'R|Number of parallel fits',
  )
  parser.add_argument('-w', '--whitelist',
    type = int, dest = 'whitelist', metavar = 'bin', required = False, nargs = '+', default = [],
    help = 'R|Whitelist bins (default: all)',
  )
  parser.add_argument('-S', '--skip-automatically',
    type = bool, dest = 'skip_automatically', required = False, default = True,
    help = 'R|Skip bins that have no data and pseudodata',
  )
  parser.add_argument('-f', '--force',
    dest = 'force', action = 'store_true', default = False, required = False,
    help = 'R|Create directory for the output file',
  )
  args = parser.parse_args()

  skip_bins = args.skip
  if args.skip_automatically:
    skip_bins = find_zero_bins(args.input_data, skip_bins)

  print("Input:              {}".format(args.input_data))
  print("Output:             {}".format(args.output))
  print("Skipping bins:      {}".format(', '.join(map(str, skip_bins))))
  print("Whitelisting bins:  {}".format(', '.join(map(str, args.whitelist))))
  print("Skipping zero bins: {}".format(args.skip_automatically))
  print("# parallel fits:    {}".format(args.jobs))

  make_fits(
    input_dir   = os.path.abspath(args.input_data),
    output_file = os.path.abspath(args.output),
    whitelist   = args.whitelist,
    skip_bins   = skip_bins,
    jobs        = args.jobs,
    force       = args.force,
  )
