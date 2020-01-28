#!/usr/bin/env python
# -*- coding: utf-8 -*-

from tthAnalysis.ChargeFlipEstimation.utils import BIN_NAMES_COMPOSITE, SmartFormatter, CHARGES

import ROOT
import numpy as np

ROOT.gROOT.SetBatch(True)

import argparse

"""@file docstring
Script for creating pseudodata and datacards, which are suitable as input to charge flip estimation

#See: <https://en.wikipedia.org/wiki/Least_squares#Weighted_least_squares>
for notation and background

@author Andres Tiko <andres.tiko@cern.ch>
@author Karl Ehatäht <karl.ehataht@cern.ch>
"""

SAMPLES = [
  "data_obs",
  "DY",
  "DY_fake",
  "WJets",
  "Singletop",
  "Diboson",
  "TTbar"
]

CATEGORIES = {
  "electron" : BIN_NAMES_COMPOSITE, # 21 categories of lepton pairs by pT and eta
  "muon"     : [ "total" ],         # For muons all in the same cateogry
}
PREFIX = 'ttH_charge_flip'

# Systematic uncertainties to add to datacard
SYSTEMATICS = {
  "electron" : [
    "CMS_ttHl_electronERUp",
    "CMS_ttHl_electronERDown",
    "CMS_ttHl_electronESEndcapUp",
    "CMS_ttHl_electronESEndcapDown",
    "CMS_ttHl_electronESBarrelUp",
    "CMS_ttHl_electronESBarrelDown",
  ],
  "muon" : [
    "CMS_ttHl_muonERUp",
    "CMS_ttHl_muonERDown",
    "CMS_ttHl_muonESEndcap1Up",
    "CMS_ttHl_muonESEndcap1Down",
    "CMS_ttHl_muonESEndcap2Up",
    "CMS_ttHl_muonESEndcap2Down",
    "CMS_ttHl_muonESBarrel1Up",
    "CMS_ttHl_muonESBarrel1Down",
    "CMS_ttHl_muonESBarrel2Up",
    "CMS_ttHl_muonESBarrel2Down",
  ],
}

def reset_bin_content(histogram, rebin, manipulate):
  histogram.Rebin(rebin)
  for bin_idx in range(1, histogram.GetNbinsX() + 1):
    content = histogram.GetBinContent(bin_idx)
    content_error = histogram.GetBinError(bin_idx)
    # If bin less than zero, set to zero with uncertainty
    if content < 0:
      histogram.SetBinContent(bin_idx, 0.)
    if content > 0. and abs(content_error) < 1e-10:
      histogram.SetBinError(bin_idx, np.sqrt(content))
      print("  Histogram: {}, bin: {}, content: {} +- {}, now {} +- {}".format(
        histogram.GetName(), bin_idx, content, content_error, histogram.GetBinContent(bin_idx),
        histogram.GetBinError(bin_idx)
      ))
  for bin_idx in [0, histogram.GetNbinsX() + 1]:
    if histogram.GetBinContent(bin_idx) > 0.:
      histogram.SetBinContent(bin_idx, 0.)
      histogram.SetBinError(bin_idx, 0.)
      print("  Histogram: {} bin {}: {} +- {}  under-/over-flow set to zero".format(
        histogram.GetName(), bin_idx, histogram.GetBinContent(bin_idx), histogram.GetBinError(bin_idx)
      ))
  if manipulate:
    print("Manipulating: {}".format(histogram.GetName()))

    eIntegral = ROOT.Double(0.0)
    integral = histogram.IntegralAndError(1, histogram.GetNbinsX(), eIntegral)
    print("Integral before: {} +- {}".format(integral, eIntegral))

    histogram.Scale(0.5 / histogram.Integral())
    for bin_idx in range(1, histogram.GetNbinsX() + 1):
      histogram.SetBinError(bin_idx, np.sqrt(max(histogram.GetBinContent(bin_idx), 0.)))
    for bin_idx in [ 0, histogram.GetNbinsX() + 1 ]:
      histogram.SetBinContent(bin_idx, 0)
      histogram.SetBinError(bin_idx, 0)

    integral = histogram.IntegralAndError(1, histogram.GetNbinsX(), eIntegral)
    print("Integral after: {} +- {}".format(integral, eIntegral))

def create_pseudodata(
      input_file_name,
      output_file_name,
      output_pseudo_file_name,
      channel,
      rebin      = 1,
      manipulate = None,
      samples    = SAMPLES,
      prefix     = PREFIX
    ):
  input_file = ROOT.TFile(input_file_name, "read")
  output_file = ROOT.TFile(output_file_name, "recreate")
  output_pseudo_file = ROOT.TFile(output_pseudo_file_name, "recreate")

  print("Input: {}".format(input_file_name))
  print("Output: {}, {}".format(output_file_name, output_pseudo_file_name))

  missing_histograms = []
  for charge in CHARGES:
    for category in CATEGORIES[channel]:
      is_first_histogram = True

      dirname_top = "{}_{}_{}".format(prefix, charge, category)
      dir_top = output_file.mkdir(dirname_top)
      dir_pseudo_top = output_pseudo_file.mkdir(dirname_top)
      dirname_rebinned = "rebinned"
      dirname_rebinned_full = "{}/{}".format(dirname_top, dirname_rebinned)
      dir_rebinned = dir_top.mkdir(dirname_rebinned)
      dir_pseudo_rebinned = dir_pseudo_top.mkdir(dirname_rebinned)

      for sample in samples:
        manipulation_enabled = ':'.join([ charge, category, sample ]) in manipulate

        # Read nominal histograms
        histo_name_nominal = "{}/{}_rebinned".format(dirname_rebinned_full, sample)
        print("  Histogram: {}".format(histo_name_nominal))
        histo_nominal = input_file.Get(histo_name_nominal)
        if not histo_nominal:
            print("\033[91m Couldn't fetch {} \033[0m".format(histo_name_nominal))
            missing_histograms.append(histo_name_nominal)
            continue

        reset_bin_content(histo_nominal, rebin = rebin, manipulate = manipulation_enabled)

        dir_rebinned.cd()
        histo_nominal.Write()

        # Don't add data to pseudodata
        if sample == "data_obs":
          continue
        dir_pseudo_rebinned.cd()
        histo_nominal.Write()

        # Add different MCs together as pseudodata
        if is_first_histogram == True:
          data_histo = histo_nominal.Clone()
          is_first_histogram = False
        else:
          data_histo.Add(histo_nominal)

        # Add systematics
        for syst in SYSTEMATICS[channel]:
          if syst.startswith(("CMS_ttHl_electronER", "CMS_ttHl_muonER")) and not sample == "DY": continue

          histo_name_sys = "{}/{}_{}_rebinned".format(dirname_rebinned_full, sample, syst)
          print("  Histogram: {}".format(histo_name_sys))
          histo_sys = input_file.Get(histo_name_sys)
          if not histo_sys:
            print("\033[91m Couldn't fetch {} \033[0m".format(histo_name_sys))
            missing_histograms.append(histo_name_sys)
            continue

          reset_bin_content(histo_sys, rebin = rebin, manipulate = manipulation_enabled)

          dir_rebinned.cd()
          histo_sys.Write()
          dir_pseudo_rebinned.cd()
          histo_sys.Write()

      data_histo.SetNameTitle("data_obs_rebinned", "data_obs_rebinned")

      # Generate poisson yields from MC expectation for pseudodata
      for bin_idx in range(1, data_histo.GetNbinsX() + 1):
        binCount1 = data_histo.GetBinContent(bin_idx)
        binCount2 = np.random.poisson(data_histo.GetBinContent(bin_idx))
        if abs(binCount2 - 0.) < 0.000005:
          while abs(binCount2 - 0.) < 0.000005:
            binCount2 = np.random.poisson(binCount1)
        data_histo.SetBinContent(bin_idx, binCount2)
        data_histo.Sumw2(ROOT.kFALSE)
      dir_pseudo_rebinned.cd()
      data_histo.Write()

  input_file.Close()
  output_file.Close()
  output_pseudo_file.Close()

  if missing_histograms:
    print("\033[91mList of Histograms could not read:\n{}\033[0m".format(",\n".join(missing_histograms)))

if __name__ == "__main__":
  parser = argparse.ArgumentParser(
    formatter_class = lambda prog: SmartFormatter(prog, max_help_position = 35)
  )
  parser.add_argument('-i', '--input',
    type = str, dest = 'input', metavar = 'file', required = True,
    help = 'R|Input datacard',
  )
  parser.add_argument('-o', '--output-data',
    type = str, dest = 'output_data', metavar = 'file', required = True,
    help = 'R|Output file name for data',
  )
  parser.add_argument('-O', '--output-pseudodata',
    type = str, dest = 'output_pseudodata', metavar = 'file', required = True,
    help = 'R|Output file name for data',
  )
  parser.add_argument('-t', '--lepton-type',
    type = str, dest = 'lepton_type', metavar = 'type', required = True, choices = CATEGORIES.keys(),
    help = 'R|Lepton type',
  )
  parser.add_argument('-p', '--prefix',
    type = str, dest = 'prefix', metavar = 'string', required = False, default = PREFIX,
    help = 'R|Prefix of input histogram names',
  )
  parser.add_argument('-r', '--rebin',
    type = str, dest = 'rebin', metavar = 'int', required = False, default = 1,
    help = 'R|Rebin',
  )
  parser.add_argument('-s', '--samples',
    type = str, dest = 'samples', metavar = 'sample', required = False, default = SAMPLES, choices = SAMPLES, nargs = '+',
    help = 'R|List of samples to consider',
  )
  parser.add_argument('-S', '--seed',
    type = int, dest = 'seed', metavar = 'int', required = False, default = 123,
    help = 'R|Seed for PRNG',
  )
  parser.add_argument('-m', '--manipulate',
    type = str, dest = 'manipulate', metavar = 'string', required = False, nargs = '+', default = [],
    help = 'R|Set bin content to 0.5 (format: charge:category:sample, eg SS:BB_LL:DY)',
  )
  args = parser.parse_args()
  np.random.seed(args.seed)

  create_pseudodata(
    input_file_name         = args.input,
    output_file_name        = args.output_data,
    output_pseudo_file_name = args.output_pseudodata,
    channel                 = args.lepton_type,
    rebin                   = args.rebin,
    manipulate              = args.manipulate,
    samples                 = args.samples,
    prefix                  = args.prefix,
  )
