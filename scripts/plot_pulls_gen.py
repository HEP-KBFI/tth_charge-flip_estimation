#!/usr/bin/env python

from tthAnalysis.ChargeFlipEstimation.utils import bin_names_composite, bin_names_single, mkdir_p
from tthAnalysis.ChargeFlipEstimation.plot_pulls import get_other_component, get_all_containers, readMisIDRatiosGen, \
                                                        readCategoryRatiosGen, make_pull_plot_21

import ROOT

"""@file docstring
Script for plotting generator-level pulls, main program superseded by plot_pulls_all.py,
but contains most of the necessary methods

@author Andres Tiko <andres.tiko@cern.ch>
"""

def get_bin_nr_composite(cat):
  return bin_names_composite.index(cat)

def get_bin_nr_single(cat):
  return bin_names_single.index(cat)

def make_pull_plot_gen(category, misIDRatios, catRatios):
  pull_plot = ROOT.TH1D(category, category, 6, 0, 6 )
  others_plot = ROOT.TH1D(category+"others", category+"others", 6, 0, 6 )
  true_plot = ROOT.TH1D(category+"true", category+"true", 6, 0, 6 )
  bin_names = get_all_containers(category)
  for b in range(1, len(bin_names)+1):
    pull_plot.GetXaxis().SetBinLabel(b,bin_names[b-1])
    others_plot.GetXaxis().SetBinLabel(b,bin_names[b-1])
    true_plot.GetXaxis().SetBinLabel(b,bin_names[b-1])
    (value, err, err_plus) = catRatios[bin_names[b-1]]
    pull_plot.SetBinContent(b, value)
    pull_plot.SetBinError(b, err)

    other = get_other_component(category, bin_names[b-1])
    (valueO, errO, err0_plus) = misIDRatios[other]
    others_plot.SetBinContent(b, valueO)
    others_plot.SetBinError(b, errO)

    true_plot.SetBinContent(b, misIDRatios[category][0])
  pull_plot.Add(others_plot, -1)
  c = ROOT.TCanvas("Plot", "Plot", 1920,1080)
  ROOT.gStyle.SetOptStat(0)
  true_plot.SetLineColor(ROOT.kRed)
  true_plot.SetLineWidth(3)
  true_plot.GetYaxis().SetRangeUser(-0.006, 0.006)
  true_plot.Draw()
  pull_plot.SetLineWidth(3)
  pull_plot.Draw("SAME")
  mydir = "pull_plots_gen/"
  mkdir_p(mydir)
  c.SaveAs("%s/%s_pulls.pdf" % (mydir, category))
  c.SaveAs("%s/%s_pulls.png" % (mydir, category))

if __name__ == "__main__":
  infile = "/hdfs/local/ttH_2tau/andres/ttHAnalysis/2016/histosCF_summer_June6/histograms/charge_flip/histograms_harvested_stage2_charge_flip_Tight.root"  
  misIDRatios = readMisIDRatiosGen(infile)
  catRatios = readCategoryRatiosGen(infile)
  make_pull_plot_21(misIDRatios, catRatios)
