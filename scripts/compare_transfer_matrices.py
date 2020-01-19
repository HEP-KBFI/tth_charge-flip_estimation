#!/usr/bin/env python

from tthAnalysis.ChargeFlipEstimation.utils import mkdir_p, bin_names_single

import ROOT

"""@file docstring
Plot transfer matrices to illustrate differences in SS and OS in reconstructing pT and eta

@author Andres Tiko <andres.tiko@cern.ch>
"""

def column_sum(histo, column):
  col_sum = 0.
  for row in range(0, histo.GetNbinsY()+2):
    col_sum += histo.GetBinContent(column, row)
  return col_sum


def column_scaler(histo):
  for col in range(0, histo.GetNbinsX()+2):
    col_sum = column_sum(histo, col)
    if col_sum == 0: continue
    for row in range(0, histo.GetNbinsY()+2):
      histo.SetBinContent(col, row, histo.GetBinContent(col, row) / col_sum)
      histo.SetBinError(col, row, histo.GetBinError(col, row) / col_sum)

def row_sum(histo, row):
  row_sum = 0.
  for col in range(0, histo.GetNbinsX()+2):
    row_sum += histo.GetBinContent(col, row)
  return row_sum

def row_scaler(histo):
  for row in range(0, histo.GetNbinsY()+2):
    r_sum = row_sum(histo, row)
    if r_sum == 0: continue
    for col in range(0, histo.GetNbinsX()+2):
      histo.SetBinContent(col, row, histo.GetBinContent(col, row) / r_sum)
      histo.SetBinError(col, row, histo.GetBinError(col, row) / r_sum)

def set_axes(histo):
  for i in range(1,7):
    histo.GetXaxis().SetBinLabel(i,bin_names_single[i-1])
    histo.GetYaxis().SetBinLabel(i,bin_names_single[i-1])

def plot_transfer_matrix(histo, name, title):
  c = ROOT.TCanvas("Plot", "Plot", 800,800)
  ROOT.gStyle.SetOptStat(0)
  ROOT.gStyle.SetPaintTextFormat("1.3f");
  set_axes(histo)
  histo.GetXaxis().SetTitle("gen")
  histo.GetYaxis().SetTitle("rec")
  histo.SetNameTitle(name, title)
  histo.Draw("colz")
  histo.Draw("text same")
  mydir = "../plots/transfer_matrices/"
  mkdir_p(mydir)
  c.SaveAs("%s/transfer_matrix_%s.pdf" % (mydir, name))
  c.SaveAs("%s/transfer_matrix_%s.png" % (mydir, name))


if __name__ == "__main__":
  ROOT.gROOT.SetBatch(True)
  infile = "/hdfs/local/ttH_2tau/andres/ttHAnalysis/2016/histosCF_summer_June6/histograms/charge_flip/histograms_harvested_stage2_charge_flip_Tight.root"
  f = ROOT.TFile(infile)
  transfer_matrix_gen_OS = f.Get("gen_ratio/transfer_matrix_noflip")
  transfer_matrix_gen_SS = f.Get("gen_ratio/transfer_matrix_flip")

  column_scaler(transfer_matrix_gen_SS)
  column_scaler(transfer_matrix_gen_OS)

  plot_transfer_matrix(transfer_matrix_gen_SS, "SS", "SS, normalized for each generated cetegory")
  plot_transfer_matrix(transfer_matrix_gen_OS, "OS", "OS, normalized for each generated cetegory")

  transfer_matrix_gen_OS = f.Get("gen_ratio/transfer_matrix_noflip")
  transfer_matrix_gen_SS = f.Get("gen_ratio/transfer_matrix_flip")
  row_scaler(transfer_matrix_gen_SS)
  row_scaler(transfer_matrix_gen_OS)

  plot_transfer_matrix(transfer_matrix_gen_SS, "SS_rec", "SS, normalized for each reconstructed cetegory")
  plot_transfer_matrix(transfer_matrix_gen_OS, "OS_rec", "OS, normalized for each reconstructed cetegory")
