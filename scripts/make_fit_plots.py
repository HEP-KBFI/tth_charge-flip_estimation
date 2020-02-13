#!/usr/bin/env python

import math
import os.path
import argparse

from tthAnalysis.ChargeFlipEstimation.utils import CHARGES, BIN_NAMES_COMPOSITE, SmartFormatter

import ROOT

COLOR_INT = 17
ALPHA = 0.40
COLOR = ROOT.gROOT.GetColor(COLOR_INT)
NEWCOLOR_INT = ROOT.gROOT.GetListOfColors().GetSize() + 1
NEWCOLOR = ROOT.TColor(NEWCOLOR_INT, COLOR.GetRed(), COLOR.GetGreen(), COLOR.GetBlue(), "", ALPHA)

def setStyle_uncertainty(histogram):
  histogram.SetLineColor(NEWCOLOR_INT)
  histogram.SetLineWidth(0)
  histogram.SetFillColor(NEWCOLOR_INT)
  histogram.SetFillStyle(1001)

def loadHistogram(inputFile, directoryName, histogramName, doHistoNorm_dNBydM):
  histogramName_full = os.path.join(directoryName, histogramName)
  print("loading histogram = {} from file = {}".format(histogramName_full, inputFile.GetName()))
  histogram = inputFile.Get(histogramName_full)
  print("histogram = {}".format(histogram))
  if not histogram:
    print("Failed to load histogram = {} from file = {} -> skipping !!".format(histogramName_full, inputFile.GetName()))
    return None
  if not histogram.GetSumw2N():
    histogram.Sumw2()

  if doHistoNorm_dNBydM: # normalized nEvents by bin width
    for iBin in range(1, histogram.GetNbinsX() + 1):
      dN  = histogram.GetBinContent(iBin)
      edN = histogram.GetBinError(iBin)
      dM  = histogram.GetBinWidth(iBin)
      histogram.SetBinContent(iBin,  dN / dM)
      histogram.SetBinError  (iBin, edN / dM)
  return histogram

def rebinHistogram(histogram_original, numBins_rebinned, xMin_rebinned, xMax_rebinned):
  histogramName_rebinned = "{}_rebinned".format(histogram_original.GetName())
  histogram_rebinned = ROOT.TH1D(
    histogramName_rebinned, histogram_original.GetTitle(), numBins_rebinned, xMin_rebinned, xMax_rebinned
  )
  histogram_rebinned.Sumw2()

  axis_original = histogram_original.GetXaxis()
  numBins_original = axis_original.GetNbins()
  for idxBin in range(1, numBins_original + 1):
    binContent_original = histogram_original.GetBinContent(idxBin)
    binError_original = histogram_original.GetBinError(idxBin)
    binCenter_original = axis_original.GetBinCenter(idxBin)
    binIndex_rebinned = histogram_rebinned.FindBin(binCenter_original)
    binContent_rebinned = histogram_rebinned.GetBinContent(binIndex_rebinned)
    binContent_rebinned += binContent_original
    histogram_rebinned.SetBinContent(binIndex_rebinned, binContent_rebinned)   
    binError_rebinned = histogram_rebinned.GetBinError(binIndex_rebinned)
    binError_rebinned = math.sqrt(binError_rebinned**2 + binError_original**2)
    histogram_rebinned.SetBinError(binIndex_rebinned, binError_rebinned)

  histogram_rebinned.SetLineColor(histogram_original.GetLineColor())
  histogram_rebinned.SetLineStyle(histogram_original.GetLineStyle())
  histogram_rebinned.SetLineWidth(histogram_original.GetLineWidth())
  histogram_rebinned.SetMarkerColor(histogram_original.GetMarkerColor())
  histogram_rebinned.SetMarkerStyle(histogram_original.GetMarkerStyle())
  histogram_rebinned.SetMarkerSize(histogram_original.GetMarkerSize())
  histogram_rebinned.SetFillColor(histogram_original.GetFillColor())
  histogram_rebinned.SetFillStyle(histogram_original.GetFillStyle())
  return histogram_rebinned

def dumpHistogram(histogram):
  print("<dumpHistogram>:")
  print(" histogram: name = {} , title = {}".format(histogram.GetName(), histogram.GetTitle()))
  xAxis = histogram.GetXaxis()
  numBins = xAxis.GetNbins()
  for iBin in range(1, numBins + 1):
    print("bin #{} (x = {}): {} +/- {}".format(
      iBin, xAxis.GetBinCenter(iBin), histogram.GetBinContent(iBin), histogram.GetBinError(iBin)
    ))


def makePlot(inputFileName_full, directoryName, xMin, xMax, xAxisTitle, yAxisTitle, yMin, yMax,
             showLegend, legendEntry_mcSignal, label1, label2, outputFileName, useLogScale, luminosity00,
             doHistoNorm_dNBydM, extensions):
  inputFile = ROOT.TFile(inputFileName_full, "read")
  if not inputFile:
    print("Failed to open input file = {}".format(inputFileName_full))
    return
  inputFile.ls()

  canvas = ROOT.TCanvas("canvas", "canvas", 550, 700)
  canvas.SetFillColor(10)
  canvas.SetBorderSize(2)
  canvas.Draw()

  topPad = ROOT.TPad("topPad", "topPad", 0.00, 0.34, 1.00, 0.995)
  topPad.SetFillColor(10)
  topPad.SetTopMargin(0.065)
  topPad.SetLeftMargin(0.20)
  topPad.SetBottomMargin(0.00)
  topPad.SetRightMargin(0.04)
  topPad.SetLogy(useLogScale)
  
  bottomPad = ROOT.TPad("bottomPad", "bottomPad", 0.00, 0.01, 1.00, 0.335)
  bottomPad.SetFillColor(10)
  bottomPad.SetTopMargin(0.085)
  bottomPad.SetLeftMargin(0.20)
  bottomPad.SetBottomMargin(0.35)
  bottomPad.SetRightMargin(0.04)
  bottomPad.SetLogy(False)

  canvas.cd()
  topPad.Draw()
  topPad.cd()

  histogram_data = loadHistogram(inputFile, directoryName, "data_obs", doHistoNorm_dNBydM)
  if not histogram_data:
    print("Failed to read histogram from {}".format(inputFileName_full))
    return

  histogram_data.SetMarkerColor(1)
  histogram_data.SetMarkerStyle(20)
  histogram_data.SetMarkerSize(1)
  histogram_data.SetLineColor(1)
  histogram_data.SetLineWidth(1)
  histogram_data.SetLineStyle(1)

  histogram_mcSignal = loadHistogram(inputFile, directoryName, "TotalSig", doHistoNorm_dNBydM)
  histogram_mcSignal.SetLineColor(ROOT.kBlue)
  histogram_mcSignal.SetLineWidth(2)
  histogram_mcSignal.SetLineStyle(1)
  histogram_mcSignal.SetFillColor(10)
  histogram_mcSignal.SetFillStyle(1001)
  
  histogram_mcBgr = loadHistogram(inputFile, directoryName, "TotalBkg", doHistoNorm_dNBydM)
  histogram_mcBgr.SetLineColor(46)
  histogram_mcBgr.SetLineWidth(0)
  histogram_mcBgr.SetLineStyle(1)
  histogram_mcBgr.SetFillColor(46)
  histogram_mcBgr.SetFillStyle(1001)

  histogramErr_mc = loadHistogram(inputFile, directoryName, "TotalProcs", doHistoNorm_dNBydM)
  setStyle_uncertainty(histogramErr_mc)

  if not (histogram_data and histogram_mcSignal and histogram_mcBgr):
    raise RuntimeError("Failed to load histograms")
  
  if histogramErr_mc:
    yMax = 10 * max(histogram_data.GetMaximum(), histogramErr_mc.GetMaximum())
    yMin1 = min(histogram_data.GetMinimum(0.), histogram_mcBgr.GetMinimum(0.))
    if yMin1 > 0.:
      yMin = 0.1 * yMin1
  
  histogramStack_mc = ROOT.THStack()
  histogramStack_mc.Add(histogram_mcBgr)
  histogramStack_mc.Add(histogram_mcSignal)
    
  histogram_ref = histogram_data
  histogram_ref.SetTitle("")
  histogram_ref.SetStats(False)
  histogram_ref.SetMaximum(yMax)
  histogram_ref.SetMinimum(yMin)

  xAxis_top = histogram_ref.GetXaxis()
  xAxis_top.SetRangeUser(xMin, xMax)
  if xAxisTitle:
    xAxis_top.SetTitle(xAxisTitle)
  xAxis_top.SetTitleOffset(1.20)
  xAxis_top.SetLabelColor(10)
  xAxis_top.SetTitleColor(10)
    
  yAxis_top = histogram_ref.GetYaxis()
  if yAxisTitle:
    yAxis_top.SetTitle(yAxisTitle)
  yAxis_top.SetTitleOffset(1.20)
  yAxis_top.SetTitleSize(0.080)
  yAxis_top.SetLabelSize(0.065)
  yAxis_top.SetTickLength(0.04)  

  histogram_ref.Draw("axis")
  histogram_mcSignal.Add(histogram_mcBgr)
  histogram_mcSignal.Draw("histsame")
  histogram_mcBgr.Draw("histsame")

  histogram_data.Draw("e1psame")
  histogram_ref.Draw("axissame")

  if showLegend:
    legend = ROOT.TLegend(0.6600, 0.7400, 0.9350, 0.9200, "", "brNDC")
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetFillColor(10)
    legend.SetTextSize(0.050)    
    histogram_data_forLegend = histogram_data.Clone()
    histogram_data_forLegend.SetMarkerSize(2)
    legend.AddEntry(histogram_data_forLegend, "Observed", "p")
    legend.AddEntry(histogram_mcSignal, legendEntry_mcSignal, "f")
    legend.AddEntry(histogram_mcBgr, "Backgrounds", "f")
    legend.Draw()

  # addLabel_CMS_preliminary
  x0 = 0.2100
  y0 = 0.9700
  x0_luminosity = 0.6350
  label_cms = ROOT.TPaveText(x0, y0 + 0.0025, x0 + 0.0950, y0 + 0.0600, "NDC")
  label_cms.AddText("CMS")
  label_cms.SetTextFont(61)
  label_cms.SetTextAlign(13)
  label_cms.SetTextSize(0.0575)
  label_cms.SetTextColor(1)
  label_cms.SetFillStyle(0)
  label_cms.SetBorderSize(0)
  label_cms.Draw()

  label_preliminary = ROOT.TPaveText(x0 + 0.1050, y0 - 0.0010, x0 + 0.2950, y0 + 0.0500, "NDC")
  label_preliminary.AddText("Preliminary")
  label_preliminary.SetTextFont(52)
  label_preliminary.SetTextAlign(13)
  label_preliminary.SetTextSize(0.050)
  label_preliminary.SetTextColor(1)
  label_preliminary.SetFillStyle(0)
  label_preliminary.SetBorderSize(0)
  label_preliminary.Draw()


  label_luminosity = ROOT.TPaveText(x0_luminosity, y0 + 0.0050, x0_luminosity + 0.1900, y0 + 0.0550, "NDC")
  if luminosity00 > 0.:
    label_luminosity.AddText("%.1f fb^{-1} (13 TeV)" % luminosity00)
  else:
    label_luminosity.AddText("13 TeV")
  label_luminosity.SetTextAlign(13)
  label_luminosity.SetTextSize(0.050)
  label_luminosity.SetTextColor(1)
  label_luminosity.SetFillStyle(0)
  label_luminosity.SetBorderSize(0)
  label_luminosity.Draw()
  # addLabel_CMS_preliminary

  label1_category = ROOT.TPaveText(0.2350, 0.8650, 0.4150, 0.9250, "NDC")
  label1_category.SetTextAlign(23)
  label1_category.AddText(label1)
  label1_category.SetTextSize(0.055)
  label1_category.SetTextColor(1)
  label1_category.SetFillStyle(0)
  label1_category.SetBorderSize(0)
  label1_category.Draw()

  label2_category = ROOT.TPaveText(0.2350, 0.7950, 0.4150, 0.8550, "NDC")
  label2_category.SetTextAlign(23)
  label2_category.AddText(label2)
  label2_category.SetTextSize(0.055)
  label2_category.SetTextColor(1)
  label2_category.SetFillStyle(0)
  label2_category.SetBorderSize(0)
  label2_category.Draw()

  canvas.cd()
  bottomPad.Draw()
  bottomPad.cd()

  histogramRatio = histogram_data.Clone("histogramRatio")
  if not histogramRatio.GetSumw2N():
    histogramRatio.Sumw2()
  histogramRatio.SetTitle("")
  histogramRatio.SetStats(False)
  histogramRatio.SetMinimum(-0.99)
  histogramRatio.SetMaximum(+0.99)
  histogramRatio.SetMarkerColor(histogram_data.GetMarkerColor())
  histogramRatio.SetMarkerStyle(histogram_data.GetMarkerStyle())
  histogramRatio.SetMarkerSize(histogram_data.GetMarkerSize())
  histogramRatio.SetLineColor(histogram_data.GetLineColor())

  histogramRatioUncertainty = histogram_data.Clone("histogramRatioUncertainty")
  if not histogramRatioUncertainty.GetSumw2N():
    histogramRatioUncertainty.Sumw2()
  histogramRatioUncertainty.SetMarkerColor(10)
  histogramRatioUncertainty.SetMarkerSize(0)
  setStyle_uncertainty(histogramRatioUncertainty)

  numBins_bottom = histogramRatio.GetNbinsX()
  for iBin in range(1, numBins_bottom + 1):
    binContent_data = histogram_data.GetBinContent(iBin)
    binError_data = histogram_data.GetBinError(iBin)
    binContent_mc = 0
    binError_mc = 0
    if histogramErr_mc:
      binContent_mc = histogramErr_mc.GetBinContent(iBin)
      binError_mc = histogramErr_mc.GetBinError(iBin)
    else:
      histograms = histogramStack_mc.GetHists()
      binError2_mc = 0.
      for histogram in histograms:
        binContent_mc += histogram.GetBinContent(iBin)
        binError2_mc += histogram.GetBinError(iBin)**2
      binError_mc = math.sqrt(binError2_mc)

    if binContent_mc > 0.:
      histogramRatio.SetBinContent(iBin, binContent_data / binContent_mc - 1.)
      histogramRatio.SetBinError  (iBin,   binError_data / binContent_mc)

      histogramRatioUncertainty.SetBinContent(iBin, 0.)
      histogramRatioUncertainty.SetBinError(iBin, binError_mc/binContent_mc)

  xAxis_bottom = histogramRatio.GetXaxis()
  xAxis_bottom.SetRangeUser(xMin, xMax)
  xAxis_bottom.SetTitle(xAxis_top.GetTitle())
  xAxis_bottom.SetLabelColor(1)
  xAxis_bottom.SetTitleColor(1)
  xAxis_bottom.SetTitleOffset(1.05)
  xAxis_bottom.SetTitleSize(0.16)
  xAxis_bottom.SetTitleFont(xAxis_top.GetTitleFont())
  xAxis_bottom.SetLabelOffset(0.02)
  xAxis_bottom.SetLabelSize(0.12)
  xAxis_bottom.SetTickLength(0.065)
  xAxis_bottom.SetNdivisions(505)
  
  yAxis_bottom = histogramRatio.GetYaxis()
  yAxis_bottom.SetTitle("#frac{Data - Expectation}{Expectation}")
  yAxis_bottom.SetLabelColor(1)
  yAxis_bottom.SetTitleColor(1)
  yAxis_bottom.SetTitleOffset(0.95)
  yAxis_bottom.SetTitleFont(yAxis_top.GetTitleFont())
  yAxis_bottom.SetNdivisions(505)
  yAxis_bottom.CenterTitle()
  yAxis_bottom.SetTitleSize(0.095)
  yAxis_bottom.SetLabelSize(0.110)
  yAxis_bottom.SetTickLength(0.04)  

  histogramRatio.Draw("ep")
  histogramRatioUncertainty.Draw("e2same")
  
  line = ROOT.TF1("line","0", xAxis_bottom.GetXmin(), xAxis_bottom.GetXmax())
  line.SetLineStyle(3)
  line.SetLineWidth(1)
  line.SetLineColor(ROOT.kBlack)
  line.Draw("same")

  
  histogramRatio.Draw("epsame")

  histogramRatio.Draw("axissame")

  canvas.Update()

  outputFileName_plot, outputFileName_full_ext = os.path.splitext(outputFileName)
  if useLogScale:
    outputFileName_plot += "_log"
  else:
    outputFileName_plot += "_linear"
  for ext in extensions:
    canvas.Print("{}.{}".format(outputFileName_plot, ext))

  del topPad
  del label1_category
  del label2_category
  del histogramRatio
  del histogramRatioUncertainty
  del line
  del bottomPad
  del canvas
  del inputFile

if __name__ == '__main__':
  parser = argparse.ArgumentParser(
    formatter_class = lambda prog: SmartFormatter(prog, max_help_position = 35)
  )
  parser.add_argument('-i', '--input',
    type = str, dest = 'input', metavar = 'directory', required = True,
    help = 'R|Input fit directory',
  )
  parser.add_argument('-n', '--normalize',
    type = bool, dest = 'normalize', metavar = 'flag', required = False, default = True,
    help = 'R|Normalize (Y-axis) bin content with bin width (dN/dM)',
  )
  parser.add_argument('-c', '--charge',
    type = str, dest = 'charge', metavar = 'charge', required = False, nargs = '+', choices = CHARGES, default = CHARGES,
    help = 'R|Charge pairs to plot',
  )
  parser.add_argument('-C', '--categories',
    type = str, dest = 'categories', metavar = 'category', required = False, nargs = '+', default = [],
    help = 'R|Selective categories to plot',
  )
  parser.add_argument('-e', '--exclude',
    type = str, dest = 'exclude', metavar = 'category', required = False, nargs = '+', default = [],
    help = 'R|Selective categories to exclude (eg SS_EE_LL)',
  )
  parser.add_argument('-x', '--xmin',
    type = float, dest = 'xmin', metavar = 'number', required = False, default = 60.,
    help = 'R|Lower pT bound',
  )
  parser.add_argument('-X', '--xmax',
    type = float, dest = 'xmax', metavar = 'number', required = False, default = 120.,
    help = 'R|Upper pT bound',
  )
  parser.add_argument('-y', '--ymin',
    type = float, dest = 'ymin', metavar = 'number', required = False, default = 1.01e-1,
    help = 'R|Minimum value on tje y-axis',
  )
  parser.add_argument('-Y', '--ymax',
    type = float, dest = 'ymax', metavar = 'number', required = False, default = 1.99e+2,
    help = 'R|Maximum value on the y-axis',
  )
  parser.add_argument('-t', '--title-x',
    type = str, dest = 'title_x', metavar = 'string', required = False, default = "m_{ee} [GeV]",
    help = 'R|Title on the x-axis',
  )
  parser.add_argument('-T', '--title-y',
    type = str, dest = 'title_y', metavar = 'string', required = False, default = "dN/dm_{ee}",
    help = 'R|Title on the y-axis',
  )
  parser.add_argument('-l', '--legend',
    type = str, dest = 'legend', metavar = 'string', required = False, default = "Z/#gamma^{*} #rightarrow ee",
    help = 'R|Legend entry',
  )
  parser.add_argument('-L', '--show-legend',
    type = bool, dest = 'show_legend', metavar = 'flag', required = False, default = True,
    help = 'R|Show legend',
  )
  parser.add_argument('-u', '--show-log',
    type = bool, dest = 'show_log', metavar = 'flag', required = False, default = True,
    help = 'R|Show in log scale',
  )
  parser.add_argument('-I', '--int-lumi',
    type = float, dest = 'int_lumi', metavar = 'figure', required = False, default = -1.,
    help = 'R|Integrated luminosity in /fb'
  )
  parser.add_argument('-E', '--extension',
    type = str, dest = 'extension', metavar = 'extension', required = False, nargs = '+',
    choices = [ 'png', 'pdf', 'root' ], default = [ 'png', 'pdf' ],
    help = 'R|Extension of output plots',
  )
  args = parser.parse_args()

  ROOT.gStyle.SetPadTickX(1)
  ROOT.gStyle.SetPadTickY(1)
  ROOT.TH1.AddDirectory(False)
  ROOT.gROOT.SetBatch(True)

  input_path = os.path.abspath(args.input)
  assert(os.path.isdir(input_path))

  if args.normalize:
    print("Normalizing (Y-axis) bin content with bin width (dN/dM)")

  categoryNum = { bin_name : bin_idx for bin_idx, bin_name in enumerate(BIN_NAMES_COMPOSITE) }
  fits = [ "prefit", "postfit" ]

  plots = []
  for charge in args.charge:
    for cat in BIN_NAMES_COMPOSITE:
      charge_cat = "{}_{}".format(charge, cat)
      if args.categories and charge_cat not in args.categories:
        continue
      if charge_cat in args.exclude:
        continue
      for fit in fits:
        plot = "{}_{}".format(charge_cat, fit)
        plots.append(plot)
        print("plotting: {}, categryBin: {}".format(plot, categoryNum[cat]))

        if charge == "OS":
          pass_or_fail = "fail"
        elif charge == "SS":
          pass_or_fail = "pass"
        else:
          assert(False)

        print(
            "Plot: %18s,  category: %s,  bin: %2i,  sCharge: %s,  fit: %s" %
            (plot, cat, categoryNum[cat], pass_or_fail, fit)
        )

        bin_dir = os.path.join(input_path, "bin{}".format(categoryNum[cat]))
        inputFileName  = os.path.join(bin_dir, "output_postfit.root")
        outputFileName = os.path.join(
          bin_dir, "FitPlots_{}_{}_{}_{}.root".format(charge, categoryNum[cat], cat, fit)
        )

        makePlot(
          inputFileName, "%s_%s" % (pass_or_fail, fit),
          args.xmin, args.xmax, args.title_x,
          args.title_y, args.ymin, args.ymax,
          args.show_legend, args.legend,
          "{} {}".format(cat, charge), fit,
          outputFileName, args.show_log,
          args.int_lumi, args.normalize, args.extension
        )
