from tthAnalysis.ChargeFlipEstimation.utils import read_category_ratios, BIN_NAMES_COMPOSITE, BIN_NAMES_COMPOSITE_NICE, \
                                                   mkdir_p, get_bin_name_single, get_bin_name, get_bin_name_nice, \
                                                   get_component_cats, make_title, BIN_NAMES_SINGLE

import ROOT
import math
import os.path
import array

def contained_in(composite_name):
  comp1 = composite_name[0]+composite_name[3]
  comp2 = composite_name[1]+composite_name[4]
  return (comp1, comp2)

def get_other_component(category, composite_category):
  comps = contained_in(composite_category)
  if comps[0] == category: return comps[1]
  else: return comps[0]

def get_all_containers(component):
  containers = []
  for c in BIN_NAMES_COMPOSITE:
    (c1, c2) = contained_in(c)
    if c1 == component or c2 == component: containers.append(c)
  return containers

def readMisIDRatiosGen(infile_name, processes = [ "DY" ], rec = False):
  ratios = {}
  ratios_num = []
  infile = ROOT.TFile.Open(infile_name, 'read')
  for process in processes:
    histogram_name = "gen_ratio/{process}/pt_eta_{process}{suffix}".format(
      process = process,
      suffix  = "_rec" if rec else "",
    )
    effs = infile.Get(histogram_name)
    totalHisto = effs.GetTotalHistogram()
    for bin_eta in range(1, totalHisto.GetNbinsY() + 1):
      for bin_pt in range(1, totalHisto.GetNbinsX() + 1):
        bin = effs.GetGlobalBin(bin_pt, bin_eta)
        eff = effs.GetEfficiency(bin)
        effErrLo = effs.GetEfficiencyErrorLow(bin)
        effErrHi = effs.GetEfficiencyErrorUp(bin)
        ratios[get_bin_name_single(bin_eta, bin_pt)] = (eff, effErrLo, effErrHi)
        ratios_num.append((eff, max(effErrLo, effErrHi)))
  infile.Close()
  return (ratios_num, ratios)

def readCategoryRatiosGen(infile_name, is_gen, exclude_bins = None):
  infile = ROOT.TFile(infile_name)
  os_err = ROOT.Double()
  ss_err = ROOT.Double()
  ratios = {}
  ratios_num = []
  gen = "gen" if is_gen else "gen_rec"
  for cat_idx, cat in enumerate(BIN_NAMES_COMPOSITE):
    histo_OS = infile.Get(os.path.join("OS", cat, "DY", gen, "mass_ll"))
    histo_SS = infile.Get(os.path.join("SS", cat, "DY", gen, "mass_ll"))
    os_count = histo_OS.IntegralAndError(0, histo_OS.GetNbinsX() + 2, os_err)
    ss_count = histo_SS.IntegralAndError(0, histo_SS.GetNbinsX() + 2, ss_err)

    if os_count > 0:
      ratio = ss_count / (ss_count + os_count)
      err = calUncertaintyR(os_count,os_err, ss_count,ss_err)
    else:
      ratio = 1.
      err = 1.
    if err == 0.:
      err = -math.log(0.32) / os_count

    bin_name_nice = get_bin_name_nice(cat_idx)
    if not (exclude_bins and bin_name_nice in exclude_bins):
      ratios[bin_name_nice] = (ratio, err, err)
      ratios_num.append((ratio, err, err))
  infile.Close()
  return (ratios_num, ratios)

def calUncertaintyR(Nos, eNos, Nss, eNss):
  term1 = Nos**2 / (Nss + Nos)**4 * eNss**2
  term2 = Nss**2 / (Nss + Nos)**4 * eNos**2
  return math.sqrt(term1 + term2)

#Creates 21 category ratios by summing the 6 and their uncertainties
def makeCatRatiosFrom6(misIDRatios, excluded = None):
  ratios = {}
  ratios_num = []
  for cat_idx, cat in enumerate(BIN_NAMES_COMPOSITE_NICE):
    if excluded and cat_idx in excluded:
      continue
    (ratio1, ratio2) = cat.split("_")
    cat_ratio = misIDRatios[ratio1][0] + misIDRatios[ratio2][0]
    err = misIDRatios[ratio1][1] + misIDRatios[ratio2][1] #TODO why not add them quadratically?
    ratios[cat] = (cat_ratio, err, err)
    ratios_num.append((cat_ratio, err, err))
  return (ratios_num, ratios)

# Makes pull plots for comparing 21 category ratios to sums obtained from 6
def make_pull_plot_21(misIDRatios, catRatios, name, output_dir, y_range = None, excluded = None):
  pull_plots = []
  sum_plots = []
  sum_plots_2 = []
  chi2s = {}

  nbins = len(BIN_NAMES_COMPOSITE_NICE)
  sum_plot = ROOT.TH1D("sum_plot", "", nbins, 0, nbins)
  gen_plot = ROOT.TH1D("gen_plot", "", nbins, 0, nbins)
  sum_plot.SetStats(0)
  gen_plot.SetStats(0)
  c = ROOT.TCanvas("Plot", "Plot", 1920, 1080)

  test1 = ROOT.TH1D("test1", "test1", 1, 0, 1)
  test2 = ROOT.TH1D("test2", "test2", 1, 0, 1)
  test1.SetStats(0)
  test2.SetStats(0)

  for bin_idx in range(1, nbins + 1):
    pull_plots.append(ROOT.TH1D("cats%d" % bin_idx, "", nbins, 0, nbins))
    sum_plots.append(ROOT.TH1D("sums%d" % bin_idx, "sums%d" % bin_idx, nbins, 0, nbins))
    sum_plots_2.append(ROOT.TH1D("sums2_%d" % bin_idx, "sums2_%d" % bin_idx, nbins, 0, nbins))
    sum_plots[-1].SetStats(0)
    sum_plots_2[-1].SetStats(0)

    gen_plot.GetXaxis().SetBinLabel(bin_idx, BIN_NAMES_COMPOSITE_NICE[bin_idx - 1])
    sum_plot.GetXaxis().SetBinLabel(bin_idx, BIN_NAMES_COMPOSITE_NICE[bin_idx - 1])

  for bin_idx in range(1, nbins + 1):
    for bin_jdx in range(1, nbins + 1):
      pull_plots [bin_idx - 1].GetXaxis().SetBinLabel(bin_jdx, BIN_NAMES_COMPOSITE_NICE[bin_jdx - 1])
      sum_plots  [bin_idx - 1].GetXaxis().SetBinLabel(bin_jdx, BIN_NAMES_COMPOSITE_NICE[bin_jdx - 1])
      sum_plots_2[bin_idx - 1].GetXaxis().SetBinLabel(bin_jdx, BIN_NAMES_COMPOSITE_NICE[bin_jdx - 1])

    if excluded and (bin_idx - 1) in excluded:
      continue

    cat1, cat2 = get_component_cats(BIN_NAMES_COMPOSITE_NICE[bin_idx - 1])
    value_gen, err, err_plus = catRatios[BIN_NAMES_COMPOSITE_NICE[bin_idx - 1]]
    pull_plots[bin_idx - 1].SetBinContent(bin_idx, value_gen)
    pull_plots[bin_idx - 1].SetBinError(bin_idx, err)

    value, err, err_plus = misIDRatios[cat1]
    sum_plots[bin_idx - 1].SetBinContent(bin_idx, value)
    sum_plots[bin_idx - 1].SetBinError(bin_idx, err)

    value, err, err_plus = misIDRatios[cat2]
    sum_plots_2[bin_idx - 1].SetBinContent(bin_idx, value)
    sum_plots_2[bin_idx - 1].SetBinError(bin_idx, err)

    sum_plots[bin_idx - 1].Add(sum_plots_2[bin_idx - 1])

    test1.SetBinContent(1, pull_plots[bin_idx - 1].GetBinContent(bin_idx))
    test1.SetBinError  (1, pull_plots[bin_idx - 1].GetBinError(bin_idx))
    test2.SetBinContent(1, sum_plots [bin_idx - 1].GetBinContent(bin_idx))
    test2.SetBinError  (1, sum_plots [bin_idx - 1].GetBinError(bin_idx))

    # Chi2 method from histogram doesn't give expected results, will calculate manually
    chi2s[BIN_NAMES_COMPOSITE_NICE[bin_idx - 1]] = abs(test1.GetBinContent(1) - test2.GetBinContent(1)) / \
                                                   (test1.GetBinError(1) + test2.GetBinError(1))

    gen_plot.Add(pull_plots[bin_idx - 1])
    sum_plot.Add(sum_plots[bin_idx - 1])

  if y_range:
    gen_plot.SetAxisRange(y_range[0], y_range[1], "Y")

  gen_plot.SetLineColor(ROOT.kRed)
  gen_plot.SetLineWidth(3)
  sum_plot.SetLineWidth(2)
  title = make_title(name)
  gen_plot.SetNameTitle(title, title)
  gen_plot.Draw("e1")
  sum_plot.Draw("e1 same")

  leg = ROOT.TLegend(0.5, 0.75, 0.9, 0.85)
  leg.SetBorderSize(0)
  leg.SetLineStyle(0)
  leg.SetTextSize(0.04)
  leg.SetFillColor(0)
  leg.AddEntry(sum_plot, "Sum of component categories", "l")
  leg.AddEntry(gen_plot, "Category for 2 electrons", "l")
  leg.Draw()

  mkdir_p(output_dir)
  c.SaveAs(os.path.join(output_dir, "pulls_{}.pdf".format(name)))
  c.SaveAs(os.path.join(output_dir, "pulls_{}.png".format(name)))
  return chi2s

def calculateError_N1byN2(N1,eN1, N2,eN2): # formala r = N1/N2
  return math.sqrt(1./(N2**2)*(eN1**2) + N1**2/(N2**4)*(eN2**2) )

def compare_misIdRatios(misIDRatiosNum_original, misIDRatiosNum_closure, misIDRatiosUncNum_closure, name, output_dir):
  hRates_original = ROOT.TH1D("Rates_original", name, len(misIDRatiosNum_original), 0, len(misIDRatiosNum_original))
  hRates_closure  = ROOT.TH1D("Rates_closure",  name, len(misIDRatiosNum_original), 0, len(misIDRatiosNum_original))

  print("Idx/bin: original vs closure")
  for i in range(len(misIDRatiosNum_original)):
    print("  %i %s: %f +- %f        %f +- %f " % (
      i, BIN_NAMES_SINGLE[i],
      misIDRatiosNum_original[i][0], misIDRatiosNum_original[i][1],
      misIDRatiosNum_closure[i],     misIDRatiosUncNum_closure[i]
    ))
    hRates_original.GetXaxis().SetBinLabel(i + 1, BIN_NAMES_SINGLE[i])
    hRates_closure.GetXaxis().SetBinLabel(i + 1, BIN_NAMES_SINGLE[i])

    hRates_original.SetBinContent(i + 1, misIDRatiosNum_original[i][0])
    hRates_original.SetBinError(i + 1, misIDRatiosNum_original[i][1])
    hRates_closure.SetBinContent(i + 1, misIDRatiosNum_closure[i])
    hRates_closure.SetBinError(i + 1, misIDRatiosUncNum_closure[i])

  hRates_original.SetStats(0)
  hRates_closure.SetStats(0)
  canvas = ROOT.TCanvas("canvas", "canvas", 650, 700)
  canvas.SetFillColor(10)
  canvas.SetBorderSize(2)
  canvas.Draw()

  topPad = ROOT.TPad("topPad", "topPad", 0.00, 0.34, 1.00, 0.995)
  topPad.SetFillColor(10)
  topPad.SetTopMargin(0.1)
  topPad.SetLeftMargin(0.15)
  topPad.SetBottomMargin(0.00)
  topPad.SetRightMargin(0.04)
  topPad.SetTicks(1, 1)

  bottomPad = ROOT.TPad("bottomPad", "bottomPad", 0.00, 0.01, 1.00, 0.335)
  bottomPad.SetFillColor(10)
  bottomPad.SetTopMargin(0.100)
  bottomPad.SetLeftMargin(0.15)
  bottomPad.SetBottomMargin(0.35)
  bottomPad.SetRightMargin(0.04)
  bottomPad.SetTicks(1, 1)

  canvas.cd()
  topPad.Draw()
  topPad.cd()

  hRates_original.SetLineColor(ROOT.kRed)
  hRates_original.SetLineWidth(3)
  hRates_closure.SetLineWidth(2)
  hRates_original.GetYaxis().SetTitle("Rates")
  hRates_original.SetAxisRange(
    0.2 * min(hRates_original.GetMinimum(), hRates_closure.GetMinimum()),
    1.3 * max(hRates_original.GetMaximum(), hRates_closure.GetMaximum()),
    "Y"
  )

  hRates_original.Draw("e1")
  hRates_closure.Draw("e1 same")

  leg = ROOT.TLegend(0.15, 0.75, 0.55, 0.85)
  leg.SetBorderSize(0)
  leg.SetLineStyle(0)
  leg.SetTextSize(0.04)
  leg.SetFillColor(0)
  leg.AddEntry(hRates_original, "Initial", "l")
  leg.AddEntry(hRates_closure, "Calculated after closure", "l")
  leg.Draw()

  canvas.cd()
  bottomPad.Draw()
  bottomPad.cd()

  hPull1 = hRates_original.Clone("pull1_misIdRates_Closure")
  if not hPull1.GetSumw2():
    hPull1.SetSumw2()
  hPull1.SetTitle("")
  hPull1.SetStats(False)
  hPull1.SetMarkerColor(hRates_original.GetMarkerColor())
  hPull1.SetMarkerStyle(hRates_original.GetMarkerStyle())
  hPull1.SetMarkerSize(hRates_original.GetMarkerSize())
  hPull1.SetLineColor(hRates_original.GetLineColor())
  hPull1.GetYaxis().SetTitle("#frac{Closure - Initial}{Initial}")

  hPull2 = hRates_original.Clone("pull2_misIdRates_Closure")
  if not hPull2.GetSumw2():
    hPull2.SetSumw2()
  hPull2.SetTitle("")
  hPull2.SetStats(False)
  hPull2.SetMarkerColor(hRates_original.GetMarkerColor())
  hPull2.SetMarkerStyle(hRates_original.GetMarkerStyle())
  hPull2.SetMarkerSize(hRates_original.GetMarkerSize())
  hPull2.SetLineColor(hRates_original.GetLineColor())
  hPull2.GetYaxis().SetTitle("#frac{Closure - Initial}{#Sigma}")

  print("(e)N0 -- original error rate; (e)N1 -- closure error rate")
  for i in range(1, hPull1.GetNbinsX() + 1):
    N0 = hRates_original.GetBinContent(i)
    eN0 = hRates_original.GetBinError(i)
    N1 = hRates_closure.GetBinContent(i)
    eN1 = hRates_closure.GetBinError(i)

    if N0 > 0.:
      err = calculateError_N1byN2(N1, eN1, N0, eN0)
      hPull1.SetBinContent(i, N1 / N0 - 1.)
      hPull1.SetBinError(i, err)

      hPull2.SetBinContent(i, (N1 - N0) / max(eN0, eN1))
      hPull2.SetBinError(i, 0)

      print("  %i: N0: %E +- %E,  N1: %E +- %E,  N1/N0: %.12f +- %E,  eN1/eN0: %f  diff/Sigma: %E" % (
        i, N0, eN0, N1, eN1, N1 / N0, err, eN1 / eN0, (N1 - N0) / max(eN0, eN1)
      ))

  hPull1.SetAxisRange(
    min(1.2 * hPull1.GetBinContent(hPull1.GetMinimumBin()), -1),
    max(1.5 * hPull1.GetBinContent(hPull1.GetMaximumBin()), 1.),
    "Y"
  )
  xAxis_bottom = hPull1.GetXaxis()
  xAxis_bottom.SetLabelSize(0.10)
  yAxis_bottom = hPull1.GetYaxis()
  yAxis_bottom.SetTitle("#frac{Closure - Initial}{Initial}")
  yAxis_bottom.CenterTitle()
  yAxis_bottom.SetNdivisions(505)
  yAxis_bottom.SetLabelSize(0.07)
  yAxis_bottom.SetTitleSize(0.07)
  yAxis_bottom.SetTitleOffset(0.9)

  hPull1.Draw("ep1")

  line0 = ROOT.TF1("line0", "0", xAxis_bottom.GetXmin(), xAxis_bottom.GetXmax())
  line0.SetLineStyle(3)
  line0.SetLineColor(1)
  line0.Draw("same")

  hPull1.Draw("ep1 same")
  hPull1.Draw("axis same")

  canvas.Update()
  mkdir_p(output_dir)
  canvas.SaveAs(os.path.join(output_dir, "pulls_MCClosure_{}_1.pdf".format(name)))
  canvas.SaveAs(os.path.join(output_dir, "pulls_MCClosure_{}_1.png".format(name)))

  canvas.cd()
  bottomPad.Draw()
  bottomPad.cd()

  hPull2.SetAxisRange(
    min(0.5 * hPull2.GetBinContent(hPull2.GetMinimumBin()), 1.3 * hPull2.GetBinContent(hPull2.GetMinimumBin())),
    1.3 * hPull2.GetBinContent(hPull2.GetMaximumBin()),
    "Y"
  )
  xAxis_bottom = hPull2.GetXaxis()
  xAxis_bottom.SetLabelSize(0.10)
  yAxis_bottom = hPull2.GetYaxis()
  yAxis_bottom.SetTitle("#frac{Closure - Initial}{#sigma}")
  yAxis_bottom.CenterTitle()
  yAxis_bottom.SetNdivisions(505)
  yAxis_bottom.SetLabelSize(0.07)
  yAxis_bottom.SetTitleSize(0.07)
  yAxis_bottom.SetTitleOffset(0.6)

  hPull2.SetMarkerStyle(20)
  hPull2.SetMarkerSize(2)
  hPull2.SetMarkerColor(hRates_original.GetLineColor())
  hPull2.SetLineColor(hRates_original.GetLineColor())
  hPull2.SetLineWidth(2)
  hPull2.Draw("ep1")

  canvas.Update()
  canvas.SaveAs(os.path.join(output_dir, "pulls_MCClosure_{}_2.pdf".format(name)))
  canvas.SaveAs(os.path.join(output_dir, "pulls_MCClosure_{}_2.png".format(name)))

def plot_ratios(ratios_data, ratios_pseudo, output_filenames):
  canvas = ROOT.TCanvas("canvas", "canvas", 1200, 700)
  canvas.cd()

  nbins = len(BIN_NAMES_COMPOSITE_NICE)
  data_plot = ROOT.TH1D("data", "", nbins, 0, nbins)
  pseudo_plot = ROOT.TH1D("pseudo", "", nbins, 0, nbins)

  data_plot.SetStats(0)
  pseudo_plot.SetStats(0)

  for bin_idx in range(1, nbins + 1):
    bin_name = BIN_NAMES_COMPOSITE_NICE[bin_idx - 1]
    data_plot.GetXaxis().SetBinLabel(bin_idx, bin_name)
    pseudo_plot.GetXaxis().SetBinLabel(bin_idx, bin_name)
    if bin_name in ratios_data:
      ratio_data = ratios_data[bin_name]
      data_plot.SetBinContent(bin_idx, 100. * ratio_data[0])
      data_plot.SetBinError(bin_idx, 100. * ratio_data[1])
    if bin_name in ratios_pseudo:
      ratio_pseudo = ratios_pseudo[bin_name]
      pseudo_plot.SetBinContent(bin_idx, 100. * ratio_pseudo[0])
      pseudo_plot.SetBinError(bin_idx, 100. * ratio_pseudo[1])

  data_plot.SetTitle("")
  data_plot.GetYaxis().SetTitle("r [%]")
  data_plot.GetXaxis().SetTitle("p_{T}-#eta bins")
  data_plot.GetYaxis().SetRangeUser(-0.03, 0.55)
  data_plot.GetXaxis().SetTitleSize(0.05)
  data_plot.GetXaxis().SetTitleOffset(0.85)
  data_plot.GetYaxis().SetTitleSize(0.05)
  data_plot.GetYaxis().SetTitleOffset(0.9)

  data_plot.SetLineWidth(1)
  data_plot.SetLineColor(2)
  data_plot.SetMarkerColor(2)
  data_plot.SetMarkerStyle(21)
  data_plot.SetMarkerSize(2)

  pseudo_plot.SetLineWidth(1)
  pseudo_plot.SetLineColor(4)
  pseudo_plot.SetMarkerColor(4)
  pseudo_plot.SetMarkerStyle(23)
  pseudo_plot.SetMarkerSize(2)

  data_plot.Draw("PE1 [] ")
  pseudo_plot.Draw("PE1 [] same")

  leg = ROOT.TLegend(0.7, 0.75, 0.95, 0.95)
  leg.SetBorderSize(0)
  leg.SetLineStyle(0)
  leg.SetTextSize(0.04)
  leg.SetFillColor(0)

  leg.AddEntry(data_plot, "Data", "lep")
  leg.AddEntry(pseudo_plot, "Pseudodata", "lep")
  leg.Draw()

  canvas.Update()
  for output_filename in output_filenames:
    canvas.SaveAs(output_filename)

def get_graph(rates, offset):
  nbins = len(BIN_NAMES_SINGLE)
  assert(nbins == len(rates))
  X  = array.array('d', [ bin_idx + offset                           for bin_idx in range(nbins) ])
  Y  = array.array('d', [ 100. * rates[BIN_NAMES_SINGLE[bin_idx]][0] for bin_idx in range(nbins) ])
  eX = array.array('d', [ 0.                                         for bin_idx in range(nbins) ])
  eY = array.array('d', [ 100. * rates[BIN_NAMES_SINGLE[bin_idx]][1] for bin_idx in range(nbins) ])
  graph = ROOT.TGraphErrors(nbins, X, Y, eX, eY)
  return graph

def plot_rates(rates_data, rates_pseudo, rates_gen, output_filenames):
  canvas = ROOT.TCanvas("canvas", "canvas", 1200, 700)
  canvas.cd()

  offset = 0.5
  nbins = len(BIN_NAMES_SINGLE)
  base = ROOT.TH1D("base", "base", nbins, offset, nbins + offset)
  for bin_idx, bin_name in enumerate(BIN_NAMES_SINGLE):
    base.GetXaxis().SetBinLabel(bin_idx + 1, bin_name)
  base.SetTitle("")
  base.GetYaxis().SetTitle("p [%]")
  base.GetXaxis().SetTitle("p_{T}-#eta bins")
  base.GetYaxis().SetRangeUser(-0.03, 0.20)
  base.GetXaxis().SetTitleSize(0.05)
  base.GetXaxis().SetTitleOffset(0.85)
  base.GetXaxis().SetLabelSize(0.05)
  base.GetYaxis().SetTitleSize(0.05)
  base.GetYaxis().SetTitleOffset(0.9)
  base.GetYaxis().SetLabelSize(0.04)
  base.GetYaxis().SetNdivisions(505)
  base.SetStats(0)

  data_graph = get_graph(rates_data, -0.15 + 2 * offset)
  pseudo_graph = get_graph(rates_pseudo, 0. + 2 * offset)
  gen_graph = get_graph(rates_gen, 0.15 + 2 * offset)

  data_graph.SetLineWidth(1)
  data_graph.SetLineColor(1)
  data_graph.SetMarkerColor(1)
  data_graph.SetMarkerStyle(21)
  data_graph.SetMarkerSize(1.3)

  pseudo_graph.SetLineWidth(1)
  pseudo_graph.SetLineColor(4)
  pseudo_graph.SetMarkerColor(4)
  pseudo_graph.SetMarkerStyle(21)
  pseudo_graph.SetMarkerSize(1.3)

  gen_graph.SetLineWidth(1)
  gen_graph.SetLineColor(2)
  gen_graph.SetMarkerColor(2)
  gen_graph.SetMarkerStyle(21)
  gen_graph.SetMarkerSize(1.3)

  base.Draw("AXIS")
  data_graph.Draw("PE1 ")
  pseudo_graph.Draw("PE1  same")
  gen_graph.Draw("PE1 same")

  leg = ROOT.TLegend(0.15, 0.70, 0.50, 0.88)
  leg.SetBorderSize(0)
  leg.SetLineStyle(0)
  leg.SetTextSize(0.04)
  leg.SetFillColor(0)
  leg.AddEntry(data_graph, "Data", "lep")
  leg.AddEntry(pseudo_graph, "Pseudodata", "lep")
  leg.AddEntry(gen_graph, "MC true", "lep")
  leg.Draw()

  canvas.Update()
  for output_filename in output_filenames:
    canvas.SaveAs(output_filename)
