from tthAnalysis.ChargeFlipEstimation.utils import read_category_ratios, BIN_NAMES_COMPOSITE, BIN_NAMES_COMPOSITE_NICE, \
                                                   mkdir_p, get_bin_name_single, get_bin_name, get_bin_name_nice, \
                                                   get_component_cats, make_title, bin_names_single

import ROOT
import math

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

def readMisIDRatiosGen(infile, processes = ["DY"], rec = ""):
  ratios = {}
  ratios_num = []
  f = ROOT.TFile(infile)
  for p in processes:
    effs = f.Get("gen_ratio/%s/pt_eta_%s%s" % (p, p, rec))
    #print "gen_ratio/pt_eta_%s" % p
    totalHisto = effs.GetTotalHistogram()
    for bin_eta in range(1, totalHisto.GetNbinsY()+1):
      for bin_pt in range(1, totalHisto.GetNbinsX()+1):
        bin = effs.GetGlobalBin(bin_pt, bin_eta)
        eff = effs.GetEfficiency(bin)
        effErrLo = effs.GetEfficiencyErrorLow(bin)
        effErrHi = effs.GetEfficiencyErrorUp(bin)
        ratios[get_bin_name_single(bin_eta, bin_pt)] = (eff, effErrLo, effErrHi)
        ratios_num.append((eff, max(effErrLo, effErrHi)))
        #print("Bin (%d, %d): Eff = %f + %f - %f" % (bin_eta, bin_pt, eff * 100, effErrHi * 100, effErrLo*100))
  return (ratios_num,ratios)

def readCategoryRatiosGen(infile, exclude_bins = [], gen = "gen"):
  f = ROOT.TFile(infile)
  i = 0
  os_err = ROOT.Double()
  ss_err = ROOT.Double()
  ratios = {}
  ratios_num = []
  for cat in BIN_NAMES_COMPOSITE:
    #histo_OS = f.Get("%s/OS/%s/mass_ll" % (gen, cat))
    #histo_SS = f.Get("%s/SS/%s/mass_ll" % (gen, cat))
    sHisto_OS0 = "%s/OS/%s/mass_ll" % (gen, cat);
    sHisto_OS = "OS/%s/DY/%s/mass_ll" % (cat,gen);
    histo_OS = f.Get(sHisto_OS)
    #print "Histo: %s, \t Histo0: %s " % (sHisto_OS,sHisto_OS0)
    sHisto_SS0 = "%s/SS/%s/mass_ll" % (gen, cat);
    sHisto_SS  = "SS/%s/DY/%s/mass_ll" % (cat, gen);
    histo_SS = f.Get(sHisto_SS)
    #print "Histo: %s, \t Histo0: %s " % (sHisto_SS,sHisto_SS0)
    os_count = histo_OS.IntegralAndError(0, histo_OS.GetNbinsX()+2, os_err)
    ss_count = histo_SS.IntegralAndError(0, histo_SS.GetNbinsX()+2, ss_err)
    if os_count > 0:
      ratio = ss_count / (ss_count + os_count)
      err  = (ss_count + ss_err) / (ss_count + ss_err + os_count - os_err) - ratio
      err1 = calUncertaintyR(os_count,os_err, ss_count,ss_err) # Siddhesh
      err = err1
    else:
      ratio = 1.
      err = 1.
    if err == 0. : err = -math.log(0.32) / os_count
    if get_bin_name_nice(i) in exclude_bins:
      i += 1
      continue
    else:
      ratios[get_bin_name_nice(i)] = (ratio, err, err)
      ratios_num.append((ratio, err, err))
      i+=1
  return (ratios_num, ratios)

def calUncertaintyR(Nos, eNos, Nss, eNss):
  term1 = (Nos**2)/((Nss+Nos)**4) * (eNss**2)
  term2 = (Nss**2)/((Nss+Nos)**4) * (eNos**2)
  return math.sqrt(term1 + term2)

#Creates 21 category ratios by summing the 6 and their uncertainties
def makeCatRatiosFrom6(misIDRatios, excluded = []):
  ratios = {}
  ratios_num = []
  for cat_idx, cat in enumerate(BIN_NAMES_COMPOSITE_NICE):
    if cat_idx in excluded: continue
    (ratio1, ratio2) = cat.split("_")
    cat_ratio = misIDRatios[ratio1][0] + misIDRatios[ratio2][0]
    #err = math.sqrt(misIDRatios[ratio1][1]**2 + misIDRatios[ratio2][1]**2)
    err = misIDRatios[ratio1][1] + misIDRatios[ratio2][1]
    ratios[cat] = (cat_ratio, err, err)
    ratios_num.append((cat_ratio, err, err))
  return (ratios_num, ratios)

# Makes pull plots for comparing 21 category ratios to sums obtained from 6
def make_pull_plot_21(misIDRatios, catRatios, name = "gen", mydir = "pull_plots_21", y_range = None, excluded = []):
  pull_plots = []
  sum_plots = []
  sum_plots_2 = []
  chi2s = {}
  sum_plot = ROOT.TH1D("sum_plot", "", 21, 0, 21)
  gen_plot = ROOT.TH1D("gen_plot", "", 21, 0, 21)
  c = ROOT.TCanvas("Plot", "Plot", 1920, 1080)
  ROOT.gStyle.SetOptStat(0)
  test1 = ROOT.TH1D("test1", "test1", 1, 0, 1)
  test2 = ROOT.TH1D("test2", "test2", 1, 0, 1)

  for b in range(1, len(BIN_NAMES_COMPOSITE_NICE) + 1):
    pull_plots.append(ROOT.TH1D("cats%d" % b, "", 21, 0, 21))
    sum_plots.append(ROOT.TH1D("sums%d" % b, "sums%d" % b, 21, 0, 21))
    sum_plots_2.append(ROOT.TH1D("sums2_%d" % b, "sums2_%d" % b, 21, 0, 21))
    gen_plot.GetXaxis().SetBinLabel(b, BIN_NAMES_COMPOSITE_NICE[b - 1])
    sum_plot.GetXaxis().SetBinLabel(b, BIN_NAMES_COMPOSITE_NICE[b - 1])

  for b in range(1, len(BIN_NAMES_COMPOSITE_NICE) + 1):
    for i in range(1, len(BIN_NAMES_COMPOSITE_NICE) + 1):
      pull_plots[b - 1].GetXaxis().SetBinLabel(i, BIN_NAMES_COMPOSITE_NICE[i - 1])
      sum_plots[b - 1].GetXaxis().SetBinLabel(i, BIN_NAMES_COMPOSITE_NICE[i - 1])
      sum_plots_2[b - 1].GetXaxis().SetBinLabel(i, BIN_NAMES_COMPOSITE_NICE[i - 1])

    if BIN_NAMES_COMPOSITE_NICE[b - 1] in excluded: continue

    (cat1, cat2) = get_component_cats(BIN_NAMES_COMPOSITE_NICE[b - 1])
    (value_gen, err, err_plus) = catRatios[BIN_NAMES_COMPOSITE_NICE[b - 1]]
    pull_plots[b - 1].SetBinContent(b, value_gen)
    pull_plots[b - 1].SetBinError(b, err)

    (value, err, err_plus) = misIDRatios[cat1]
    sum_plots[b - 1].SetBinContent(b, value)
    sum_plots[b - 1].SetBinError(b, err)

    (value, err, err_plus) = misIDRatios[cat2]
    sum_plots_2[b - 1].SetBinContent(b, value)
    sum_plots_2[b - 1].SetBinError(b, err)

    sum_plots[b - 1].Add(sum_plots_2[b - 1])

    test1.SetBinContent(1, pull_plots[b - 1].GetBinContent(b))
    test1.SetBinError(1, pull_plots[b - 1].GetBinError(b))
    test2.SetBinContent(1, sum_plots[b - 1].GetBinContent(b))
    test2.SetBinError(1, sum_plots[b - 1].GetBinError(b))
    # chi2s.append(test1.Chi2Test(test2, "WW"))
    # Chi2 method from histogram doesn't give expected results, will calculate manually
    chi2s[BIN_NAMES_COMPOSITE_NICE[b - 1]] = abs(test1.GetBinContent(1) - test2.GetBinContent(1)) / (
          test1.GetBinError(1) + test2.GetBinError(1))
    # print b, test1.GetBinContent(1), test2.GetBinContent(1)

    gen_plot.Add(pull_plots[b - 1])
    sum_plot.Add(sum_plots[b - 1])

  if y_range:
    gen_plot.SetAxisRange(y_range[0], y_range[1], "Y")

  gen_plot.SetLineColor(ROOT.kRed)
  gen_plot.SetLineWidth(3)
  sum_plot.SetLineWidth(2)
  title = make_title(name, excluded)
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

  mkdir_p(mydir)
  c.SaveAs("%s/pulls_%s.pdf" % (mydir, name))
  c.SaveAs("%s/pulls_%s.png" % (mydir, name))
  return chi2s

def calculateError_N1byN2(N1,eN1, N2,eN2): # formala r = N1/N2
  return math.sqrt(1./(N2**2)*(eN1**2) + N1**2/(N2**4)*(eN2**2) )

def compare_misIdRatios(misIDRatiosNum_original, misIDRatiosNum_closure, misIDRatiosUncNum_closure, name = "",
                        mydir = "pull_plots_21", outFileName = ""):
  print("\n\ncompare_misIdRatios():: \n Bin \t rates origianl \t\t rates closure")

  hRates_original = ROOT.TH1D("Rates_original", name, len(misIDRatiosNum_original), 0, len(misIDRatiosNum_original))
  hRates_closure = ROOT.TH1D("Rates_closure", name, len(misIDRatiosNum_original), 0, len(misIDRatiosNum_original))
  for i in range(len(misIDRatiosNum_original)):
    print("  %i %s: %f +- %f \t %f +- %f " % (
    i, bin_names_single[i], misIDRatiosNum_original[i][0], misIDRatiosNum_original[i][1], misIDRatiosNum_closure[i],
    misIDRatiosUncNum_closure[i]))
    hRates_original.GetXaxis().SetBinLabel(i + 1, bin_names_single[i])
    hRates_closure.GetXaxis().SetBinLabel(i + 1, bin_names_single[i])

    hRates_original.SetBinContent(i + 1, misIDRatiosNum_original[i][0])
    hRates_original.SetBinError(i + 1, misIDRatiosNum_original[i][1])
    hRates_closure.SetBinContent(i + 1, misIDRatiosNum_closure[i])
    hRates_closure.SetBinError(i + 1, misIDRatiosUncNum_closure[i])

  ROOT.gStyle.SetOptStat(0)
  canvas = ROOT.TCanvas("canvas", "canvas", 650, 700)
  canvas.SetFillColor(10);
  canvas.SetBorderSize(2);
  canvas.Draw();

  topPad = ROOT.TPad("topPad", "topPad", 0.00, 0.34, 1.00, 0.995);
  topPad.SetFillColor(10);
  topPad.SetTopMargin(0.1);  # 0.065
  topPad.SetLeftMargin(0.15);
  topPad.SetBottomMargin(0.00);
  topPad.SetRightMargin(0.04);
  topPad.SetTicks(1, 1)

  bottomPad = ROOT.TPad("bottomPad", "bottomPad", 0.00, 0.01, 1.00, 0.335);
  bottomPad.SetFillColor(10);
  bottomPad.SetTopMargin(0.100);
  bottomPad.SetLeftMargin(0.15);
  bottomPad.SetBottomMargin(0.35);
  bottomPad.SetRightMargin(0.04);
  bottomPad.SetTicks(1, 1)

  canvas.cd();
  topPad.Draw();
  topPad.cd();

  print("hRates_original_1:: axis range: {},  {}".format(hRates_original.GetMinimum(), hRates_original.GetMaximum()))
  hRates_original.SetLineColor(ROOT.kRed)
  hRates_original.SetLineWidth(3)
  hRates_closure.SetLineWidth(2)
  hRates_original.GetYaxis().SetTitle("Rates")
  hRates_original.SetAxisRange(0.2 * min(hRates_original.GetMinimum(), hRates_closure.GetMinimum()),
                               1.3 * max(hRates_original.GetMaximum(), hRates_closure.GetMaximum()), "Y")
  print("hRates_original_2:: axis range: {},  {}".format(hRates_original.GetMinimum(), hRates_original.GetMaximum()))

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

  canvas.cd();
  bottomPad.Draw();
  bottomPad.cd();

  hPull1 = hRates_original.Clone("pull1_misIdRates_Closure")
  if (not hPull1.GetSumw2()):
    hPull1.SetSumw2()
  hPull1.SetTitle("");
  hPull1.SetStats(False);
  # hPull1.SetMinimum(-1.2);
  # hPull1.SetMaximum(+1.2);
  hPull1.SetMarkerColor(hRates_original.GetMarkerColor());
  hPull1.SetMarkerStyle(hRates_original.GetMarkerStyle());
  hPull1.SetMarkerSize(hRates_original.GetMarkerSize());
  hPull1.SetLineColor(hRates_original.GetLineColor());
  hPull1.GetYaxis().SetTitle("#frac{Closure - Initial}{Initial}")

  hPull2 = hRates_original.Clone("pull2_misIdRates_Closure")
  if (not hPull2.GetSumw2()):
    hPull2.SetSumw2()
  hPull2.SetTitle("");
  hPull2.SetStats(False);
  # hPull2.SetMinimum(-1.2);
  # hPull2.SetMaximum(+1.2);
  hPull2.SetMarkerColor(hRates_original.GetMarkerColor());
  hPull2.SetMarkerStyle(hRates_original.GetMarkerStyle());
  hPull2.SetMarkerSize(hRates_original.GetMarkerSize());
  hPull2.SetLineColor(hRates_original.GetLineColor());
  hPull2.GetYaxis().SetTitle("#frac{Closure - Initial}{#Sigma}")

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

      print(" \t bin %i: N0: %E +- %E,  N1: %E +- %E,  N1/N0: %.12f +- %E,  eN1/eN0: %f  diff/Sigma: %E" % (
      i, N0, eN0, N1, eN1, N1 / N0, err, eN1 / eN0, (N1 - N0) / max(eN0, eN1)))

  print("hPull1_1:: axis range: {},  {}".format(hPull1.GetMinimum(), hPull1.GetMaximum()))
  # hPull1.SetAxisRange(min(1.2*hPull1.GetMinimum(), -1), 1.3*hPull1.GetMaximum(), "Y")
  hPull1.SetAxisRange(min(1.2 * hPull1.GetBinContent(hPull1.GetMinimumBin()), -1),
                      max(1.5 * hPull1.GetBinContent(hPull1.GetMaximumBin()), 1.), "Y")
  # hPull1.SetAxisRange(-1, 3, "Y")
  xAxis_bottom = hPull1.GetXaxis();
  xAxis_bottom.SetLabelSize(0.10)
  yAxis_bottom = hPull1.GetYaxis();
  yAxis_bottom.SetTitle("#frac{Closure - Initial}{Initial}")
  yAxis_bottom.CenterTitle()
  yAxis_bottom.SetNdivisions(505)
  yAxis_bottom.SetLabelSize(0.07)
  yAxis_bottom.SetTitleSize(0.07)
  yAxis_bottom.SetTitleOffset(0.9)
  print("hPull1_2:: axis range: {},  {}".format(hPull1.GetMinimum(), hPull1.GetMaximum()))
  for i in range(1, hPull1.GetNbinsX() + 1):
    print("\thPull1:: %i, %.12f %E" % (i, hPull1.GetBinContent(i), hPull1.GetBinError(i)))

  hPull1.Draw("ep1")

  line0 = ROOT.TF1("line0", "0", xAxis_bottom.GetXmin(), xAxis_bottom.GetXmax())
  line0.SetLineStyle(3);
  line0.SetLineColor(1);
  line0.Draw("same");

  hPull1.Draw("ep1 same")
  hPull1.Draw("axis same")

  canvas.Update()
  mkdir_p(mydir)
  canvas.SaveAs("%s/pulls_MCClosure_%s_1.pdf" % (mydir, name))
  canvas.SaveAs("%s/pulls_MCClosure_%s_1.png" % (mydir, name))

  canvas.cd();
  bottomPad.Draw();
  bottomPad.cd();

  print("hPull2_1:: axis range: {},  {}".format(hPull2.GetMinimum(), hPull2.GetMaximum()))
  # hPull2.SetAxisRange(0.8*hPull2.GetMinimum(), 1.3*hPull2.GetMaximum(), "Y")
  hPull2.SetAxisRange(
    min(0.5 * hPull2.GetBinContent(hPull2.GetMinimumBin()), 1.3 * hPull2.GetBinContent(hPull2.GetMinimumBin())),
    1.3 * hPull2.GetBinContent(hPull2.GetMaximumBin()), "Y")
  # hPull2.SetAxisRange(-2,2, "Y")
  xAxis_bottom = hPull2.GetXaxis();
  xAxis_bottom.SetLabelSize(0.10)
  yAxis_bottom = hPull2.GetYaxis();
  yAxis_bottom.SetTitle("#frac{Closure - Initial}{#sigma}")
  yAxis_bottom.CenterTitle()
  yAxis_bottom.SetNdivisions(505)
  yAxis_bottom.SetLabelSize(0.07)
  yAxis_bottom.SetTitleSize(0.07)
  yAxis_bottom.SetTitleOffset(0.6)
  print("hPull2_2:: axis range: {},  {}".format(hPull2.GetMinimum(), hPull2.GetMaximum()))
  for i in range(1, hPull2.GetNbinsX() + 1):
    print("\thPull2:: %i, %.12f %E" % (i, hPull2.GetBinContent(i), hPull2.GetBinError(i)))
  hPull2.SetMarkerStyle(20)
  hPull2.SetMarkerSize(2)
  hPull2.SetMarkerColor(hRates_original.GetLineColor())
  hPull2.SetLineColor(hRates_original.GetLineColor())
  hPull2.SetLineWidth(2)
  print("hRates_original: GetMarkerStyle: {}, GetMarkerSize: {}, GetMarkerColor: {}, GetLineColor: {}".format(
    hRates_original.GetMarkerStyle(), hRates_original.GetMarkerSize(), hRates_original.GetMarkerColor(),
    hRates_original.GetLineColor()))

  hPull2.Draw("ep1")

  canvas.Update()
  mkdir_p(mydir)
  canvas.SaveAs("%s/pulls_MCClosure_%s_2.pdf" % (mydir, name))
  canvas.SaveAs("%s/pulls_MCClosure_%s_2.png" % (mydir, name))
