from tthAnalysis.ChargeFlipEstimation.utils import read_category_ratios, bin_names_composite, bin_names_composite_nice, \
                                                   mkdir_p, get_bin_name_single, get_bin_name, get_bin_name_nice, \
                                                   get_component_cats, make_title

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
  for c in bin_names_composite:
    (c1, c2) = contained_in(c)
    if c1 == component or c2 == component: containers.append(c)
  return containers

def readMisIDRatiosGen(infile, processes = ["DY"], rec = ""):
  ratios = {}
  ratios_num = []
  f = ROOT.TFile(infile)
  for p in processes:
    effs = f.Get("gen_ratio/pt_eta_%s%s" % (p, rec))
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
        #print "Bin (%d, %d): Eff = %f + %f - %f" % (bin_eta, bin_pt, eff * 100, effErrHi * 100, effErrLo*100)
  return (ratios_num,ratios)

def readCategoryRatiosGen(infile, exclude_bins = [], gen = "gen"):
  f = ROOT.TFile(infile)
  i = 0
  os_err = ROOT.Double()
  ss_err = ROOT.Double()
  ratios = {}
  ratios_num = []
  for cat in bin_names_composite:
    histo_OS = f.Get("%s/OS/%s/mass_ll" % (gen, cat))
    histo_SS = f.Get("%s/SS/%s/mass_ll" % (gen, cat))
    os_count = histo_OS.IntegralAndError(0, histo_OS.GetNbinsX()+2, os_err)
    ss_count = histo_SS.IntegralAndError(0, histo_SS.GetNbinsX()+2, ss_err)
    if os_count > 0:
      ratio = ss_count / (ss_count + os_count)
      err = (ss_count + ss_err) / (ss_count + ss_err + os_count - os_err) - ratio
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

#Creates 21 category ratios by summing the 6 and their uncertainties
def makeCatRatiosFrom6(misIDRatios, excluded = []):
  ratios = {}
  ratios_num = []
  for cat in bin_names_composite_nice:
    if cat in excluded: continue
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

  for b in range(1, len(bin_names_composite_nice) + 1):
    pull_plots.append(ROOT.TH1D("cats%d" % b, "", 21, 0, 21))
    sum_plots.append(ROOT.TH1D("sums%d" % b, "sums%d" % b, 21, 0, 21))
    sum_plots_2.append(ROOT.TH1D("sums2_%d" % b, "sums2_%d" % b, 21, 0, 21))
    gen_plot.GetXaxis().SetBinLabel(b, bin_names_composite_nice[b - 1])
    sum_plot.GetXaxis().SetBinLabel(b, bin_names_composite_nice[b - 1])

  for b in range(1, len(bin_names_composite_nice) + 1):
    for i in range(1, len(bin_names_composite_nice) + 1):
      pull_plots[b - 1].GetXaxis().SetBinLabel(i, bin_names_composite_nice[i - 1])
      sum_plots[b - 1].GetXaxis().SetBinLabel(i, bin_names_composite_nice[i - 1])
      sum_plots_2[b - 1].GetXaxis().SetBinLabel(i, bin_names_composite_nice[i - 1])

    if bin_names_composite_nice[b - 1] in excluded: continue

    (cat1, cat2) = get_component_cats(bin_names_composite_nice[b - 1])
    (value_gen, err, err_plus) = catRatios[bin_names_composite_nice[b - 1]]
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
    chi2s[bin_names_composite_nice[b - 1]] = abs(test1.GetBinContent(1) - test2.GetBinContent(1)) / (
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
