

#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>
#include <TH1.h>
#include <THStack.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TPaveText.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TList.h>
#include <TF1.h>
#include <TColor.h>

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <assert.h>

bool doHistoNorm_dNBydM = false;

double luminosity00 = -1.;

//void addLabel_CMS_luminosity(double x0_cms, double y0, double x0_luminosity)
//{
//  TPaveText* label_cms = new TPaveText(x0_cms, y0 + 0.050, x0_cms + 0.1900, y0 + 0.100, "NDC");
//  label_cms->AddText("CMS");
//  label_cms->SetTextFont(61);
//  label_cms->SetTextAlign(23);
//  label_cms->SetTextSize(0.055);
//  label_cms->SetTextColor(1);
//  label_cms->SetFillStyle(0);
//  label_cms->SetBorderSize(0);
//  label_cms->Draw();
// 
//  TPaveText* label_luminosity = new TPaveText(x0_luminosity, y0 + 0.050, x0_luminosity + 0.1900, y0 + 0.100, "NDC");
//  label_luminosity->AddText("12.9 fb^{-1} (13 TeV)");
//  label_luminosity->SetTextAlign(13);
//  label_luminosity->SetTextSize(0.050);
//  label_luminosity->SetTextColor(1);
//  label_luminosity->SetFillStyle(0);
//  label_luminosity->SetBorderSize(0);
//  label_luminosity->Draw();
//}

 
void addLabel_CMS_preliminary(double x0, double y0, double x0_luminosity)
{
  TPaveText* label_cms = new TPaveText(x0, y0 + 0.0025, x0 + 0.0950, y0 + 0.0600, "NDC");
  label_cms->AddText("CMS");
  label_cms->SetTextFont(61);
  label_cms->SetTextAlign(13);
  label_cms->SetTextSize(0.0575);
  label_cms->SetTextColor(1);
  label_cms->SetFillStyle(0);
  label_cms->SetBorderSize(0);
  label_cms->Draw();
  
  TPaveText* label_preliminary = new TPaveText(x0 + 0.1050, y0 - 0.0010, x0 + 0.2950, y0 + 0.0500, "NDC");
  label_preliminary->AddText("Preliminary");
  label_preliminary->SetTextFont(52);
  label_preliminary->SetTextAlign(13);
  label_preliminary->SetTextSize(0.050);
  label_preliminary->SetTextColor(1);
  label_preliminary->SetFillStyle(0);
  label_preliminary->SetBorderSize(0);
  label_preliminary->Draw();


  TPaveText* label_luminosity = new TPaveText(x0_luminosity, y0 + 0.0050, x0_luminosity + 0.1900, y0 + 0.0550, "NDC");
  label_luminosity->AddText(Form("%.1f fb^{-1} (13 TeV)",luminosity00));
  label_luminosity->SetTextAlign(13);
  label_luminosity->SetTextSize(0.050);
  label_luminosity->SetTextColor(1);
  label_luminosity->SetFillStyle(0);
  label_luminosity->SetBorderSize(0);
  label_luminosity->Draw();
}

void setStyle_uncertainty(TH1* histogram)
{
  const int color_int = 17;
  const double alpha = 0.40;
  TColor* color = gROOT->GetColor(color_int);
  static int newColor_int = -1;
  static TColor* newColor = 0;
  if ( !newColor ) {
    newColor_int = gROOT->GetListOfColors()->GetSize() + 1;
    newColor = new TColor(newColor_int, color->GetRed(), color->GetGreen(), color->GetBlue(), "", alpha);
  }
  histogram->SetLineColor(newColor_int);
  histogram->SetLineWidth(0);
  histogram->SetFillColor(newColor_int);
  histogram->SetFillStyle(1001);
}

TH1* loadHistogram(TFile* inputFile, const std::string& directoryName, const std::string& histogramName)
{
  std::string histogramName_full = directoryName;
  if ( histogramName_full.length() > 0 && histogramName_full.back() != '/' ) histogramName_full.append("/");
  histogramName_full.append(histogramName);
  std::cout << "loading histogram = " << histogramName_full << " from file = " << inputFile->GetName() << std::endl;
  TH1* histogram = dynamic_cast<TH1*>(inputFile->Get(histogramName_full.data()));
  std::cout << "histogram = " << histogram << std::endl;
  if ( !histogram ) {
    std::cerr << "Failed to load histogram = " << histogramName_full << " from file = " << inputFile->GetName() << " --> skipping !!" << std::endl;
    return 0;
  }
  if ( !histogram->GetSumw2N() ) histogram->Sumw2();

  if (doHistoNorm_dNBydM) { // normalized nEvents by bin width
    for (int iBin=1; iBin<=histogram->GetNbinsX(); iBin++) {
      double dN  = histogram->GetBinContent(iBin);
      double edN = histogram->GetBinError(iBin);
      double dM  = histogram->GetBinWidth(iBin);
      histogram->SetBinContent(iBin,  dN/dM);
      histogram->SetBinError(iBin,   edN/dM);      
    }
  }
    
  return histogram;
}

TH1* rebinHistogram(const TH1* histogram_original, unsigned numBins_rebinned, double xMin_rebinned, double xMax_rebinned)
{
  std::string histogramName_rebinned = std::string(histogram_original->GetName()).append("_rebinned");
  TH1* histogram_rebinned = new TH1D(histogramName_rebinned.data(), histogram_original->GetTitle(), numBins_rebinned, xMin_rebinned, xMax_rebinned);
  histogram_rebinned->Sumw2();
  const TAxis* axis_original = histogram_original->GetXaxis();
  int numBins_original = axis_original->GetNbins();
  for ( int idxBin = 1; idxBin <= numBins_original; ++idxBin ) {
    double binContent_original = histogram_original->GetBinContent(idxBin);
    double binError_original = histogram_original->GetBinError(idxBin);
    double binCenter_original = axis_original->GetBinCenter(idxBin);
    int binIndex_rebinned = histogram_rebinned->FindBin(binCenter_original);
    double binContent_rebinned = histogram_rebinned->GetBinContent(binIndex_rebinned);
    binContent_rebinned += binContent_original;
    histogram_rebinned->SetBinContent(binIndex_rebinned, binContent_rebinned);   
    double binError_rebinned = histogram_rebinned->GetBinError(binIndex_rebinned);
    binError_rebinned = TMath::Sqrt(binError_rebinned*binError_rebinned + binError_original*binError_original);
    histogram_rebinned->SetBinError(binIndex_rebinned, binError_rebinned);
  }
  histogram_rebinned->SetLineColor(histogram_original->GetLineColor());
  histogram_rebinned->SetLineStyle(histogram_original->GetLineStyle());
  histogram_rebinned->SetLineWidth(histogram_original->GetLineWidth());
  histogram_rebinned->SetMarkerColor(histogram_original->GetMarkerColor());
  histogram_rebinned->SetMarkerStyle(histogram_original->GetMarkerStyle());
  histogram_rebinned->SetMarkerSize(histogram_original->GetMarkerSize());
  histogram_rebinned->SetFillColor(histogram_original->GetFillColor());
  histogram_rebinned->SetFillStyle(histogram_original->GetFillStyle());
  return histogram_rebinned;
}

void dumpHistogram(TH1* histogram)
{
  std::cout << "<dumpHistogram>:" << std::endl;
  std::cout << " histogram: name = " << histogram->GetName() << ", title = " << histogram->GetTitle() << std::endl;
  TAxis* xAxis = histogram->GetXaxis();
  int numBins = xAxis->GetNbins();
  for ( int iBin = 1; iBin <= numBins; ++iBin ) {
    std::cout << "bin #" << iBin << " (x = " << xAxis->GetBinCenter(iBin) << "): " << histogram->GetBinContent(iBin) << " +/- " << histogram->GetBinError(iBin) << std::endl;
  }
}

double square(double x)
{
  return x*x;
}


void makePlot(
       const std::string& inputFilePath, const std::string& inputFileName, const std::string& directoryName, 
       unsigned numBins, double xMin, double xMax, const std::string& xAxisTitle, 
       const std::string& yAxisTitle, double yMin, double yMax, 
       bool showLegend, const std::string& legendEntry_mcSignal,
       const std::string& label1, const std::string& label2, 
       const std::string& outputFilePath, const std::string& outputFileName,
       bool useLogScale)
{
  std::string inputFileName_full = Form("%s%s", inputFilePath.data(), inputFileName.data());
  TFile* inputFile = new TFile(inputFileName_full.data());
  if ( !inputFile->IsOpen() ) {
    std::cerr << "\n\n\nFailed to open input file = " << inputFileName_full << " !!" << std::endl;
    //assert(0);
    return;
  }

  inputFile->ls();

  //TCanvas* canvas = new TCanvas("canvas", "canvas", 950, 1100);
  TCanvas* canvas = new TCanvas("canvas", "canvas", 550, 700);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->Draw();

  TPad* topPad = new TPad("topPad", "topPad", 0.00, 0.34, 1.00, 0.995);
  topPad->SetFillColor(10);
  topPad->SetTopMargin(0.065);
  topPad->SetLeftMargin(0.20);
  topPad->SetBottomMargin(0.00);
  topPad->SetRightMargin(0.04);
  topPad->SetLogy(useLogScale);
  
  TPad* bottomPad = new TPad("bottomPad", "bottomPad", 0.00, 0.01, 1.00, 0.335);
  bottomPad->SetFillColor(10);
  bottomPad->SetTopMargin(0.085);
  bottomPad->SetLeftMargin(0.20);
  bottomPad->SetBottomMargin(0.35);
  bottomPad->SetRightMargin(0.04);
  bottomPad->SetLogy(false);

  canvas->cd();
  topPad->Draw();
  topPad->cd();

  TH1* histogram_data = loadHistogram(inputFile, directoryName, "data_obs");
  if ( !histogram_data) {
    std::cerr << "\n\n\nFailed to read histogram from " << inputFileName_full << " !!" << std::endl;
    //assert(0);
    return;
  }
  histogram_data->SetMarkerColor(1);
  histogram_data->SetMarkerStyle(20);
  histogram_data->SetMarkerSize(1);
  histogram_data->SetLineColor(1);
  histogram_data->SetLineWidth(1);
  histogram_data->SetLineStyle(1);

  TH1* histogram_mcSignal = loadHistogram(inputFile, directoryName, "TotalSig");
  histogram_mcSignal->SetLineColor(kBlue);
  histogram_mcSignal->SetLineWidth(2);
  histogram_mcSignal->SetLineStyle(1);
  histogram_mcSignal->SetFillColor(10);
  histogram_mcSignal->SetFillStyle(1001);
  
  TH1* histogram_mcBgr = loadHistogram(inputFile, directoryName, "TotalBkg");
  histogram_mcBgr->SetLineColor(46);
  histogram_mcBgr->SetLineWidth(0);
  histogram_mcBgr->SetLineStyle(1);
  histogram_mcBgr->SetFillColor(46);
  histogram_mcBgr->SetFillStyle(1001);

  TH1* histogramErr_mc = loadHistogram(inputFile, directoryName, "TotalProcs");
  setStyle_uncertainty(histogramErr_mc);

  if ( !(histogram_data && histogram_mcSignal && histogram_mcBgr) ) {
    std::cerr << "Failed to load histograms !!" << std::endl;
    assert(0);
  }
  
  if (histogramErr_mc) {
    yMax = 10*TMath::Max(histogram_data->GetMaximum(), histogramErr_mc->GetMaximum());
    //double yMin1 = TMath::Min(histogram_data->GetMinimum(0.), histogramErr_mc->GetMinimum(0.));
    double yMin1 = TMath::Min(histogram_data->GetMinimum(0.), histogram_mcBgr->GetMinimum(0.));
    //std::cout << "yMin1: " << yMin1 << std::endl;
    if (yMin1 > 0.) yMin = 0.1*yMin1;
  }
  
  THStack* histogramStack_mc = new THStack();
  histogramStack_mc->Add(histogram_mcBgr);
  histogramStack_mc->Add(histogram_mcSignal);
    
  TH1* histogram_ref = histogram_data;
  histogram_ref->SetTitle("");
  histogram_ref->SetStats(false);
  histogram_ref->SetMaximum(yMax);
  histogram_ref->SetMinimum(yMin);

  TAxis* xAxis_top = histogram_ref->GetXaxis();
  xAxis_top->SetRangeUser(xMin, xMax);
  if ( xAxisTitle != "" ) xAxis_top->SetTitle(xAxisTitle.data());
  xAxis_top->SetTitleOffset(1.20);
  xAxis_top->SetLabelColor(10);
  xAxis_top->SetTitleColor(10);
    
  TAxis* yAxis_top = histogram_ref->GetYaxis();
  if ( yAxisTitle != "" ) yAxis_top->SetTitle(yAxisTitle.data());
  yAxis_top->SetTitleOffset(1.20);
  yAxis_top->SetTitleSize(0.080);
  yAxis_top->SetLabelSize(0.065);
  yAxis_top->SetTickLength(0.04);  

  histogram_ref->Draw("axis");

  // CV: avoid using THStack, as it causes segmentation violation 
  //histogramStack_mc->Draw("histsame");

  histogram_mcSignal->Add(histogram_mcBgr);
  histogram_mcSignal->Draw("histsame");
  
  histogram_mcBgr->Draw("histsame");

  if ( histogramErr_mc ) {
    //setStyle_uncertainty(histogramErr_mc);
    //histogramErr_mc->Draw("e2same");
    //histogramErr_mc->Draw("HIST same");
  }

  histogram_data->Draw("e1psame");

  histogram_ref->Draw("axissame");

  if ( showLegend ) {
    TLegend* legend = new TLegend(0.6600, 0.7400, 0.9350, 0.9200, NULL, "brNDC");
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetFillColor(10);
    legend->SetTextSize(0.050);    
    TH1* histogram_data_forLegend = (TH1*)histogram_data->Clone();
    histogram_data_forLegend->SetMarkerSize(2);
    legend->AddEntry(histogram_data_forLegend, "Observed", "p");
    legend->AddEntry(histogram_mcSignal, legendEntry_mcSignal.data(), "f");
    legend->AddEntry(histogram_mcBgr, "Backgrounds", "f");
    legend->Draw();
  }

  //addLabel_CMS_luminosity(0.2100, 0.9700, 0.6350);
  addLabel_CMS_preliminary(0.2100, 0.9700, 0.6350);

  TPaveText* label1_category = new TPaveText(0.2350, 0.8650, 0.4150, 0.9250, "NDC");
  label1_category->SetTextAlign(23);
  label1_category->AddText(label1.data());
  label1_category->SetTextSize(0.055);
  label1_category->SetTextColor(1);
  label1_category->SetFillStyle(0);
  label1_category->SetBorderSize(0);
  label1_category->Draw();

  TPaveText* label2_category = new TPaveText(0.2350, 0.7950, 0.4150, 0.8550, "NDC");
  label2_category->SetTextAlign(23);
  label2_category->AddText(label2.data());
  label2_category->SetTextSize(0.055);
  label2_category->SetTextColor(1);
  label2_category->SetFillStyle(0);
  label2_category->SetBorderSize(0);
  label2_category->Draw();

  canvas->cd();
  bottomPad->Draw();
  bottomPad->cd();

  TH1* histogramRatio = (TH1*)histogram_data->Clone("histogramRatio");
  if ( !histogramRatio->GetSumw2N() ) histogramRatio->Sumw2();
  histogramRatio->SetTitle("");
  histogramRatio->SetStats(false);
  histogramRatio->SetMinimum(-0.99);
  histogramRatio->SetMaximum(+0.99);
  histogramRatio->SetMarkerColor(histogram_data->GetMarkerColor());
  histogramRatio->SetMarkerStyle(histogram_data->GetMarkerStyle());
  histogramRatio->SetMarkerSize(histogram_data->GetMarkerSize());
  histogramRatio->SetLineColor(histogram_data->GetLineColor());

  TH1* histogramRatioUncertainty = (TH1*)histogram_data->Clone("histogramRatioUncertainty");
  if ( !histogramRatioUncertainty->GetSumw2N() ) histogramRatioUncertainty->Sumw2();
  histogramRatioUncertainty->SetMarkerColor(10);
  histogramRatioUncertainty->SetMarkerSize(0);
  setStyle_uncertainty(histogramRatioUncertainty);

  int numBins_bottom = histogramRatio->GetNbinsX();
  for ( int iBin = 1; iBin <= numBins_bottom; ++iBin ) {
    double binContent_data = histogram_data->GetBinContent(iBin);
    double binError_data = histogram_data->GetBinError(iBin);
    double binContent_mc = 0;
    double binError_mc = 0;
    if ( histogramErr_mc ) {
      binContent_mc = histogramErr_mc->GetBinContent(iBin);
      binError_mc = histogramErr_mc->GetBinError(iBin);
    } else {
      TList* histograms = histogramStack_mc->GetHists();
      TIter nextHistogram(histograms);
      double binError2_mc = 0.;
      while ( TH1* histogram = dynamic_cast<TH1*>(nextHistogram()) ) {
	binContent_mc += histogram->GetBinContent(iBin);
	binError2_mc += square(histogram->GetBinError(iBin));
      }
      binError_mc = TMath::Sqrt(binError2_mc);
    }

    if ( binContent_mc > 0. ) {
      histogramRatio->SetBinContent(iBin, binContent_data/binContent_mc - 1.0);
      histogramRatio->SetBinError(iBin, binError_data/binContent_mc);

      histogramRatioUncertainty->SetBinContent(iBin, 0.);
      histogramRatioUncertainty->SetBinError(iBin, binError_mc/binContent_mc);
    }
  }

  TAxis* xAxis_bottom = histogramRatio->GetXaxis();
  xAxis_bottom->SetRangeUser(xMin, xMax);
  xAxis_bottom->SetTitle(xAxis_top->GetTitle());
  xAxis_bottom->SetLabelColor(1);
  xAxis_bottom->SetTitleColor(1);
  xAxis_bottom->SetTitleOffset(1.05);
  xAxis_bottom->SetTitleSize(0.16);
  xAxis_bottom->SetTitleFont(xAxis_top->GetTitleFont());
  xAxis_bottom->SetLabelOffset(0.02);
  xAxis_bottom->SetLabelSize(0.12);
  xAxis_bottom->SetTickLength(0.065);
  xAxis_bottom->SetNdivisions(505);
  
  TAxis* yAxis_bottom = histogramRatio->GetYaxis();
  yAxis_bottom->SetTitle("#frac{Data - Expectation}{Expectation}");
  yAxis_bottom->SetLabelColor(1);
  yAxis_bottom->SetTitleColor(1);
  yAxis_bottom->SetTitleOffset(0.95);
  yAxis_bottom->SetTitleFont(yAxis_top->GetTitleFont());
  yAxis_bottom->SetNdivisions(505);
  yAxis_bottom->CenterTitle();
  yAxis_bottom->SetTitleSize(0.095);
  yAxis_bottom->SetLabelSize(0.110);
  yAxis_bottom->SetTickLength(0.04);  

  histogramRatio->Draw("ep");
  histogramRatioUncertainty->Draw("e2same");
  
  TF1* line = new TF1("line","0", xAxis_bottom->GetXmin(), xAxis_bottom->GetXmax());
  line->SetLineStyle(3);
  line->SetLineWidth(1);
  line->SetLineColor(kBlack);
  line->Draw("same");

  
  histogramRatio->Draw("epsame");

  histogramRatio->Draw("axissame");

  canvas->Update();

  std::string outputFileName_full = Form("%s%s", outputFilePath.data(), outputFileName.data());
  size_t idx = outputFileName_full.find_last_of('.');
  std::string outputFileName_plot = std::string(outputFileName_full, 0, idx);
  if ( useLogScale ) outputFileName_plot.append("_log");
  else outputFileName_plot.append("_linear");
  //canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  //canvas->Print(std::string(outputFileName_plot).append(".root").data());

  //delete label_cms;
  delete topPad;
  delete label1_category;
  delete label2_category;
  delete histogramRatio;
  delete histogramRatioUncertainty;
  delete line;
  delete bottomPad;    
  delete canvas;

  delete inputFile;
}

// to be used to read histograms from prepareDatacards.root
void makePlot1(
       const std::string& inputFilePath, const std::string& inputFileName, const std::string& directoryName, 
       unsigned numBins, double xMin, double xMax, const std::string& xAxisTitle, 
       const std::string& yAxisTitle, double yMin, double yMax, 
       bool showLegend, const std::string& legendEntry_mcSignal,
       const std::string& label1, const std::string& label2, 
       const std::string& outputFilePath, const std::string& outputFileName,
       bool useLogScale)
{
  std::string inputFileName_full = Form("%s%s", inputFilePath.data(), inputFileName.data());
  TFile* inputFile = new TFile(inputFileName_full.data());
  if ( !inputFile ) {
    std::cerr << "\n\n\nFailed to open input file = " << inputFileName_full << " !!" << std::endl;
    //assert(0);
    return; 
  }

  inputFile->ls();
  
  //TCanvas* canvas = new TCanvas("canvas", "canvas", 950, 1100);
  TCanvas* canvas = new TCanvas("canvas", "canvas", 550, 700);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->Draw();

  TPad* topPad = new TPad("topPad", "topPad", 0.00, 0.34, 1.00, 0.995);
  topPad->SetFillColor(10);
  topPad->SetTopMargin(0.065);
  topPad->SetLeftMargin(0.20);
  topPad->SetBottomMargin(0.00);
  topPad->SetRightMargin(0.04);
  topPad->SetLogy(useLogScale);
  
  TPad* bottomPad = new TPad("bottomPad", "bottomPad", 0.00, 0.01, 1.00, 0.335);
  bottomPad->SetFillColor(10);
  bottomPad->SetTopMargin(0.085);
  bottomPad->SetLeftMargin(0.20);
  bottomPad->SetBottomMargin(0.35);
  bottomPad->SetRightMargin(0.04);
  bottomPad->SetLogy(false);

  canvas->cd();
  topPad->Draw();
  topPad->cd();

  TH1* histogram_data = loadHistogram(inputFile, directoryName, "data_obs_rebinned");
  histogram_data->SetMarkerColor(1);
  histogram_data->SetMarkerStyle(20);
  histogram_data->SetMarkerSize(2);
  histogram_data->SetLineColor(1);
  histogram_data->SetLineWidth(1);
  histogram_data->SetLineStyle(1);


  TH1* histogram_mcSignal = loadHistogram(inputFile, directoryName, "DY_rebinned");
  histogram_mcSignal->SetLineColor(1);
  histogram_mcSignal->SetLineWidth(2);
  histogram_mcSignal->SetLineStyle(1);
  histogram_mcSignal->SetFillColor(10);
  histogram_mcSignal->SetFillStyle(1001);

  TH1* histogram_mcBgr1 = loadHistogram(inputFile, directoryName, "DY_fake_rebinned");
  TH1* histogram_mcBgr2 = loadHistogram(inputFile, directoryName, "Diboson_rebinned");
  TH1* histogram_mcBgr3 = loadHistogram(inputFile, directoryName, "Singletop_rebinned");
  TH1* histogram_mcBgr4 = loadHistogram(inputFile, directoryName, "TTbar_rebinned");
  
  TH1* histogram_mcBgr = (TH1*)histogram_mcBgr1->Clone("TotalBk");
  histogram_mcBgr->Add(histogram_mcBgr2);
  histogram_mcBgr->Add(histogram_mcBgr3);
  histogram_mcBgr->Add(histogram_mcBgr4);
  
  histogram_mcBgr->SetLineColor(46);
  histogram_mcBgr->SetLineWidth(0);
  histogram_mcBgr->SetLineStyle(1);
  histogram_mcBgr->SetFillColor(46);
  histogram_mcBgr->SetFillStyle(1001);

  TH1* histogramErr_mc = (TH1*)histogram_mcSignal->Clone("TotalSignal");
  histogramErr_mc->Add(histogram_mcBgr);
  setStyle_uncertainty(histogramErr_mc);

  if ( !(histogram_data && histogram_mcSignal && histogram_mcBgr) ) {
    std::cerr << "Failed to load histograms !!" << std::endl;
    assert(0);
  }

  if (histogramErr_mc) {
    yMax = 1.3*TMath::Max(histogram_data->GetMaximum(), histogramErr_mc->GetMaximum());
  }

  THStack* histogramStack_mc = new THStack();
  histogramStack_mc->Add(histogram_mcBgr);
  histogramStack_mc->Add(histogram_mcSignal);
    
  TH1* histogram_ref = histogram_data;
  histogram_ref->SetTitle("");
  histogram_ref->SetStats(false);
  histogram_ref->SetMaximum(yMax);
  histogram_ref->SetMinimum(yMin);

  TAxis* xAxis_top = histogram_ref->GetXaxis();
  xAxis_top->SetRangeUser(xMin, xMax);
  if ( xAxisTitle != "" ) xAxis_top->SetTitle(xAxisTitle.data());
  xAxis_top->SetTitleOffset(1.20);
  xAxis_top->SetLabelColor(10);
  xAxis_top->SetTitleColor(10);
    
  TAxis* yAxis_top = histogram_ref->GetYaxis();
  if ( yAxisTitle != "" ) yAxis_top->SetTitle(yAxisTitle.data());
  yAxis_top->SetTitleOffset(1.20);
  yAxis_top->SetTitleSize(0.080);
  yAxis_top->SetLabelSize(0.065);
  yAxis_top->SetTickLength(0.04);  

  histogram_ref->Draw("axis");

  // CV: avoid using THStack, as it causes segmentation violation 
  //histogramStack_mc->Draw("histsame");

  histogram_mcSignal->Add(histogram_mcBgr);
  histogram_mcSignal->Draw("histsame");
  
  histogram_mcBgr->Draw("histsame");

  if ( histogramErr_mc ) {
    setStyle_uncertainty(histogramErr_mc);
    histogramErr_mc->Draw("e2same");
  }
    
  histogram_data->Draw("e1psame");

  histogram_ref->Draw("axissame");

  if ( showLegend ) {
    TLegend* legend = new TLegend(0.6600, 0.7400, 0.9350, 0.9200, NULL, "brNDC");
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetFillColor(10);
    legend->SetTextSize(0.050);    
    TH1* histogram_data_forLegend = (TH1*)histogram_data->Clone();
    histogram_data_forLegend->SetMarkerSize(2);
    legend->AddEntry(histogram_data_forLegend, "Observed", "p");
    legend->AddEntry(histogram_mcSignal, legendEntry_mcSignal.data(), "f");
    legend->AddEntry(histogram_mcBgr, "Backgrounds", "f");
    legend->Draw();
  }

  //addLabel_CMS_luminosity(0.2100, 0.9700, 0.6350);
  //addLabel_CMS_preliminary(0.2100, 0.9700, 0.6350);

  TPaveText* label1_category = new TPaveText(0.2350, 0.8650, 0.4150, 0.9250, "NDC");
  label1_category->SetTextAlign(23);
  label1_category->AddText(label1.data());
  label1_category->SetTextSize(0.055);
  label1_category->SetTextColor(1);
  label1_category->SetFillStyle(0);
  label1_category->SetBorderSize(0);
  label1_category->Draw();

  TPaveText* label2_category = new TPaveText(0.2350, 0.7950, 0.4150, 0.8550, "NDC");
  label2_category->SetTextAlign(23);
  label2_category->AddText(label2.data());
  label2_category->SetTextSize(0.055);
  label2_category->SetTextColor(1);
  label2_category->SetFillStyle(0);
  label2_category->SetBorderSize(0);
  label2_category->Draw();

  canvas->cd();
  bottomPad->Draw();
  bottomPad->cd();
  
  TH1* histogramRatio = (TH1*)histogram_data->Clone("histogramRatio");
  if ( !histogramRatio->GetSumw2N() ) histogramRatio->Sumw2();
  histogramRatio->SetTitle("");
  histogramRatio->SetStats(false);
  histogramRatio->SetMinimum(-0.99);
  histogramRatio->SetMaximum(+0.99);
  histogramRatio->SetMarkerColor(histogram_data->GetMarkerColor());
  histogramRatio->SetMarkerStyle(histogram_data->GetMarkerStyle());
  histogramRatio->SetMarkerSize(histogram_data->GetMarkerSize());
  histogramRatio->SetLineColor(histogram_data->GetLineColor());

  TH1* histogramRatioUncertainty = (TH1*)histogram_data->Clone("histogramRatioUncertainty");
  if ( !histogramRatioUncertainty->GetSumw2N() ) histogramRatioUncertainty->Sumw2();
  histogramRatioUncertainty->SetMarkerColor(10);
  histogramRatioUncertainty->SetMarkerSize(0);
  setStyle_uncertainty(histogramRatioUncertainty);

  int numBins_bottom = histogramRatio->GetNbinsX();
  for ( int iBin = 1; iBin <= numBins_bottom; ++iBin ) {
    double binContent_data = histogram_data->GetBinContent(iBin);
    double binError_data = histogram_data->GetBinError(iBin);
    double binContent_mc = 0;
    double binError_mc = 0;
    if ( histogramErr_mc ) {
      binContent_mc = histogramErr_mc->GetBinContent(iBin);
      binError_mc = histogramErr_mc->GetBinError(iBin);
    } else {
      TList* histograms = histogramStack_mc->GetHists();
      TIter nextHistogram(histograms);
      double binError2_mc = 0.;
      while ( TH1* histogram = dynamic_cast<TH1*>(nextHistogram()) ) {
	binContent_mc += histogram->GetBinContent(iBin);
	binError2_mc += square(histogram->GetBinError(iBin));
      }
      binError_mc = TMath::Sqrt(binError2_mc);
    }

    if ( binContent_mc > 0. ) {
      histogramRatio->SetBinContent(iBin, binContent_data/binContent_mc - 1.0);
      histogramRatio->SetBinError(iBin, binError_data/binContent_mc);

      histogramRatioUncertainty->SetBinContent(iBin, 0.);
      histogramRatioUncertainty->SetBinError(iBin, binError_mc/binContent_mc);
    }
  }

  TAxis* xAxis_bottom = histogramRatio->GetXaxis();
  xAxis_bottom->SetRangeUser(xMin, xMax);
  xAxis_bottom->SetTitle(xAxis_top->GetTitle());
  xAxis_bottom->SetLabelColor(1);
  xAxis_bottom->SetTitleColor(1);
  xAxis_bottom->SetTitleOffset(1.05);
  xAxis_bottom->SetTitleSize(0.16);
  xAxis_bottom->SetTitleFont(xAxis_top->GetTitleFont());
  xAxis_bottom->SetLabelOffset(0.02);
  xAxis_bottom->SetLabelSize(0.12);
  xAxis_bottom->SetTickLength(0.065);
  xAxis_bottom->SetNdivisions(505);
  
  TAxis* yAxis_bottom = histogramRatio->GetYaxis();
  yAxis_bottom->SetTitle("#frac{Data - Expectation}{Expectation}");
  yAxis_bottom->SetLabelColor(1);
  yAxis_bottom->SetTitleColor(1);
  yAxis_bottom->SetTitleOffset(0.95);
  yAxis_bottom->SetTitleFont(yAxis_top->GetTitleFont());
  yAxis_bottom->SetNdivisions(505);
  yAxis_bottom->CenterTitle();
  yAxis_bottom->SetTitleSize(0.095);
  yAxis_bottom->SetLabelSize(0.110);
  yAxis_bottom->SetTickLength(0.04);  

  histogramRatio->Draw("ep");
  
  TF1* line = new TF1("line","0", xAxis_bottom->GetXmin(), xAxis_bottom->GetXmax());
  line->SetLineStyle(3);
  //line->SetLineWidth(1.5);
  line->SetLineColor(kBlack);
  line->Draw("same");

  histogramRatioUncertainty->Draw("e2same");  
  histogramRatio->Draw("epsame");

  histogramRatio->Draw("axissame");

  canvas->Update();

  std::string outputFileName_full = Form("%s%s", outputFilePath.data(), outputFileName.data());
  size_t idx = outputFileName_full.find_last_of('.');
  std::string outputFileName_plot = std::string(outputFileName_full, 0, idx);
  if ( useLogScale ) outputFileName_plot.append("_log");
  else outputFileName_plot.append("_linear");
  //canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  //canvas->Print(std::string(outputFileName_plot).append(".root").data());

  //delete label_cms;
  delete topPad;
  delete label1_category;
  delete label2_category;
  delete histogramRatio;
  delete histogramRatioUncertainty;
  delete line;
  delete bottomPad;    
  delete canvas;

  delete inputFile;
}


void makeFitPlots()
{
  //gROOT->SetBatch(true);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
 
  TH1::AddDirectory(false);

  doHistoNorm_dNBydM = true; // normalize (Y-axis) bin content with bin width (dN/dM)
  if (doHistoNorm_dNBydM) std::cout << "Option doHistoNorm_dNBydM set: normalize (Y-axis) bin content with bin width (dN/dM)" << std::endl;
  
  std::vector<std::string> categories = {"BB_LL", "BB_ML", "BB_MM", "BB_HL", "BB_HM", "BB_HH", "EE_LL",
					 "EE_ML", "EE_MM", "EE_HL", "EE_HM", "EE_HH", "BE_LL", "BE_ML",
					 "EB_ML", "BE_MM", "BE_HL", "EB_HL", "BE_HM", "EB_HM", "BE_HH"};
  std::vector<std::string> charge = {"OS", "SS"};  
  std::map<std::string, int> categoryNum = {{"BB_LL", 0}, {"BB_ML", 1}, {"BB_MM", 2}, {"BB_HL", 3}, {"BB_HM", 4}, {"BB_HH", 5}, {"EE_LL", 6},
					    {"EE_ML", 7}, {"EE_MM", 8}, {"EE_HL", 9}, {"EE_HM",10}, {"EE_HH",11}, {"BE_LL",12}, {"BE_ML",13},
					    {"EB_ML",14}, {"BE_MM",15}, {"BE_HL",16}, {"EB_HL",17}, {"BE_HM",18}, {"EB_HM",19}, {"BE_HH",20}};

  std::vector<std::string> fits = {"prefit", "postfit"};
  
  bool isPlotSelectiveCategories = false;
  std::vector<std::string> plotsPrint;
  plotsPrint.push_back("OS_BB_LL");
  plotsPrint.push_back("SS_BB_LL");
  
  if ( !isPlotSelectiveCategories) {
    for (std::vector<std::string>::const_iterator cat = categories.begin();
	 cat != categories.end(); ++cat) {
      TString sCat = (*cat);
      if (categoryNum[sCat.Data()] == 6) continue;
      
      for (std::vector<std::string>::const_iterator chr = charge.begin();
	   chr != charge.end(); ++chr) {
	TString sCharge = (*chr);
	TString sChrCat = Form("%s_%s",sCharge.Data(),sCat.Data());
	
	plotsPrint.push_back(sChrCat.Data());     
	
      }
    }
  }
 
  
  std::vector<std::string> plots;
  for (std::vector<std::string>::const_iterator chr = charge.begin();
       chr != charge.end(); ++chr) {
    TString sCharge = (*chr);
    for (std::vector<std::string>::const_iterator cat = categories.begin();
	 cat != categories.end(); ++cat) {
      TString sCat = (*cat);
      TString sChrCat = Form("%s_%s",sCharge.Data(),sCat.Data());
      if (isPlotSelectiveCategories &&
	  !(std::find(plotsPrint.begin(), plotsPrint.end(), sChrCat.Data()) != plotsPrint.end())) continue;
      
      for (std::vector<std::string>::const_iterator fit = fits.begin();
	   fit != fits.end(); ++fit) {
	TString sFit = (*fit);
	TString sPlot = Form("%s_%s_%s",sCharge.Data(),sCat.Data(),sFit.Data());
	plots.push_back(sPlot.Data());
	printf("plot: %s, \t categryBin: %i\n",sPlot.Data(),categoryNum[sCat.Data()]);
      }
    }
  }  
  
  std::string inputFilePath = "./";
  std::map<std::string, std::string> inputFileNames; // key = plot
  //inputFileNames["e_EE_MM_SS"] = "prepareDatacards_charge_flip_mass_ll_v2.root";

  std::map<std::string, std::string> directoryNames; // key = plot
  //directoryNames["e_BB_MM_SS"] = "ttH_charge_flip_SS_BB_MM/rebinned";

  
  std::string outputFilePath = "";
  std::map<std::string, std::string> outputFileNames; // key = plot
  //outputFileNames["e_BB_MM_SS"] = "eTest1.root";

  std::map<std::string, double> xMin; // key = plot
  std::map<std::string, double> xMax; // key = plot
  std::map<std::string, std::string> xAxisTitles; // key = plot
  std::map<std::string, std::string> yAxisTitles; // key = plot
  std::map<std::string, unsigned> numBins; // key = plot
  std::map<std::string, double> yMin; // key = plot
  std::map<std::string, double> yMax; // key = plot
  std::map<std::string, std::string> legendEntries_mcSignal; // key = plot
  std::map<std::string, std::string> labels1; // key = plot
  std::map<std::string, std::string> labels2; // key = plot
  
  for ( std::vector<std::string>::const_iterator plot = plots.begin();
	plot != plots.end(); ++plot ) {

    TString sPlot = *plot;
    TString sCharge(sPlot(0,2));
    TString sCategory(sPlot(3,5));
    TString sCharge1,sFit;
    if (sPlot.Contains("OS"))      sCharge1 = "fail";
    if (sPlot.Contains("SS"))      sCharge1 = "pass"; 
    if (sPlot.Contains("prefit"))  sFit     = "prefit";
    if (sPlot.Contains("postfit")) sFit     = "postfit";
    printf("Plot: %18s,  category: %s,  bin: %2i,  sCharge: %s,  sFit: %s\n",
	   sPlot.Data(), sCategory.Data(), categoryNum[sCategory.Data()],
	   sCharge1.Data(), sFit.Data());

    TString myDir121 = "2017_v3_wSyst_20200115";
    luminosity00 = 41.5; // 2016: 35.9,  2017: 41.5,  2018: 51..
    
    inputFileNames[*plot]  = Form("fit_output_data_%s/bin%i/output_postfit.root",myDir121.Data(),categoryNum[sCategory.Data()]);
    outputFileNames[*plot] = Form("fit_output_data_%s/bin%i/FitPlots_%s_%i_%s_%s.root",myDir121.Data(),categoryNum[sCategory.Data()],sCharge.Data(),categoryNum[sCategory.Data()],sCategory.Data(),sFit.Data()); //Form("%s.root",sPlot.Data());
    //inputFileNames[*plot]  = Form("fit_output_pseudodata_%s/bin%i/output_postfit.root",myDir121.Data(),categoryNum[sCategory.Data()]);
    //outputFileNames[*plot] = Form("fit_output_pseudodata_%s/bin%i/FitPlots_%s_%i_%s_%s.root",myDir121.Data(),categoryNum[sCategory.Data()],sCharge.Data(),categoryNum[sCategory.Data()],sCategory.Data(),sFit.Data()); //Form("%s.root",sPlot.Data());     
    directoryNames[*plot]  = Form("%s_%s",sCharge1.Data(),sFit.Data());
    

    /*
    // tmp
    inputFileNames[*plot]  = Form("tmp/bin%i/output_postfit.root",categoryNum[sCategory.Data()]);
    outputFileNames[*plot] = Form("tmp/bin%i/FitPlots_%s_%s_%s.root",categoryNum[sCategory.Data()],sCharge.Data(),sCategory.Data(),sFit.Data()); //Form("%s.root",sPlot.Data());
    directoryNames[*plot]  = Form("%s_%s",sCharge1.Data(),sFit.Data());
    */

    xMin[*plot] =  60.;
    xMax[*plot] = 120.;
    xAxisTitles[*plot] = "m_{ee} [GeV]";
    yAxisTitles[*plot] = "dN/dm_{ee} ";
    numBins[*plot] = 60;
    yMin[*plot] = 1.01e-1;
    yMax[*plot] = 1.99e+2;
    legendEntries_mcSignal[*plot] = "Z/#gamma^{*} #rightarrow ee";
    labels1[*plot] = Form("%s %s",sCategory.Data(), sCharge.Data());
    labels2[*plot] = sFit.Data(); 
  }
  

  
  for ( std::vector<std::string>::const_iterator plot = plots.begin();
	plot != plots.end(); ++plot ) {
    cout << "\n\nPlot making for " << *plot << endl;
    if ( inputFileNames.find(*plot) == inputFileNames.end() ) {
      std::cerr << "No input file defined for plot = " << plot->data() << " !!" << std::endl;
      assert(0);
    }
    if ( outputFileNames.find(*plot) == outputFileNames.end() ) {
      std::cerr << "No output file defined for input file = " << plot->data() << " !!" << std::endl;
      assert(0);
    }
    bool showLegend = true;
    makePlot(
      inputFilePath, inputFileNames[*plot], directoryNames[*plot],
      numBins[*plot], xMin[*plot], xMax[*plot], xAxisTitles[*plot], 
      yAxisTitles[*plot], yMin[*plot], yMax[*plot], 
      showLegend, legendEntries_mcSignal[*plot], 
      labels1[*plot], labels2[*plot], 
      outputFilePath, outputFileNames[*plot], true);
    //makePlot(
    //  inputFilePath, inputFileNames[*plot], directoryNames[*plot],
    //  numBins[*plot], xMin[*plot], xMax[*plot], xAxisTitles[*plot], 
    //  yAxisTitles[*plot], yMin[*plot], yMax[*plot], 
    //  showLegend, legendEntries_mcSignal[*plot], 
    //  labels1[*plot], labels2[*plot], 
    //  outputFilePath, outputFileNames[*plot], false);
  } 
  
  
}
