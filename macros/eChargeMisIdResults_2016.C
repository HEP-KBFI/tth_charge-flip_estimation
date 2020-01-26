

void eChargeMisIdResults_2016() {

  TFile *tFIn=new TFile("../fit_output_data_2016_v3_wSyst_20200115/fit_res_exclusions.root");
  TH2D *h=(TH2D*)tFIn->Get("chargeMisId_data_exclusions");
  h->SetName("eChargeMisIdRates");
  h->SetTitle("eChargeMisIdRates");
  
  TCanvas *c1=new TCanvas("c1","c1",500,400);
  c1->cd();
  gPad->SetLogx();
  
  h->Draw("colz text");

  printf("nBinsX: pT %i [",h->GetNbinsX());
  for (int i=1; i <= h->GetNbinsX(); i++)
    printf("%g, ",h->GetXaxis()->GetBinLowEdge(i));
  printf("%g)\n",h->GetXaxis()->GetXmax());

  printf("nBinsY: abs(eta) %i [",h->GetNbinsY());
  for (int i=1; i <= h->GetNbinsY(); i++)
    printf("%g, ",h->GetYaxis()->GetBinLowEdge(i));
  printf("%g)\n",h->GetYaxis()->GetXmax());

  printf("eChergeMisId rates in percentage::\n");
  for (int iY=1; iY <= h->GetNbinsY(); iY++) {
    printf("eta bin: %i :: ",iY);
    for (int iX=1; iX <= h->GetNbinsX(); iX++) {
      printf(" \t %.4f +- %.4f ",h->GetBinContent(iX,iY)*100,h->GetBinError(iX,iY)*100);
    }
    printf("\n");
  }
  
  TFile *tOut=new TFile("ElectronChargeMisIdRates_era2016_2020Jan15.root","RECREATE");
  tOut->cd();
  h->Write();
  tOut->Close();
  
}
