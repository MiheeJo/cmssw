#include <TSystem.h>
#include <TROOT.h>
#include <TTree.h>
#include <TKey.h>
#include <TH1.h>
#include <TMath.h>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TPave.h>
#include <TPaveStats.h>
#include <TStyle.h>
#include <TLatex.h>
#include <iostream>
#include <JpsiFunc.h>

using namespace std;


double getAvgEffInRapPt(TH1D *h, double xmin, double xmax) {
  double avgEff = 0;
  int nbins = 0;
  int xaxismin = h->FindBin(xmin);
  int xaxismax = h->FindBin(xmax);
  for (int i=xaxismin; i<=xaxismax; i++) {
    double bincont = h->GetBinContent(i+1);
    if (bincont) {
      avgEff += bincont;
      nbins++;
    }
  }
  avgEff /= nbins;
  return avgEff;
}

void getCorrectedEffErr(const int nbins, TH1 *hrec, TH1 *hgen, TH1 *heff) {
  for (int a=0; a<nbins; a++) {
    double genInt = hgen->GetBinContent(a+1);
    double genErr = hgen->GetBinError(a+1);
    double recInt = hrec->GetBinContent(a+1);
    double recErr = hrec->GetBinError(a+1);
    double eff = recInt / genInt;

    double tmpErrGen1 = TMath::Power(eff,2) / TMath::Power(genInt,2);
    double tmpErrRec1 = TMath::Power(recErr,2);
    double tmpErr1 = tmpErrGen1 * tmpErrRec1;

    double tmpErrGen2 = TMath::Power(1-eff,2) / TMath::Power(genInt,2);
    double tmpErrRec2 = TMath::Abs(TMath::Power(genErr,2) - TMath::Power(recErr,2));
    double tmpErr2 = tmpErrGen2 * tmpErrRec2;
    double effErr = TMath::Sqrt(tmpErr1 + tmpErr2);

    if (genInt == 0) {
      heff->SetBinContent(a+1, 0);
      heff->SetBinError(a+1, 0);
    } else {
      heff->SetBinContent(a+1, eff);
//      heff->SetBinError(a+1, effErr);
      heff->SetBinError(a+1, 0); // This study doesn't have correct errors
    }
  }
}


////////////////
/// * MAIN * ///
////////////////
int getEff(void)
{
  gROOT->Macro("~/JpsiStyle.C");
  gStyle->SetOptStat(1);

  TH1::SetDefaultSumw2();
  
//  TFile *froot = new TFile("/tmp/miheejo/validation_jpsiMuMu_JpsiPt69_v5.root");
  TFile *froot = new TFile("/afs/cern.ch/user/m/miheejo/public/HIN-14-005/TnP/validation_jpsiMuMu_AllSet_v5.root");
//  TFile *froot = new TFile("/afs/cern.ch/user/m/miheejo/public/HIN-14-005/TnP/validationEDM.root");
  if (froot->IsZombie()) return 1;
  froot->cd("DQMData/Muons/RecoMuonV/MultiTrack/standAloneMuons_UpdatedAtVtx_tpToStaUpdMuonAssociation/");
  TDirectory *root_dir = gDirectory;
  root_dir->cd();
 
  double rapArrRes[] = {0, 0.8, 1.6, 2.4};
  double rapArr[] = {0, 0.8, 1.6, 2.4};
  double ptArrRes[] = {0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20};
//  double ptArrRes[] = {0, 3, 6.5, 10, 20};
  const int nintPtArrRes = sizeof(ptArrRes)/sizeof(double) -1;
  const int nintRapArrRes = sizeof(rapArrRes)/sizeof(double) -1;
  const int nintRapArr = sizeof(rapArr)/sizeof(double) -1;
  TH1F *etaSimRecRes[50], *ptSimRecRes[50];
  TH1F *etaSimRecResAll[50], *ptSimRecResAll[50];
  TH1F *num_assoc_simToReco_pT[50], *num_simul_pT[50];
  TH1F *num_assoc_simToReco_cent[50], *num_simul_cent[50];
  TH1F *num_assoc_simToReco_eta, *num_simul_eta;
  TH1F *etaSIM;
  TH1F *pTSIM;

  const int fitSize = nintPtArrRes*nintRapArrRes;
  double ptSimRecResFitResults[fitSize]={0}, etaSimRecResFitResults[fitSize]={0};
  double ptSimRecResFitErrors[fitSize]={0}, etaSimRecResFitErrors[fitSize]={0};

  root_dir->GetObject("ptSimRecRes",ptSimRecRes[0]);
  root_dir->GetObject("etaSimRecRes",etaSimRecRes[0]);

  // Resolution summary (x-axis: pt)
  // 1/ For same eta bin, read each pt plot and fit
  // 2/ Resolution from fit == sigma of fit func
  // 3/ 1 eta resol hist is filled with # of pT bins
  for (int irap=0; irap<nintRapArrRes; irap++) {
    ptSimRecResAll[irap] = new TH1F(Form("ptSimRecResAll_Rap%.0f-%.0f_Pt%.0f-%.0f",rapArrRes[irap]*10,rapArrRes[irap+1]*10,ptArrRes[0]*10,ptArrRes[nintPtArrRes]*10),";p_{T} (GeV/c);#sigma of 1 Gaussian fit",nintPtArrRes,ptArrRes);
    for (int ipt=0; ipt<nintPtArrRes; ipt++) {
      int idx = nintPtArrRes*irap + ipt+1;
      root_dir->GetObject(Form("ptSimRecRes_Rap%.0f-%.0f_Pt%.0f-%.0f",rapArrRes[irap]*10,rapArrRes[irap+1]*10,ptArrRes[ipt]*10,ptArrRes[ipt+1]*10),ptSimRecRes[idx]);

      int maxbin = ptSimRecRes[idx]->GetMaximumBin();
      double maxbinCenter = ptSimRecRes[idx]->GetBinCenter(maxbin);
      TF1 *fpt = new TF1("fpt","gaus",maxbinCenter-0.23,maxbinCenter+0.23);
      fpt->SetLineColor(kRed);
      fpt->SetParameters(0,1E-3,4E-2);
      TFitResultPtr fitres = ptSimRecRes[idx]->Fit("fpt","RS");
      int fitr = fitres;
      if (fitr != 0) {
        ptSimRecRes[idx]->Fit("fpt","RS");
      }
      ptSimRecResFitResults[idx] = fpt->GetParameter(2);
      ptSimRecResFitErrors[idx] = fpt->GetParError(2);
      cout << ptSimRecRes[idx]->GetName() << " " << ptSimRecResFitResults[idx] << " " << ptSimRecResFitErrors[idx] << endl;

      TCanvas *ctmp = new TCanvas("ctmp","ctmp",600,600);
      double min = ptSimRecRes[idx]->GetMinimum() ? ptSimRecRes[idx]->GetMinimum() : 0.1;
      double max = ptSimRecRes[idx]->GetMaximum()*1.5;
      ptSimRecRes[idx]->GetYaxis()->SetRangeUser(min,max);
      ptSimRecRes[idx]->Draw("pe");
      ctmp->SaveAs(Form("%s_Lin.png",ptSimRecRes[idx]->GetName()));
      ctmp->SaveAs(Form("%s_Lin.pdf",ptSimRecRes[idx]->GetName()));

      delete fpt;
      delete ctmp;

      double cont = ptSimRecResFitResults[idx];
      double err = 0; // ptSimRecResFitErrors[idx];
      if (idx!=1) {
        ptSimRecResAll[irap]->SetBinContent(ipt+1,cont);
        ptSimRecResAll[irap]->SetBinError(ipt+1,err);
      }
    }
  }

  // Resolution summary plot drawing (x-axis: eta)
  gStyle->SetOptStat(0);
  TCanvas *ctmp = new TCanvas("ctmp","ctmp",600,600);
  ctmp->SetLogy(1);
  SetHistStyle(ptSimRecResAll[0],0,0,1E-3,0.6);
  ptSimRecResAll[0]->Draw("p");
  
  TLegend *legpt = new TLegend(0.5,0.15,0.9,0.35);
  SetLegendStyle(legpt);
  legpt->AddEntry(ptSimRecResAll[0],"|#eta|<0.8","lp");
  for (int irap=1; irap<nintRapArrRes; irap++) {
    SetHistStyle(ptSimRecResAll[irap],irap,irap,1E-3,0.6);
    ptSimRecResAll[irap]->Draw("p same");
  }
  legpt->AddEntry(ptSimRecResAll[1],"0.8<|#eta|<1.6","lp");
  legpt->AddEntry(ptSimRecResAll[2],"1.6<|#eta|<2.4","lp");
  legpt->Draw("same");
  ctmp->SaveAs("./ptSimRecResFitResults.pdf");
  ctmp->SaveAs("./ptSimRecResFitResults.png");
  delete ctmp;
  delete legpt;
  gStyle->SetOptStat(1);


  // Resolution summary (x-axis: eta)
  // 1/ For same pT bin, read each eta plots and fit
  // 2/ Resolution from fit == sigma of fit func
  // 3/ 1 pT resol hist is filled with # of eta bins
  for (int ipt=0; ipt<nintPtArrRes; ipt++) {
    etaSimRecResAll[ipt] = new TH1F(Form("etaSimRecResAll_Rap%.0f-%.0f_Pt%.0f-%.0f",rapArrRes[0]*10,rapArrRes[nintRapArrRes]*10,ptArrRes[ipt]*10,ptArrRes[ipt+1]*10),";|#eta|;#sigma of 1 Gaussian fit",nintRapArrRes,rapArrRes);
    for (int irap=0; irap<nintRapArrRes; irap++) {
      int idx = nintRapArrRes*ipt + irap+1;
      root_dir->GetObject(Form("etaSimRecRes_Rap%.0f-%.0f_Pt%.0f-%.0f",rapArrRes[irap]*10,rapArrRes[irap+1]*10,ptArrRes[ipt]*10,ptArrRes[ipt+1]*10),etaSimRecRes[idx]);

      int maxbin = etaSimRecRes[idx]->GetMaximumBin();
      double maxbinCenter = etaSimRecRes[idx]->GetBinCenter(maxbin);
      TF1 *frap = new TF1("frap","gaus",maxbinCenter-0.17,maxbinCenter+0.17);
      frap->SetLineColor(kRed);
      frap->SetParameters(0,1E-3,4E-2);
      TFitResultPtr fitres = etaSimRecRes[idx]->Fit("frap","RS");
      int fitr = fitres;
      if (fitr != 0) {
        etaSimRecRes[idx]->Fit("frap","RS");
      }
      etaSimRecResFitResults[idx] = frap->GetParameter(2);
      etaSimRecResFitErrors[idx] = frap->GetParError(2);
      cout << etaSimRecRes[idx]->GetName() << " " << etaSimRecResFitResults[idx] << " " << etaSimRecResFitErrors[idx] << endl;

      TCanvas *ctmp = new TCanvas("ctmp","ctmp",600,600);
      double min = etaSimRecRes[idx]->GetMinimum() ? etaSimRecRes[idx]->GetMinimum() : 0.1;
      double max = etaSimRecRes[idx]->GetMaximum()*1.5;
      etaSimRecRes[idx]->GetYaxis()->SetRangeUser(min,max);
      etaSimRecRes[idx]->Draw("pe");
      ctmp->SaveAs(Form("%s_Lin.png",etaSimRecRes[idx]->GetName()));
      ctmp->SaveAs(Form("%s_Lin.pdf",etaSimRecRes[idx]->GetName()));

      delete frap;
      delete ctmp;

      double cont = etaSimRecResFitResults[idx];
      double err = 0; //etaSimRecResFitErrors[idx];
      if (idx!=1) {
        etaSimRecResAll[ipt]->SetBinContent(irap+1,cont);
        etaSimRecResAll[ipt]->SetBinError(irap+1,err);
      }
    }
  }

  // Resolution summary plot drawing (x-axis: eta)
  gStyle->SetOptStat(0);
  ctmp = new TCanvas("ctmp","ctmp",600,600);
  ctmp->SetLogy(1);
  SetHistStyle(etaSimRecResAll[0],0,0,1E-3,1.0);
  etaSimRecResAll[0]->Draw("p");

  TLegend *legeta = new TLegend(0.15,0.15,0.9,0.35);
  SetLegendStyle(legeta);
  legeta->SetNColumns(2);
  legeta->AddEntry(etaSimRecResAll[0],Form("%.1f-%.1f GeV/c",ptArrRes[0],ptArrRes[1]),"lp");
  for (int ipt=1; ipt<nintPtArrRes/2; ipt++) {
    SetHistStyle(etaSimRecResAll[ipt],ipt,ipt,1E-3,0.6);
    legeta->AddEntry(etaSimRecResAll[ipt],Form("%.1f-%.1f GeV/c",ptArrRes[ipt],ptArrRes[ipt+1]),"lp");
    etaSimRecResAll[ipt]->Draw("p same");
  }
  legeta->Draw("same");
  ctmp->SaveAs("./etaSimRecResFitResults_1.pdf");
  ctmp->SaveAs("./etaSimRecResFitResults_1.png");
  ctmp->Clear();
  delete legeta;
  legeta = new TLegend(0.15,0.15,0.9,0.35);
  SetLegendStyle(legeta);
  legeta->SetNColumns(2);
  SetHistStyle(etaSimRecResAll[nintPtArrRes/2],nintPtArrRes/2,nintPtArrRes/2,1E-3,0.6);
  etaSimRecResAll[nintPtArrRes/2]->Draw("p");
  for (int ipt=nintPtArrRes/2+1; ipt<nintPtArrRes; ipt++) {
    SetHistStyle(etaSimRecResAll[ipt],ipt,ipt,1E-3,0.6);
    legeta->AddEntry(etaSimRecResAll[ipt],Form("%.1f-%.1f GeV/c",ptArrRes[ipt],ptArrRes[ipt+1]),"lp");
    etaSimRecResAll[ipt]->Draw("p same");
  }
  legeta->Draw("same");
  ctmp->SaveAs("./etaSimRecResFitResults_2.pdf");
  ctmp->SaveAs("./etaSimRecResFitResults_2.png");
  delete ctmp;
  delete legeta;
  gStyle->SetOptStat(1);


  // Get numerator and denominator histos for efficiency
  root_dir->GetObject("num_assoc(simToReco)_pT",num_assoc_simToReco_pT[0]);
  root_dir->GetObject("num_simul_pT",num_simul_pT[0]);
  for (int irap=0; irap<nintRapArr; irap++) {
    for (int ipt=0; ipt<nintPtArrRes; ipt++) {
      int idx = nintPtArrRes*irap + ipt+1;
      cout << "eff of pT: " << idx << " " << irap << " " << ipt+1 << endl;
      root_dir->GetObject(Form("num_assoc(simToReco)_pT_Rap%.0f-%.0f_Pt%.0f-%.0f",rapArr[irap]*10,rapArr[irap+1]*10,ptArrRes[ipt]*10,ptArrRes[ipt+1]*10),num_assoc_simToReco_pT[idx]);
      root_dir->GetObject(Form("num_simul_pT_Rap%.0f-%.0f_Pt%.0f-%.0f",rapArr[irap]*10,rapArr[irap+1]*10,ptArrRes[ipt]*10,ptArrRes[ipt+1]*10),num_simul_pT[idx]);
    }
    // In a same rap region, all pT denominator and numerator histos will be merged to have 1 complete pT x-axis
    // Hence ipt should start from 1 to avoid double adding of [0] histo
    for (int ipt=1; ipt<nintPtArrRes; ipt++) {
      int idx = nintPtArrRes*irap + ipt+1;
      // Because the first histo(0-3 pT bin) will have all other bin's content, ipt==0
      int idx2 = nintPtArrRes*irap + 0+1;
      cout << "merging pT histos: " << idx << " " << idx2 << " " << irap << " " << ipt+1 << endl;
      num_assoc_simToReco_pT[idx2]->Add(num_assoc_simToReco_pT[idx]);
      num_simul_pT[idx2]->Add(num_simul_pT[idx]);
    }
  }
  root_dir->GetObject("num_assoc(simToReco)_cent",num_assoc_simToReco_cent[0]);
  root_dir->GetObject("num_simul_cent",num_simul_cent[0]);
  for (int irap=0; irap<nintRapArr; irap++) {
    for (int ipt=0; ipt<nintPtArrRes; ipt++) {
      int idx = nintPtArrRes*irap + ipt+1;
      root_dir->GetObject(Form("num_assoc(simToReco)_cent_Rap%.0f-%.0f_Pt%.0f-%.0f",rapArr[irap]*10,rapArr[irap+1]*10,ptArrRes[ipt]*10,ptArrRes[ipt+1]*10),num_assoc_simToReco_cent[idx]);
      root_dir->GetObject(Form("num_simul_cent_Rap%.0f-%.0f_Pt%.0f-%.0f",rapArr[irap]*10,rapArr[irap+1]*10,ptArrRes[ipt]*10,ptArrRes[ipt+1]*10),num_simul_cent[idx]);
    }
  }
  root_dir->GetObject("num_assoc(simToReco)_eta",num_assoc_simToReco_eta);
  root_dir->GetObject("num_simul_eta",num_simul_eta);
  root_dir->GetObject("simulation/etaSIM",etaSIM);
  root_dir->GetObject("simulation/ptSIM",pTSIM);

  // Check binnng in numerator and denominator histos for efficiency
  if (num_assoc_simToReco_eta->GetNbinsX() != num_simul_eta->GetNbinsX()) {
    cout << "num_assoc_simToReco_eta and num_simul_eta have different number of bins in X-axis" << endl;
  }
  if (num_assoc_simToReco_pT[0]->GetNbinsX() != num_simul_pT[0]->GetNbinsX()) {
    cout << "num_assoc_simToReco_pT and num_simul_pT have different number of bins in X-axis" << endl;
  }
  if (num_assoc_simToReco_cent[0]->GetNbinsX() != num_simul_cent[0]->GetNbinsX()) {
    cout << "num_assoc_simToReco_cent and num_simul_cent have different number of bins in X-axis" << endl;
  }

  const int nbins_eta = num_assoc_simToReco_eta->GetNbinsX();
  const int nbins_pT = num_assoc_simToReco_pT[0]->GetNbinsX();
  const int nbins_cent = num_assoc_simToReco_cent[0]->GetNbinsX();
  double *bins_eta = new double[nbins_eta+1];
  double *bins_pT = new double[nbins_pT+1];
  double *bins_cent = new double[nbins_cent+1];
  cout << "nbins_eta: " << nbins_eta << " nbins_pT: " << nbins_pT << " nbins_cent: " << nbins_cent << endl;
  cout << "bins_eta: " <<endl;
  for (int i=0; i<=nbins_eta; i++) {
    bins_eta[i] = num_assoc_simToReco_eta->GetXaxis()->GetBinUpEdge(i);
    cout << bins_eta[i] << "\t";
  }
  cout << endl << "bins_pT: " << endl;
  for (int i=0; i<=nbins_pT; i++) {
    bins_pT[i] = num_assoc_simToReco_pT[0]->GetXaxis()->GetBinUpEdge(i);
    cout << num_simul_pT[0]->GetXaxis()->GetBinUpEdge(i) << "\t";
    cout << bins_pT[i] << "\t";
  }
  cout << endl << "bins_cent: " << endl;
  for (int i=0; i<=nbins_cent; i++) {
    bins_cent[i] = num_assoc_simToReco_cent[0]->GetXaxis()->GetBinUpEdge(i);
    cout << num_simul_cent[0]->GetXaxis()->GetBinUpEdge(i) << "\t";
    cout << bins_cent[i] << "\t";
  }
  cout << endl;

  // Create Efficiency plots
  TH1F *eff_simToReco_eta = new TH1F("eff_simToReco_eta",";#eta;Efficiency",nbins_eta,bins_eta);
  TH1F *eff_simToReco_pT[50];
  TH1F *eff_simToReco_cent[50];

  eff_simToReco_pT[0] = new TH1F("eff_simToReco_pT",";p_{T} (GeV/c);Efficiency",nbins_pT,bins_pT);
  eff_simToReco_cent[0] = new TH1F("eff_simToReco_cent",";Centrality;Efficiency",nbins_cent,bins_cent);
  for (int irap=0; irap<nintRapArr; irap++) {
    for (int ipt=0; ipt<nintPtArrRes; ipt++) {
      int idx = nintPtArrRes*irap + ipt+1;
      eff_simToReco_pT[idx] = new TH1F(Form("eff_simToReco_pT_Rap%.0f-%.0f_Pt%.0f-%.0f",rapArr[irap]*10,rapArr[irap+1]*10,ptArrRes[ipt]*10,ptArrRes[ipt+1]*10),";p_{T} (GeV/c);Efficiency",nbins_pT,bins_pT);
      eff_simToReco_cent[idx] = new TH1F(Form("eff_simToReco_cent_Rap%.0f-%.0f_Pt%.0f-%.0f",rapArr[irap]*10,rapArr[irap+1]*10,ptArrRes[ipt]*10,ptArrRes[ipt+1]*10),";Centrality;Efficiency",nbins_cent,bins_cent);
    }
  }
  cout << eff_simToReco_pT[0]->GetName() << endl;
  cout << eff_simToReco_pT[1]->GetName() << endl;

  // rap, pT and cent efficiency plots, integrated
  getCorrectedEffErr(nbins_pT, num_assoc_simToReco_pT[0], num_simul_pT[0], eff_simToReco_pT[0]);
  getCorrectedEffErr(nbins_cent, num_assoc_simToReco_cent[0], num_simul_cent[0], eff_simToReco_cent[0]);
  getCorrectedEffErr(nbins_eta, num_assoc_simToReco_eta, num_simul_eta, eff_simToReco_eta);

  // Get pT and cent efficiency plots, differential
  for (int irap=0; irap<nintRapArr; irap++) {
    for (int ipt=0; ipt<nintPtArrRes; ipt++) {
      int idx = nintPtArrRes*irap + ipt+1;
      cout << "pT:" << nbins_pT << " " << num_assoc_simToReco_pT[idx] << " " << num_simul_pT[idx] << " " << eff_simToReco_pT[idx] << endl;
      cout << "cent:" << nbins_cent << " " << num_assoc_simToReco_cent[idx] << " " << num_simul_cent[idx] << " " << eff_simToReco_cent[idx] << endl;
      getCorrectedEffErr(nbins_pT, num_assoc_simToReco_pT[idx], num_simul_pT[idx], eff_simToReco_pT[idx]);
      getCorrectedEffErr(nbins_cent, num_assoc_simToReco_cent[idx], num_simul_cent[idx], eff_simToReco_cent[idx]);
    }
  }

  TLatex *lat = new TLatex(); lat->SetTextSize(0.035); lat->SetNDC();

  // Draw resolution plots for each pt, eta
  TCanvas *canv = new TCanvas("canv","canv",600,600);
  canv->SetLogy(1);
  for (int irap=0; irap<nintRapArrRes; irap++) {
    for (int ipt=0; ipt<nintPtArrRes; ipt++) {
      int idx = nintPtArrRes*irap + ipt+1;
      cout << idx << " " << irap << " " << ipt << endl;
      ptSimRecRes[idx]->GetXaxis()->SetTitle("(Reco-Sim)/Sim of p_{T}");
      ptSimRecRes[idx]->GetXaxis()->SetRangeUser(-1,1);
      ptSimRecRes[idx]->Draw("pe");
      lat->DrawLatex(0.15,0.9,Form("%.1f<|#eta|<%.1f, %.1f-%.1f GeV/c",rapArrRes[irap],rapArrRes[irap+1],ptArrRes[ipt],ptArrRes[ipt+1]));
      lat->DrawLatex(0.15,0.85,Form("#sigma: %.2f #pm %.2f",ptSimRecResFitResults[idx],ptSimRecResFitErrors[idx]));
      canv->SaveAs(Form("./%s.png",ptSimRecRes[idx]->GetName()));
      canv->SaveAs(Form("./%s.pdf",ptSimRecRes[idx]->GetName()));
    }
  }
  for (int ipt=0; ipt<nintPtArrRes; ipt++) {
    for (int irap=0; irap<nintRapArrRes; irap++) {
      int idx = nintRapArrRes*ipt + irap+1;
      etaSimRecRes[idx]->GetXaxis()->SetTitle("(Reco-Sim)/Sim of #eta");
      etaSimRecRes[idx]->GetXaxis()->SetRangeUser(-1,1);
      etaSimRecRes[idx]->Draw("pe");
      lat->DrawLatex(0.15,0.9,Form("%.1f<|#eta|<%.1f, %.1f-%.1f GeV/c",rapArrRes[irap],rapArrRes[irap+1],ptArrRes[ipt],ptArrRes[ipt+1]));
      lat->DrawLatex(0.15,0.85,Form("#sigma: %.2f #pm %.2f",etaSimRecResFitResults[idx],etaSimRecResFitErrors[idx]));
      canv->SaveAs(Form("./%s.png",etaSimRecRes[idx]->GetName()));
      canv->SaveAs(Form("./%s.pdf",etaSimRecRes[idx]->GetName()));
    }
  }
  canv->SetLogy(0);

  // Draw efficiency plots integrated over regions
  eff_simToReco_eta->GetXaxis()->SetRangeUser(0,2.5);
  eff_simToReco_eta->GetYaxis()->SetRangeUser(0,1.2);
  SetHistStyleDefault(eff_simToReco_eta,0,0);
  eff_simToReco_eta->Draw("p");
  canv->SaveAs("./eff_simToReco_eta.png");
  canv->SaveAs("./eff_simToReco_eta.pdf");
  SetHistStyleDefault(eff_simToReco_pT[0],0,0);
  eff_simToReco_pT[0]->GetYaxis()->SetRangeUser(0,1.2);
  eff_simToReco_pT[0]->Draw("p");
  canv->SaveAs(Form("./%s.png",eff_simToReco_pT[0]->GetName()));
  canv->SaveAs(Form("./%s.pdf",eff_simToReco_pT[0]->GetName()));
  SetHistStyleDefault(eff_simToReco_cent[0],0,0);
  eff_simToReco_cent[0]->GetYaxis()->SetRangeUser(0,1.2);
  eff_simToReco_cent[0]->Draw("p");
  canv->SaveAs(Form("./%s.png",eff_simToReco_cent[0]->GetName()));
  canv->SaveAs(Form("./%s.pdf",eff_simToReco_cent[0]->GetName()));

  // Draw efficiency plots differential
  gStyle->SetOptStat(0);
  TCanvas *canv3 = new TCanvas("canv3","canv3",600,600);
  legpt = new TLegend(0.5,0.22,0.9,0.42);
  SetLegendStyle(legpt);
  for (int irap=0; irap<nintRapArrRes; irap++) {
    TCanvas *canv2 = new TCanvas("canv2","canv2",600,600);
    TLegend *legcent = new TLegend(0.5,0.65,0.93,0.95);
//    TLegend *legcent = new TLegend(0.53,0.15,0.93,0.42);
    SetLegendStyle(legcent);
    legcent->SetNColumns(2);
    legcent->SetHeader(Form("%.1f<|#eta|<%.1f",rapArr[irap],rapArr[irap+1]));
    for (int ipt=0; ipt<nintPtArrRes; ipt++) {
      int idx = nintPtArrRes*irap + ipt+1;
      
      canv->cd();
      SetHistStyleDefault(eff_simToReco_pT[idx],irap,irap);
      eff_simToReco_pT[idx]->GetYaxis()->SetRangeUser(0,1.2);
      eff_simToReco_pT[idx]->Draw("p");
      // Because the first histo(0-3 pT bin) will have all other bin's content, ipt==0
      int idx2 = nintPtArrRes*irap + 0+1;
      if (idx!=idx2) {
        lat->DrawLatex(0.15,0.9,Form("%.1f<|#eta|<%.1f, %.1f-%.1f GeV/c",rapArr[irap],rapArr[irap+1],ptArrRes[ipt],ptArrRes[ipt+1]));
      } else {
        // This is merged pT case!
        lat->DrawLatex(0.15,0.9,Form("%.1f<|#eta|<%.1f, %.1f-%.1f GeV/c",rapArr[irap],rapArr[irap+1],ptArrRes[0],ptArrRes[nintPtArrRes]));
      }
      canv->SaveAs(Form("./%s.png",eff_simToReco_pT[idx]->GetName()));
      canv->SaveAs(Form("./%s.pdf",eff_simToReco_pT[idx]->GetName()));

      canv2->cd();
      SetHistStyleDefault(eff_simToReco_cent[idx],ipt,ipt);
      eff_simToReco_cent[idx]->GetYaxis()->SetRangeUser(0,1.2);
      if (ipt==0) {
        eff_simToReco_cent[idx]->Draw("p");
      } else {
        eff_simToReco_cent[idx]->Draw("p same");
      }
      legcent->AddEntry(eff_simToReco_cent[idx],Form("%.1f-%.1f GeV/c",ptArrRes[ipt],ptArrRes[ipt+1]),"pl");
    }
    int idx3 = nintPtArrRes*irap + 0+1;
    canv2->cd();
    legcent->Draw("same");
    canv2->SaveAs(Form("./%s.png",eff_simToReco_cent[idx3]->GetName()));
    canv2->SaveAs(Form("./%s.pdf",eff_simToReco_cent[idx3]->GetName()));
    delete legcent;
    delete canv2;

    canv3->cd();
    if (irap==0) {
      eff_simToReco_pT[idx3]->Draw("p");
      legpt->AddEntry(eff_simToReco_pT[idx3],Form("%.1f<|#eta|<%.1f",rapArr[irap],rapArr[irap+1]),"pl");
    } else if (irap!=0) {
      eff_simToReco_pT[idx3]->Draw("p same");
      legpt->AddEntry(eff_simToReco_pT[idx3],Form("%.1f<|#eta|<%.1f",rapArr[irap],rapArr[irap+1]),"pl");
    }
  }
  legpt->Draw("same");
  canv3->SaveAs("./eff_SimToReco_pT_same.png");
  canv3->SaveAs("./eff_SimToReco_pT_same.pdf");
  delete legpt;
  delete canv3;

  froot->Close();

  return 0;
}

