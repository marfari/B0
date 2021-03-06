#include <iomanip>
#include <sstream>
#include <vector>
#include <TStyle.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TNtupleD.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <TSystem.h>
#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooGlobalFunc.h>
#include <RooGenericPdf.h>
#include <RooConstVar.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooChebychev.h>
#include <RooBernstein.h>
#include <RooExponential.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooPlotable.h>
#include <RooBifurGauss.h>
#include <RooCBShape.h>
#include "TMath.h"
#include <RooGenericPdf.h>
#include "TRatioPlot.h"
#include <RooBifurGauss.h>
#include <RooProduct.h>
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooMCStudy.h"
#include <iostream>
#include <TF1.h>
#include <RooPolynomial.h>
#include <fstream>
#include <TGraph.h>
#include "TMultiGraph.h"
#include <TEfficiency.h>


void set_up_workspace_variables(RooWorkspace& w);
void read_data(RooWorkspace& w, TString f_input);
void build_pdf (RooWorkspace& w, std::string choice = "nominal");
void plot_mass_fit(RooWorkspace& w);
void cut_tree(TString f_input);

#define apply_cut 0

// 1 = creates new tree with cuts applied
// 0 = does not create new tree

void b0(){

  TString input_file_data = "/lstore/cms/nuno/ppdata2017/001127_trkqual_b0/BZData.root";
  //TString input_file_data = "./results/cutted_tree.root";

  if(apply_cut == 1){cut_tree(input_file_data);}
  
  RooWorkspace* ws = new RooWorkspace("ws");
  set_up_workspace_variables(*ws);
  read_data(*ws,input_file_data);

  build_pdf(*ws, "nominal");
  plot_mass_fit(*ws);

}


void plot_mass_fit(RooWorkspace& w){

  RooAbsPdf* model = w.pdf("model");
  RooDataSet* data = (RooDataSet*) w.data("data"); 

  RooRealVar Bmass = *(w.var("Bmass"));
  RooRealVar* lambda   = w.var("lambda");

  Bmass.setRange("all", Bmass.getMin(),Bmass.getMax());
  model->fitTo(*data, RooFit::Range("all"));

  RooPlot* massframe = Bmass.frame();

  data->plotOn(massframe, RooFit::Name("Data"));
  model->plotOn(massframe, RooFit::Name("Fit"), RooFit::Range("all"), RooFit::LineColor(kRed), RooFit::LineStyle(1), RooFit::LineWidth(2));
  model->plotOn(massframe, RooFit::Name("Total Signal"), RooFit::Components("pdf_t_signal"), RooFit::Range("all"), RooFit::LineColor(kOrange), RooFit::LineStyle(kDashed));
  model->plotOn(massframe, RooFit::Name("Corr Tag"), RooFit::Components("pdf_c_signal"), RooFit::Range("all"), RooFit::LineColor(kGreen), RooFit::LineStyle(kDashed));
  model->plotOn(massframe, RooFit::Name("Mis Tag"), RooFit::Components("k_pi_swap"), RooFit::Range("all"), RooFit::LineColor(kGray), RooFit::LineStyle(kDashed));
  model->plotOn(massframe, RooFit::Name("Combinatorial"), RooFit::Components("pdf_m_combinatorial"), RooFit::Range("all"), RooFit::LineColor(kBlue), RooFit::LineStyle(kDashed));

  model->paramOn(massframe, RooFit::Layout(0.65,0.95,0.95));

  TCanvas d;
  d.SetTitle("");

  TPad *p1 = new TPad("p1","p1",0.0,0.27,0.82,0.99);
  p1->SetTitle("");
  p1->SetBorderMode(1);
  p1->SetFrameBorderMode(0);
  p1->SetBorderSize(2);

  p1->SetBottomMargin(0.10);

  p1->Draw();

  TPad *p2 = new TPad("p2","p2",0.0,0.065,0.82,0.24);
  p2->SetTitle("");
  p2->SetTopMargin(0.);
  p2->SetBottomMargin(0.2);

  p2->SetBorderMode(1);

  p2->Draw();

  p1->cd();

  massframe->Draw();
  TLatex* tex11 = new TLatex(0.6,0.8,"302.3 pb^{-1} (pp) 5.02 TeV");

  tex11->SetNDC(kTRUE);
  tex11->SetLineWidth(2);
  tex11->SetTextSize(0.04);
  //tex11->Draw();
  tex11 = new TLatex(0.6,0.85,"CMS Preliminary");
  tex11->SetNDC(kTRUE);
  tex11->SetTextFont(42);
  tex11->SetTextSize(0.04);
  tex11->SetLineWidth(2);
  //tex11->Draw();
  
  double lambda_str = lambda->getVal();
  double lambda_err = lambda->getError();

  double chis = massframe->chiSquare();

  TLatex* tex12 = new TLatex(0.15, 0.85, Form("#lambda_{exp} = %.3lf #pm %.3lf",lambda_str,lambda_err));
  tex12->SetNDC(kTRUE);
  tex12->SetTextFont(42);
  tex12->SetTextSize(0.04);

  TLatex* tex13 = new TLatex(0.15, 0.8, Form("#chi/DOF = %.3lf",chis));
  tex13->SetNDC(kTRUE);
  tex13->SetTextFont(42);
  tex13->SetTextSize(0.04);

  TLegend *leg = new TLegend (0.7, 0.2, 0.9, 0.4);

  leg->SetTextSize(0.03);
  leg->AddEntry(massframe->findObject("Data"), "Data", "l");
  leg->AddEntry(massframe->findObject("Fit"),"Fit","l");
  leg->AddEntry(massframe->findObject("Total Signal"), "Total sig", "l");
  leg->AddEntry(massframe->findObject("Corr Tag"), "Corr. sig", "l");
  leg->AddEntry(massframe->findObject("Mis Tag"), "Mis-tag sig", "l");
  leg->AddEntry(massframe->findObject("Combinatorial"), "Comb. bkg", "l");

  //leg->Draw("same");

  //pull dists
  RooHist* pull_hist = massframe->pullHist("Data","Fit");
  RooPlot *pull_plot = Bmass.frame();

  pull_plot->addPlotable(static_cast<RooPlotable*>(pull_hist),"P");
  pull_plot->SetTitle("");
  pull_plot->GetXaxis()->SetTitle("");
  pull_plot->GetXaxis()->SetTitleFont(42);
  pull_plot->GetXaxis()->SetTitleSize(0.17);
  pull_plot->GetXaxis()->SetTitleOffset(1.09);

  pull_plot->GetXaxis()->SetLabelFont(42);
  pull_plot->GetXaxis()->SetLabelSize(0.15);
  pull_plot->GetXaxis()->SetLabelOffset(0.01);
  pull_plot->GetXaxis()->SetTickLength(0.13);

  pull_plot->GetYaxis()->SetTitle("Pull hist");
  pull_plot->GetYaxis()->SetTitleFont(42);
  pull_plot->GetYaxis()->SetTitleSize(0.10);
  pull_plot->GetYaxis()->SetTitleOffset(1.09);

  pull_plot->GetYaxis()->SetLabelFont(42);
  pull_plot->GetYaxis()->SetLabelSize(0.13);
  pull_plot->GetYaxis()->SetLabelOffset(0.005);

  pull_plot->GetYaxis()->SetNdivisions(305);

  p2->cd();
  pull_plot->Draw();

  d.SaveAs("./results/mass_fit_B0.pdf");
  d.SaveAs("./results/mass_fit_B0.gif");

}


void build_pdf(RooWorkspace& w, std::string choice){

  cout << "choice " << choice << endl;

  RooRealVar Bmass = *(w.var("Bmass"));
  RooDataSet* data = (RooDataSet*) w.data("data");

  double mass_peak = 5.280;
 
  //NORMALIZATIONS	
  double n_signal_initial = data->sumEntries(TString::Format("abs(Bmass-%g)<0.05",mass_peak)) - data->sumEntries(TString::Format("abs(Bmass-%g)<0.10&&abs(Bmass-%g)>0.05",mass_peak,mass_peak));
  double n_combinatorial_initial = data->sumEntries() - n_signal_initial;
  RooRealVar n_signal("n_signal","n_signal",n_signal_initial,0.,data->sumEntries());
  RooRealVar n_combinatorial("n_combinatorial","n_combinatorial",n_combinatorial_initial,0.,data->sumEntries());

  //SIGNAL (correctly tagged)
  //sum of two gaussians
  RooRealVar mean("mean","mean",mass_peak,mass_peak-0.1,mass_peak+0.1);
  RooRealVar sigma1("sigma1","sigma1",0.02,0.005,0.05);
  RooGaussian signal1("signal1","signal_gauss1",Bmass,mean,sigma1);
  RooRealVar sigma2("sigma2","sigma2",0.01,0.005,0.05);
  RooGaussian signal2("signal2","signal_gauss2",Bmass,mean,sigma2);
  RooRealVar cofs("cofs", "cofs", 0.3, 0., 1.);
  RooAddPdf pdf_c_signal("pdf_c_signal", "pdf_c_signal", RooArgList(signal1,signal2),cofs);
 
  //BACKGROUND
  //exponential
  RooRealVar lambda("lambda","lambda",-2.,-5.,1.0);
  RooExponential pdf_m_combinatorial("pdf_m_combinatorial", "pdf_m_combinatorial", Bmass, lambda);

  //K Pi SWAP COMPONENT
  RooRealVar sigma_swapped1("sigma_swapped1","sigma_swapped1", 0.1133, 0.010, 0.150);
  RooRealVar sigma_swapped2("sigma_swapped2","sigma_swapped2", 0.01529, 0.010, 0.150);
  RooRealVar sigma_swapped3("sigma_swapped3","sigma_swapped3", 0.0424, 0.010, 0.150);
  RooRealVar alpha1("alpha1","alpha1", 1.78, -20., 20.);
  RooRealVar alpha2("alpha2","alpha2", 0.150, -20., 20.);
  RooRealVar alpha3("alpha3","alpha3", -6.802, -20., 20.);
  RooRealVar n1_parameter("n1_parameter", "n1_parameter", 32., 0., 300.);
  RooRealVar n2_parameter("n2_parameter", "n2_parameter", 98., 0., 300.);
  RooRealVar n3_parameter("n3_parameter", "n3_parameter", 179., 0., 300.);
  RooRealVar r1("r1","r1", 0.249, 0.0, 1.0); 
  RooRealVar r2("r2","r2", 0.3922, 0.0, 1.0);  

  sigma_swapped1.setConstant(kTRUE);
  sigma_swapped2.setConstant(kTRUE);
  sigma_swapped3.setConstant(kTRUE);
  alpha1.setConstant(kTRUE);
  alpha2.setConstant(kTRUE);
  alpha3.setConstant(kTRUE);
  n1_parameter.setConstant(kTRUE);
  n2_parameter.setConstant(kTRUE);
  n3_parameter.setConstant(kTRUE);
  r1.setConstant(kTRUE);
  r2.setConstant(kTRUE);

  RooCBShape swapped1("swapped1","swapped1", Bmass, mean, sigma_swapped1, alpha1, n1_parameter);
  RooCBShape swapped2("swapped2","swapped2", Bmass, mean, sigma_swapped2, alpha2, n2_parameter);
  RooCBShape swapped3("swapped3","swapped3", Bmass, mean, sigma_swapped3, alpha3, n3_parameter);

  RooAddPdf k_pi_swap("k_pi_swap","k_pi_swap", RooArgSet(swapped1,swapped2,swapped3), RooArgSet(r1,r2));

  // FULL MODEL
  RooRealVar f_swap("f_swap","f_swap", 0.5);
  RooAddPdf pdf_t_signal("pdf_t_signal","pdf_t_signal",RooArgList(k_pi_swap,pdf_c_signal),RooArgList(f_swap));

  if (choice == "nominal"){
    RooAddPdf model("model","model", RooArgList(pdf_t_signal, pdf_m_combinatorial), RooArgList(n_signal, n_combinatorial));
    w.import(model);}
}


void cut_tree(TString f_input){
  TFile* fin_data = new TFile(f_input);
  TTree* t1_data = (TTree*)fin_data->Get("ntKstar");

  float Bmass;
  float By;
  float Bpt;
  float Btrk1pt;
  float Btrk2pt;
  float Btrk1eta;
  float Btrk2eta;
  float Btrk1pterr;
  float Btrk2pterr;
  float Bchi2cl;
  float Bsvpvdistance;
  float Bsvpvdiserr;
  float Bsvpvdistance_2D;
  float Bsvpvdiserr_2D;
  float Bmumumass;
  float Bmu1eta;
  float Bmu2eta;
  float Bmu1pt;
  float Bmu2pt;
  float Bmu1dxyPV;
  float Bmu2dxyPV;
  float Bmu1dzPV;
  float Bmu2dzPV;
  float Bd0;
  float Bd0Err;
  float Bdtheta;
  float Balpha;
  float Btrk1Dz1;
  float Btrk2Dz1;
  float Btrk1DzError1;
  float Btrk2DzError1;
  float Btrk1Dxy1;
  float Btrk2Dxy1;
  float Btrk1DxyError1;
  float Btrk2DxyError1;
  float Bmumueta;
  float Bmumuphi;
  float Bmumupt;

  t1_data->SetBranchAddress("Bmass", &Bmass);
  t1_data->SetBranchAddress("By", &By);
  t1_data->SetBranchAddress("Bpt", &Bpt);
  t1_data->SetBranchAddress("Btrk1Pt", &Btrk1pt);
  t1_data->SetBranchAddress("Btrk2Pt", &Btrk2pt);
  t1_data->SetBranchAddress("Btrk1Eta", &Btrk1eta);
  t1_data->SetBranchAddress("Btrk2Eta", &Btrk2eta);
  t1_data->SetBranchAddress("Btrk1PtErr", &Btrk1pterr);
  t1_data->SetBranchAddress("Btrk2PtErr", &Btrk2pterr);
  t1_data->SetBranchAddress("Bchi2cl", &Bchi2cl);
  t1_data->SetBranchAddress("BsvpvDistance", &Bsvpvdistance);
  t1_data->SetBranchAddress("BsvpvDisErr", &Bsvpvdiserr);
  t1_data->SetBranchAddress("BsvpvDistance_2D", &Bsvpvdistance_2D);
  t1_data->SetBranchAddress("BsvpvDisErr_2D", &Bsvpvdiserr_2D);
  t1_data->SetBranchAddress("Bmumumass", &Bmumumass);
  t1_data->SetBranchAddress("Bmu1eta", &Bmu1eta);
  t1_data->SetBranchAddress("Bmu2eta", &Bmu2eta);
  t1_data->SetBranchAddress("Bmu1pt", &Bmu1pt);
  t1_data->SetBranchAddress("Bmu2pt", &Bmu2pt);
  t1_data->SetBranchAddress("Bmu1dxyPV", &Bmu1dxyPV);
  t1_data->SetBranchAddress("Bmu2dxyPV", &Bmu2dxyPV);
  t1_data->SetBranchAddress("Bmu1dzPV", &Bmu1dzPV);
  t1_data->SetBranchAddress("Bmu2dzPV", &Bmu2dzPV);
  t1_data->SetBranchAddress("Bd0", &Bd0);
  t1_data->SetBranchAddress("Bd0Err", &Bmu2dzPV);
  t1_data->SetBranchAddress("Bdtheta", &Bdtheta);
  t1_data->SetBranchAddress("Balpha", &Balpha);
  t1_data->SetBranchAddress("Btrk1Dz1", &Btrk1Dz1);
  t1_data->SetBranchAddress("Btrk2Dz1", &Btrk2Dz1);
  t1_data->SetBranchAddress("Btrk1DzError1", &Btrk1DzError1);
  t1_data->SetBranchAddress("Btrk2DzError1", &Btrk2DzError1);
  t1_data->SetBranchAddress("Btrk1Dxy1", &Btrk1Dxy1);
  t1_data->SetBranchAddress("Btrk2Dxy1", &Btrk2Dxy1);
  t1_data->SetBranchAddress("Btrk1DxyError1", &Btrk1DxyError1);
  t1_data->SetBranchAddress("Btrk2DxyError1", &Btrk2DxyError1);
  t1_data->SetBranchAddress("Bmumueta", &Bmumueta);
  t1_data->SetBranchAddress("Bmumuphi", &Bmumuphi);
  t1_data->SetBranchAddress("Bmumupt", &Bmumupt);

  TFile* f_tree = new TFile("./results/cutted_tree.root", "recreate");
  f_tree->cd();

  TTree *t_new = new TTree("ntKstar", "ntKstar");

  float bmass;
  float by;
  float bpt;
  float btrk1pt;
  float btrk2pt;
  float btrk1eta;
  float btrk2eta;
  float btrk1pterr;
  float btrk2pterr;
  float bchi2cl;
  float bsvpvdistance;
  float bsvpvdiserr;
  float bsvpvdistance_2D;
  float bsvpvdiserr_2D;
  float bmumumass;
  float bmu1eta;
  float bmu2eta;
  float bmu1pt;
  float bmu2pt;
  float bmu1dxyPV;
  float bmu2dxyPV;
  float bmu1dzPV;
  float bmu2dzPV;
  float bd0;
  float bd0Err;
  float bdtheta;
  float balpha;
  float btrk1Dz1;
  float btrk2Dz1;
  float btrk1DzError1;
  float btrk2DzError1;
  float btrk1Dxy1;
  float btrk2Dxy1;
  float btrk1DxyError1;
  float btrk2DxyError1;
  float bmumueta;
  float bmumuphi;
  float bmumupt;

  t_new->Branch("Bmass", &bmass);
  t_new->Branch("By", &by);
  t_new->Branch("Bpt", &bpt);
  t_new->Branch("Btrk1Pt", &btrk1pt);
  t_new->Branch("Btrk2Pt", &btrk2pt);
  t_new->Branch("Btrk1Eta", &btrk1eta);
  t_new->Branch("Btrk2Eta", &btrk2eta);
  t_new->Branch("Btrk1PtErr", &btrk1pterr);
  t_new->Branch("Btrk2PtErr", &btrk2pterr);
  t_new->Branch("Bchi2cl", &bchi2cl);
  t_new->Branch("BsvpvDistance", &bsvpvdistance);
  t_new->Branch("BsvpvDisErr", &bsvpvdiserr);
  t_new->Branch("BsvpvDistance_2D", &bsvpvdistance_2D);
  t_new->Branch("BsvpvDisErr_2D", &bsvpvdiserr_2D);
  t_new->Branch("Bmumumass", &bmumumass);
  t_new->Branch("Bmu1eta", &bmu1eta);
  t_new->Branch("Bmu2eta", &bmu2eta);
  t_new->Branch("Bmu1pt", &bmu1pt);
  t_new->Branch("Bmu2pt", &bmu2pt);
  t_new->Branch("Bmu1dxyPV", &bmu1dxyPV);
  t_new->Branch("Bmu2dxyPV", &bmu2dxyPV);
  t_new->Branch("Bmu1dzPV", &bmu1dzPV);
  t_new->Branch("Bmu2dzPV", &bmu2dzPV);
  t_new->Branch("Bd0", &bd0);
  t_new->Branch("Bd0Err", &bmu2dzPV);
  t_new->Branch("Bdtheta", &bdtheta);
  t_new->Branch("Balpha", &balpha);
  t_new->Branch("Btrk1Dz1", &btrk1Dz1);
  t_new->Branch("Btrk2Dz1", &btrk2Dz1);
  t_new->Branch("Btrk1DzError1", &btrk1DzError1);
  t_new->Branch("Btrk2DzError1", &btrk2DzError1);
  t_new->Branch("Btrk1Dxy1", &btrk1Dxy1);
  t_new->Branch("Btrk2Dxy1", &btrk2Dxy1);
  t_new->Branch("Btrk1DxyError1", &btrk1DxyError1);
  t_new->Branch("Btrk2DxyError1", &btrk2DxyError1);
  t_new->Branch("Bmumueta", &bmumueta);
  t_new->Branch("Bmumuphi", &bmumuphi);
  t_new->Branch("Bmumupt", &bmumupt);

  float real_bmass = 5.280;

  for(int i = 0; i < t1_data->GetEntries(); i++){
    t1_data->GetEntry(i);

    if(abs(Bmass-real_bmass) < 0.280){
      bmass = Bmass;
      by = By;
      bpt = Bpt;
      btrk1pt = Btrk1pt;
      btrk2pt = Btrk2pt;
      btrk1eta = Btrk1eta;
      btrk2eta = Btrk2eta;
      btrk1pterr = Btrk1pterr;
      btrk2pterr = Btrk2pterr;
      bchi2cl = Bchi2cl;
      bsvpvdistance = Bsvpvdistance;
      bsvpvdiserr = Bsvpvdiserr;
      bsvpvdistance_2D = Bsvpvdistance_2D;
      bsvpvdiserr_2D = Bsvpvdiserr_2D;
      bmumumass = Bmumumass;
      bmu1eta = Bmu1eta;
      bmu2eta = Bmu2eta;
      bmu1pt = Bmu1pt;
      bmu2pt = Bmu2pt;
      bmu1dxyPV = Bmu1dxyPV;
      bmu2dxyPV = Bmu2dxyPV;
      bmu1dzPV = Bmu1dzPV;
      bmu2dzPV = Bmu2dzPV;
      bd0 = Bd0;
      bd0Err = Bd0Err;
      bdtheta = Bdtheta;
      balpha = Balpha;
      btrk1Dz1 = Btrk1Dz1;
      btrk2Dz1 = Btrk1Dz1;
      btrk1DzError1 = Btrk1DzError1;
      btrk2DzError1 = Btrk2DzError1;
      btrk1Dxy1 = Btrk1Dxy1;
      btrk2Dxy1 = Btrk2Dxy1;
      btrk1DxyError1 = Btrk1DxyError1;
      btrk2DxyError1 = Btrk2DxyError1;
      bmumueta = Bmumueta;
      bmumuphi = Bmumuphi;
      bmumupt = Bmumupt;
      t_new->Fill();
    }
  }
 f_tree->Write();
 f_tree->Close();
}

void read_data(RooWorkspace& w, TString f_input){
  TFile* fin_data = new TFile(f_input);
  TTree* t1_data = (TTree*)fin_data->Get("ntKstar");

  RooArgList arg_list ("arg_list");

  arg_list.add(*(w.var("Bmass")));
  arg_list.add(*(w.var("By")));
  arg_list.add(*(w.var("Bpt")));
  arg_list.add(*(w.var("Btrk1Pt")));
  arg_list.add(*(w.var("Btrk2Pt")));
  arg_list.add(*(w.var("Btrk1Eta")));
  arg_list.add(*(w.var("Btrk2Eta")));
  arg_list.add(*(w.var("Btrk1PtErr")));
  arg_list.add(*(w.var("Btrk2PtErr")));
  arg_list.add(*(w.var("Bchi2cl")));
  arg_list.add(*(w.var("BsvpvDistance")));
  arg_list.add(*(w.var("BsvpvDisErr")));
  arg_list.add(*(w.var("BsvpvDistance_2D")));
  arg_list.add(*(w.var("BsvpvDisErr_2D")));
  arg_list.add(*(w.var("Bmumumass")));
  arg_list.add(*(w.var("Bmu1eta")));
  arg_list.add(*(w.var("Bmu2eta")));
  arg_list.add(*(w.var("Bmu1pt")));
  arg_list.add(*(w.var("Bmu2pt")));
  arg_list.add(*(w.var("Bmu1dxyPV")));
  arg_list.add(*(w.var("Bmu2dxyPV")));
  arg_list.add(*(w.var("Bmu1dzPV")));
  arg_list.add(*(w.var("Bmu2dzPV")));
  arg_list.add(*(w.var("Bd0")));
  arg_list.add(*(w.var("Bd0Err")));
  arg_list.add(*(w.var("Bdtheta")));
  arg_list.add(*(w.var("Balpha")));
  arg_list.add(*(w.var("Btrk1Dz1")));
  arg_list.add(*(w.var("Btrk2Dz1")));
  arg_list.add(*(w.var("Btrk1DzError1")));
  arg_list.add(*(w.var("Btrk2DzError1")));
  arg_list.add(*(w.var("Btrk1Dxy1")));
  arg_list.add(*(w.var("Btrk2Dxy1")));
  arg_list.add(*(w.var("Btrk1DxyError1")));
  arg_list.add(*(w.var("Btrk2DxyError1")));
  arg_list.add(*(w.var("Bmumueta")));
  arg_list.add(*(w.var("Bmumuphi")));
  arg_list.add(*(w.var("Bmumupt")));


  RooDataSet* data = new RooDataSet("data","data",t1_data,arg_list);
  w.import(*data);
 
}


void set_up_workspace_variables(RooWorkspace& w){

  double mass_min, mass_max;
  double y_min, y_max;
  double pt_min, pt_max;
  double trk1pt_min, trk1pt_max;
  double trk2pt_min, trk2pt_max;
  double trk1eta_min, trk1eta_max;
  double trk2eta_min, trk2eta_max;
  double trk1pterr_min, trk1pterr_max;
  double trk2pterr_min, trk2pterr_max;
  double chi2cl_min, chi2cl_max;
  double svpvDistance_min, svpvDistance_max;
  double svpvDistanceErr_min, svpvDistanceErr_max;
  double svpvDistance2D_min, svpvDistance2D_max;
  double svpvDistanceErr2D_min, svpvDistanceErr2D_max;
  double mumumass_min, mumumass_max;
  double mu1eta_min, mu1eta_max;
  double mu2eta_min, mu2eta_max;
  double mu1pt_min, mu1pt_max;
  double mu2pt_min, mu2pt_max;
  double mu1dxyPV_min, mu1dxyPV_max;
  double mu2dxyPV_min, mu2dxyPV_max;
  double mu1dzPV_min, mu1dzPV_max;
  double mu2dzPV_min, mu2dzPV_max;
  double d0_min, d0_max;
  double d0err_min, d0err_max;
  double dtheta_min, dtheta_max;
  double alpha_min, alpha_max;
  double trk1Dz1_min, trk1Dz1_max;
  double trk2Dz1_min, trk2Dz1_max;
  double trk1DzError1_min, trk1DzError1_max;
  double trk2DzError1_min, trk2DzError1_max;
  double trk1Dxy1_min, trk1Dxy1_max;
  double trk2Dxy1_min, trk2Dxy1_max;
  double trk1DxyError1_min, trk1DxyError1_max;
  double trk2DxyError1_min, trk2DxyError1_max;
  double mumueta_min, mumueta_max;
  double mumuphi_min, mumuphi_max;
  double mumupt_min, mumupt_max;


  mass_min = 5.0;
  mass_max = 5.6;

  y_min = -2.4;
  y_max = 2.4;

  pt_min = 0.;
  pt_max = 100.;

  trk1pt_min = 0.;
  trk1pt_max = 25.;

  trk2pt_min = 0.;
  trk2pt_max = 25.;

  trk1eta_min = -2.4;
  trk1eta_max = 2.4;

  trk2eta_min = -2.4;
  trk2eta_max = 2.4;

  trk1pterr_min = 0.;
  trk1pterr_max = 0.4;

  trk2pterr_min = 0.;
  trk2pterr_max = 0.4;

  chi2cl_min = 0.;
  chi2cl_max = 1.05;

  svpvDistance_min = 0.;
  svpvDistance_max = 25.;

  svpvDistanceErr_min = 0.;
  svpvDistanceErr_max = 0.2;

  svpvDistance2D_min = 0.;
  svpvDistance2D_max = 1. ;

  svpvDistanceErr2D_min = 0.;
  svpvDistanceErr2D_max = 0.025;

  mumumass_min = 2.95;
  mumumass_max = 3.25;

  mu1eta_min = -2.4;
  mu1eta_max = 2.4;

  mu2eta_min = -2.4;
  mu2eta_max = 2.4;

  mu1pt_min = 0.;
  mu1pt_max = 50.;

  mu2pt_min = 0.;
  mu2pt_max = 50.;

  mu1dxyPV_min = -0.15;
  mu1dxyPV_max = 0.15;

  mu2dxyPV_min = -0.15;
  mu2dxyPV_max = 0.15;

  mu1dzPV_min = -20.;
  mu1dzPV_max = 20.;

  mu2dzPV_min = -20.;
  mu2dzPV_max = 20.;

  d0_min = 0.;
  d0_max = 0.5;

  d0err_min = 0.;
  d0err_max = 0.0002;

  dtheta_min = 0.;
  dtheta_max = 3.1;

  alpha_min = 0.;
  alpha_max = 3.1;

  trk1Dz1_min = -20.;
  trk1Dz1_max = 20.;

  trk2Dz1_min = -20.;
  trk2Dz1_max = 20.;

  trk1DzError1_min = 0.;
  trk1DzError1_max = 0.5;

  trk2DzError1_min = 0.;
  trk2DzError1_max = 0.5;

  trk1Dxy1_min = -0.3;
  trk1Dxy1_max = 0.3;

  trk2Dxy1_min = -0.3;
  trk2Dxy1_max = 0.3;

  trk1DxyError1_min = 0.;
  trk1DxyError1_max = 0.1;

  trk2DxyError1_min = 0.;
  trk2DxyError1_max = 0.1;

  mumueta_min = -6.;
  mumueta_max = 6.;

  mumuphi_min = -3.1;
  mumuphi_max = 3.1;

  mumupt_min = 0.;
  mumupt_max = 60.;


  RooRealVar Bmass("Bmass","Bmass",mass_min,mass_max);
  RooRealVar By("By","By",y_min,y_max);
  RooRealVar Bpt("Bpt","Bpt",pt_min,pt_max);
  RooRealVar Btrk1Pt("Btrk1Pt","Btrk1Pt",trk1pt_min,trk1pt_max);
  RooRealVar Btrk2Pt("Btrk2Pt","Btrk2Pt",trk2pt_min,trk2pt_max);
  RooRealVar Btrk1Eta("Btrk1Eta","Btrk1Eta",trk1eta_min,trk1eta_max);
  RooRealVar Btrk2Eta("Btrk2Eta","Btrk2Eta",trk2eta_min,trk2eta_max);
  RooRealVar Btrk1PtErr("Btrk1PtErr","Btrk1PtErr",trk1pterr_min,trk1pterr_max);
  RooRealVar Btrk2PtErr("Btrk2PtErr","Btrk2PtErr",trk2pterr_min,trk2pterr_max);
  RooRealVar Bchi2cl("Bchi2cl","Bchi2cl",chi2cl_min,chi2cl_max);
  RooRealVar BsvpvDistance("BsvpvDistance", "BsvpvDistance", svpvDistance_min, svpvDistance_max);
  RooRealVar BsvpvDistance_2D("BsvpvDistance_2D", "BsvpvDistance_2D", svpvDistance2D_min, svpvDistance2D_max);
  RooRealVar BsvpvDisErr("BsvpvDisErr", "BsvpvDisErr", svpvDistanceErr_min, svpvDistanceErr_max);
  RooRealVar BsvpvDisErr_2D("BsvpvDisErr_2D", "BsvpvDisErr_2D", svpvDistanceErr2D_min, svpvDistanceErr2D_max);
  RooRealVar Bmumumass("Bmumumass", "Bmumumass", mumumass_min, mumumass_max);
  RooRealVar Bmu1eta("Bmu1eta","Bmu1eta",mu1eta_min,mu1eta_max);
  RooRealVar Bmu2eta("Bmu2eta","Bmu2eta",mu2eta_min,mu2eta_max);
  RooRealVar Bmu1pt("Bmu1pt","Bmu1pt",mu1pt_min,mu1pt_max);
  RooRealVar Bmu2pt("Bmu2pt","Bmu2pt",mu2pt_min,mu2pt_max);
  RooRealVar Bmu1dxyPV("Bmu1dxyPV", "Bmu1dxyPV", mu1dxyPV_min, mu1dxyPV_max);
  RooRealVar Bmu2dxyPV("Bmu2dxyPV", "Bmu2dxyPV", mu2dxyPV_min, mu2dxyPV_max);
  RooRealVar Bmu1dzPV("Bmu1dzPV", "Bmu1dzPV", mu1dzPV_min, mu1dzPV_max);
  RooRealVar Bmu2dzPV("Bmu2dzPV", "Bmu2dzPV", mu2dzPV_min, mu2dzPV_max);
  RooRealVar Bd0("Bd0", "Bd0", d0_min, d0_max);
  RooRealVar Bd0Err("Bd0Err", "Bd0Err", d0err_min, d0err_max);
  RooRealVar Bdtheta("Bdtheta", "Bdtheta", dtheta_min, dtheta_max);
  RooRealVar Balpha("Balpha", "Balpha", alpha_min, alpha_max);
  RooRealVar Btrk1Dz1("Btrk1Dz1","Btrk1Dz1",trk1Dz1_min,trk1Dz1_max);
  RooRealVar Btrk2Dz1("Btrk2Dz1","Btrk2Dz1",trk2Dz1_min,trk2Dz1_max);
  RooRealVar Btrk1DzError1("Btrk1DzError1","Btrk1DzError1",trk1DzError1_min,trk1DzError1_max);
  RooRealVar Btrk2DzError1("Btrk2DzError1","Btrk2DzError1",trk2DzError1_min,trk2DzError1_max);
  RooRealVar Btrk1Dxy1("Btrk1Dxy1","Btrk1Dxy1",trk1Dxy1_min,trk1Dxy1_max);
  RooRealVar Btrk2Dxy1("Btrk2Dxy1","Btrk2Dxy1",trk2Dxy1_min,trk2Dxy1_max);
  RooRealVar Btrk1DxyError1("Btrk1DxyError1","Btrk1DxyError1",trk1DxyError1_min,trk1DxyError1_max);
  RooRealVar Btrk2DxyError1("Btrk2DxyError1","Btrk2DxyError1",trk2DxyError1_min,trk2DxyError1_max);
  RooRealVar Bmumueta("Bmumueta", "Bmumueta", mumueta_min, mumueta_max);
  RooRealVar Bmumuphi("Bmumuphi", "Bmumuphi", mumuphi_min, mumuphi_max);
  RooRealVar Bmumupt("Bmumupt", "Bmumupt", mumupt_min, mumupt_max);


  w.import(Bmass);
  w.import(By);
  w.import(Bpt);
  w.import(Btrk1Pt);
  w.import(Btrk2Pt);
  w.import(Btrk1Eta);
  w.import(Btrk2Eta);
  w.import(Btrk1PtErr);
  w.import(Btrk2PtErr);
  w.import(Bchi2cl);
  w.import(BsvpvDistance);
  w.import(BsvpvDistance_2D);
  w.import(BsvpvDisErr);
  w.import(BsvpvDisErr_2D);
  w.import(Bmumumass);
  w.import(Bmu1eta);
  w.import(Bmu2eta);
  w.import(Bmu1pt);
  w.import(Bmu2pt);
  w.import(Bmu1dxyPV);
  w.import(Bmu2dxyPV);
  w.import(Bmu1dzPV);
  w.import(Bmu2dzPV);
  w.import(Bd0);
  w.import(Bd0Err);
  w.import(Bdtheta);
  w.import(Balpha);
  w.import(Btrk1Dz1);
  w.import(Btrk2Dz1);
  w.import(Btrk1DzError1);
  w.import(Btrk2DzError1);
  w.import(Btrk1Dxy1);
  w.import(Btrk2Dxy1);
  w.import(Btrk1DxyError1);
  w.import(Btrk2DxyError1);
  w.import(Bmumueta);
  w.import(Bmumuphi);
  w.import(Bmumupt);
}










































