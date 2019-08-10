//root
#include <TLine.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TString.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TMath.h>
#include <TLegend.h>
#include <TCut.h>
#include <THStack.h>
#include <TGaxis.h>
#include <TF1.h>
#include <TError.h> // root verbosity level
#include <TApplication.h>

//C, C++
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

void plot_ph(){
 	
	// ROOT GUI
	bool is_interactive = 0;
	// store results
	bool store_result = 1;

	TApplication * ROOTapp;
	if (is_interactive){ROOTapp = new TApplication("ROOT interactive", &argc, argv);}

	/***** 
	INITIALIZE 
	*****/

	gStyle->SetErrorX(0);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	// gStyle->SetOptStat("e");
	// gStyle->SetOptFit(1);
	gStyle->SetStatW(0.15); gStyle->SetStatH(0.1); // stats box size
	gStyle->SetStatX(0.9); gStyle->SetStatY(0.9); // stats box position
	gStyle->SetGridStyle(3);
	gStyle->SetGridWidth(1);
	gStyle->SetGridColor(16);
	// gStyle->SetImageScaling(30.);
	gStyle->SetLineScalePS(1);

	// parse input arguments
	// string run_name = "/Users/julianschliwinski/Studium/CERN_Summer_2019/analysis/cosmics_piTrigger/comb_root/cosmics_7to16_32to50_1325CS_100ns_isCalib0_newBL.root";
	// int run_nr = atoi(argv[2]);

	TString dir = "./comb_root/";
	TString f = "cosmics_7to16_32to50_1325CS_100ns_isCalib0_newBL";

	TFile* file0 = new TFile(dir+f+".root");


	TTree* tree;
	file0->GetObject("ntuple",tree);

	TCanvas * C1 = new TCanvas("tree overview","",700,500);
	TH1D * h = new TH1D("h","h",300,-0.1,1.5);
	h->SetTitle("SiPM pulse-height spectrum");
	h->GetXaxis()->SetTitle("pulse-height [arb. unit]");
	h->GetYaxis()->SetTitle("entries");
	h->SetLineWidth(2);
	h->SetLineColorAlpha(kBlack,0.7);
    h->SetMarkerStyle(7);

	tree->Draw("charge_alt>>h","bl_value<-0.001 && run_nr>=34 && run_nr<=39","goff");
	// h->Scale(1/h->Integral(),"");
	// h->Scale(1/h->Integral(),"width");
	// h->Scale(1/h->GetEntries());
	h->Draw("HISTE");

	if (is_interactive){ROOTapp->Run();}
}