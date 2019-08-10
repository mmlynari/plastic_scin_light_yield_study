//root
#include <TLine.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2D.h>
#include "TFitResult.h"
#include <THStack.h>
#include <TStyle.h>
#include <TString.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <Math/VavilovAccuratePdf.h>
#include <TLegend.h>
#include <TPaveStats.h>
#include <TCut.h>
#include <THStack.h>
#include <TGaxis.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TError.h> // root verbosity level

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

/*******************
__ FUNCTIONS ______
*******************/

// Vavilov fit function
struct Vavilov_Func { 
   Vavilov_Func() {}

   double operator() (const double *x, const double *p) { 
      double kappa = p[0]; 
      double beta2 = p[1];
      return p[4]*( pdf.Pdf( (x[0]-p[2])/p[3], kappa,beta2) );
   }

   ROOT::Math::VavilovAccurate pdf; 
};


void vav_fit_nsetup(int n_setup, vector<string> run_name_list, vector<double> bl_cut_list, vector< vector<int> > run_nr_list, vector<string> label_list, TH1D * h_vec[n_setup],TF1 * v_func_vec[n_setup], double* v_max, double* v_mean, double* v_rchi2, double* h_max_vec, int* h_entries)
{
  // fit light yield distribution for each data set (number of data sets -> n_setup)
  // data set definded by run range
  // two pre-fits (Gaussian, Landau) to get starting values for Vavilov fit
  for (int i = 0; i < n_setup; ++i)
  {
    /***** 
    __ SETTINGS ___________________________
    *****/

    string filename = run_name_list[i];     // root file name
    double bl_cut = bl_cut_list[i];         // baseline cut to filter electronic noice
    int run_nr_l = run_nr_list[i][0];       // run range, lower bound
    int run_nr_u = run_nr_list[i][1];       // run range, upper bound
    string run_label = label_list[i];       // run description, legend label


    string path = "./comb_root/"+filename+".root";  // root file directory
    string target_path = "./comb_root/vavilov/";    // target directory for exports
    
    // histogram settings
    double xmin = -9.5;
    double xmax = 100.5;
    int nBins = (xmax-xmin)*1; 

    /***** 
    __ LIGHT YIELD DISTRIBUTION in run range ___________________________
    *****/

    // open tree
    TFile* file = new TFile(path.c_str());
    TTree* tree;
    file->GetObject("ntuple",tree);

    // draw histogram
    TCanvas* C1 = new TCanvas("C1","GP fit",900,700);
    TH1D* h = new TH1D("h","",nBins,xmin,xmax);

    TString h_title;
    h_title.Form("Vavilov fit: %s",run_label.c_str());
    h->SetTitle(h_title);
    h->SetTitleSize(0.3);
    h->SetLineColorAlpha(kBlack,0.7);
    h->SetMarkerStyle(7);
    h->SetMarkerColorAlpha(kBlack,0.6);
    h->GetYaxis()->SetTitle("entries");
    h->GetXaxis()->SetTitle("number-of-photoelectrons [N_{pe}]");
    h->GetXaxis()->SetTitleSize(0.035);
    h->GetXaxis()->SetTitleOffset(1.3);

    // cuts:
    // exclude pedestal region < 0.5 N_pe
    // exclude high baseline values
    // select run range
    TString cut;
    cut.Form("charge_alt>0.5 && bl_value<%f && run_nr>=%d && run_nr<=%d ",bl_cut,run_nr_l,run_nr_u);
    tree->Draw(Form("charge_alt>>h"),cut,"HISTE");

    // save histogram outside loop
    h_vec[i] = (TH1D*)h->Clone(Form("h%d",i));

    /***** 
    __ PRE FIT - Single Gauss fits___________________________
    *****/

    Double_t par_single[3]; // to store fit results

    float range = 4.5; // symmetric fit range around hist maximum
    float h_max = h->GetBinCenter( h->GetMaximumBin() );
    TF1 * peak_single = new TF1("peak","gaus",h_max-range,h_max+range);
    h->Fit("peak","RMQ+");
    peak_single->GetParameters(&par_single[0]);

    /***** 
    __ LANDAU FIT ___________________________
    *****/

    TF1 *f = new TF1("fit_l","landau",-1,30);
    f->SetLineColor(3);
    f->SetNpx(1000); // draw function with high resolution 
    f->SetParameters(&par_single[0]);

    h->Fit("fit_l","RMQ");
    // h->Draw("sameFUNC");
    
     // store fit results
    float lan_mpv = f->GetParameter(1);
    float lan_mpv_err = f->GetParError(1);
    float lan_rchi2 = f->GetChisquare()/f->GetNDF();

    /***** 
    __ Vavilov FIT ___________________________
    *****/

    Vavilov_Func * func = new Vavilov_Func();
    TF1 * f1 = new TF1("f1",func, 1,60,5,"Vavilov_Func");
    f1->SetParNames("kappa","beta2","mean","sigma","Amp"); 
    f1->SetParameters(0.01,0.03,f->GetParameter(1),f->GetParameter(2),h->GetEntries());
    f1->SetParLimits(0,0.01,10); // Vavilov model valid only in this kappa range
    // f1->SetParLimits(1,0,1);
    f1->FixParameter(1,1); // fix beta2 to "light speed", no dependence this parameter observed
    f1->SetNpx(1000); // draw function with high resolution 

    h->Fit("f1","RQM");

    h->Draw("sameFUNC");

    double h_mean =  h->GetMean();
    double h_max_y = h->GetMaximum();
    // store fit results
    double vav_mean = f1->GetParameter(2);
    double vav_max = f1->GetMaximumX(vav_mean-5,vav_mean);
    double vav_rchi2 =  f1->GetChisquare()/f1->GetNDF();

    // draw lines for MPV and mean
    TLine * ln_vav_mean = new TLine(vav_mean,0,vav_mean,h_max_y);
    ln_vav_mean->SetLineColor(8);
    ln_vav_mean->SetLineStyle(2);
    ln_vav_mean->Draw("same");

    TLine * ln_vav_max = new TLine(vav_max,0,vav_max,h_max_y);
    ln_vav_max->SetLineColor(9);
    ln_vav_max->SetLineStyle(9);
    ln_vav_max->Draw("same");

    h->Draw("sameFUNC");

    // custom histogram legend with fit results
    TLegend *h_vav_leg = new TLegend(0.50,0.60,0.9,0.9);
    h_vav_leg->AddEntry(h,Form("#bf{data}"),"elp");
    h_vav_leg->AddEntry((TObject*)0,Form("entries = %1.f",h->GetEntries()),"");
    h_vav_leg->AddEntry(f1,Form("Vavilov fit"),"l");
    h_vav_leg->AddEntry((TObject*)0,Form("#chi^{2}/ndf = %1.1f",vav_rchi2),"");
    h_vav_leg->AddEntry((TObject*)0,Form("#kappa = %1.3f #pm %1.3f",f1->GetParameter(0),f1->GetParError(0)),"");
    h_vav_leg->AddEntry((TObject*)0,Form("#beta^{2} = %1.f",f1->GetParameter(1)),"");
    h_vav_leg->AddEntry(ln_vav_mean,Form("mean = %1.2f #pm %1.2f",vav_mean,f1->GetParError(2)),"l");
    h_vav_leg->AddEntry((TObject*)0,Form("#sigma = %1.3f #pm %1.3f",f1->GetParameter(3),f1->GetParError(3)),"");
    h_vav_leg->AddEntry((TObject*)0,Form("amplitude = %1.1f #pm %1.1f",f1->GetParameter(4),f1->GetParError(4)),"");
    h_vav_leg->AddEntry(ln_vav_max,Form("MPV = %1.2f",vav_max),"l");

    h_vav_leg->Draw();


    /***** 
    __ Save Result outside LOOP ___________________________
    *****/    
    v_func_vec[i] = f1;
    v_max[i] = vav_max;
    v_mean[i] = vav_mean;
    v_rchi2[i] = vav_rchi2;
    h_max_vec[i] = h_max;
    h_entries[i] = h->GetEntries();

    /***** 
    __ Export Result ___________________________
    *****/

    printf("Data file: %s\n",filename.c_str() );
    printf("runs %d to %d\n",run_nr_l,run_nr_u );
    printf("Landau Fit: MPV: %1.4f ± %1.4f | red. chi2 = %1.2f \n",lan_mpv,lan_mpv_err,lan_rchi2);
    printf("Vavilov Fit: MPV: %1.4f | red. chi2 = %1.2f \n",vav_max,vav_rchi2);
    printf("\n");

    string pdf_filename1 = target_path+filename+Form("_vavilov_fit_%dto_%d.pdf",run_nr_l,run_nr_u);

    gErrorIgnoreLevel = kError; // suppress root terminal output 
    C1->Print(pdf_filename1.c_str());
    gErrorIgnoreLevel = kUnset ; // allow root terminal output
    delete C1;
  }
  // end fit loop
}

// plot multiple data sets (number of data sets = n_select) including Vavilov fits
void plot_multi_dist(int n_select, TString vav_title, string plot_id, double plt_ymax, vector<double> leg_xy, double leg_fsize, vector<string> label_list, vector<int> c_list, vector<int> l_list, vector<int> m_list, TH1D * h_vec[n_select], TF1 * v_func_vec[n_select], double* v_max, double* v_mean, double* v_rchi2, int* h_entries )
{
  // __ INITIALIZE CANVAS & HIST __
  double c_alpha = 0.07;

  double range_l = -2.5;
  double range_u = 55.5;
  double nbins_new = (range_u-range_l)*1;

  TCanvas * C_comb = new TCanvas("Vavilov Overview","", 1500,1000);
  C_comb->SetTopMargin(0.90);
  C_comb->SetBottomMargin(0.12);
  C_comb->SetLeftMargin(0.13);

  TH2D * h_frame = new TH2D("","",nbins_new,range_l,range_u,10,0,plt_ymax);
  h_frame->SetTitle(vav_title);
  h_frame->GetXaxis()->SetTitle("number-of-photoelectrons [N_{pe}] ");
  h_frame->GetXaxis()->SetTitleSize(0.035);
  h_frame->GetXaxis()->SetTitleOffset(1.5);
  h_frame->GetYaxis()->SetTitle("relative probability");
  h_frame->Draw();

  TLegend *leg = new TLegend(leg_xy[0],leg_xy[1],0.9,0.9);
  leg->SetTextSize(leg_fsize);

  // __ DRAW __
  for (int i = 0; i < n_select; ++i)
  {
    // draw normalized histogram
    h_vec[i]->Scale(1/h_vec[i]->GetEntries());
    // h_vec[i]->Scale(1/h_vec[i]->Integral(),"width");
    h_vec[i]->SetLineColor(c_list[i]);
    h_vec[i]->SetLineStyle(l_list[i]);
    h_vec[i]->SetMarkerStyle(m_list[i]);
    h_vec[i]->SetMarkerColor(c_list[i]);
    h_vec[i]->SetFillColorAlpha(c_list[i],c_alpha-0.005*i);
    h_vec[i]->SetTitle(Form("#bf{%s}",label_list[i].c_str()));
    h_vec[i]->Draw("sameHISTEA");
    leg->AddEntry(Form("h%d",i),"","fple");

    // draw normalized Vavilov function
    v_func_vec[i]->SetRange(-5,80);
    v_func_vec[i]->SetNormalized(1);
    v_func_vec[i]->SetLineColor(c_list[i]);
    v_func_vec[i]->Draw("same");
    leg->AddEntry(v_func_vec[i],Form("#splitline{Vavilov fit:}{MPV = %1.2f N_{pe}   mean = %1.2f N_{pe}  #chi^{2} /ndf = %1.2f  entries: %d}",v_max[i],v_mean[i],v_rchi2[i],h_entries[i]),"l");
  }
  leg->Draw();

  // export
  gErrorIgnoreLevel = kError; // suppress root terminal output 
  C_comb->Print(Form("./comb_root/vavilov/vavilov_%s.pdf",plot_id.c_str()));
  gErrorIgnoreLevel = kUnset ; // allow root terminal output

}



/********************
______ MAIN   ______
********************/
int main(int argc, char *argv[])
{
  // ____ style options ____
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetStatW(0.15); gStyle->SetStatH(0.1); // stats box size
  gStyle->SetStatX(0.9); gStyle->SetStatY(0.9); // stats box position
  // gStyle->SetTitleSize(0.3); // title font size
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);
  gStyle->SetGridColor(16);

  /********************
  __ INITIALIZE ______
  ********************/

  // ____ INPUT objects ____
  // ____ DATA initialization ____

  vector<string> run_name_list, label_list;
  vector< vector<int> > run_nr_list;
  vector<int> c_list, m_list, l_list;
  vector<int> c_list_mpv, m_list_mpv, l_list_mpv, w_id_list_mpv;
  vector<double> bl_cut_list, u_gain;

  // root file names
  run_name_list.push_back("cosmics_1to6_015C_20ns_isCalib1_newBL");
  run_name_list.push_back("cosmics_1to6_015C_20ns_isCalib1_newBL");
  run_name_list.push_back("cosmics_7to16_32to50_1325CS_100ns_isCalib1_newBL");
  run_name_list.push_back("cosmics_1to6_015C_20ns_isCalib1_newBL");
  run_name_list.push_back( "cosmics_7to16_32to50_1325CS_100ns_isCalib1_newBL");
  run_name_list.push_back( "cosmics_7to16_32to50_1325CS_100ns_isCalib1_newBL");
  run_name_list.push_back( "cosmics_7to16_32to50_1325CS_100ns_isCalib1_newBL");
  run_name_list.push_back( "cosmics_7to16_32to50_1325CS_100ns_isCalib1_newBL");
  run_name_list.push_back( "cosmics_7to16_32to50_1325CS_100ns_isCalib1_newBL");
  run_name_list.push_back( "cosmics_7to16_32to50_1325CS_100ns_isCalib1_newBL");
  run_name_list.push_back( "cosmics_7to16_32to50_1325CS_100ns_isCalib1_newBL");
  run_name_list.push_back( "cosmics_7to16_32to50_1325CS_100ns_isCalib1_newBL");
  run_name_list.push_back( "cosmics_7to16_32to50_1325CS_100ns_isCalib1_newBL");
  
  // group runs to corresponding setups
  run_nr_list.push_back(vector<int>{4,4});
  run_nr_list.push_back(vector<int>{5,5});
  run_nr_list.push_back(vector<int>{9,9});
  run_nr_list.push_back(vector<int>{6,6});
  run_nr_list.push_back(vector<int>{40,41});
  run_nr_list.push_back(vector<int>{7,8});
  run_nr_list.push_back(vector<int>{32,33});
  run_nr_list.push_back(vector<int>{34,39});
  run_nr_list.push_back(vector<int>{10,11});
  run_nr_list.push_back(vector<int>{47,50});
  run_nr_list.push_back(vector<int>{12,16});
  run_nr_list.push_back(vector<int>{42,43});
  run_nr_list.push_back(vector<int>{44,46});
  
  // legend entries/plot titles
  label_list.push_back("setup 0 b), tile: short 10, fiber: 1.80 m, wrapping: Tyvek, SiPM: 015C");
  label_list.push_back("setup 0 b), tile: short 10, fiber: 0.45 m, wrapping: Tyvek, SiPM: 015C");
  label_list.push_back("setup 0 b), tile: short 10, fiber: 0.45 m, wrapping: Tyvek, SiPM: 1325CS");
  label_list.push_back("setup 3 b), tile: 10, fiber: 0.45 m, wrapping: Tyvek, SiPM: 015C");
  label_list.push_back("setup 3 c), tile: 10, fiber: 0.45 m, wrapping: Aluminum, SiPM: 1325CS");
  label_list.push_back("setup 3 b), tile: 10, fiber: 0.45 m, wrapping: Tyvek, SiPM: 1325CS");
  label_list.push_back("setup 3 a), tile: 10, fiber: 0.45 m, wrapping: Mylar, SiPM: 1325CS");
  label_list.push_back("setup 2 a), tile: 10, fiber: 1.80 m, wrapping: Mylar, SiPM: 1325CS");
  label_list.push_back("setup 4 a), tile: 2#times 1, fiber: 0.45 m, wrapping: Mylar, SiPM: 1325CS");
  label_list.push_back("setup 1 c), tile: 2#times 1, fiber: 1.80 m, wrapping: Aluminum, SiPM: 1325CS");
  label_list.push_back("setup 1 b), tile: 2#times 1, fiber: 1.80 m, wrapping: Tyvek, SiPM: 1325CS");
  label_list.push_back("setup 1 b)*, tile: 2#times 1, fiber: 1.80 m, wrapping: 2#times Tyvek, SiPM: 1325CS");
  label_list.push_back("setup 1 a), tile: 2#times 1, fiber: 1.80 m, wrapping: Mylar, SiPM: 1325CS");
  
  // marker style, line style, color
  m_list.push_back(30); l_list.push_back(1); c_list.push_back(820+9);  //kSpring
  m_list.push_back(31); l_list.push_back(2); c_list.push_back(840+4);  //kTeal
  m_list.push_back(32); l_list.push_back(3); c_list.push_back(416+1);  //kGreen 
  m_list.push_back(29); l_list.push_back(1); c_list.push_back(432-3);  //kCyan  
  m_list.push_back(20); l_list.push_back(2); c_list.push_back(860-5);  //kAzure 
  m_list.push_back(21); l_list.push_back(3); c_list.push_back(860+7);  //kAzure 
  m_list.push_back(22); l_list.push_back(4); c_list.push_back(600+1);  //kBlue 
  m_list.push_back(23); l_list.push_back(5); c_list.push_back(860+3);  //kAzure 
  m_list.push_back(24); l_list.push_back(6); c_list.push_back(632-4);  //kRed 
  m_list.push_back(32); l_list.push_back(7); c_list.push_back(800);  //kOrange
  m_list.push_back(25); l_list.push_back(8); c_list.push_back(800+7);  //kOrange 
  m_list.push_back(26); l_list.push_back(9); c_list.push_back(632+2);  //kRed 
  m_list.push_back(27); l_list.push_back(10); c_list.push_back(900+7);  //kPink 
  
  // wrapping vs. mvp plot: marker style, line style, color
  m_list_mpv.push_back(30); l_list_mpv.push_back(9); c_list_mpv.push_back(840+4);  //kBlue 
  m_list_mpv.push_back(29); l_list_mpv.push_back(9); c_list_mpv.push_back(840+4);  //kBlue 
  m_list_mpv.push_back(29); l_list_mpv.push_back(9); c_list_mpv.push_back(840+4);  //kBlue 
  m_list_mpv.push_back(29); l_list_mpv.push_back(6); c_list_mpv.push_back(600+1);  //kBlue 
  m_list_mpv.push_back(20); l_list_mpv.push_back(1); c_list_mpv.push_back(600+1);  //kBlue 
  m_list_mpv.push_back(20); l_list_mpv.push_back(1); c_list_mpv.push_back(600+1);  //kBlue 
  m_list_mpv.push_back(20); l_list_mpv.push_back(1); c_list_mpv.push_back(600+1);  //kBlue 
  m_list_mpv.push_back(24); l_list_mpv.push_back(1); c_list_mpv.push_back(600+1);  //kBlue 
  m_list_mpv.push_back(21); l_list_mpv.push_back(5); c_list_mpv.push_back(632-4);  //kRed 
  m_list_mpv.push_back(25); l_list_mpv.push_back(5); c_list_mpv.push_back(632-4);  //kRed 
  m_list_mpv.push_back(25); l_list_mpv.push_back(5); c_list_mpv.push_back(632-4);  //kRed 
  m_list_mpv.push_back(30); l_list_mpv.push_back(5); c_list_mpv.push_back(632-4);  //kRed 
  m_list_mpv.push_back(25); l_list_mpv.push_back(5); c_list_mpv.push_back(632-4);  //kRed 
  
  // wrapping vs. mvp plot: wrapping type 1=a,2=b,3=c
  w_id_list_mpv.push_back(2);
  w_id_list_mpv.push_back(2);
  w_id_list_mpv.push_back(2);
  w_id_list_mpv.push_back(2);
  w_id_list_mpv.push_back(3);
  w_id_list_mpv.push_back(2);
  w_id_list_mpv.push_back(1);
  w_id_list_mpv.push_back(1);
  w_id_list_mpv.push_back(1);
  w_id_list_mpv.push_back(3);
  w_id_list_mpv.push_back(2);
  w_id_list_mpv.push_back(2);
  w_id_list_mpv.push_back(1);
  
  // baseline cut, gain uncertainty from calibration 
  bl_cut_list.push_back(0);  u_gain.push_back(0.041);  
  bl_cut_list.push_back(0);  u_gain.push_back(0.041);  
  bl_cut_list.push_back(-0.001);  u_gain.push_back(0.0008);
  bl_cut_list.push_back(0);  u_gain.push_back(0.041);  
  bl_cut_list.push_back(-0.001);  u_gain.push_back(0.0008);
  bl_cut_list.push_back(-0.001);  u_gain.push_back(0.0008);
  bl_cut_list.push_back(-0.001);  u_gain.push_back(0.0008);
  bl_cut_list.push_back(-0.001);  u_gain.push_back(0.0008);
  bl_cut_list.push_back(-0.001);  u_gain.push_back(0.0008);
  bl_cut_list.push_back(-0.001);  u_gain.push_back(0.0008);
  bl_cut_list.push_back(-0.001);  u_gain.push_back(0.0008);
  bl_cut_list.push_back(-0.001);  u_gain.push_back(0.0008);
  bl_cut_list.push_back(-0.001);  u_gain.push_back(0.0008);   

  int n_setup;
  n_setup = run_name_list.size();
  
  string plot_id;
  TString vav_title;

  // ____ OUTPUT objects ____

  TH1D * h_vec[n_setup];
  TF1 * v_func_vec[n_setup];
  double v_max[n_setup], v_max_err[n_setup], v_mean[n_setup], v_rchi2[n_setup], h_max_vec[n_setup];
  int h_entries[n_setup];


  /********************
  __ FIT ROUTINE______
  ********************/

  vav_fit_nsetup(n_setup, run_name_list, bl_cut_list, run_nr_list, label_list, h_vec, v_func_vec, v_max, v_mean, v_rchi2, h_max_vec, h_entries);

  /***********************
  __ COMBINED PLOTS ______
  ***********************/

  // select data sets & plotting parameters
  // comment in/out desired data set selection

  // // __ all data __
  // vector <int> data_selector = {0,1,2,3,4,5,6,7,8,9,10,11,12};
  // plot_id = "overview_allSetups";
  // vav_title.Form("MIP light yield - overview all setups");
  // double plt_ymax = 0.25; // plot y-limit 
  // vector<double> leg_xy = {0.53, 0.15}; // point coordinates of lower left legend corner
  // double leg_fsize = 0.013; // legend font size

  // __ all data w/o setup 0 __
  vector <int> data_selector = {3,4,5,6,7,8,9,10,11,12};
  plot_id = "overview_allSetups_woShort";
  vav_title.Form("MIP light yield - overview all setups");
  double plt_ymax = 0.12; // plot y-limit 
  vector<double> leg_xy = {0.50, 0.15}; // point coordinates of lower left legend corner
  double leg_fsize = 0.014; // legend font size

  // // __ compare fiber length -> setup 1-4, Mylar __
  // vector <int> data_selector = {6,7,8,12};
  // plot_id = "tile_1_10_mylar";
  // vav_title.Form("MIP light yield  - FCC tile: layer 1 & 10 - compare fiber length");
  // double plt_ymax = 0.11; // plot y-limit 
  // vector<double> leg_xy = {0.35, 0.52}; // point coordinates of lower left legend corner
  // double leg_fsize = 0.02; // legend font size

  // // __ compare wrapping layer 1 tile -> all setup 1 __
  // vector <int> data_selector = {9,10,11,12};
  // plot_id = "setup_1";
  // vav_title.Form("MIP light yield  - FCC tile: layer 1 - setup 1 - compare wrapping");
  // double plt_ymax = 0.095; // plot y-limit 
  // vector<double> leg_xy = {0.35, 0.52}; // point coordinates of lower left legend corner
  // double leg_fsize = 0.02; // legend font size

  // __ compare wrapping layer 10 tile -> all setup 3 __
  // vector <int> data_selector = {4,5,6};
  // plot_id = "setup_3";
  // vav_title.Form("MIP light yield - FCC tile: layer 10 - setup 3 - compare wrapping");
  // double plt_ymax = 0.115; // plot y-limit 
  // vector<double> leg_xy = {0.35, 0.55}; // point coordinates of lower left legend corner
  // double leg_fsize = 0.02; // legend font size

  // // __ compare SiPMs -> setup 3b) __
  // vector <int> data_selector = {3,5};
  // plot_id = "tile_10_SiPM_old_new";
  // vav_title.Form("MIP light yield - FCC tile: layer 10 - setup 3 b) - compare SiPMs");
  // double plt_ymax = 0.12; // plot y-limit 
  // vector<double> leg_xy = {0.35, 0.65}; // point coordinates of lower left legend corner
  // double leg_fsize = 0.02; // legend font size
  

  // initialize objects to copy selected data to
  int n_select = data_selector.size();
  TH1D * h_vec_select[n_select];
  TF1 * v_func_vec_select[n_select];
  double v_max_select[n_select], v_mean_select[n_select], v_rchi2_select[n_select], h_max_vec_select[n_select];
  int h_entries_select[n_select];
  vector<string> label_list_select(n_select);
  vector<int> c_list_select(n_select);
  vector<int> l_list_select(n_select);
  vector<int> m_list_select(n_select);

  for (int i = 0; i < n_select; ++i)
  {
    // map selected 
    label_list_select[i] = label_list[ data_selector[i] ];
    c_list_select[i] = c_list[ data_selector[i] ];
    l_list_select[i] = l_list[ data_selector[i] ];
    m_list_select[i] = m_list[ data_selector[i] ];

    // h_vec_select[i] = h_vec[ data_selector[i] ];
    h_vec_select[i] = (TH1D*)h_vec[ data_selector[i] ]->Clone(Form("h%d",i));
    v_func_vec_select[i]  = v_func_vec[ data_selector[i] ];
    v_max_select[i] = v_max[ data_selector[i] ];
    v_mean_select[i]  = v_mean[ data_selector[i] ];
    v_rchi2_select[i] = v_rchi2[ data_selector[i] ];
    h_max_vec_select[i] = h_max_vec[ data_selector[i] ];
    h_entries_select[i] = h_entries[ data_selector[i] ];
  }

  plot_multi_dist( n_select, vav_title, plot_id, plt_ymax, leg_xy, leg_fsize, label_list_select, c_list_select, l_list_select, m_list_select, h_vec_select, v_func_vec_select, v_max_select, v_mean_select, v_rchi2_select, h_entries_select );

  /*************************
  __ MPV VS. WRAPPING ______
  *************************/

  // __ all data w/o setup 0 __
  vector <int> data_selector_mpv = {3,4,5,6,7,8,9,10,11,12};
  string plot_id_mpv = "overview_allSetups_woShort";
  double plt_ymax_mpv = 12; // plot y-limit 
  vector<double> leg_xy_mpv = {0.25,0.5}; // point coordinates of lower left legend corner
  double leg_fsize_mpv = 0.02; // legend font size

  // initialize objects to copy selected data to
  int n_select_mpv = data_selector_mpv.size();
  double v_max_select_mpv[n_select_mpv];
  vector<double> u_gain_select_mpv(n_select_mpv);
  vector<string> label_list_select_mpv(n_select_mpv);
  vector<int> c_list_select_mpv(n_select_mpv);
  vector<int> l_list_select_mpv(n_select_mpv);
  vector<int> m_list_select_mpv(n_select_mpv);
  vector<int> w_id_list_select_mpv(n_select_mpv);

  // copy selected data
  for (int i = 0; i < n_select_mpv; ++i)
  {
    label_list_select_mpv[i] = label_list[ data_selector_mpv[i] ];
    c_list_select_mpv[i] = c_list_mpv[ data_selector_mpv[i] ];
    l_list_select_mpv[i] = l_list_mpv[ data_selector_mpv[i] ];
    m_list_select_mpv[i] = m_list_mpv[ data_selector_mpv[i] ];
    u_gain_select_mpv[i] = u_gain[ data_selector_mpv[i] ];
    w_id_list_select_mpv[i] = w_id_list_mpv[ data_selector_mpv[i] ];

    v_max_select_mpv[i] = v_max[ data_selector_mpv[i] ];
  }

  TCanvas *C_mpv = new TCanvas("mpv overview","",1500,1000);
  C_mpv->SetGrid();
  TLegend *gr_mpv_leg = new TLegend(leg_xy_mpv[0],leg_xy_mpv[1],0.9,0.9);
  gr_mpv_leg->SetTextSize(leg_fsize_mpv);

  // asymmetric error values from binning dependend fluctuation analyis, hard-coded
  // list of graphs, each one mpv point
  TGraphAsymmErrors *gr_mpv_list[n_select];
  double x,y,xel,yel,xeh,yeh;
  for (int i = 0; i < n_select; ++i)
  {
    x = w_id_list_select_mpv[i];
    y = v_max_select_mpv[i];
    xel=0; xeh=0;
    yel = sqrt( u_gain_select_mpv[i]*u_gain_select_mpv[i] + 0.03*0.03)*v_max_select_mpv[i];
    yeh = sqrt( u_gain_select_mpv[i]*u_gain_select_mpv[i] + 0.05*0.05)*v_max_select_mpv[i];
    gr_mpv_list[i] = new TGraphAsymmErrors(1,&x,&y,&xel,&xeh,&yel,&yeh);
    gr_mpv_list[i]->SetMarkerColor(c_list_select_mpv[i]);
    gr_mpv_list[i]->SetLineColor(c_list_select_mpv[i]);
    gr_mpv_list[i]->SetMarkerStyle(m_list_select_mpv[i]);
    gr_mpv_list[i]->SetMarkerSize(1.3);
    // gr_mpv_list[i]->SetName(Form("%s",label_list[i].c_str()));
    gr_mpv_leg->AddEntry(gr_mpv_list[i],Form("%s, MVP = %1.2f{}_{-%1.2f}^{+%1.2f}",label_list_select_mpv[i].c_str(),v_max_select_mpv[i],yel,yeh),"pe1");
  }
  gr_mpv_list[0]->GetYaxis()->SetRangeUser(3,plt_ymax_mpv);
  gr_mpv_list[0]->GetXaxis()->SetLimits(0.5,3.5);
  gr_mpv_list[0]->GetXaxis()->SetNdivisions(3);
  gr_mpv_list[0]->GetYaxis()->SetNdivisions(9);
  gr_mpv_list[0]->GetXaxis()->SetLabelOffset(-1.5); // get rid of numberd axis labels
  gr_mpv_list[0]->SetTitle("Vavilov fit: most probable value (MPV), compare wrapping materials");
  gr_mpv_list[0]->GetXaxis()->SetTitle("wrapping material");
  gr_mpv_list[0]->GetXaxis()->SetTitleOffset(1.15);
  gr_mpv_list[0]->GetYaxis()->SetTitle("MPV [N_{pe}]");
  gr_mpv_list[0]->Draw("ape1");
  gr_mpv_leg->Draw();

  // Draw labels on the x axis
  TText *t = new TText();
  t->SetTextAlign(40);
  t->SetTextSize(0.035);
  t->SetTextFont(42);
  vector<string> x_labels = {"a) Mylar","b) Tyvek","c) Aluminum"};
  int xpos_label[3] = {1,2,3}; // x,y position in absolute values
  for (Int_t i=0;i<3;i++){ t->DrawText(xpos_label[i]+0.2,2.65,x_labels[i].c_str());}
  // draw mpv points
  for (int i = 0; i < n_select; ++i){gr_mpv_list[i]->Draw("pe1");}
  
  // export
  gErrorIgnoreLevel = kError; // suppress root terminal output 
  C_mpv->Print(Form("./comb_root/vavilov/mvp_vs_wrapping_%s.pdf",plot_id_mpv.c_str()));
  gErrorIgnoreLevel = kUnset ; // allow root terminal output
  return 0;

}


