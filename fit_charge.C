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
#include <TGraphErrors.h>
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

const int n_peaks = 20;

// define a fit function, with single parameter corresponding to ptp distance
Double_t fitf(Double_t *x,Double_t *p)
{
  //0 - N0
  //1 - mu for poison
  //2 - muXT for crosstalk probability
  //3,4 -sogma0, sigma1
  //5,6 - G,B
  double sum = 0;
  for(int k=0;k<=30;k++)
  {
  	Double_t sigma0 = p[3];
  	Double_t sigma1 = p[4];
  	Double_t sigmaK = sqrt(sigma0*sigma0+k*sigma1*sigma1);
  	Double_t mu = p[1];
  	Double_t muXT = p[2];
  	Double_t G= p[5];
  	Double_t B= p[6];
  	
  	sum=sum+p[0]*mu*TMath::Power((mu+k*muXT),k-1)*TMath::Exp(-(mu+k*muXT))/TMath::Factorial(k) * (1./sqrt(2.*TMath::Pi())/sigmaK)*TMath::Exp(-1./sqrt(2)*TMath::Power( ( (x[0]-(k*G+B))/sigmaK ),2 ));
  }
  return sum;
}

// fit function, sum of n_peaks Gaussians 
Double_t alt_f(Double_t *x,Double_t *p)
{
  Double_t gaus_single[n_peaks];
  for (int i = 0; i < n_peaks; ++i)
  {
    gaus_single[i] = p[0+3*i]*TMath::Exp(-TMath::Power( p[1+3*i]-x[0],2) / (2*TMath::Power(p[2+3*i],2)) );
  }

  double gaus_comb = 0;
  for (int i = 0; i < n_peaks; ++i)
    {
      gaus_comb = gaus_comb + gaus_single[i];
    }  
  return gaus_comb;
}

/********************
__ FIT ROUTINE______
********************/

int main(int argc, char *argv[])
{
  // style options
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetStatW(0.1); gStyle->SetStatH(0.1); // stats box size
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);
  gStyle->SetGridColor(16);

  /***** 
  __ INITIALIZE ___________________________
  *****/

  // __SETTINGS_____

  string path = "./comb_root/"+(string)argv[1]+".root";
  string target_path = "./comb_root/pe_spectrum/";
  string out_list = "calib_factor_1325CS_runs_7to50";

  // fit range continuous fits
  float l_range = 0.03; 
  float u_range = 0.03;
 
   // histogram settings
  double xmin = -0.1;
  double xmax = 1.8;
  int nBins = (xmax-xmin)*500;

  // y-range of calib factor plot
  double l_values = 0.071, u_values =0.077;

  bool print_verb = 1;
  bool print_cf_array = 0;
  bool print_cf_err_array = 0;

  // __FIT RANGE__
  // read off spectrum, ranges for individual peaks
  // (minima between peaks)
  std::vector<double> ranges = {-0.04, 0.04, 0.12, 0.19, 0.26, 0.34, 0.42, 0.49, 0.56, 0.63, 0.71, 0.78, 0.85, 0.93, 1.0, 1.07, 1.15, 1.22, 1.29, 1.36, 1.44, 1.51, 1.59, 1.66}; // run 7-46/50
  // std::vector<double> ranges = {-0.04, 0.02, 0.09, 0.16, 0.24, 0.32, 0.39, 0.47, 0.53, 0.61, 0.68, 0.76, 0.83, 0.91, 0.98, 1.05, 1.12, 1.19, 1.27, 1.34, 1.42, 1.49, 1.56, 1.63}; // run 7-41
  // std::vector<double> ranges = {-0.04, 0.02, 0.09, 0.16, 0.24, 0.32, 0.39, 0.47, 0.53, 0.61, 0.68, 0.76, 0.83, 0.91, 0.98, 1.05, 1.12, 1.19, 1.27, 1.34, 1.42, 1.49, 1.56, 1.63}; // run 7-33


  // open tree
  TFile* file = new TFile(path.c_str());
  TTree* tree;
  file->GetObject("ntuple",tree);

  // to show generalized poisson fit
  TCanvas* C1 = new TCanvas("C1","GP fit",1000,800);
  TH1D* h = new TH1D("h","pulse-height spectrum - Hamamatsu S13360-1325CS",nBins,xmin,xmax);
  h->SetLineColorAlpha(kBlack,0.7);
  h->SetMarkerStyle(7);
  h->SetMarkerColorAlpha(kBlack,0.6);

  TString cut("bl_value<-0.001");
  tree->Draw(Form("charge_alt>>h"),cut,"HISTE");

  // to show alternative fit
  TCanvas *C2 = new TCanvas("C2","alt fit",1000,800);
  TH1D* h2 = new TH1D("h2","",nBins,xmin,xmax);
  h2 = (TH1D*)h->Clone();
  h2->GetXaxis()->SetRangeUser(xmin,xmax);

  /***** 
  __ PRE FIT - Single Gauss fits___________________________
  *****/
  
  // get starting value for gauss mean parameter
  float pos_peak[n_peaks];
  for (int i = 0; i < n_peaks; ++i)
  {
    h->GetXaxis()->SetRange( h->GetXaxis()->FindBin(ranges[i]) ,h->GetXaxis()->FindBin(ranges[i+1]));
    pos_peak[i] = h->GetXaxis()->GetBinCenter(h->GetMaximumBin());
  } 
  h->GetXaxis()->SetRange(0,nBins);

  // fit
  TF1 * peak_single[n_peaks];
  Double_t par_single[3*n_peaks];
  for (int i = 0; i < n_peaks; ++i)
  { 
    float range = 0.026; // symmetric fit range for individual peaks
    peak_single[i] = new TF1("peak","gaus",pos_peak[i]-range,pos_peak[i]+range);
    peak_single[i]->SetParameter(2,0.15);
    peak_single[i]->SetParameter(1,pos_peak[i]);
    h->Fit("peak","RQ+");
    peak_single[i]->GetParameters(&par_single[0+3*i]);
  }

  /***** 
  __ GENERALIZED POISSON * GAUSS FIT ___________________________
  *****/

  TF1 *f = new TF1("fitf",fitf,-0.08,1.7,7);
  h->GetYaxis()->SetTitle("entries");
  h->GetXaxis()->SetTitle("charge [V #times ns]");
  f->SetLineColor(3);
  f->SetNpx(1000);

  f->SetParName(0,"N0"); f->SetParameter(0,300 );
  f->SetParName(1,"#mu"); f->SetParameter(1,5.5);
  f->SetParName(2,"#mu_{XT}"); f->SetParameter(2,0.6); f->SetParLimits(2,0,1);
  f->SetParName(3,"#sigma_{0 p.e.}"); f->SetParameter(3,par_single[2]); f->SetParLimits(3,0,10);
  f->SetParName(4,"#sigma_{1 p.e.}"); f->SetParameter(4,par_single[5]); f->SetParLimits(4,0,10);
  f->SetParName(5,"Gain"); f->SetParameter(5,par_single[4]-par_single[1]);
  f->SetParName(6,"Base line"); f->SetParameter(6,par_single[1]);

  h->Fit("fitf","RMQ");
  // h->Fit("fitf","RQ0");
  // h->Fit("fitf","RM0");
  // h->Fit("fitf","RE");

  C1->cd(); 
  h->Draw("sameFUNC");
  


  double par_GP[7];
  f->GetParameters(&par_GP[0]);

  double calib_factor = f->GetParameter(5);
  double calib_factor_err = f->GetParError(5);
  double GP_chi2_ndof = f->GetChisquare()/f->GetNDF();
  double baseline = f->GetParameter(6);
  double baseline_err = f->GetParError(6);
  double norm = f->GetParameter(0);
  double norm_err = f->GetParError(0);
  double mu = f->GetParameter(1);
  double mu_err = f->GetParError(1);
  double mu_xt = f->GetParameter(2);
  double mu_xt_err = f->GetParError(2);
  double sig0 = f->GetParameter(3);
  double sig0_err = f->GetParError(3);
  double sig1 = f->GetParameter(4);
  double sig1_err = f->GetParError(4);

  /***** 
  __ ALTERNATIVE MULTI-GAUSS FIT ___________________________
  *****/
  C2->cd();
    // C2->SetTopMargin(0.90);
  // C2->SetBottomMargin(0.12);
  C2->SetLeftMargin(0.125);
  C2->SetRightMargin(0.05);

  TF1 *alt = new TF1("alt",alt_f, pos_peak[0]-l_range, pos_peak[n_peaks-1]+u_range, 3*n_peaks);
  alt->SetParameters(&par_single[0]);
  alt->SetNpx(1000);
  alt->SetLineStyle(9);
  h2->Fit("alt","RQM");

  h2->GetYaxis()->SetTitle("entries");
  h2->GetXaxis()->SetTitle("charge [V #times ns]");
  h2->Draw("HISTE");

  Double_t par_alt[3+n_peaks];
  Double_t par_alt_err[3+n_peaks];
  alt->GetParameters(&par_alt[0]);

  // draw individual gaussians
  TF1* f_s[n_peaks];
  for (int i = 0; i < n_peaks; ++i)
  {
    f_s[i] = new TF1("s_peak","gaus",0,1);
    f_s[i]->SetRange(par_alt[1+3*i]-0.07,par_alt[1+3*i]+0.07);
    f_s[i]->SetParameters(&par_alt[0+3*i]);
    f_s[i]->SetLineColorAlpha(4,0.4);
    f_s[i]->SetLineStyle(2);
    f_s[i]->SetNpx(1000);
    f_s[i]->Draw("same");
  }

  alt->Draw("same"); // draw fit above individual fit graphs
  h->Draw("FUNCsame"); // also draw GP fit
  C2->cd();

  double gauss_chi2_ndof = alt->GetChisquare()/alt->GetNDF();

  // custom legend
  TLegend *h_leg = new TLegend(0.65,0.45,0.95,0.9);
  h_leg->SetTextSize(0.02);
  // h_leg->SetHeader("integration window 20 ns","C");
  h_leg->AddEntry(h2,Form("#bf{data}"),"lpe");
  h_leg->AddEntry((TObject*)0, Form("entries: %1.0f",h2->GetEntries()), "");
  h_leg->AddEntry((TObject*)0, Form("integration window 100 ns"), "");
  h_leg->AddEntry(f,Form("Generalized Poisson + Gaussian fit"),"l");
  h_leg->AddEntry((TObject*)0, Form("#chi^{2}/ndf = %1.2f",GP_chi2_ndof), "");
  h_leg->AddEntry((TObject*)0, Form("N = %1.f #pm %1.f",norm, norm_err), "");
  h_leg->AddEntry((TObject*)0, Form("#mu = %1.2f #pm %1.2f",mu, mu_err), "");
  h_leg->AddEntry((TObject*)0, Form("#mu_{XT} = %1.3f #pm %1.3f",mu_xt, mu_xt_err), "");
  h_leg->AddEntry((TObject*)0, Form("#sigma_{0 p.e.} = %1.5f #pm %1.5f",sig0, sig0_err), "");
  h_leg->AddEntry((TObject*)0, Form("#sigma_{1 p.e.} = %1.5f #pm %1.5f",sig1, sig1_err), "");
  h_leg->AddEntry((TObject*)0, Form("G_{fit} = %1.5f #pm %1.5f",calib_factor, calib_factor_err), "");
  h_leg->AddEntry((TObject*)0, Form("B = %1.5f #pm %1.5f",baseline, baseline_err), "");
  h_leg->AddEntry(alt,Form("multi-Gaussian fit"),"l");
  h_leg->AddEntry((TObject*)0, Form("#chi^{2}/ndf = %1.2f",gauss_chi2_ndof), "");
  h_leg->AddEntry(f_s[1],Form("multi-Gaussian fit, individual peaks"),"l");

  h_leg->Draw();


  double pos_alt[n_peaks], u_pos_alt[n_peaks];
  for (int i = 0; i < n_peaks; ++i)
  {
    pos_alt[i] = alt->GetParameter(1+3*i);
    u_pos_alt[i] = alt->GetParError(1+3*i);
  }

  // peak-to-peak distance
  double diff[n_peaks-2], u_diff[n_peaks-2], wght[n_peaks-2];
  double w_mean, u_w_mean, sum1, sum2;
  for (int i = 0; i < n_peaks-2; ++i)
  {
    diff[i] = (pos_alt[i+2]-pos_alt[i])/2;
    u_diff[i] = 1./2*sqrt( u_pos_alt[i+2]*u_pos_alt[i+2] + u_pos_alt[i]*u_pos_alt[i] );
    // printf("%f\n",u_diff[i] );
    wght[i] = 1./( u_diff[i] * u_diff[i] );
    sum1 += diff[i]*wght[i];
    sum2 += wght[i];
  }
  // weighted mean, neglecting correlation
  // (currently not used)
  w_mean = sum1/sum2;
  u_w_mean = sqrt (1./sum2);

  /***** 
  __ SYSTEMATICS ___________________________
  *****/

  double var_sum = 0, calib_factor_sys, calib_factor_err_t;

  for (int i = 0; i < n_peaks-2; ++i)
  {
    var_sum = TMath::Power(diff[i]-calib_factor,2);
  } 
  calib_factor_sys = sqrt(var_sum/(n_peaks-2) );
  calib_factor_err_t = sqrt(TMath::Power(calib_factor_sys,2)+TMath::Power(calib_factor_err,2));
  
  double rel_err = calib_factor_err_t/calib_factor;

  /***** 
  __ COMPARE MULTI-GAUSS <-> GEN.POISS*GAUSS ___________________________
  *****/

  double x_values[n_peaks-1];
  double x_values_err[n_peaks-1];
  for (int i = 0; i < n_peaks-1; ++i)
  {
    x_values[i] = i+1;
    x_values_err[i] = 0;
  }

  TCanvas *C3 = new TCanvas("C3","Compare Fit Results",800,500);
  C3->SetGrid();
  gPad->SetGridx(); gPad->SetGridy();
  TGraphErrors * gr_alt = new TGraphErrors(n_peaks-2,x_values,diff,x_values_err,u_diff);
  gr_alt->GetYaxis()->SetRangeUser(l_values,u_values);
  gr_alt->SetName("gr_alt");
  gr_alt->SetTitle("pulse-height spectrum: peak-to-peak distance - Hamamatsu S13360-1325CS");
  gr_alt->SetMarkerColor(4);
  gr_alt->SetMarkerSize(0.5);
  gr_alt->SetMarkerStyle(21);
  gr_alt->GetXaxis()->SetTitle("peak number");
  gr_alt->GetYaxis()->SetTitle("peak-to-peak distance #Delta_{ptp} [V #times ns]");
  gr_alt->Draw("AP");
  

  double x_points[2]; x_points[0] = 1; x_points[1] = n_peaks-2;
  double x_points_err[2]; x_points_err[0] = 0; x_points_err[1] = 0;
  double y_points[2];
  y_points[0]=calib_factor; y_points[1]=calib_factor;
  double y_points_err1[2];
  y_points_err1[0] = calib_factor_err_t; y_points_err1[1] = calib_factor_err_t;
  double y_points_err2[2];
  y_points_err2[0] = calib_factor_err; y_points_err2[1] = calib_factor_err;

   TGraphErrors * gr_comb = new TGraphErrors(2,x_points,y_points,x_points_err,y_points_err1);

  gr_comb->SetName("gr_comb");
  gr_comb->SetLineColor(2);
  gr_comb->SetMarkerColor(2);
  gr_comb->SetMarkerStyle(11);
  gr_comb->SetFillColor(2);
  gr_comb->SetFillColorAlpha(2, 0.2);
  gr_comb->Draw("l3");

  TGraphErrors * gr_comb2 = new TGraphErrors(2,x_points,y_points,x_points_err,y_points_err2);
  gr_comb2->SetName("gr_comb2");
  gr_comb2->SetLineColor(2);
  gr_comb2->SetMarkerColor(2);
  gr_comb2->SetMarkerStyle(11);
  gr_comb2->SetFillColor(2);
  gr_comb2->SetFillColorAlpha(2, 0.45);
  gr_comb2->Draw("l3");
  
  TLegend *leg = new TLegend(0.37,0.6,0.9,0.9);
  leg->SetTextSize(0.03);
  leg->AddEntry("gr_alt","#splitline{individual peak distances: #Delta_{ptp,i}}{multi-Gaussian fit g(x)}","lep");
  leg->AddEntry("gr_comb2",Form("gain G_{fit} = %1.5f",calib_factor),"l");
  leg->AddEntry("gr_comb2",Form("stat: #sigma_{G_{fit}} = %1.5f",calib_factor_err),"lf");
  leg->AddEntry("gr_comb",Form("stat+sys: #sqrt{#sigma^{2}_{G_{fit}} + #sigma^{2}_{#Delta_{ptp}}} = %1.5f #rightarrow rel. uncert.: %1.2f %%",calib_factor_err_t,rel_err*100),"lf");

  leg->Draw();

  /***** 
  __ Export Result ___________________________
  *****/

  if (print_verb)
  {
    // print results
    printf("SiPM : 1pe dist: %1.4f ± %1.4f ± %1.4f | red. chi2 = %1.2f | %1.2f %%\n",calib_factor,calib_factor_err,calib_factor_sys,GP_chi2_ndof,rel_err*100 );
    printf("Baseline (0 p.e. peak) = %1.8f ± %1.8f\n",baseline,baseline_err );
  }
  if (print_cf_array)
  {
    printf("%f, ",calib_factor );
  }
  if (print_cf_err_array)
  {
    printf("%f, ",calib_factor_err );
  }

  string pdf_filename1 = target_path+(string)argv[1]+"_GP_fit.pdf";
  string pdf_filename2 = target_path+(string)argv[1]+"_contG_fit.pdf";
  string pdf_filename3 = target_path+(string)argv[1]+"_values.pdf";
  string list_filename = target_path+out_list+".txt";

  FILE * factor_list;
  factor_list = fopen(list_filename.c_str(),"a");
  fprintf(factor_list, "SiPM : 1pe dist: %f ± %f ± %f ± %f  | red. chi2 = %f | %1.2f %%\n",calib_factor,calib_factor_err,calib_factor_sys,calib_factor_err_t,GP_chi2_ndof,rel_err*100 );
  fclose(factor_list);

  gErrorIgnoreLevel = kError; // suppress root terminal output 
    C1->Print(pdf_filename1.c_str());
    C2->Print(pdf_filename2.c_str());
    C3->Print(pdf_filename3.c_str());
  
  return 0;

}


