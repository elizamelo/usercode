//STANDARD ROOT INCLUDES
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TLegend.h>
#include <TDirectory.h>
#include <TStyle.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TMatrixT.h>
#include <TVectorT.h>

//STANDARD C++ INCLUDES
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
using namespace std;

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <TH1F.h>
#include <TH2F.h>
#include <THStack.h>
//#include "Styles.C"
#include <TCanvas.h>
#include <TLegend.h>
#include <TRandom.h>
#include <TMath.h>
#include "TLatex.h"

using std::endl;
using std::cerr;
using std::cout;

#include "TLine.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TLegend.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TH2F.h"

#include "TF1.h"
#include "TGraphErrors.h"
//#include "tdrstyle.C"
Bool_t debug = true;
Bool_t RPRight(true);
Bool_t RPLeft(!RPRight);
void PlotterDimuonsMCZB(){
    StyleNice();
    cout<<"RPRight: "<<RPRight<<endl;
    cout<<"RPLeft: "<<RPLeft<<endl;
    //______________________________________________________
    //Selecting Jpsi , rp right  (with t cut)
    if (RPRight){         
        // MakePlot("proton_right_xi_kint","xi_proton_signal_right_kint","xi_proton_backg_right_kint","xi_proton_both_right_kint", 1,-0.049,1.0,0.0001,1000.0,0.75, 0.5, 0.9, 0.75,0.1,1E0,0.1,.5E0,0.1,1E0,0,1,"#xi_{Right}","N Events",0);
        MakePlot("proton_right_xi_sel","xi_proton_signal_right_kin_cut","xi_proton_backg_right_kin_cut","xi_proton_both_right_kin_cut", 1,-0.049,1.0,0.0001,1000.0,0.75, 0.5, 0.9, 0.75,0.1,1E0,0.1,.5E0,0.1,1E0,0,1,"#xi_{Right}","N Events",0);
        MakePlot("proton_right_t_sel","t_proton_signal_right_kin_cut","t_proton_backg_right_kin_cut", "t_proton_backg_both_kin_cut", 1,-0.049,1.0,0.0001,1000.0,0.75, 0.5, 0.9, 0.75,0.1,1E0,0.1,.5E0,0.1,1E0,0,1,"|t|_{Right} [GeV^{2}]","dN/d|t|", 1);
        //MakePlot("proton_right_xi_halosel","xi_proton_signal_right_kin_cut_halo","xi_proton_backg_right_kin_cut_halo","xi_proton_both_right_kin_cut_halo", 1,-0.049,1.0,0.0001,1000.0,0.75, 0.5, 0.9, 0.75,0.1,1E0,0.1,.5E0,0.1,1E0,0,1,"#xi_{Right}","N Events", 0);
        //MakePlot("proton_right_t_halosel","t_proton_signal_right_kin_cut_halo","t_proton_backg_right_kin_cut_halo", "t_proton_backg_both_kin_cut_halo", 1,-0.049,1.0,0.0001,1000.0,0.75, 0.5, 0.9, 0.75,0.1,1E0,0.1,.5E0,0.1,1E0,0,1,"|t|_{Right} [GeV^{2}]","dN/d|t|", 1);
        //MakePlot("xitotem_xicms_rightRPs_kin","xi_cms_totem_right_signal_kin","xi_cms_totem_right_zb_kin","xi_cms_totem_right_both_kin",  1,-0.29,0.43,0.0001,10000.0,0.85, 0.7, 0.9, 0.95,0.1,1E0,0.1,.5E0,0.1,1E0,0,0,"#xi^{CMS-}-#xi^{TOTEM-Right}","N Events",0);

    }
    //__________________________________________________ 
    //Selecting Jpsi , rp left  (with t cut)
    if(RPLeft){
        //MakePlot("proton_Left_xi_kint","xi_proton_signal_left_kint","xi_proton_backg_left_kint","xi_proton_both_left_kint", 1,-0.049,1.0,0.0001,1000.0,0.75, 0.5, 0.9, 0.75,0.1,1E0,0.1,.5E0,0.1,1E0,1,0,"#xi_{left}","N Events",0);
        //MakePlot("proton_Left_xi_sel","xi_proton_signal_left_kin_cut","xi_proton_backg_left_kin_cut","xi_proton_both_left_kin_cut",  1,-0.049,1.0,0.0001,1000.0,0.75, 0.5, 0.9, 0.75,0.1,1E0,0.1,.5E0,0.1,1E0,0,1,"#xi_{Left}","N Events",0);
        //MakePlot("proton_Left_t_sel","t_proton_signal_left_kin_cut","t_proton_backg_left_kin_cut", "t_proton_left_backg_both_kin_cut", 1,- 0.049,1.0,0.0001,1000.0,0.75, 0.5, 0.9, 0.75,0.1,1E0,0.1,.5E0,0.1,1E0,0,1,"|t|_{Left} [GeV^{2}]","dN/d|t|", 1);
        //MakePlot("proton_Left_xi_halosel","xi_proton_signal_left_kin_cut_halo","xi_proton_backg_left_kin_cut_halo", "xi_proton_both_left_kin_cut_halo", 1,-0.049,1.0,0.0001,1000.0,0.75, 0.5, 0.9, 0.75,0.1,1E0,0.1,.5E0,0.1,1E0,0,1,"#xi_{Left}","N Events",  0);
        //MakePlot("proton_Left_t_halosel","t_proton_signal_left_kin_cut_halo","t_proton_backg_left_kin_cut_halo",  "t_proton_left_backg_both_kin_cut_halo", 1,-0.049,1.0,0.0001,1000.0,0.75, 0.5, 0.9, 0.75,0.1,1E0,0.1,.5E0,0.1,1E0,0,1,"|t|_{Left} [GeV^{2}]","dN/d|t|", 1);
        //MakePlot("xitotem_xicms_leftRPs_kin","xi_cms_totem_left_signal_kin","xi_cms_totem_left_zb_kin","xi_cms_totem_left_both_kin",  1,-0.29,0.43,0.0001,10000.0,0.85, 0.7, 0.9, 0.95,0.1,1E0,0.1,.5E0,0.1,1E0,0,0,"#xi^{CMS-}-#xi^{TOTEM-Left}","N Events",0);
    }       

    //-------------------------------

}

void MakePlot(TString name,TString name2, TString name3,TString name4,bool logscale,
        double xmin,double xmax,double ymin,double ymax,
        double xminleg, double yminleg, double xmaxleg, double ymaxleg,
        double xtex1, double ytex1, double xtex2, double ytex2, double
        xtex3, double ytex3, bool t_xi_cut_proton_left, bool t_xi_cut_proton_right,
        TString xname,TString yname, bool width){
    //double bsize = 0.0001;
    TCanvas *can = new TCanvas(name+TString("_AllMCvsData_12May2014"), name+TString("_AllMCvsData_12May2014"), 0, 0, 600, 700);
    //stack = new THStack("MC", "MC");
    if(logscale) can->SetLogy(1);
    TLegend* leg = new TLegend(xminleg,yminleg, xmaxleg, ymaxleg);
    TLegend* leg2 = new TLegend(xminleg,yminleg, xmaxleg, ymaxleg);
    leg->SetBorderSize(0);
    leg->SetLineStyle(0);
    leg->SetTextSize(.05);
    leg->SetFillStyle(0);
    leg->SetFillColor(0);

    leg2->SetBorderSize(0);
    leg->SetLineStyle(0);
    leg2->SetTextSize(.05);
    leg2->SetFillStyle(0);
    leg2->SetFillColor(0);
    //============================================
    can -> Draw();

    //============================================ 
    //TFile *data = TFile::Open("RootFiles/data_histos_12May.root");
    TFile *data = TFile::Open("RootFiles/histos_DATASDJpsi_25March2014.root"); 
    //TFile *data = TFile::Open("RootFiles/histos_DATASDJpsi_12March2014.root");
    //============================================
    ///MC SD MINUS
    if(RPRight){
        //TFile *pythia = TFile::Open("RootFiles/pythia_11May2014fzb_pu.root");
        //TFile *pompyt= TFile::Open("RootFiles/pompyt_SDMinus_11May2014fzb_pu.root");
        TFile *pythia = TFile::Open("RootFiles/pythia_24MarchMod_pu.root");
        TFile *pompyt= TFile::Open("RootFiles/pompyt_SDMinus_24MarchMod_pu.root");
    }

    /**/
    //============================================
    ///MC SD PLUS
    if(RPLeft){
        TFile *pythia = TFile::Open("RootFiles/pythia_24MarchMod_pu.root");
        TFile *pompyt= TFile::Open("RootFiles/pompyt_SDPlus_24MarchMod_pu.root");
    }
    //============================================
    TH1F* h_data = (TH1F*)data->Get(name); 
    TH1F* h_pompyt_sig = (TH1F*)pompyt->Get(name2); 
    TH1F* h_pythia_sig = (TH1F*)pythia->Get(name2);
    TH1F* h_pompyt_bckg = (TH1F*)pompyt->Get(name3); 
    TH1F* h_pythia_bckg = (TH1F*)pythia->Get(name3);
    TH1F* h_pompyt_both = (TH1F*)pompyt->Get(name4); 
    TH1F* h_pythia_both = (TH1F*)pythia->Get(name4);

    //TH1F* h_pomwig = (TH1F*)pomwig->Get(name2);

    //Label
    TString labelData = "Data ";  
    TString labelMC = "All MC + ZB";
    TString labelBackg = "Background";

    TString labelMCPOMPYTSig = "p from MC SD";  
    TString labelMCPOMPYTBckg = "p from ZB";  
    TString labelMCPOMPYTBoth = "p from MC and from ZB";  

    TString labelMCPYTHIASig = "p from MC ND";  
    TString labelMCPYTHIABckg = "p from ZB in MC ND";  
    TString labelMCPYTHIABoth = "p from MC ND and from ZB";  

    //(#xi_{CMS}-#xi_{TOTEM}<0.009)
    TH1F *hist_data = (TH1F*)h_data->Clone("hist_data");
    TH1F *hist_pompytSig = (TH1F*)h_pompyt_sig->Clone("hist_pompytSig");
    TH1F *hist_pompytBckg = (TH1F*)h_pompyt_bckg->Clone("hist_pompytBckg");
    TH1F *hist_pompytBoth = (TH1F*)h_pompyt_both->Clone("hist_pompytBoth");
    TH1F *hist_pythiaSig = (TH1F*)h_pythia_sig ->Clone("hist_pythiaSig");
    TH1F *hist_pythiaBckg = (TH1F*)h_pythia_bckg->Clone("hist_pythiaBckg");
    TH1F *hist_pythiaBoth = (TH1F*)h_pythia_both->Clone("hist_pythiaBoth");


    //THStack *stack;

    hist_data -> SetDefaultSumw2();
    hist_data -> SetStats(0); 

    hist_data->GetYaxis()->SetTitle(yname);
    hist_data->GetXaxis()->SetTitle(xname);
    hist_data  -> GetYaxis()->SetRangeUser(ymin,ymax);
    hist_data -> GetXaxis()->SetRangeUser(xmin,xmax);

    hist_pompytSig->GetYaxis()->SetTitle(yname);
    hist_pompytSig->GetXaxis()->SetTitle(xname);
    hist_pompytSig->GetYaxis()->SetRangeUser(ymin,ymax);
    hist_pompytSig->GetXaxis()->SetRangeUser(xmin,xmax);

    hist_pompytBckg->GetYaxis()->SetTitle(yname);
    hist_pompytBckg->GetXaxis()->SetTitle(xname);
    hist_pompytBckg -> GetYaxis()->SetRangeUser(ymin,ymax);
    hist_pompytBckg-> GetXaxis()->SetRangeUser(xmin,xmax);

    hist_pompytBoth->GetYaxis()->SetTitle(yname);
    hist_pompytBoth->GetXaxis()->SetTitle(xname);
    hist_pompytBoth->GetYaxis()->SetRangeUser(ymin,ymax);
    hist_pompytBoth->GetXaxis()->SetRangeUser(xmin,xmax);

    hist_pythiaSig->GetYaxis()->SetTitle(yname);
    hist_pythiaSig->GetXaxis()->SetTitle(xname);
    hist_pythiaSig->GetYaxis()->SetRangeUser(ymin,ymax);
    hist_pythiaSig->GetXaxis()->SetRangeUser(xmin,xmax);

    hist_pythiaBckg->GetYaxis()->SetTitle(yname);
    hist_pythiaBckg->GetXaxis()->SetTitle(xname);
    hist_pythiaBckg -> GetYaxis()->SetRangeUser(ymin,ymax);
    hist_pythiaBckg-> GetXaxis()->SetRangeUser(xmin,xmax);

    hist_pythiaBoth->GetYaxis()->SetTitle(yname);
    hist_pythiaBoth->GetXaxis()->SetTitle(xname);
    hist_pythiaBoth->GetYaxis()->SetRangeUser(ymin,ymax);
    hist_pythiaBoth->GetXaxis()->SetRangeUser(xmin,xmax);

    hist_data->SetMarkerStyle(20);
    hist_data->SetMarkerColor(kBlack);
    hist_data->SetLineColor(kBlack);
    hist_data->SetMarkerSize(1.3);


    leg->AddEntry(hist_data,labelData,"EP");
    leg->AddEntry(hist_pompytSig,labelMCPOMPYTSig,"LF");
    leg->AddEntry(hist_pompytBckg,labelMCPOMPYTBckg,"LF");
    leg->AddEntry(hist_pompytBoth,labelMCPOMPYTBoth,"LF");
    leg->AddEntry(hist_pythiaSig,labelMCPYTHIASig,"LF");
    leg->AddEntry(hist_pythiaBckg,labelMCPYTHIABckg,"LF");
    leg->AddEntry(hist_pythiaBoth,labelMCPYTHIABoth,"LF");



    hist_pompytSig->SetLineColor(214);
    hist_pompytSig->SetLineWidth(4);

    hist_pompytBckg->SetLineColor(217);
    hist_pompytBckg->SetLineWidth(4);
    hist_pompytBckg->SetLineStyle(2);

    hist_pompytBoth->SetLineColor(209);
    hist_pompytBoth->SetLineWidth(4);
    //hist_pompytBoth->SetLineStyle(3);

    hist_pythiaSig->SetLineColor(205);
    hist_pythiaSig->SetLineWidth(4);

    hist_pythiaBckg->SetLineColor(95);
    hist_pythiaBckg->SetLineWidth(4);
    hist_pythiaBckg->SetLineStyle(2);

    hist_pythiaBoth->SetLineColor(223);
    hist_pythiaBoth->SetLineWidth(4);
    hist_pythiaBoth->SetLineStyle(3);
    //___________________________________________________
    //***********Norm to Lumi*********************
    //****** MC SD MINUS
    if(RPRight){
        //float cross_section_pw = 12.4; //nb cross section for pomwig minus
        float cross_section_po = 1.9635e+03*0.06;//258000.0*0.06; //nb cross section for pompyt ;
        float cross_section_py = 232000.0*0.06*0.0463;//1.378e+04*0.06; //nb cross section for pythia 
        float luminosity = 49.15;//nb-1
        //float n_data_pw = luminosity*cross_section_pw;
        float n_data_po = luminosity*cross_section_po;
        float n_data_py = luminosity*cross_section_py;
        //float n_MC_pw = 796446.0;
        float n_MC_po = 920788.0;
        float n_MC_py =  1.3298e+06;/*Minus*/  
        //******* MC SD PLUS
}else if(RPLeft){
    //float cross_section_pw = 12.7; //nb cross section for pomwig minus
    float cross_section_po = 1.9635e+03*0.06;//257000.0*0.06; //nb cross section for pompyt 
    float cross_section_py = 1.378e+04*0.06;//232000.0227000.0*0.06; //nb cross section for pythia 
    float luminosity = 49.15;//nb-1
    //float n_data_pw = luminosity*cross_section_pw;
    float n_data_po = luminosity*cross_section_po;
    float n_data_py = luminosity*cross_section_py;
    //float n_MC_pw = 786263.0;
    float n_MC_po = 927662.0;
    float n_MC_py = 1.3298e+06;} /*Plus */

    //double scalePomwig = n_data_pw/n_MC_pw;
    double scalePompyt = n_data_po/n_MC_po; 
    double scalePythia = n_data_py/n_MC_py;//L*xsec/Nr
    if(debug){
        cout << "scalePompyt: "<< scalePompyt <<endl;
        cout << "scalePythia: "<< scalePythia <<endl;
    }
/*hist_pythiaSig ->Scale(scalePythia);
hist_pythiaBckg ->Scale(scalePythia);
hist_pythiaBoth ->Scale(scalePythia);
hist_pompytSig ->Scale(scalePompyt);
hist_pompytBckg ->Scale(scalePompyt);
hist_pompytBoth ->Scale(scalePompyt);/**/
    
    //POMPYT+PYTHIA List
    TList *listMCbackg = new TList;
    listMCbackg->Add(hist_pompytBckg);                                                                                                  
    listMCbackg->Add(hist_pythiaBckg);
    TH1F *hist_MCbackg = (TH1F*)hist_pompytBckg->Clone("hist_MCbackg");
    hist_MCbackg->Reset();
    hist_MCbackg->Merge(listMCbackg);
    cout<<"Integral bckg before scale ="<<hist_MCbackg->Integral()<<endl;
    /*Double_t xmin_b = hist_MCbackg->GetXaxis()->GetXmin();
    Double_t xmax_b = hist_MCbackg->GetXaxis()->GetXmax();
    TAxis *axis_simu_b =  hist_MCbackg->GetXaxis();
    int bmin_simu_b = axis_simu_b->FindBin(xmin_b);
    int bmax_simu_b = axis_simu_b->FindBin(xmax_b);                                                                                
    int b_exato_b = axis_simu_b->FindBin(0.009);
    cout<<"xmin: "<<xmin_b<<" "<<"xmax: "<<xmax_b<<endl;
    cout<<"bmin: "<<bmin_simu_b<<" "<<"bmax: "<<bmax_simu_b<<endl;
    cout<<"binexato="<<b_exato_b<<endl;
    double area_bckg_B = hist_MCbackg->Integral(b_exato_b+1,bmax_simu_b);
    double area_signal_B = hist_MCbackg->Integral(bmin_simu_b,b_exato_b-1);
*/
if(RPRight){
    cout<<"Integral proton from ZB POMPYT ="<< hist_pompytBckg->Integral()<<endl;                                                       
    cout<<"Integral proton from ZB PYTHIA ="<< hist_pythiaBckg->Integral()<<endl; 
    
    //cout<<"bckg in signal area= "<<area_signal_B <<endl;
    //cout<<"bckg in bckg area= "<<area_bckg_B <<endl;
hist_pythiaSig ->Scale(0.005267);
hist_pythiaBckg ->Scale(0.005267);
hist_pythiaBoth ->Scale(0.005267);
/*hist_pompytSig ->Scale(scalePompyt*0.143);
hist_pompytBckg ->Scale(scalePompyt*0.143);
hist_pompytBoth ->Scale(scalePompyt*0.143);  */
hist_pompytSig ->Scale(scalePompyt);
hist_pompytBckg ->Scale(scalePompyt);
hist_pompytBoth ->Scale(scalePompyt);   
    
/*hist_pythiaSig ->Scale(0.005267);
hist_pythiaBckg ->Scale(0.005267);
hist_pythiaBoth ->Scale(0.005267);
hist_pompytSig ->Scale(9.89168603788515780e-05);
hist_pompytBckg ->Scale(9.89168603788515780e-05);
hist_pompytBoth ->Scale(9.89168603788515780e-05);*/

    cout<<"Af.Scale: Integral proton from ZB POMPYT ="<< hist_pompytBckg->Integral()<<endl;
    cout<<"Af.Scale: Integral proton from ZB PYTHIA ="<< hist_pythiaBckg->Integral()<<endl;

}else if(RPLeft){
    cout<<"Integral proton from ZB POMPYT ="<< hist_pompytBckg->Integral()<<endl;
    cout<<"Integral proton from ZB PYTHIA ="<< hist_pythiaBckg->Integral()<<endl;
    //cout<<"bckg in signal area= "<<area_signal_B <<endl;
    //cout<<"bckg in bckg area= "<<area_bckg_B <<endl;
    hist_pythiaSig ->Scale(6.51041666666666696e-03);
    hist_pythiaBckg ->Scale(6.51041666666666696e-03);
    hist_pythiaBoth ->Scale(6.51041666666666696e-03);
    hist_pompytSig ->Scale(6.63118980123008554e-05);
    hist_pompytBckg ->Scale(6.63118980123008554e-05);
    hist_pompytBoth ->Scale(6.63118980123008554e-05);
    cout<<"Af.Scale: Integral proton from ZB POMPYT ="<< hist_pompytBckg->Integral()<<endl;
    cout<<"Af.Scale: Integral proton from ZB PYTHIA ="<< hist_pythiaBckg->Integral()<<endl;
}
//===========================================
TList *listMCsignal = new TList;
listMCsignal->Add(hist_pompytSig);
listMCsignal->Add(hist_pythiaSig);
TH1F *hist_MCsignal = (TH1F*)hist_pompytSig->Clone("hist_MCsignal");
hist_MCsignal->Reset();
hist_MCsignal->Merge(listMCsignal);
hist_MCsignal -> SetStats(0);  
//*********Sum of MC+ZB*********
TList *listMCZB = new TList;
listMCZB->Add(hist_pompytBoth);
listMCZB->Add(hist_pythiaBoth);
TH1F *hist_MCZB = (TH1F*)hist_pompytBoth->Clone("hist_MCZB");
hist_MCZB->Reset();
hist_MCZB->Merge(listMCZB);
hist_MCZB->SetStats(0);
//*********SUM OF BCKG**********
TList *listBackg = new TList;
listBackg->Add(hist_pompytBckg);
listBackg->Add(hist_pythiaBckg);
TH1F *hist_Backg = (TH1F*)hist_pompytBckg->Clone("hist_Backg");
hist_Backg->Reset();
hist_Backg->Merge(listBackg);
hist_Backg->SetStats(0);
//*********SUM OF 3 Contributions*******
TList *listAll3 = new TList;
listAll3->Add(hist_MCsignal);
listAll3->Add(hist_MCZB);
listAll3->Add(hist_Backg);
TH1F *hist_All3 = (TH1F*)hist_MCsignal->Clone("hist_All3");
hist_All3->Reset();
hist_All3->Merge(listAll3);
hist_All3->SetStats(0);
//============================================
leg2->AddEntry(hist_MCsignal,"proton from MC(Pompyt+Pythia)","LF");
leg2->AddEntry(hist_Backg,"proton from ZB","LF");
leg2->AddEntry(hist_data,labelData,"EP");
leg2->AddEntry(hist_MCZB,"protons from both MC and ZB","LF");
leg2->AddEntry(hist_All3,"MC+(MC&ZB)+ZB","LF");
//============================================
hist_MCZB->SetLineColor(209);
hist_MCZB->SetLineWidth(4);

hist_MCsignal->SetLineColor(214);
hist_MCsignal->SetLineWidth(4);

hist_Backg->SetFillColor(220);
//hist_Backg->SetLineWidth(4);
hist_Backg->SetFillStyle(3002);

hist_All3->SetLineColor(kRed);
hist_All3->SetLineWidth(4);
hist_All3->SetLineStyle(4);
//============================================
if (width){
    hist_data->Scale(1,"width");
    hist_Backg->Scale(1,"width");
    hist_MCZB->Scale(1,"width");
    hist_MCsignal->Scale(1,"width");
    hist_All3->Scale(1,"width");

    hist_pompytSig ->Scale(1,"width");
    hist_pompytBckg ->Scale(1,"width");
    hist_pompytBoth ->Scale(1,"width");
    hist_pythiaSig ->Scale(1,"width");
    hist_pythiaBckg ->Scale(1,"width");
    hist_pythiaBoth ->Scale(1,"width");

}
//=============================================
cout<<"Nr of BCKG from ZB after scale: "<<hist_Backg->Integral()<<endl;
/*Double_t xmin = hist_MCZB->GetXaxis()->GetXmin();
Double_t xmax = hist_MCZB->GetXaxis()->GetXmax();
TAxis *axis_simu =  hist_MCZB->GetXaxis();
int bmin_simu = axis_simu->FindBin(xmin);
int bmax_simu = axis_simu->FindBin(xmax);                                                                                
int b_exato = axis_simu->FindBin(0.009);
cout<<"xmin: "<<xmin<<" "<<"xmax: "<<xmax<<endl;
cout<<"bmin: "<<bmin_simu<<" "<<"bmax: "<<bmax_simu<<endl;
cout<<"binexato="<<b_exato<<endl;
double area_bckg = hist_MCZB->Integral(b_exato+1,bmax_simu);
double area_signal = hist_MCZB->Integral(bmin_simu,b_exato-1);
cout<<"bckg in signal area= "<<area_signal <<endl;
cout<<"bckg in bckg area= "<<area_bckg <<endl;
*/
//=============================================
  /*hist_pompytSig ->Draw("histsame");
  hist_pompytBckg ->Draw("histsame");
  hist_pompytBoth ->Draw("histsame");
  hist_pythiaSig ->Draw("histsame");
  hist_pythiaBckg ->Draw("histsame");
  hist_pythiaBoth ->Draw("histsame");/* */

hist_Backg->Draw("hist");
hist_MCsignal->Draw("histsame");
hist_MCZB->Draw("histsame");
hist_All3->Draw("histsame");/* */
hist_data ->Draw("EPsame");
//hist_MC->Draw();
//============================================


//leg->Draw();
//leg->Draw("same");
leg2->Draw("same");

//============================================
gPad->RedrawAxis();  

if(!(t_xi_cut_proton_left ||t_xi_cut_proton_right) ){
    TLatex* CMStext1 = new TLatex(xtex1,ytex1,"CMSTOTEM - Work in progress"); CMStext1->SetTextFont(42); //CMStext1 -> Draw("sames");
    TLatex* CMStext2 = new TLatex(xtex2,ytex2,"HLT_L1_DoubleMu0 Data - p+p #sqrt{s} = 8 Tev"); CMStext2->SetTextFont(42); //CMStext2 -> Draw("sames");
    TLatex* CMStext3 = new TLatex(xtex3,ytex3,"Runs 198902-198903"); CMStext3->SetTextFont(42);}  //CMStext3 -> Draw("sames");}

    if(t_xi_cut_proton_left){
        TLatex* CMStext1 = new TLatex(xtex1,ytex1,"HLT_L1_DoubleMu0 Data - p+p #sqrt{s} = 8 Tev"); CMStext1->SetTextFont(42);// CMStext1 -> Draw("sames");
        TLatex* CMStext2 = new TLatex(xtex2,ytex2,"RPLeft - 0.0 < |#xi| < 0.15"); CMStext2->SetTextFont(42); //CMStext2 -> Draw("sames");
        TLatex* CMStext3 = new TLatex(xtex3,ytex3,"RPLeft: 0.03<-t< 1.0"); CMStext3->SetTextFont(42); CMStext3->SetTextColor(52);}//CMStext3 -> Draw("sames");}
        if(t_xi_cut_proton_right){
            TLatex* CMStext1 = new TLatex(xtex1,ytex1,"HLT_L1_DoubleMu0 Data - p+p #sqrt{s} = 8 Tev"); CMStext1->SetTextFont(42);// CMStext1 -> Draw("sames");
            TLatex* CMStext2 = new TLatex(xtex2,ytex2,"RPRight - 0.0 < |#xi| < 0.15"); CMStext2->SetTextFont(42);// CMStext2 -> Draw("sames");
            TLatex* CMStext3 = new TLatex(xtex3,ytex3,"RPRight: 0.03<-t< 1.0"); CMStext3->SetTextFont(42);CMStext3->SetTextColor(52);} //CMStext3 -> Draw("sames");}
            //============================================    


            }

// C O S M E T I C S

void StyleNice(){

    //gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    gStyle->SetPadLeftMargin(0.2);
    gStyle->SetPadBottomMargin(0.2);
    gStyle->SetLegendBorderSize(0);

    gStyle->SetMarkerStyle(20);
    gStyle->SetMarkerSize(1.);

    gStyle->SetTitleSize(0.07,"X");
    gStyle->SetTitleSize(0.07,"Y");

    gStyle->SetLabelSize(0.05,"X");
    gStyle->SetLabelSize(0.05,"Y");

    gStyle->SetLabelOffset(0.01,"Y");
    gStyle->SetLabelOffset(0.01,"X");

    gStyle->SetTitleOffset(1.06,"Y");
    gStyle->SetTitleOffset(1.06,"X");

    gStyle->SetLineWidth(3);
    gStyle->SetHistLineWidth(2);

    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetTickLength(0.01,"Y");
    gStyle->SetTickLength(0.02,"X");

    gROOT->ForceStyle();

    /*cout<<hist_data->GetNbinsX()<<endl;
      for (int i=1; i<= hist_data->GetNbinsX(); i++){
      double scaleFactor;
      if(hist_data->GetBinContent(i) > 0.){
      scaleFactor = hist_data->GetBinContent(i)/hist_pompytBoth->GetBinContent(i); cout<<"scaleFactor = " << scaleFactor<<endl;
      hist_MC->SetBinContent(i,hist_MC->GetBinContent(i)*scaleFactor);
      hist_MC->SetBinError(i,hist_MC->GetBinError(i)*scaleFactor);
      hist_Backg->SetBinContent(i,hist_Backg->GetBinContent(i)*scaleFactor);
      hist_Backg->SetBinError(i,hist_Backg->GetBinError(i)*scaleFactor);
      }else{
      hist_MC->SetBinContent(i,hist_MC->GetBinContent(i)*factornon);
      hist_MC->SetBinError(i,hist_MC->GetBinError(i)*factornon);
      hist_Backg->SetBinContent(i,hist_Backg->GetBinContent(i)*factornon);
      hist_Backg->SetBinError(i,hist_Backg->GetBinError(i)*factornon);
      }
      }*/

}
