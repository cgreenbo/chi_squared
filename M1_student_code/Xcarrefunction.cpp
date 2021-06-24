#include <iostream>
#include "TStyle.h"
#include "../include/Xcarre.hpp"


#define Xcarrefunction_cpp



double StatisticAndSystematic_function(double syst_ttbar, double syst_multijets, double syst_t_channel, double syst_tbar_channel, double syst_ElectroWeak, double lumi,bool Graph_Drawing)
{
    
    Config config;
    read (&config);

    string file_hist_top = "results/histo/hist_top"+config.version+".root";   //genered by sensibility
    string file_hist_antitop = "results/histo/hist_antitop"+config.version+".root";



    TFile* Filet = TFile::Open(file_hist_top.c_str());      //file top
    TFile* Filetbar = TFile::Open(file_hist_antitop.c_str());     //file antitop

    
    int Nbin = 1e5;
    double pas =  2.0*config.range/(double) Nbin;
    
    
  //signal = modulation *b
  //on va chercher les valeurs des bruits de fonds dans fichiers histo genere par sensibility
    TH1F* datt = (TH1F*)Filet->Get("Asimov");
    TH1F* Sigt = (TH1F*)Filet->Get("Modulation_t");
    TH1F* ttbart = (TH1F*)Filet->Get("ttbar");
    TH1F* electroweakt = (TH1F*)Filet->Get("Electroweak");
    TH1F* multijett = (TH1F*)Filet->Get("Multijet");
    TH1F* SingleTopSMt = (TH1F*)Filet->Get("t_channel");

    TH1F* chiCarre = new TH1F ("Chi Carre", "" , Nbin, -config.range, config.range);
    TF1* fithist = new TF1 ("pol2","pol2", -config.range,config.range);
    fithist->SetParameters(0,100.0);
    fithist->SetParameters(1,0.);
    fithist->SetParameters(2,0.);

    TH1F* dattbar = (TH1F*)Filetbar->Get("Asimov");
    TH1F* Sigtbar = (TH1F*)Filetbar->Get("Modulation_tbar");
    TH1F* ttbartbar = (TH1F*)Filetbar->Get("ttbar"); //comes from hist_tchannel
    TH1F* electroweaktbar = (TH1F*)Filetbar->Get("Electroweak");
    TH1F* SingleTopSMtbar = (TH1F*)Filetbar->Get("tbar_channel");
    TH1F* multijettbar = (TH1F*)Filetbar->Get("Multijet");

    double asimovt = datt->GetBinContent(1);        //access the bin content of a given bin, here bin=1
    double ttbartbkg = ttbart->GetBinContent(1)*syst_ttbar;
    double multijettbkg = multijett->GetBinContent(1)*syst_multijets;
    double SingleTopSMtbkg = SingleTopSMt->GetBinContent(1)*syst_t_channel;
    double electroweaktbkg = electroweakt->GetBinContent(1)*syst_ElectroWeak; //multipile par 1

    double asimovtbar = dattbar->GetBinContent(1);
    double ttbartbarbkg = ttbartbar->GetBinContent(1)*syst_ttbar;
    double multijettbarbkg = multijettbar->GetBinContent(1)*syst_multijets;
    double electroweaktbarbkg = electroweaktbar->GetBinContent(1)*syst_ElectroWeak;
    double SingleTopSMtbarbkg = SingleTopSMtbar->GetBinContent(1)*syst_tbar_channel;
    double sigmatott = sqrt(lumi*(ttbartbkg+multijettbkg+electroweaktbkg+SingleTopSMtbkg));
    double sigmatottbar = sqrt(lumi*(ttbartbarbkg+multijettbarbkg+electroweaktbarbkg+SingleTopSMtbarbkg));


    double dividet = 1.0/(pow(sigmatott,2)*24);      //barres d'erreurs top
    double dividetbar = 1.0/(pow(sigmatottbar,2)*24);     //barres d'erreurs antitop
    int bin = 1;
    
 //tracage par Bin
    for (double i = 0 ; i < Nbin; i++)
    {
      double b = -config.range + (double) i * pas;
      double tchannel = 0;
      double tbarchannel = 0;

      for (int p=1; p<=24; p++)   //tracage horaire 
      {//X^2
        tchannel += pow(asimovt-lumi*((b*Sigt->GetBinContent(p)+1)*SingleTopSMtbkg+ttbartbkg+electroweaktbkg+multijettbkg),2);
        tbarchannel += pow(asimovtbar-lumi*((b*Sigtbar->GetBinContent(p)+1)*SingleTopSMtbarbkg+ttbartbarbkg+electroweaktbarbkg+multijettbarbkg),2);
      }

      chiCarre->SetBinContent(bin,tchannel*dividet+tbarchannel*dividetbar);
      bin++;
    }
    //cout<<"Avant"<<endl;
    chiCarre->Fit(fithist,"R"); //minimisation
    //cout<<"Apres"<<endl;
    //exit(1);
    if(Graph_Drawing==true)
    {   
    gStyle->SetOptFit(1); 
    auto legend = new TLegend(0.65,0.75,0.90,0.85);
    chiCarre->GetXaxis()->SetTitle("b_{#mu} (GeV)");
    chiCarre->GetYaxis()->SetTitle("#delta#chi^{2}");
    chiCarre->GetXaxis()->CenterTitle(kTRUE);
    chiCarre->GetYaxis()->CenterTitle(kTRUE);
    chiCarre->SetStats(kFALSE);
    string nameXroot = "XCarre"+config.version+".root";
    TFile* Xroot = new TFile(nameXroot.c_str(),"RECREATE");
    Xroot->cd();
    chiCarre->Write("chiCarre");
    chiCarre->SetLineColor(kRed);
    chiCarre->Draw();
    legend->AddEntry(chiCarre, "#chi_{tot}^{2} = #chi_{t}^{2} +#chi_{#bar{t}}^{2}", "l");
    legend->Draw();
    TLine* deltaChi = new TLine (-config.range, 1, config.range, 1);
    int min=chiCarre->GetMinimumBin();
    double content=chiCarre->GetBinContent(min)+1;
    TLine* deltaChiPlus = new TLine (-abs(fithist->GetX(content,-config.range,config.range)), 0, -abs(fithist->GetX(content,-config.range,config.range)) , 1);
    TLine* deltaChiMoins = new TLine (abs(fithist->GetX(content,-config.range,config.range)), 0, abs(fithist->GetX(content,-config.range,config.range)), 1);
    //deltaChi->SetLineColor(kBlue);
    deltaChi->Draw();
    //legend->AddEntry(deltaChi, "#delta #chi^{2}", "l");
    deltaChiPlus->Draw();
    deltaChiMoins->Draw();
    //Xroot->Close();
    TCanvas* feuille = new TCanvas ("", "", 450, 450);
    feuille->SetLeftMargin(0.14);
    feuille->SetRightMargin(0.03);
    feuille->SaveAs(("results/graph/XCarre"+config.version+".pdf").c_str());
    Xroot->Close();
    }
    int min=chiCarre->GetMinimumBin();
    double content=chiCarre->GetBinContent(min)+1;
    //chiCarre->Fit(fithist,"R"); //minimisation
    //cout<<"Delta Chi at 1 GeV = "<<fithist->GetX(chiCarre->GetBinContent(chiCarre->GetMinimumBin())+1)<<endl;
  //**** pour bu (stat+syst) ^^^^^^^^si on met les syst != 1.0
return config.time_modulation*fithist->GetX(content,-config.range,config.range);       //statisical error
}




/********************************FONCTION CHARME******************************/

double StatisticAndSystematic_Charm(double syst_TTbar, double syst_Multijets, double syst_singleTop, double syst_Diboson, double syst_Zjets, double luminosity, double syst_Wlight, double sysyt_WCharm,double syst_WAntiCharm,bool Graph_Drawing)
{
    
    Config config;
    read (&config);

    string file_hist_charm = "results/histo/hist_charme"+config.version+".root";   //genered by sensibility
    string file_hist_anticharm = "results/histo/hist_anticharme"+config.version+".root";



    TFile* Filec = TFile::Open(file_hist_charm.c_str());      //file top
    TFile* Filecbar = TFile::Open(file_hist_anticharm.c_str());     //file antitop

    
    int Nbin = 1e5;
    double pas =  2.0*config.range/(double) Nbin;
    
    
  //signal = modulation *b
  //on va chercher les valeurs des bruits de fonds dans fichiers histo genere par sensibility
    TH1F* Data_C = (TH1F*)Filec->Get("Asimov");
    TH1F* Signal_c = (TH1F*)Filec->Get("Modulation_c");

    TH1F* ttbar_c = (TH1F*)Filec->Get("ttbar");
    TH1F* Zjets_c = (TH1F*)Filec->Get("Zjets");
    TH1F* multijet_c = (TH1F*)Filec->Get("Multijets");
    TH1F* Diboson_c = (TH1F*)Filec->Get("Diboson");
    TH1F* W_light_c = (TH1F*)Filec->Get("W+light");
    TH1F* SingleTop_c = (TH1F*)Filec->Get("single_top");
    TH1F* W_c = (TH1F*)Filec->Get("W+c");

    TH1F* chiCarre = new TH1F ("Chi Carre", "" , Nbin, -config.range, config.range);
    TF1* fithist = new TF1 ("pol2","pol2", -config.range,config.range);
    fithist->SetParameters(0,100.0);
    fithist->SetParameters(1,0.);
    fithist->SetParameters(2,0.);

    TH1F* Data_Cbar = (TH1F*)Filecbar->Get("Asimov");
    TH1F* Signal_cbar = (TH1F*)Filecbar->Get("Modulation_cbar");

    TH1F* ttbar_cbar = (TH1F*)Filecbar->Get("ttbar");
    TH1F* Zjets_cbar = (TH1F*)Filecbar->Get("Zjets");
    TH1F* multijet_cbar = (TH1F*)Filecbar->Get("Multijet");
    TH1F* Diboson_cbar = (TH1F*)Filecbar->Get("Diboson");
    TH1F* W_light_cbar = (TH1F*)Filecbar->Get("W+light");
    TH1F* SingleTop_cbar = (TH1F*)Filecbar->Get("singletop");
    TH1F* W_cbar = (TH1F*)Filecbar->Get("W+cbar");

   

    double asimovc = Data_C->GetBinContent(1);        //access the bin content of a given bin, here bin=1
    double ttbar_c_bkg = ttbar_c->GetBinContent(1)*syst_TTbar;
    double multijet_c_bkg = multijet_c->GetBinContent(1)*syst_Multijets;
    double SingleTop_c_bkg = SingleTop_c->GetBinContent(1)*syst_singleTop;
    double Zjets_c_bkg = Zjets_c->GetBinContent(1)*syst_Zjets; 
    double Diboson_c_bkg = Diboson_c->GetBinContent(1)*syst_Diboson; 
    double W_light_c_bkg = W_light_c->GetBinContent(1)*syst_Wlight; 
    double W_c_data = W_c->GetBinContent(1)*sysyt_WCharm;


    
    double asimovcbar = Data_Cbar->GetBinContent(1);        //access the bin content of a given bin, here bin=1
    double ttbar_cbar_bkg = ttbar_cbar->GetBinContent(1)*syst_TTbar;
    double multijet_cbar_bkg = multijet_cbar->GetBinContent(1)*syst_Multijets;
    double SingleTop_cbar_bkg = SingleTop_cbar->GetBinContent(1)*syst_singleTop;
    double Zjets_cbar_bkg = Zjets_cbar->GetBinContent(1)*syst_Zjets; 
    double Diboson_cbar_bkg = Diboson_cbar->GetBinContent(1)*syst_Diboson; 
    double W_light_cbar_bkg = W_light_cbar->GetBinContent(1)*syst_Wlight;
    double W_cbar_data = W_cbar->GetBinContent(1)*syst_WAntiCharm;
    

    double bkg_SUM_c=ttbar_c_bkg+multijet_c_bkg+SingleTop_c_bkg+Zjets_c_bkg+Diboson_c_bkg+W_light_c_bkg+W_c_data;
    double bkg_SUM_cbar=ttbar_cbar_bkg+multijet_cbar_bkg+SingleTop_cbar_bkg+Zjets_cbar_bkg+Diboson_cbar_bkg+W_light_cbar_bkg+W_cbar_data;
    double sigmatott = sqrt(luminosity*bkg_SUM_c);
    double sigmatottbar = sqrt(luminosity*bkg_SUM_cbar);

    double dividet = 1.0/(pow(sigmatott,2)*24);      //barres d'erreurs top
    double dividetbar = 1.0/(pow(sigmatottbar,2)*24);     //barres d'erreurs antitop
    int bin = 1;
    
 //tracage par Bin
    for (double i = 0 ; i < Nbin; i++)
    {
      double b = -config.range + (double) i * pas;
      double c_channel = 0;
      double cbar_channel = 0;

      for (int p=1; p<=24; p++)   //tracage horaire 
      {//X^2
        c_channel += pow(asimovc-luminosity*((b*Signal_c->GetBinContent(p)+1)*bkg_SUM_c),2);
        cbar_channel += pow(asimovcbar-luminosity*((b*Signal_cbar->GetBinContent(p)+1)*bkg_SUM_cbar),2);
      }

      chiCarre->SetBinContent(bin,c_channel*dividet+cbar_channel*dividetbar);
      bin++;
    }

    chiCarre->Fit(fithist,"R"); //minimisation

    if(Graph_Drawing==true)
    {    
      gStyle->SetOptFit(1);
      auto legend = new TLegend(0.65,0.75,0.90,0.85);
      chiCarre->GetXaxis()->SetTitle("b_{#mu} (GeV)");
      chiCarre->GetYaxis()->SetTitle("#delta#chi^{2}");
      chiCarre->GetXaxis()->CenterTitle(kTRUE);
      chiCarre->GetYaxis()->CenterTitle(kTRUE);
      chiCarre->SetStats(kFALSE);

    //creation fichier root
      string nameXroot = "XCarreCharm"+config.version+".root";
      TFile* Xroot = new TFile(nameXroot.c_str(),"RECREATE");
      Xroot->cd();
      chiCarre->Write("chiCarre");
      chiCarre->SetLineColor(kRed);
      chiCarre->Draw();
      legend->AddEntry(chiCarre, "#chi_{tot}^{2} = #chi_{c}^{2} +#chi_{#bar{c}}^{2}", "l");
      legend->Draw();
      int min=chiCarre->GetMinimumBin();
      double content=chiCarre->GetBinContent(min)+1;
      TLine* deltaChi = new TLine (-config.range, 1, config.range, 1);
      TLine* deltaChiPlus = new TLine (-abs(fithist->GetX(content)), 0, -abs(fithist->GetX(content)) , 1);
      TLine* deltaChiMoins = new TLine (abs(fithist->GetX(content)), 0, abs(fithist->GetX(content)), 1);
      //deltaChi->SetLineColor(kBlue);
      deltaChi->Draw();
      //legend->AddEntry(deltaChi, "#delta #chi^{2}", "l");
      deltaChiPlus->Draw();
      deltaChiMoins->Draw();
      //Xroot->Close();
      TCanvas* feuille = new TCanvas ("", "", 450, 450);
      feuille->SetLeftMargin(0.14);
      feuille->SetRightMargin(0.03);
      feuille->SaveAs(("results/graph/XCarreCharm"+config.version+".pdf").c_str());
      Xroot->Close();
    }

    int min=chiCarre->GetMinimumBin();
    double content=chiCarre->GetBinContent(min)+1;
    //chiCarre->Fit(fithist,"R"); //minimisation
    //cout<<"Delta Chi at 1 GeV = "<<fithist->GetX(chiCarre->GetBinContent(chiCarre->GetMinimumBin())+1)<<endl;
    //**** pour bu (stat+syst) ^^^^^^^^si on met les syst != 1.0
  return config.time_modulation*fithist->GetX(content,-config.range,config.range);       //statisical error
}
