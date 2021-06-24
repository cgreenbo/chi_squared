#include <iostream>
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#include <sstream>
#include <TH1.h>
#include <TGraph.h>
#include <TF1.h>
#include <TFormula.h>

using namespace std;



double Stat_and_Syst(double syst_tt, double syst_tw, double syst_wzjets, double syst_QCD_EC, double syst_QCD_Barrel)
{
    string EFT = "cbwi"; //Switch betwwen ctwi and cbwi
    bool Draw_Graph = false;
    string KV = "phi"; //Switch between phi and cos

    int nbfiles = 6;
    int nb_bins = 20;
    string suffix[nbfiles];
    
    suffix[0] = "hist_reco_phiStar.root";
    suffix[1] = "hist_reco_cosThetaStar.root";
    suffix[2] = "signal_proc_cosThetaStar_cbwi_elmu_20bins_Rwgt_cbwi2p5.root";
    suffix[3] = "signal_proc_cosThetaStar_ctwi_elmu_20bins_Rwgt_cbwi2p5.root";
    suffix[4] = "signal_proc_PhiStar_cbwi_elmu_20bins_Rwgt_cbwi2p5.root";
    suffix[5] = "signal_proc_PhiStar_ctwi_elmu_20bins_Rwgt_cbwi2p5.root";

    TFile* fInput[nbfiles]; //Import .root files
    string inputDirectory; //pwd to .root files directory
    TH1D* tInput[16]; //Import signals and background noise
    //Import TF1 fits and formulas for Wilson coefficients
    TF1* Fit_phiS_cbwi_Input[nb_bins];
    TF1* Fit_phiS_ctwi_Input[nb_bins];
    TF1* Fit_cosS_cbwi_Input[nb_bins];
    TF1* Fit_cosS_ctwi_Input[nb_bins];
    TFormula* Wilson_phiS_cbwi_Input[nb_bins];
    TFormula* Wilson_phiS_ctwi_Input[nb_bins];
    TFormula* Wilson_cosS_cbwi_Input[nb_bins];
    TFormula* Wilson_cosS_ctwi_Input[nb_bins];

    for(int i=0 ; i<nbfiles ; i++)
    {
        inputDirectory = "data/" + suffix[i];
        fInput[i] = new TFile(inputDirectory.c_str(),"READ");
    }

    tInput[0] = (TH1D*) fInput[0]->Get("elmu__top__nominal");
    tInput[1] = (TH1D*) fInput[0]->Get("elmu__antitop__nominal");
    tInput[2] = (TH1D*) fInput[0]->Get("elmu__tt__nominal");
    tInput[3] = (TH1D*) fInput[0]->Get("elmu__tw__nominal");
    tInput[4] = (TH1D*) fInput[0]->Get("elmu__wzjets__nominal");
    tInput[5] = (TH1D*) fInput[0]->Get("elmu__data_obs__nominal");
    tInput[6] = (TH1D*) fInput[0]->Get("elmu__QCD_DD_EC__nominal");
    tInput[7] = (TH1D*) fInput[0]->Get("elmu__QCD_DD_Barrel__nominal");

    tInput[8] = (TH1D*) fInput[1]->Get("elmu__top__nominal");
    tInput[9] = (TH1D*) fInput[1]->Get("elmu__antitop__nominal");
    tInput[10] = (TH1D*) fInput[1]->Get("elmu__tt__nominal");
    tInput[11] = (TH1D*) fInput[1]->Get("elmu__tw__nominal");
    tInput[12] = (TH1D*) fInput[1]->Get("elmu__wzjets__nominal");
    tInput[13] = (TH1D*) fInput[1]->Get("elmu__data_obs__nominal");
    tInput[14] = (TH1D*) fInput[1]->Get("elmu__QCD_DD_EC__nominal");
    tInput[15] = (TH1D*) fInput[1]->Get("elmu__QCD_DD_Barrel__nominal");


    //Fill Tables with the EFT/SM values. [0]->Wilson=-2; [1]->Wilson=-1; [2]->Wilson=0; [3]->Wilson=1; [4]->Wilson=2;
    for(int i=0 ; i<nb_bins ; i++)
    {
        string bin_name = "bin_content_";
        bin_name += to_string(i+1);
        Fit_cosS_cbwi_Input[i] = (TF1*) fInput[2]->Get(bin_name.c_str());
        Wilson_cosS_cbwi_Input[i] = Fit_cosS_cbwi_Input[i]->GetFormula();
    }

    for(int i=0 ; i<nb_bins ; i++)
    {
        string bin_name = "bin_content_";
        bin_name += to_string(i+1);
        Fit_cosS_ctwi_Input[i] = (TF1*) fInput[3]->Get(bin_name.c_str());
        Wilson_cosS_ctwi_Input[i] = Fit_cosS_ctwi_Input[i]->GetFormula();
    }

    for(int i=0 ; i<nb_bins ; i++)
    {
        string bin_name = "bin_content_";
        bin_name += to_string(i+1);
        Fit_phiS_cbwi_Input[i] = (TF1*) fInput[4]->Get(bin_name.c_str());
        Wilson_phiS_cbwi_Input[i] = Fit_phiS_cbwi_Input[i]->GetFormula();
    }

    for(int i=0 ; i<nb_bins ; i++)
    {
        string bin_name = "bin_content_";
        bin_name += to_string(i+1);
        Fit_phiS_ctwi_Input[i] = (TF1*) fInput[5]->Get(bin_name.c_str());
        Wilson_phiS_ctwi_Input[i] = Fit_phiS_ctwi_Input[i]->GetFormula();
    }


    int ChiSquare_bins = 20;
    int bin = 1;
    double c = -2;
    double step = 0.001 , range = 4;
    double nb_iterations = range/step;

    TH1F* ChiSquare_phiS = new TH1F("ChiSquare_phiS", "" , nb_iterations, -range, range);


    for(int j=0 ; j<=nb_iterations ; j++)
    {
        //cout<<"c= "<<c<<endl;
        double EFT_vs_SM = 0;
        double sigma = 0;
        double chiSquare = 0;
        //cout<<"Wilson before for= "<<Wilson<<endl;

        for(int i=1 ; i<=20 ; i++)
        {
            if(EFT == "cbwi")
            {
                EFT_vs_SM = 0;
                EFT_vs_SM = Wilson_phiS_cbwi_Input[i-1]->Eval(c); 
                //cout<<"WIlson = "<<Wilson<<" in "<<i<<endl;         
            }
            else
            {
                EFT_vs_SM = 0;
                EFT_vs_SM = Wilson_phiS_ctwi_Input[i-1]->Eval(c); 
                //cout<<"WIlson = "<<Wilson<<" in "<<i<<endl;         
            }


            sigma = tInput[0]->GetBinContent(i) + tInput[1]->GetBinContent(i) + tInput[2]->GetBinContent(i) + 
                tInput[3]->GetBinContent(i) + tInput[4]->GetBinContent(i) + tInput[6]->GetBinContent(i) + tInput[7]->GetBinContent(i);
                
            chiSquare += (pow(sigma - (EFT_vs_SM*(tInput[0]->GetBinContent(i) + tInput[1]->GetBinContent(i)) + tInput[2]->GetBinContent(i)*syst_tt + 
                tInput[3]->GetBinContent(i)*syst_tw + tInput[4]->GetBinContent(i)*syst_wzjets + tInput[6]->GetBinContent(i)*syst_QCD_EC + tInput[7]->GetBinContent(i)*syst_QCD_Barrel),2))/sigma;
            
        }
 
        //cout<<"tchannel= "<<tchannel<<"\n"<<endl;
        ChiSquare_phiS->SetBinContent(bin,chiSquare/ChiSquare_bins);
        bin++;
        c+=step;
    }



    
    bin = 1;
    c = -2;
    TH1F* ChiSquare_cosThetaS = new TH1F("ChiSquare_cosS", "" , nb_iterations, -range, range);

    for(int j=0 ; j<nb_iterations ; j++)
    {
        //cout<<"Iteration: "<<c<<endl;
        double EFT_vs_SM = 0;
        double sigma = 0;
        double chiSquare = 0;
        //cout<<"Wilson before for= "<<Wilson<<endl;

        for(int i=1; i<=20 ; i++)
        {

            if(EFT == "cbwi")
            {
                EFT_vs_SM = 0;
                EFT_vs_SM = Wilson_cosS_cbwi_Input[i-1]->Eval(c); 
                //cout<<"WIlson = "<<Wilson<<" in "<<i<<endl;         
            }
            else
            {
                EFT_vs_SM = 0;
                EFT_vs_SM = Wilson_cosS_ctwi_Input[i-1]->Eval(c); 
                //cout<<"WIlson = "<<Wilson<<" in "<<i<<endl;         
            }


            sigma = tInput[8]->GetBinContent(i) + tInput[9]->GetBinContent(i) + tInput[10]->GetBinContent(i) + 
                tInput[11]->GetBinContent(i) + tInput[12]->GetBinContent(i) + tInput[14]->GetBinContent(i) + tInput[15]->GetBinContent(i);
                
            chiSquare += (pow(sigma - (EFT_vs_SM*(tInput[8]->GetBinContent(i) + tInput[9]->GetBinContent(i)) + tInput[10]->GetBinContent(i)*syst_tt + 
                tInput[11]->GetBinContent(i)*syst_tw + tInput[12]->GetBinContent(i)*syst_wzjets + tInput[14]->GetBinContent(i)*syst_QCD_EC + tInput[15]->GetBinContent(i)*syst_QCD_Barrel),2))/sigma;
            
        }
        //cout<<"tchannel= "<<tchannel<<endl;
        ChiSquare_cosThetaS->SetBinContent(bin,chiSquare/ChiSquare_bins);
        bin++;
        c+=step;
    }

    int range_fit = 5;

    TFile* file = new TFile("results/Chi2_test_cbtwi.root","RECREATE");
    file->cd();
    ChiSquare_phiS->Write();
    ChiSquare_cosThetaS->Write();

    TF1* fithist_phiS = new TF1 ("pol4","pol4",-range_fit, range_fit);
    fithist_phiS->SetParameters(0,100.0);
    fithist_phiS->SetParameters(1,0.);
    fithist_phiS->SetParameters(2,0.);

    ChiSquare_phiS->Fit(fithist_phiS);
    int min_bin_phiS = ChiSquare_phiS->GetMinimumBin();
    double min_value_phiS = ChiSquare_phiS->GetBinContent(min_bin_phiS) + 1;
    double low_x = fithist_phiS->GetX(min_value_phiS,-5,0);
    double up_x = fithist_phiS->GetX(min_value_phiS,0,5);
    double delta_C_phiS = (up_x - low_x)/2;
    cout<<"Delta "<<EFT<<" for phistar = "<<delta_C_phiS<<endl;

    TF1* fithist_cosS = new TF1 ("pol4","pol4",-range_fit, range_fit);
    fithist_cosS->SetParameters(0,100.0);
    fithist_cosS->SetParameters(1,0.);
    fithist_cosS->SetParameters(2,0.);

    ChiSquare_cosThetaS->Fit(fithist_cosS);
    int min_bin_cosS = ChiSquare_cosThetaS->GetMinimumBin();
    double min_value_cosS = ChiSquare_cosThetaS->GetBinContent(min_bin_cosS) + 1;
    double low_x2 = fithist_cosS->GetX(min_value_cosS,-5,0);
    double up_x2 = fithist_cosS->GetX(min_value_cosS,0,5);
    double delta_C_cos = (up_x2 - low_x2)/2;
    cout<<"Delta "<<EFT<<" for cosThetaStar = "<<delta_C_cos<<endl;


    if(Draw_Graph == true)
    {
        if(EFT == "cbwi") ChiSquare_phiS->GetXaxis()->SetTitle("C_{bW}^{I}");
        if(EFT == "ctwi") ChiSquare_phiS->GetXaxis()->SetTitle("C_{tW}^{I}");
        ChiSquare_phiS->GetYaxis()->SetTitle("#chi^{2}");
        TCanvas* c1 = new TCanvas("Canvas","Canvas");
        ChiSquare_phiS->Draw("");
        fithist_phiS->Draw("SAME");
        string name = "results/ChiSquare_phiS_" + EFT + ".pdf";
        c1->Print(name.c_str());
    }

    if(Draw_Graph == true)
    {
        if(EFT == "cbwi") ChiSquare_cosThetaS->GetXaxis()->SetTitle("C_{bW}^{I}");
        if(EFT == "ctwi") ChiSquare_cosThetaS->GetXaxis()->SetTitle("C_{tW}^{I}");
        ChiSquare_cosThetaS->GetYaxis()->SetTitle("#chi^{2}");
        TCanvas* c2 = new TCanvas("Canvas","Canvas");
        ChiSquare_cosThetaS->Draw("");
        fithist_cosS->Draw("SAME");
        string name2 = "results/ChiSquare_cosThetaS_" + EFT + ".pdf";
        c2->Print(name2.c_str());
    }

    file->Close();

    if(KV == "phi") return delta_C_phiS;
    if(KV == "cos") return delta_C_cos;
}

int main()
{

    double Lumi_2017_syst = 0.0023;

    double tt_normalization_syst = 1.06;
    double tw_normalization_syst = 1.11;
    double wzjets_normalization = 1.1;

    double qcd_normalization_2017_syst = 1.5;

    double syst_tt=1, syst_tw=1, syst_wzjets=1, syst_QCD_EC=1, syst_QCD_Barrel=1;
    //Syst + Stat for Lumi 2017:
    //Up:
    syst_tt +=  Lumi_2017_syst ;
    syst_tw += Lumi_2017_syst ;
    syst_wzjets += Lumi_2017_syst ;
    syst_QCD_EC += Lumi_2017_syst ;
    syst_QCD_Barrel += Lumi_2017_syst ;
    double stat_syst_Lumi_2017_up = Stat_and_Syst(syst_tt, syst_tw, syst_wzjets, syst_QCD_EC, syst_QCD_Barrel);

    syst_tt=1; syst_tw=1; syst_wzjets=1; syst_QCD_EC=1; syst_QCD_Barrel=1;
    //Down:
    syst_tt -=  Lumi_2017_syst ;
    syst_tw -= Lumi_2017_syst ;
    syst_wzjets -= Lumi_2017_syst ;
    syst_QCD_EC -= Lumi_2017_syst ;
    syst_QCD_Barrel -= Lumi_2017_syst ;
    double stat_syst_Lumi_2017_down = Stat_and_Syst(syst_tt, syst_tw, syst_wzjets, syst_QCD_EC, syst_QCD_Barrel);

    cout<<stat_syst_Lumi_2017_up<<";"<<stat_syst_Lumi_2017_down<<endl;

    double max_Lumi_2017 = max(abs(stat_syst_Lumi_2017_up),abs(stat_syst_Lumi_2017_down));
    cout<<"max_Lumi_2017="<<max_Lumi_2017<<endl;

    double syst_Lumi_2017 = sqrt( abs( pow(max_Lumi_2017,2) - pow (2.23841,2) ) );
    cout<<"syst_Lumi_2017="<<syst_Lumi_2017<<endl;
   
    return 0;    
}



/*

    TCanvas* c1 = new TCanvas("Canvas","Canvas");
    tInput[4]->Draw();
    string name = "results/top2.pdf";
    c1->Print(name.c_str());


    double Asimov_phiS = 0;
    double Asimov_cosS = 0;
    for(int i=0 ; i<5 ; i++)
    {
        Asimov_phiS+= phiS_top[i] + phiS_tt[i] + phiS_antitop[i] +  phiS_tW[i] + phiS_wzjets[i] + phiS_data[i] + phiS_QCD_EC[i] + phiS_QCD_Barrel[i];
        Asimov_cosS+= cosS_top[i] + cosS_tt[i] + cosS_antitop[i] +  cosS_tW[i] + cosS_wzjets[i] + cosS_data[i] + cosS_QCD_EC[i] + cosS_QCD_Barrel[i];
    }



    int k = 2;
    int j = 1;
    for(int i=0 ; i<20 ; i++)
    {
        string bin_name = "bin_content_part_";
        bin_name += to_string(j);
        gInput[i] = (TGraph*) fInput[k]->Get(bin_name.c_str());
        //cout<<"bin_name="<<bin_name<<endl;
        //cout<<"k="<<k<<endl;
        if(i==4) k++;
        if(i==9) k++;
        if(i==14) k++;
        j++;
        if(j==6) j=1;
    }

    for(int i=0 ; i<16 ; i++) //Rebin Phistar and CosThetaStar Histos to 5 bins 
    {
        tInput[i]->Rebin(4);
    }
    
    int bin1 = 1; int bin2 = 2; int bin3 = 3; int bin4 = 4; int bin5 = 5;

    //Fill tables with phiS and cosS Histograms values per bin per channel. Channel={signal,noise}    
    double phiS_top[5] = {tInput[0]->GetBinContent(bin1), tInput[0]->GetBinContent(bin2), tInput[0]->GetBinContent(bin3), tInput[0]->GetBinContent(bin4), tInput[0]->GetBinContent(bin5)};
    double phiS_tt[5] = {tInput[1]->GetBinContent(bin1), tInput[1]->GetBinContent(bin2), tInput[1]->GetBinContent(bin3), tInput[1]->GetBinContent(bin4), tInput[1]->GetBinContent(bin5)};
    double phiS_antitop[5] = {tInput[2]->GetBinContent(bin1), tInput[2]->GetBinContent(bin2), tInput[2]->GetBinContent(bin3), tInput[2]->GetBinContent(bin4), tInput[2]->GetBinContent(bin5)};
    double phiS_tW[5] = {tInput[3]->GetBinContent(bin1), tInput[3]->GetBinContent(bin2), tInput[3]->GetBinContent(bin3), tInput[3]->GetBinContent(bin4), tInput[3]->GetBinContent(bin5)};
    double phiS_wzjets[5] = {tInput[4]->GetBinContent(bin1), tInput[4]->GetBinContent(bin2), tInput[4]->GetBinContent(bin3), tInput[4]->GetBinContent(bin4), tInput[4]->GetBinContent(bin5)};
    double phiS_data[5] = {tInput[5]->GetBinContent(bin1), tInput[5]->GetBinContent(bin2), tInput[5]->GetBinContent(bin3), tInput[5]->GetBinContent(bin4), tInput[5]->GetBinContent(bin5)};
    double phiS_QCD_EC[5] = {tInput[6]->GetBinContent(bin1), tInput[6]->GetBinContent(bin2), tInput[6]->GetBinContent(bin3), tInput[6]->GetBinContent(bin4), tInput[6]->GetBinContent(bin5)};
    double phiS_QCD_Barrel[5] = {tInput[7]->GetBinContent(bin1), tInput[7]->GetBinContent(bin2), tInput[7]->GetBinContent(bin3), tInput[7]->GetBinContent(bin4), tInput[7]->GetBinContent(bin5)};

    double cosS_top[5] = {tInput[8]->GetBinContent(bin1), tInput[8]->GetBinContent(bin2), tInput[8]->GetBinContent(bin3), tInput[8]->GetBinContent(bin4), tInput[8]->GetBinContent(bin5)};
    double cosS_tt[5] = {tInput[9]->GetBinContent(bin1), tInput[9]->GetBinContent(bin2), tInput[9]->GetBinContent(bin3), tInput[9]->GetBinContent(bin4), tInput[9]->GetBinContent(bin5)};
    double cosS_antitop[5] = {tInput[10]->GetBinContent(bin1), tInput[10]->GetBinContent(bin2), tInput[10]->GetBinContent(bin3), tInput[10]->GetBinContent(bin4), tInput[10]->GetBinContent(bin5)};
    double cosS_tW[5] = {tInput[11]->GetBinContent(bin1), tInput[11]->GetBinContent(bin2), tInput[11]->GetBinContent(bin3), tInput[11]->GetBinContent(bin4), tInput[11]->GetBinContent(bin5)};
    double cosS_wzjets[5] = {tInput[12]->GetBinContent(bin1), tInput[12]->GetBinContent(bin2), tInput[12]->GetBinContent(bin3), tInput[12]->GetBinContent(bin4), tInput[12]->GetBinContent(bin5)};
    double cosS_data[5] = {tInput[13]->GetBinContent(bin1), tInput[13]->GetBinContent(bin2), tInput[13]->GetBinContent(bin3), tInput[13]->GetBinContent(bin4), tInput[13]->GetBinContent(bin5)};
    double cosS_QCD_EC[5] = {tInput[14]->GetBinContent(bin1), tInput[14]->GetBinContent(bin2), tInput[14]->GetBinContent(bin3), tInput[14]->GetBinContent(bin4), tInput[14]->GetBinContent(bin5)};
    double cosS_QCD_Barrel[5] = {tInput[15]->GetBinContent(bin1), tInput[15]->GetBinContent(bin2), tInput[15]->GetBinContent(bin3), tInput[15]->GetBinContent(bin4), tInput[15]->GetBinContent(bin5)};

                if(EFT=="cbwi")  //Slection of the correct ROOT file and CbWi value
            {
                Wilson = EFT_phiS_cbwi_Input[i]->GetPointY(c);
                //cout<<"Wilson in "<<i<<" = "<<Wilson<<endl;
            }
            else  //Slection of the correct ROOT file and CtWi value
            {
                Wilson = EFT_phiS_ctwi_Input[i]->GetPointY(c);
                //cout<<"Wilson in "<<i<<" = "<<Wilson<<endl;
            }

    double x = T->Eval(-1.5);

*/