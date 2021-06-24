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
#include <iostream>
#include <fstream>

using namespace std;



double Stat_and_Syst(string EFT, string KV, TH1D* phi_Input[7], TH1D* cos_Input[7], double syst_tt, double syst_tw, double syst_wzjets, double syst_QCD_EC, double syst_QCD_Barrel, bool Draw_Graph, TFormula* Wilson_phiS_cbwi_Input[20], TFormula* Wilson_phiS_ctwi_Input[20], TFormula* Wilson_cosS_cbwi_Input[20], TFormula* Wilson_cosS_ctwi_Input[20])
{
    int nb_bins = 20;
    int bin = 1;
    double c = -2;
    double step = 0.001 , range = 4;
    double nb_iterations = range/step;
    int range_fit = 5;


    if(KV == "phi")
    {
        cout<<"PhiStar_"<<EFT<<endl;
        TH1F* ChiSquare_phiS = new TH1F("ChiSquare_phiS", "" , nb_iterations, -range, range);

        //PhiStar
        for(int j=0 ; j<=nb_iterations ; j++)
        {
            //cout<<"c= "<<c<<endl;
            double EFT_vs_SM = 0;
            double sigma = 0;
            double chiSquare = 0;
            //cout<<"Wilson before for= "<<Wilson<<endl;

            for(int i=1 ; i<=nb_bins ; i++)
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


                sigma = phi_Input[0]->GetBinContent(i) + phi_Input[1]->GetBinContent(i) + phi_Input[2]->GetBinContent(i) + 
                    phi_Input[3]->GetBinContent(i) + phi_Input[4]->GetBinContent(i) + phi_Input[5]->GetBinContent(i) + phi_Input[6]->GetBinContent(i);
                
                chiSquare += (pow(sigma - (EFT_vs_SM*(phi_Input[0]->GetBinContent(i) + phi_Input[1]->GetBinContent(i)) + phi_Input[2]->GetBinContent(i)*syst_tt + 
                    phi_Input[3]->GetBinContent(i)*syst_tw + phi_Input[4]->GetBinContent(i)*syst_wzjets + phi_Input[5]->GetBinContent(i)*syst_QCD_EC + phi_Input[6]->GetBinContent(i)*syst_QCD_Barrel),2))/sigma;
            
            }
 
            //cout<<"tchannel= "<<tchannel<<"\n"<<endl;
            ChiSquare_phiS->SetBinContent(bin,chiSquare/nb_bins);
            bin++;
            c+=step;
        }

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

        return delta_C_phiS;
    }

    if(KV == "cos")
    {
        cout<<"CosThetaStar_"<<EFT<<endl;
        TH1F* ChiSquare_cosThetaS = new TH1F("ChiSquare_cosS", "" , nb_iterations, -range, range);

        for(int j=0 ; j<nb_iterations ; j++)
        {
            //cout<<"Iteration: "<<c<<endl;
            double EFT_vs_SM = 0;
            double sigma = 0;
            double chiSquare = 0;
            //cout<<"Wilson before for= "<<Wilson<<endl;

            for(int i=1; i<=nb_bins ; i++)
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


                sigma = cos_Input[0]->GetBinContent(i) + cos_Input[1]->GetBinContent(i) + cos_Input[2]->GetBinContent(i) + 
                    cos_Input[3]->GetBinContent(i) + cos_Input[4]->GetBinContent(i) + cos_Input[5]->GetBinContent(i) + cos_Input[6]->GetBinContent(i);
                
                chiSquare += (pow(sigma - (EFT_vs_SM*(cos_Input[0]->GetBinContent(i) + cos_Input[1]->GetBinContent(i)) + cos_Input[2]->GetBinContent(i)*syst_tt + 
                    cos_Input[3]->GetBinContent(i)*syst_tw + cos_Input[4]->GetBinContent(i)*syst_wzjets + cos_Input[5]->GetBinContent(i)*syst_QCD_EC + cos_Input[6]->GetBinContent(i)*syst_QCD_Barrel),2))/sigma;
            
            }
            //cout<<"tchannel= "<<tchannel<<endl;
            ChiSquare_cosThetaS->SetBinContent(bin,chiSquare/nb_bins);
            bin++;
            c+=step;
        }

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

        return delta_C_cos;
    } 

    return 0;
}


int main()
{

    ////////////////////////////Import all TFormulas to get EFT/SM ratio per value of Wilson Coeff////////////////////////////
    int nb_bins = 20;
    int EFT_nb_files = 4;

    string EFT_file[EFT_nb_files];
    EFT_file[0] = "signal_proc_cosThetaStar_cbwi_elmu_20bins_Rwgt_cbwi2p5.root";
    EFT_file[1] = "signal_proc_cosThetaStar_ctwi_elmu_20bins_Rwgt_cbwi2p5.root";
    EFT_file[2] = "signal_proc_PhiStar_cbwi_elmu_20bins_Rwgt_cbwi2p5.root";
    EFT_file[3] = "signal_proc_PhiStar_ctwi_elmu_20bins_Rwgt_cbwi2p5.root";



    //Import TF1 fits and formulas for Wilson coefficients:
    TF1* Fit_phiS_cbwi_Input[nb_bins];
    TF1* Fit_phiS_ctwi_Input[nb_bins];
    TF1* Fit_cosS_cbwi_Input[nb_bins];
    TF1* Fit_cosS_ctwi_Input[nb_bins];
    TFormula* Wilson_phiS_cbwi_Input[nb_bins];
    TFormula* Wilson_phiS_ctwi_Input[nb_bins];
    TFormula* Wilson_cosS_cbwi_Input[nb_bins];
    TFormula* Wilson_cosS_ctwi_Input[nb_bins];
    string inputDirectory; //pwd to .root files directory


    TFile* EFT_fInput[EFT_nb_files]; //Import .root files
    
    for(int i=0 ; i<EFT_nb_files ; i++)
    {
        inputDirectory = "data/" + EFT_file[i];
        EFT_fInput[i] = new TFile(inputDirectory.c_str(),"READ");
    }

    //Fill Tables with the EFT/SM values. [0]->Wilson=-2; [1]->Wilson=-1; [2]->Wilson=0; [3]->Wilson=1; [4]->Wilson=2;
    string bin_name;
    for(int i=0 ; i<nb_bins ; i++)
    {
        bin_name = "bin_content_";
        bin_name += to_string(i+1);
        Fit_cosS_cbwi_Input[i] = (TF1*) EFT_fInput[0]->Get(bin_name.c_str());
        Wilson_cosS_cbwi_Input[i] = Fit_cosS_cbwi_Input[i]->GetFormula();
    }

    for(int i=0 ; i<nb_bins ; i++)
    {
        bin_name = "bin_content_";
        bin_name += to_string(i+1);
        Fit_cosS_ctwi_Input[i] = (TF1*) EFT_fInput[1]->Get(bin_name.c_str());
        Wilson_cosS_ctwi_Input[i] = Fit_cosS_ctwi_Input[i]->GetFormula();
    }

    for(int i=0 ; i<nb_bins ; i++)
    {
        bin_name = "bin_content_";
        bin_name += to_string(i+1);
        Fit_phiS_cbwi_Input[i] = (TF1*) EFT_fInput[2]->Get(bin_name.c_str());
        Wilson_phiS_cbwi_Input[i] = Fit_phiS_cbwi_Input[i]->GetFormula();
    }

    for(int i=0 ; i<nb_bins ; i++)
    {
        bin_name = "bin_content_";
        bin_name += to_string(i+1);
        Fit_phiS_ctwi_Input[i] = (TF1*) EFT_fInput[3]->Get(bin_name.c_str());
        Wilson_phiS_ctwi_Input[i] = Fit_phiS_ctwi_Input[i]->GetFormula();
    }


    ////////////////////////////Import STreco data////////////////////////////
    int nbfiles = 2;
    string suffix[nbfiles];
    
    suffix[0] = "hist_reco_phiStar.root";
    suffix[1] = "hist_reco_cosThetaStar.root";

    TFile* fInput[nbfiles]; //Import .root files
    TH1D* nominal_Input_phi[8]; //Import signals and background noise nominals for PhiStar
    TH1D* nominal_Input_cos[8]; //Import signals and background noise nominals for CosThetaStar
    //TH1D* syst_Input_phi[40]; //Import signals and background noise Systematic Uncertainties

    for(int i=0 ; i<nbfiles ; i++)
    {
        inputDirectory = "data/" + suffix[i];
        fInput[i] = new TFile(inputDirectory.c_str(),"READ");
    }
    //Fill data from STreco
    //nominals = stats uncertainties
    nominal_Input_phi[0] = (TH1D*) fInput[0]->Get("elmu__top__nominal");
    nominal_Input_phi[1] = (TH1D*) fInput[0]->Get("elmu__antitop__nominal");
    nominal_Input_phi[2] = (TH1D*) fInput[0]->Get("elmu__tt__nominal");
    nominal_Input_phi[3] = (TH1D*) fInput[0]->Get("elmu__tw__nominal");
    nominal_Input_phi[4] = (TH1D*) fInput[0]->Get("elmu__wzjets__nominal");
    nominal_Input_phi[5] = (TH1D*) fInput[0]->Get("elmu__QCD_DD_EC__nominal");
    nominal_Input_phi[6] = (TH1D*) fInput[0]->Get("elmu__QCD_DD_Barrel__nominal");

    nominal_Input_cos[0] = (TH1D*) fInput[1]->Get("elmu__top__nominal");
    nominal_Input_cos[1] = (TH1D*) fInput[1]->Get("elmu__antitop__nominal");
    nominal_Input_cos[2] = (TH1D*) fInput[1]->Get("elmu__tt__nominal");
    nominal_Input_cos[3] = (TH1D*) fInput[1]->Get("elmu__tw__nominal");
    nominal_Input_cos[4] = (TH1D*) fInput[1]->Get("elmu__wzjets__nominal");
    nominal_Input_cos[5] = (TH1D*) fInput[1]->Get("elmu__QCD_DD_EC__nominal");
    nominal_Input_cos[6] = (TH1D*) fInput[1]->Get("elmu__QCD_DD_Barrel__nominal");



    ofstream out_file;
    out_file.open("Uncerts.txt");
    out_file << "Only Systematic Uncertainties:\n";
    
    //Compute Statistical uncert:
    double delta_cbwi_phi = Stat_and_Syst("cbwi", "phi", nominal_Input_phi, nominal_Input_cos, 1, 1, 1, 1, 1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);
    cout<<"delta_cbwi_phiStar = "<<delta_cbwi_phi<<endl;

    double delta_ctwi_phi = Stat_and_Syst("ctwi", "phi", nominal_Input_phi, nominal_Input_cos, 1, 1, 1, 1, 1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);
    cout<<"delta_ctwi_phiStar = "<<delta_ctwi_phi<<endl;

    double delta_cbwi_cos = Stat_and_Syst("cbwi", "cos", nominal_Input_phi, nominal_Input_cos, 1, 1, 1, 1, 1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);
    cout<<"delta_cbwi_cosThetaStar = "<<delta_cbwi_cos<<endl;

    double delta_ctwi_cos = Stat_and_Syst("ctwi", "cos", nominal_Input_phi, nominal_Input_cos, 1, 1, 1, 1, 1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);
    cout<<"delta_ctwi_cosThetaStar = "<<delta_ctwi_cos<<endl;


    out_file << "delta_cbwi_phi = " << delta_cbwi_phi << "\n";
    out_file << "delta_ctwi_phi = " << delta_ctwi_phi << "\n";
    out_file << "delta_cbwi_cos = " << delta_cbwi_cos << "\n";
    out_file << "delta_ctwi_cos = " << delta_ctwi_cos << "\n";

/*

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
    double stat_syst_Lumi_2017_up = Stat_and_Syst(EFT, KV, syst_tt, syst_tw, syst_wzjets, syst_QCD_EC, syst_QCD_Barrel, Draw_Graph);

    syst_tt=1; syst_tw=1; syst_wzjets=1; syst_QCD_EC=1; syst_QCD_Barrel=1;
    //Down:
    syst_tt -=  Lumi_2017_syst ;
    syst_tw -= Lumi_2017_syst ;
    syst_wzjets -= Lumi_2017_syst ;
    syst_QCD_EC -= Lumi_2017_syst ;
    syst_QCD_Barrel -= Lumi_2017_syst ;
    double stat_syst_Lumi_2017_down = Stat_and_Syst(EFT, KV, syst_tt, syst_tw, syst_wzjets, syst_QCD_EC, syst_QCD_Barrel, Draw_Graph);

    cout<<stat_syst_Lumi_2017_up<<";"<<stat_syst_Lumi_2017_down<<endl;

    double max_Lumi_2017 = max(abs(stat_syst_Lumi_2017_up),abs(stat_syst_Lumi_2017_down));
    cout<<"max_Lumi_2017="<<max_Lumi_2017<<endl;

    double syst_Lumi_2017 = sqrt( abs( pow(max_Lumi_2017,2) - pow (2.23841,2) ) );
    cout<<"syst_Lumi_2017="<<syst_Lumi_2017<<endl;


    myfile.close();

    */

    return 0;    
}


