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
            ChiSquare_cosThetaS->SetBinContent(bin, chiSquare/nb_bins);
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

        return delta_C_cos;
    } 

    return 0;
}


int main()
{

    ////////////////////////////Import all TFormulas to get EFT/SM ratio per value of Wilson Coeff////////////////////////////
    int nb_bins = 20;
    int EFT_nb_files = 4;
    int chiffres_signif_s= 2;
    //int chiffres_signif_ss = 3;

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
    //Import .root files
    TFile* phi_Input;
    TFile* cos_Input; 
    TH1D* nominal_Input_phi[7]; //Import signals and background noise nominals for PhiStar
    TH1D* nominal_Input_cos[7]; //Import signals and background noise nominals for CosThetaStar
    
    phi_Input = new TFile(("data/" + suffix[0]).c_str(),"READ");
    cos_Input = new TFile(("data/" + suffix[1]).c_str(),"READ");

    //Fill data from STreco
    //nominals = stats uncertainties
    nominal_Input_phi[0] = (TH1D*) phi_Input->Get("elmu__top__nominal");
    nominal_Input_phi[1] = (TH1D*) phi_Input->Get("elmu__antitop__nominal");
    nominal_Input_phi[2] = (TH1D*) phi_Input->Get("elmu__tt__nominal");
    nominal_Input_phi[3] = (TH1D*) phi_Input->Get("elmu__tw__nominal");
    nominal_Input_phi[4] = (TH1D*) phi_Input->Get("elmu__wzjets__nominal");
    nominal_Input_phi[5] = (TH1D*) phi_Input->Get("elmu__QCD_DD_EC__nominal");
    nominal_Input_phi[6] = (TH1D*) phi_Input->Get("elmu__QCD_DD_Barrel__nominal");

    nominal_Input_cos[0] = (TH1D*) cos_Input->Get("elmu__top__nominal");
    nominal_Input_cos[1] = (TH1D*) cos_Input->Get("elmu__antitop__nominal");
    nominal_Input_cos[2] = (TH1D*) cos_Input->Get("elmu__tt__nominal");
    nominal_Input_cos[3] = (TH1D*) cos_Input->Get("elmu__tw__nominal");
    nominal_Input_cos[4] = (TH1D*) cos_Input->Get("elmu__wzjets__nominal");
    nominal_Input_cos[5] = (TH1D*) cos_Input->Get("elmu__QCD_DD_EC__nominal");
    nominal_Input_cos[6] = (TH1D*) cos_Input->Get("elmu__QCD_DD_Barrel__nominal");



    ofstream out_file;
    out_file.open("Uncerts.txt");
    out_file << "LaTeX table:\n\n";

    out_file << " " << "   &   " << "C_{tW}^{I}" << "   &   " << "C_{bW}^{I}" << "\\" << "\\" << "\n";
    out_file << "\\hline\n";

    
    //Compute Statistical uncert:
    double delta_cbwi_phi_stat = Stat_and_Syst("cbwi", "phi", nominal_Input_phi, nominal_Input_cos, 1, 1, 1, 1, 1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);
    cout<<"delta_cbwi_phiStar_stat = "<<delta_cbwi_phi_stat<<endl;

    double delta_ctwi_phi_stat = Stat_and_Syst("ctwi", "phi", nominal_Input_phi, nominal_Input_cos, 1, 1, 1, 1, 1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);
    cout<<"delta_ctwi_phiStar_stat = "<<delta_ctwi_phi_stat<<endl;

    double delta_cbwi_cos_stat = Stat_and_Syst("cbwi", "cos", nominal_Input_phi, nominal_Input_cos, 1, 1, 1, 1, 1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);
    cout<<"delta_cbwi_cosThetaStar_stat = "<<delta_cbwi_cos_stat<<endl;

    double delta_ctwi_cos_stat = Stat_and_Syst("ctwi", "cos", nominal_Input_phi, nominal_Input_cos, 1, 1, 1, 1, 1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);
    cout<<"delta_ctwi_cosThetaStar_stat = "<<delta_ctwi_cos_stat<<endl;


    //out_file << "delta_cbwi_phi_stat = " << delta_cbwi_phi_stat << "\n";
    //out_file << "delta_ctwi_phi_stat = " << delta_ctwi_phi_stat << "\n";
    //out_file << "delta_cbwi_cos_stat = " << delta_cbwi_cos_stat << "\n";
    //out_file << "delta_ctwi_cos_stat = " << delta_ctwi_cos_stat << "\n";

    out_file << "$\\sigma_{stat}$:" << "   &   " << std::setprecision(chiffres_signif_s) << min(delta_ctwi_phi_stat,delta_ctwi_cos_stat) << "   &   " << std::setprecision(chiffres_signif_s) << min(delta_cbwi_phi_stat,delta_ctwi_cos_stat) << "\\" << "\\" << "\n"; //"   &   " << std::setprecision(chiffres_signif_ss) << max_value_between_cosT_and_phiS_ctwi_ss << "   &   " << std::setprecision(chiffres_signif_ss) << max_value_between_cosT_and_phiS_cbwi_ss << "\\" << "\\" << "\n";
    out_file << "\\hline\n";
    

    ////////////////////////////////////Systematics///////////////////////////////////
    int nb_of_systs = 31;
    out_file << "$\\sigma_{syst}$:" << "   &   " << " " << "   &   " << " " << "\\" << "\\" << "\n";
    string prefix = "elmu";
    string physics_process[7] = {"__top__", "__antitop__", "__tt__", "__tw__", "__wzjets__", "__QCD_DD_EC__","__QCD_DD_Barrel__"};
    
    string syst_unc[nb_of_systs] = {"bsfHFcorr", "bsfHFuncorr", "bsfLF", "bsf", "elsf", "jer", "jesAbsolute", "jesAbsolute_2017", "jesBBEC1", "jesBBEC1_2017",
        "jesEC2", "jesEC2_2017", "jesFlavorMerged", "jesHF", "jesHF_2017", "jesRelativeBal", "jesRelativeSample_2017", "jesTotal", "musf", "pdf", 
            "prefire", "pu", "ewk_q2", "single_top_q2", "tchannel_isr", "top_fsr", "top_pt", "tt_isr", "tt_q2", "tw_isr", "unclustEn"};
    
    string Up_Down[2] = {"Up", "Down"};

    TH1D* syst_Input_phi[7];
    TH1D* syst_Input_cos[7];

    double stat_syst_cbwi_phi[2];
    double stat_syst_ctwi_phi[2];
    double stat_syst_cbwi_cos[2];
    double stat_syst_ctwi_cos[2];

    double max_cbwi_phi;
    double max_ctwi_phi;
    double max_cbwi_cos;
    double max_ctwi_cos;

    double Systematic_cbwi_phi;
    double Systematic_ctwi_phi;
    double Systematic_cbwi_cos;
    double Systematic_ctwi_cos;

    double max_value_between_cosT_and_phiS_cbwi_ss;
    double max_value_between_cosT_and_phiS_ctwi_ss;
    double max_value_between_cosT_and_phiS_cbwi_s;
    double max_value_between_cosT_and_phiS_ctwi_s;


    for(int i=0 ; i<nb_of_systs ; i++) //iteration over systematics
    {
        for(int j=0 ; j<2 ; j++)
        {
            syst_Input_phi[0] = (TH1D*) phi_Input->Get((prefix + physics_process[0] + syst_unc[i] + Up_Down[j]).c_str());
            syst_Input_phi[1] = (TH1D*) phi_Input->Get((prefix + physics_process[1] + syst_unc[i] + Up_Down[j]).c_str());
            syst_Input_phi[2] = (TH1D*) phi_Input->Get((prefix + physics_process[2] + syst_unc[i] + Up_Down[j]).c_str());
            syst_Input_phi[3] = (TH1D*) phi_Input->Get((prefix + physics_process[3] + syst_unc[i] + Up_Down[j]).c_str());
            syst_Input_phi[4] = (TH1D*) phi_Input->Get((prefix + physics_process[4] + syst_unc[i] + Up_Down[j]).c_str());
            syst_Input_phi[5] = (TH1D*) phi_Input->Get("elmu__QCD_DD_EC__nominal");
            syst_Input_phi[6] = (TH1D*) phi_Input->Get("elmu__QCD_DD_Barrel__nominal");

            syst_Input_cos[0] = (TH1D*) cos_Input->Get((prefix + physics_process[0] + syst_unc[i] + Up_Down[j]).c_str());
            syst_Input_cos[1] = (TH1D*) cos_Input->Get((prefix + physics_process[1] + syst_unc[i] + Up_Down[j]).c_str());
            syst_Input_cos[2] = (TH1D*) cos_Input->Get((prefix + physics_process[2] + syst_unc[i] + Up_Down[j]).c_str());
            syst_Input_cos[3] = (TH1D*) cos_Input->Get((prefix + physics_process[3] + syst_unc[i] + Up_Down[j]).c_str());
            syst_Input_cos[4] = (TH1D*) cos_Input->Get((prefix + physics_process[4] + syst_unc[i] + Up_Down[j]).c_str());
            syst_Input_cos[5] = (TH1D*) cos_Input->Get("elmu__QCD_DD_EC__nominal");
            syst_Input_cos[6] = (TH1D*) cos_Input->Get("elmu__QCD_DD_Barrel__nominal");
            
            ///////Listing Exceptions//////////
            if(i == 22)
            {
                syst_Input_phi[2] = (TH1D*) phi_Input->Get("elmu__tt__nominal");
                syst_Input_cos[2] = (TH1D*) cos_Input->Get("elmu__tt__nominal");
            }

            if(i == 23 || i == 24)
            {
                syst_Input_phi[2] = (TH1D*) phi_Input->Get("elmu__tt__nominal");
                syst_Input_phi[3] = (TH1D*) phi_Input->Get("elmu__tw__nominal");
                syst_Input_phi[4] = (TH1D*) phi_Input->Get("elmu__wzjets__nominal");

                syst_Input_cos[2] = (TH1D*) cos_Input->Get("elmu__tt__nominal");
                syst_Input_cos[3] = (TH1D*) cos_Input->Get("elmu__tw__nominal");
                syst_Input_cos[4] = (TH1D*) cos_Input->Get("elmu__wzjets__nominal");
            }

            if(i == 25)
            {
                syst_Input_phi[4] = (TH1D*) phi_Input->Get("elmu__wzjets__nominal");
                syst_Input_cos[4] = (TH1D*) cos_Input->Get("elmu__wzjets__nominal");
            }

            if(i == 26 || i == 27 || i == 28)
            {
                syst_Input_phi[0] = (TH1D*) phi_Input->Get("elmu__top__nominal");
                syst_Input_phi[1] = (TH1D*) phi_Input->Get("elmu__antitop__nominal");
                syst_Input_phi[3] = (TH1D*) phi_Input->Get("elmu__tw__nominal");
                syst_Input_phi[4] = (TH1D*) phi_Input->Get("elmu__wzjets__nominal");

                syst_Input_cos[0] = (TH1D*) cos_Input->Get("elmu__top__nominal");
                syst_Input_cos[1] = (TH1D*) cos_Input->Get("elmu__antitop__nominal");
                syst_Input_cos[3] = (TH1D*) cos_Input->Get("elmu__tw__nominal");
                syst_Input_cos[4] = (TH1D*) cos_Input->Get("elmu__wzjets__nominal");
            }

            if(i == 29)
            {
                syst_Input_phi[0] = (TH1D*) phi_Input->Get("elmu__top__nominal");
                syst_Input_phi[1] = (TH1D*) phi_Input->Get("elmu__antitop__nominal");
                syst_Input_phi[2] = (TH1D*) phi_Input->Get("elmu__tt__nominal");
                syst_Input_phi[4] = (TH1D*) phi_Input->Get("elmu__wzjets__nominal");

                syst_Input_cos[0] = (TH1D*) cos_Input->Get("elmu__top__nominal");
                syst_Input_cos[1] = (TH1D*) cos_Input->Get("elmu__antitop__nominal");
                syst_Input_cos[2] = (TH1D*) cos_Input->Get("elmu__tt__nominal");
                syst_Input_cos[4] = (TH1D*) cos_Input->Get("elmu__wzjets__nominal");     
            }
            
            stat_syst_cbwi_phi[j] = Stat_and_Syst("cbwi", "phi", syst_Input_phi, syst_Input_cos,1,1,1,1,1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);
            stat_syst_ctwi_phi[j] = Stat_and_Syst("ctwi", "phi", syst_Input_phi, syst_Input_cos,1,1,1,1,1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);
            stat_syst_cbwi_cos[j] = Stat_and_Syst("cbwi", "cos", syst_Input_phi, syst_Input_cos,1,1,1,1,1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);
            stat_syst_ctwi_cos[j] = Stat_and_Syst("ctwi", "cos", syst_Input_phi, syst_Input_cos,1,1,1,1,1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);

        }

        //syst_unc[i].replace(syst_unc[i].find("_"),1,"\\_");

        max_cbwi_phi = max(abs(stat_syst_cbwi_phi[0]),abs(stat_syst_cbwi_phi[1]));
        max_ctwi_phi = max(abs(stat_syst_ctwi_phi[0]),abs(stat_syst_ctwi_phi[1]));
        max_cbwi_cos = max(abs(stat_syst_cbwi_cos[0]),abs(stat_syst_cbwi_cos[1]));
        max_ctwi_cos = max(abs(stat_syst_ctwi_cos[0]),abs(stat_syst_ctwi_cos[1]));

        //out_file <<"(Stat + Syst) delta_cbwi_phiStar for "<<syst_unc[i]<< " = " << max_cbwi_phi << "\n";
        //out_file <<"(Stat + Syst) delta_ctwi_phiStar for "<<syst_unc[i]<< " = " << max_ctwi_phi << "\n";
        //out_file <<"(Stat + Syst) delta_cbwi_cosThetaStar for "<<syst_unc[i]<< " = " << max_cbwi_cos << "\n";
        //out_file <<"(Stat + Syst) delta_ctwi_cosThetaStar for "<<syst_unc[i]<< " = " << max_ctwi_cos << "\n" << "\n";

        max_value_between_cosT_and_phiS_cbwi_ss = min(max_cbwi_phi,max_cbwi_cos);
        //out_file<<"\n (Stat + Syst) Min Value for CbWi in " << syst_unc[i] << " = " << max_value_between_cosT_and_phiS_cbwi_ss << "\n";
        max_value_between_cosT_and_phiS_ctwi_ss = min(max_ctwi_phi,max_ctwi_cos);
        //out_file<<"\n (Stat + Syst) Min Value for CtWi in " << syst_unc[i] << " = " << max_value_between_cosT_and_phiS_ctwi_ss << "\n";

        Systematic_cbwi_phi = sqrt( abs( pow(max_cbwi_phi,2) - pow (delta_cbwi_phi_stat,2) ) );
        Systematic_ctwi_phi = sqrt( abs( pow(max_ctwi_phi,2) - pow (delta_ctwi_phi_stat,2) ) );
        Systematic_cbwi_cos = sqrt( abs( pow(max_cbwi_cos,2) - pow (delta_cbwi_cos_stat,2) ) );
        Systematic_ctwi_cos = sqrt( abs( pow(max_ctwi_cos,2) - pow (delta_ctwi_cos_stat,2) ) );

        //out_file <<"(Syst) delta_cbwi_phiStar for "<<syst_unc[i]<< " = " << Systematic_cbwi_phi << "\n";
        //out_file <<"(Syst) delta_ctwi_phiStar for "<<syst_unc[i]<< " = " << Systematic_ctwi_phi << "\n";
        //out_file <<"(Syst) delta_cbwi_cosThetaStar for "<<syst_unc[i]<< " = " << Systematic_cbwi_cos << "\n";
        //out_file <<"(Syst) delta_ctwi_cosThetaStar for "<<syst_unc[i]<< " = " << Systematic_ctwi_cos << "\n" << "\n";

        max_value_between_cosT_and_phiS_cbwi_s = min(Systematic_cbwi_phi,Systematic_cbwi_cos);
        //out_file<<"\n (Syst) Min Value for CbWi in " << syst_unc[i] << " = " << max_value_between_cosT_and_phiS_cbwi << "\n";
        max_value_between_cosT_and_phiS_ctwi_s = min(Systematic_ctwi_phi,Systematic_ctwi_cos);
        //out_file<<"\n (Syst) Min Value for CtWi in " << syst_unc[i] << " = " << max_value_between_cosT_and_phiS_ctwi << "\n \n \n";

        out_file << syst_unc[i] << "   &   " << std::setprecision(chiffres_signif_s) << max_value_between_cosT_and_phiS_ctwi_s << "   &   " << std::setprecision(chiffres_signif_s) << max_value_between_cosT_and_phiS_cbwi_s << "\\" << "\\" << "\n"; //"   &   " << std::setprecision(chiffres_signif_ss) << max_value_between_cosT_and_phiS_ctwi_ss << "   &   " << std::setprecision(chiffres_signif_ss) << max_value_between_cosT_and_phiS_cbwi_ss << "\\" << "\\" << "\n";

    }

    double Lumi_2017_syst = 0.0023;
    double tt_normalization_syst = 0.06;
    double tw_normalization_syst = 0.11;
    double wzjets_normalization_syst = 0.1;
    double qcd_normalization_2017_syst = 0.15;

    double Errors[5] = {Lumi_2017_syst, tt_normalization_syst, tw_normalization_syst, wzjets_normalization_syst, qcd_normalization_2017_syst};

    double up_error_cbwi_phi;
    double down_error_cbwi_phi;

    double up_error_ctwi_phi;
    double down_error_ctwi_phi;

    double up_error_cbwi_cos;
    double down_error_cbwi_cos;

    double up_error_ctwi_cos;
    double down_error_ctwi_cos;

    double min_cbwi;
    double min_ctwi;

    string source[5] = {"Lumi_2017_syst", "tt_normalization_syst", "tw_normalization_syst", "wzjets_normalization_syst", "qcd_normalization_2017_syst"};

    for(int i=0 ; i<5 ; i++)
    { 
        if(i==0)
        {
            cout<<"treating i = "<<i<<" ; for "<<source[i]<<"Error= "<<Errors[i]<<endl;
            up_error_cbwi_phi = Stat_and_Syst("cbwi", "phi", nominal_Input_phi, nominal_Input_cos, 1+Errors[i], 1+Errors[i], 1+Errors[i], 1+Errors[i], 1+Errors[i], true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);
            down_error_cbwi_phi = Stat_and_Syst("cbwi", "phi", nominal_Input_phi, nominal_Input_cos, 1-Errors[i], 1-Errors[i], 1-Errors[i], 1-Errors[i], 1-Errors[i], true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);

            up_error_ctwi_phi = Stat_and_Syst("ctwi", "phi", nominal_Input_phi, nominal_Input_cos, 1+Errors[i], 1+Errors[i], 1+Errors[i], 1+Errors[i], 1+Errors[i], true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);
            down_error_ctwi_phi = Stat_and_Syst("ctwi", "phi", nominal_Input_phi, nominal_Input_cos, 1-Errors[i], 1-Errors[i], 1-Errors[i], 1-Errors[i], 1-Errors[i], true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);

            up_error_cbwi_cos = Stat_and_Syst("cbwi", "cos", nominal_Input_phi, nominal_Input_cos, 1+Errors[i], 1+Errors[i], 1+Errors[i], 1+Errors[i], 1+Errors[i], true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);
            down_error_cbwi_cos = Stat_and_Syst("cbwi", "cos", nominal_Input_phi, nominal_Input_cos, 1-Errors[i], 1-Errors[i], 1-Errors[i], 1-Errors[i], 1-Errors[i], true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);

            up_error_ctwi_cos = Stat_and_Syst("ctwi", "cos", nominal_Input_phi, nominal_Input_cos, 1+Errors[i], 1+Errors[i], 1+Errors[i], 1+Errors[i], 1+Errors[i], true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);
            down_error_ctwi_cos = Stat_and_Syst("ctwi", "cos", nominal_Input_phi, nominal_Input_cos, 1-Errors[i], 1-Errors[i], 1-Errors[i], 1-Errors[i], 1-Errors[i], true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);
        }

        if(i==1)
        {
            cout<<"treating i = "<<i<<" ; for "<<source[i]<<"Error= "<<Errors[i]<<endl;            
            up_error_cbwi_phi = Stat_and_Syst("cbwi", "phi", nominal_Input_phi, nominal_Input_cos, 1+Errors[i], 1, 1, 1, 1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);
            down_error_cbwi_phi = Stat_and_Syst("cbwi", "phi", nominal_Input_phi, nominal_Input_cos, 1-Errors[i], 1, 1, 1, 1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);

            up_error_ctwi_phi = Stat_and_Syst("ctwi", "phi", nominal_Input_phi, nominal_Input_cos, 1+Errors[i], 1, 1, 1, 1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);
            down_error_ctwi_phi = Stat_and_Syst("ctwi", "phi", nominal_Input_phi, nominal_Input_cos, 1-Errors[i], 1, 1, 1, 1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);

            up_error_cbwi_cos = Stat_and_Syst("cbwi", "cos", nominal_Input_phi, nominal_Input_cos, 1+Errors[i], 1, 1, 1, 1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);
            down_error_cbwi_cos = Stat_and_Syst("cbwi", "cos", nominal_Input_phi, nominal_Input_cos, 1-Errors[i], 1, 1, 1, 1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);

            up_error_ctwi_cos = Stat_and_Syst("ctwi", "cos", nominal_Input_phi, nominal_Input_cos, 1+Errors[i], 1, 1, 1, 1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);
            down_error_ctwi_cos = Stat_and_Syst("ctwi", "cos", nominal_Input_phi, nominal_Input_cos, 1-Errors[i], 1, 1, 1, 1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);

        }

        if(i==2)
        {
            cout<<"treating i = "<<i<<" ; for "<<source[i]<<"Error= "<<Errors[i]<<endl;
            up_error_cbwi_phi = Stat_and_Syst("cbwi", "phi", nominal_Input_phi, nominal_Input_cos, 1, 1+Errors[i], 1, 1, 1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);
            down_error_cbwi_phi = Stat_and_Syst("cbwi", "phi", nominal_Input_phi, nominal_Input_cos, 1, 1-Errors[i], 1, 1, 1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);

            up_error_ctwi_phi = Stat_and_Syst("ctwi", "phi", nominal_Input_phi, nominal_Input_cos, 1, 1+Errors[i], 1, 1, 1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);
            down_error_ctwi_phi = Stat_and_Syst("ctwi", "phi", nominal_Input_phi, nominal_Input_cos, 1, 1-Errors[i], 1, 1, 1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);

            up_error_cbwi_cos = Stat_and_Syst("cbwi", "cos", nominal_Input_phi, nominal_Input_cos, 1, 1+Errors[i], 1, 1, 1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);
            down_error_cbwi_cos = Stat_and_Syst("cbwi", "cos", nominal_Input_phi, nominal_Input_cos, 1, 1-Errors[i], 1, 1, 1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);

            up_error_ctwi_cos = Stat_and_Syst("ctwi", "cos", nominal_Input_phi, nominal_Input_cos, 1, 1+Errors[i], 1, 1, 1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);
            down_error_ctwi_cos = Stat_and_Syst("ctwi", "cos", nominal_Input_phi, nominal_Input_cos, 1, 1-Errors[i], 1, 1, 1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);
         
        }

        if(i==3)
        {
            cout<<"treating i = "<<i<<" ; for "<<source[i]<<"Error= "<<Errors[i]<<endl;
            up_error_cbwi_phi = Stat_and_Syst("cbwi", "phi", nominal_Input_phi, nominal_Input_cos, 1, 1, 1+Errors[i], 1, 1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);
            down_error_cbwi_phi = Stat_and_Syst("cbwi", "phi", nominal_Input_phi, nominal_Input_cos, 1, 1, 1-Errors[i], 1, 1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);

            up_error_ctwi_phi = Stat_and_Syst("ctwi", "phi", nominal_Input_phi, nominal_Input_cos, 1, 1, 1+Errors[i], 1, 1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);
            down_error_ctwi_phi = Stat_and_Syst("ctwi", "phi", nominal_Input_phi, nominal_Input_cos, 1, 1, 1-Errors[i], 1, 1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);

            up_error_cbwi_cos = Stat_and_Syst("cbwi", "cos", nominal_Input_phi, nominal_Input_cos, 1, 1, 1+Errors[i], 1, 1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);
            down_error_cbwi_cos = Stat_and_Syst("cbwi", "cos", nominal_Input_phi, nominal_Input_cos, 1, 1, 1-Errors[i], 1, 1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);

            up_error_ctwi_cos = Stat_and_Syst("ctwi", "cos", nominal_Input_phi, nominal_Input_cos, 1, 1, 1+Errors[i], 1, 1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);
            down_error_ctwi_cos = Stat_and_Syst("ctwi", "cos", nominal_Input_phi, nominal_Input_cos, 1, 1, 1-Errors[i], 1, 1, true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);
         
        }

        if(i==4)
        {
            cout<<"treating i = "<<i<<" ; for "<<source[i]<<"Error= "<<Errors[i]<<endl;
            up_error_cbwi_phi = Stat_and_Syst("cbwi", "phi", nominal_Input_phi, nominal_Input_cos, 1, 1, 1, 1+Errors[i], 1+Errors[i], true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);
            down_error_cbwi_phi = Stat_and_Syst("cbwi", "phi", nominal_Input_phi, nominal_Input_cos, 1, 1, 1, 1-Errors[i], 1-Errors[i], true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);

            up_error_ctwi_phi = Stat_and_Syst("ctwi", "phi", nominal_Input_phi, nominal_Input_cos, 1, 1, 1, 1+Errors[i], 1+Errors[i], true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);
            down_error_ctwi_phi = Stat_and_Syst("ctwi", "phi", nominal_Input_phi, nominal_Input_cos, 1, 1, 1, 1-Errors[i], 1-Errors[i], true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);

            up_error_cbwi_cos = Stat_and_Syst("cbwi", "cos", nominal_Input_phi, nominal_Input_cos, 1, 1, 1, 1+Errors[i], 1+Errors[i], true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);
            down_error_cbwi_cos = Stat_and_Syst("cbwi", "cos", nominal_Input_phi, nominal_Input_cos, 1, 1, 1, 1-Errors[i], 1-Errors[i], true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);

            up_error_ctwi_cos = Stat_and_Syst("ctwi", "cos", nominal_Input_phi, nominal_Input_cos, 1, 1, 1, 1+Errors[i], 1+Errors[i], true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);
            down_error_ctwi_cos = Stat_and_Syst("ctwi", "cos", nominal_Input_phi, nominal_Input_cos, 1, 1, 1, 1-Errors[i], 1-Errors[i], true, Wilson_phiS_cbwi_Input, Wilson_phiS_ctwi_Input, Wilson_cosS_cbwi_Input, Wilson_cosS_ctwi_Input);
         
        }



        //source[i].replace(source[i].find("_"),1,"\\_");


        max_cbwi_phi = max(abs(up_error_cbwi_phi),abs(down_error_cbwi_phi)); 
        max_ctwi_phi = max(abs(up_error_ctwi_phi),abs(down_error_ctwi_phi));
        max_cbwi_cos = max(abs(up_error_cbwi_cos),abs(down_error_cbwi_cos));
        max_ctwi_cos = max(abs(up_error_ctwi_cos),abs(down_error_ctwi_cos));

        //out_file <<"(Stat + Syst) delta_cbwi_phiStar for "<<source[i]<<" and " << " = " << max_cbwi_phi << "\n";
        //out_file <<"(Stat + Syst) delta_ctwi_phiStar for "<<source[i]<<" and " << " = " << max_ctwi_phi << "\n";
        //out_file <<"(Stat + Syst) delta_cbwi_cosThetaStar for "<<source[i]<<" and " << " = " << max_cbwi_cos << "\n";
        //out_file <<"(Stat + Syst) delta_ctwi_cosThetaStar for "<<source[i]<<" and "<< " = " << max_ctwi_cos << "\n" << "\n";

        max_value_between_cosT_and_phiS_cbwi_ss = min(max_cbwi_phi,max_cbwi_cos);
        //out_file<<"\n (Stat + Syst) Min Value for CbWi in " << source[i] << " = " << max_value_between_cosT_and_phiS_cbwi_ss << "\n";
        max_value_between_cosT_and_phiS_ctwi_ss = min(max_ctwi_phi,max_ctwi_cos);
        //out_file<<"\n (Stat + Syst) Min Value for CtWi in " << source[i] << " = " << max_value_between_cosT_and_phiS_ctwi_ss << "\n";


        Systematic_cbwi_phi = sqrt( abs( pow(max_cbwi_phi,2) - pow (delta_cbwi_phi_stat,2) ) );
        Systematic_ctwi_phi = sqrt( abs( pow(max_ctwi_phi,2) - pow (delta_ctwi_phi_stat,2) ) );
        Systematic_cbwi_cos = sqrt( abs( pow(max_cbwi_cos,2) - pow (delta_cbwi_cos_stat,2) ) );
        Systematic_ctwi_cos = sqrt( abs( pow(max_ctwi_cos,2) - pow (delta_ctwi_cos_stat,2) ) );

        //out_file <<"(Syst) delta_cbwi_phiStar for "<<source[i]<<" and " << " = " << Systematic_cbwi_phi << "\n";
        //out_file <<"(Syst) delta_ctwi_phiStar for "<<source[i]<<" and " << " = " << Systematic_ctwi_phi << "\n";
        //out_file <<"(Syst) delta_cbwi_cosThetaStar for "<<source[i]<<" and " << " = " << Systematic_cbwi_cos << "\n";
        //out_file <<"(Syst) delta_ctwi_cosThetaStar for "<<source[i]<<" and "<< " = " << Systematic_ctwi_cos << "\n" << "\n";

        min_cbwi = min(Systematic_cbwi_phi,Systematic_cbwi_cos);
        //out_file<<"\n (Syst) Min Value for CbWi in " << source[i] << " = " << min_cbwi << "\n";
        min_ctwi = min(Systematic_ctwi_phi,Systematic_ctwi_cos);
        //out_file<<"\n (Syst) Min Value for CtWi in " << source[i] << " = " << min_ctwi << "\n \n \n";

        out_file << source[i] << "   &   " << std::setprecision(chiffres_signif_s) << min_ctwi << "   &   " << std::setprecision(chiffres_signif_s) << min_cbwi << "\\" << "\\" << "\n"; //"   &   " << std::setprecision(chiffres_signif_ss) << max_value_between_cosT_and_phiS_ctwi_ss << "   &   " << std::setprecision(chiffres_signif_ss) << max_value_between_cosT_and_phiS_cbwi_ss << "\\" << "\\" << "\n";


        up_error_cbwi_phi=0;
        down_error_cbwi_phi=0;

        up_error_ctwi_phi=0;
        down_error_ctwi_phi=0;

        up_error_cbwi_cos=0;
        down_error_cbwi_cos=0;

        up_error_ctwi_cos=0;
        down_error_ctwi_cos=0;
    }

    //out_file << "\\hline\n";


    out_file.close();

    return 0;    
}



/*
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



*/

