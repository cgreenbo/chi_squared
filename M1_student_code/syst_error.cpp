#include  "../include/root_simu.hpp"
#include  "../include/matrixSingleTop.hpp"
#include <iostream>
#include <vector>
#include <cmath>
//#include <ofstream>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TF1.h>
#include <TH1F.h>
#include <TLine.h>
#include <sstream>
#include "../include/Xcarre.hpp"

using namespace std;
 

int main ()
{
  Config config;
  read (&config);
  
  //string nameFile3 = "results/uncertainty/syst_stat_error"+version+".txt";
  ofstream file("results/uncertainty/syst_stat_error"+config.b_mu+"_"+config.version+".txt");

  double statsyst1[6];
  double statsyst2[6];
  string process[6] = {"ttbar", "multi jet", "t channel", "tbar channel", "electroweak", "L"};
                    
                //ttbar, multi jet, t channel, tbar channel, electro, L (Sane report)

  //double Systematic_error[6] = {0.009, 0.017, 0.033, 0.055, 0.05, 0.025};


  double Statistic_errror = StatisticAndSystematic_function(syst_ttbar,syst_multijets,syst_t_channel,syst_tbar_channel,syst_ElectroWeak, lumi,true);

    //stat + syst error
  statsyst1[0] = StatisticAndSystematic_function(syst_ttbar+Systematic_error[0], syst_multijets, syst_t_channel, syst_tbar_channel, syst_ElectroWeak, lumi,false);//ttbar
  statsyst1[1] = StatisticAndSystematic_function(syst_ttbar, syst_multijets+Systematic_error[1], syst_t_channel, syst_tbar_channel, syst_ElectroWeak, lumi,false); //multi jet
  statsyst1[2] = StatisticAndSystematic_function(syst_ttbar, syst_multijets, syst_t_channel+Systematic_error[2], syst_tbar_channel+Systematic_error[3], syst_ElectroWeak, lumi,false);//single t
  statsyst1[3] = StatisticAndSystematic_function(syst_ttbar, syst_multijets, syst_t_channel+Systematic_error[2], syst_tbar_channel-Systematic_error[3], syst_ElectroWeak, lumi,false);//single tbar
  statsyst1[4] = StatisticAndSystematic_function(syst_ttbar, syst_multijets, syst_t_channel, syst_tbar_channel, syst_ElectroWeak+Systematic_error[4], lumi,false);//electro
  statsyst1[5] = StatisticAndSystematic_function(syst_ttbar, syst_multijets, syst_t_channel, syst_tbar_channel, syst_ElectroWeak, lumi+Systematic_error[5],false);//L


  statsyst2[0] = StatisticAndSystematic_function(syst_ttbar-Systematic_error[0], syst_multijets, syst_t_channel, syst_tbar_channel,syst_ElectroWeak, lumi,false); //ttbar
  statsyst2[1] = StatisticAndSystematic_function(syst_ttbar, syst_multijets-Systematic_error[1], syst_t_channel, syst_tbar_channel, syst_ElectroWeak, lumi,false); //multi jet
  statsyst2[2] = StatisticAndSystematic_function(syst_ttbar, syst_multijets, syst_t_channel-Systematic_error[2], syst_tbar_channel-Systematic_error[3], syst_ElectroWeak, lumi,false);//single t dd
  statsyst2[3] = StatisticAndSystematic_function(syst_ttbar, syst_multijets, syst_t_channel-Systematic_error[2], syst_tbar_channel+Systematic_error[3], syst_ElectroWeak, lumi,false);//single tbar
  statsyst2[4] = StatisticAndSystematic_function(syst_ttbar, syst_multijets, syst_t_channel, syst_tbar_channel, syst_ElectroWeak-Systematic_error[4], lumi,false);//electro
  statsyst2[5] = StatisticAndSystematic_function(syst_ttbar, syst_multijets, syst_t_channel, syst_tbar_channel, syst_ElectroWeak, lumi-Systematic_error[5],false);//L

  double syst[6];
  double max_statsyst[6];
  double squaredSum=0;


//suivant le range la valeur des max change 
    for(int i=0; i<6; i++)
    {
        max_statsyst[i] = max(abs(statsyst1[i]),abs(statsyst2[i]));       //we choose the max value btw positiv or negativ

      if(i==3)    //tbar
      {
        process[i] = process[i] + " " + process[i-1];         //tbar and t have the same value so we choose the max btw them
        max_statsyst[i] = max(abs(max_statsyst[i]),abs(max_statsyst[i-1]));
      }

      syst[i] = sqrt( abs( pow(max_statsyst[i],2) - pow (Statistic_errror,2) ) );   
      squaredSum += pow(syst[i],2);      //systematic error

      if(i==2)      //t
      {
      } 
      else
      {
        cout << process[i] << " : " << syst[i] << endl;     //affichage
        file << process[i] << " : " << syst[i] << endl;
      }

}
cout<<"statical error of b_mu when delta Chi^2 at 1 GeV = "<<Statistic_errror/config.time_modulation<<endl;
file << "statical error of b_mu when delta Chi^2 at 1 GeV = "<< Statistic_errror/config.time_modulation <<endl;

cout << "total systematic error of b mu : " << sqrt(squaredSum)/config.time_modulation<<endl; 
file << "total systematic error of b mu : " << sqrt(squaredSum)/config.time_modulation<<endl;

file.close();

return 0;

}
