#include  "../include/root_simu.hpp"
#include  "../include/matrixSingleTop.hpp"
#include <iostream>
#include <vector>
#include <cmath>
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
  
  double StatError = StatisticAndSystematic_function(syst_ttbar, syst_multijets, syst_t_channel, syst_tbar_channel, syst_ElectroWeak, lumi,true);
  cout<<"statistical error of b_mu when delta Chi^2 at 1 GeV = "<<StatError/config.time_modulation<<endl;

  return 0;
}
