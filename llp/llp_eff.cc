//llp_eff.cc is used to generate ~g -> g ~G events for long-lived particle searches.
//1000 events are generated in each of the bins from m_~g = 1000 GeV to 3000GeV and ctau_~g = 10^2.5 mm/c to 10^6.0 mm/c.
//Only the hard_scattered gluino particle informations and missing transverse energies are saved.

#include "Pythia8/Pythia.h"

#include "TTree.h"
#include "TFile.h"

#include <cstdlib>
#include <math.h>

using namespace Pythia8;

int main()
{

TFile *file = TFile::Open("tree_llp_met.root","recreate");
TTree *tree = new TTree("tree","tree");


  int evt;
  float met;
  std::vector<int> id;
  std::vector<int> status;
  std::vector<int> mom1;
  std::vector<int> mom2;
  std::vector<int> dtr1;
  std::vector<int> dtr2;
  std::vector<float> px;
  std::vector<float> py;
  std::vector<float> pz;
  std::vector<float> e;
  std::vector<float> m;
  std::vector<float> tau; // mm/c
  std::vector<float> tau0; // mm/c
  std::vector<float> vtx_x; // mm
  std::vector<float> vtx_y;
  std::vector<float> vtx_z;
  std::vector<float> vtx_t;
  std::vector<float> dvtx_x;
  std::vector<float> dvtx_y;
  std::vector<float> dvtx_z;
  std::vector<float> dvtx_t;


  tree->Branch("evt",     &evt);
  tree->Branch("id",      &id);
  tree->Branch("status",  &status);
  tree->Branch("mom1",    &mom1);
  tree->Branch("mom2",    &mom2);
  tree->Branch("dtr1",    &dtr1);
  tree->Branch("dtr2",    &dtr2);
  tree->Branch("px",      &px);
  tree->Branch("py",      &py);
  tree->Branch("pz",      &pz);
  tree->Branch("e",       &e);
  tree->Branch("m",       &m);
  tree->Branch("tau",     &tau);
  tree->Branch("tau0",     &tau0);
  tree->Branch("vtx_x",   &vtx_x);
  tree->Branch("vtx_y",   &vtx_y);
  tree->Branch("vtx_z",   &vtx_z);
  tree->Branch("vtx_t",   &vtx_t);
  tree->Branch("dvtx_x",   &dvtx_x);
  tree->Branch("dvtx_y",   &dvtx_y);
  tree->Branch("dvtx_z",   &dvtx_z);
  tree->Branch("dvtx_t",   &dvtx_t);
  tree->Branch("met",     &met);

for (float m_0 = 1000; m_0 < 3001; m_0 += 100){

  for (float tau_0 = 2.5; tau_0 <6.01; tau_0 += 0.1){


  int NUM_EVENTS=1000;

  Pythia pythia;
  pythia.readString("Beams:eCM = 13000");
  pythia.readString("Beams:frameType = 1");
  pythia.readString("Tune:pp  = 5");
  pythia.readString("SUSY:gg2gluinogluino  = on");
  pythia.readString("SUSY:qqbar2gluinogluino  = on");
  pythia.readString("RHadrons:allow  = on");
  pythia.readString("RHadrons:allowDecay = on");
  pythia.readString("RHadrons:setMasses = on");
  pythia.readString("RHadrons:probGluinoball = 0.1");
  pythia.readString("PhaseSpace:pTHatMin = 0.5");
  pythia.readString(Form("SLHA:file = slha/tau%.1f.spc", tau_0));   //SLHA files should be kept as it is.
  pythia.readString("SLHA:allowUserOverride = on");
  pythia.readString(Form("1000021:m0 = %f", m_0));    //Only mass can be overwritten
  pythia.readString("Next:showScaleAndVertex = on");
  pythia.init();




  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < NUM_EVENTS; ++iEvent)
  {

    float metx = 0.0;
    float mety = 0.0;


    id.clear();
    status.clear();
    mom1.clear();
    mom2.clear();
    dtr1.clear();
    dtr2.clear();
    px.clear();
    py.clear();
    pz.clear();
    e.clear();
    m.clear();
    tau.clear();
    tau0.clear();
    vtx_x.clear();
    vtx_y.clear();
    vtx_z.clear();
    vtx_t.clear();
    dvtx_x.clear();
    dvtx_y.clear();
    dvtx_z.clear();
    dvtx_t.clear();


    if (!pythia.next()) continue;

    evt = iEvent;

    //Missing E_T calculation
    for(int iptl=0; iptl<pythia.event.size(); iptl++)
    {    

    float vtx_x_f = pythia.event.at(iptl).xProd();
    float vtx_y_f = pythia.event.at(iptl).yProd();
    int status_i = pythia.event.at(iptl).status();
    int id_i = pythia.event.at(iptl).id();
    
    if (abs(status_i) == 23 && abs(id_i) == 21){

     float dist;
     
     dist = sqrt( vtx_x_f * vtx_x_f + vtx_y_f * vtx_y_f);
      
     if ( dist < 2850 )
        {
         metx = metx + pythia.event.at(iptl).px(); 
         mety = mety + pythia.event.at(iptl).py();
        }
    
      }




      bool store = false;
      int abs_status = abs(pythia.event.at(iptl).status());
      int abs_id = abs(pythia.event.at(iptl).id());
      if(abs_status>20 && abs_status<30 && abs_id == 1000021)   store = true;
      if(!store) continue;
      if(0)
      {
        cout <<  pythia.event.at(iptl).id() << "\t";
        cout <<  pythia.event.at(iptl).status() << "\t";
        cout <<  pythia.event.at(iptl).mother1() << "\t";
        cout <<  pythia.event.at(iptl).mother2() << "\t";
        cout <<  pythia.event.at(iptl).daughter1() << "\t";
        cout <<  pythia.event.at(iptl).daughter2() << "\t";
        cout <<  pythia.event.at(iptl).px() << "\t";
        cout <<  pythia.event.at(iptl).py() << "\t";
        cout <<  pythia.event.at(iptl).pz() << "\t";
        cout <<  pythia.event.at(iptl).e() << "\t";
        cout <<  pythia.event.at(iptl).m() << "\t";
        cout <<  pythia.event.at(iptl).tau() << "\t";
        cout <<  pythia.event.at(iptl).tau0() << "\t";
        cout <<  pythia.event.at(iptl).xProd() << "\t";
        cout <<  pythia.event.at(iptl).yProd() << "\t";
        cout <<  pythia.event.at(iptl).zProd() << "\t";
        cout <<  pythia.event.at(iptl).tProd() << "\t";
        cout <<  pythia.event.at(iptl).xDec() << "\t";
        cout <<  pythia.event.at(iptl).yDec() << "\t";
        cout <<  pythia.event.at(iptl).zDec() << "\t";
        cout <<  pythia.event.at(iptl).tDec() << "\t";
        cout << endl;
      }

      id.push_back(pythia.event.at(iptl).id());
      status.push_back(pythia.event.at(iptl).status());
      mom1.push_back(pythia.event.at(iptl).mother1());
      mom2.push_back(pythia.event.at(iptl).mother2());
      dtr1.push_back(pythia.event.at(iptl).daughter1());
      dtr2.push_back(pythia.event.at(iptl).daughter2());
      px.push_back(pythia.event.at(iptl).px());
      py.push_back(pythia.event.at(iptl).py());
      pz.push_back(pythia.event.at(iptl).pz());
      e.push_back(pythia.event.at(iptl).e());
      m.push_back(pythia.event.at(iptl).m());
      tau.push_back(pythia.event.at(iptl).tau());
      tau0.push_back(pythia.event.at(iptl).tau0());
      vtx_x.push_back(pythia.event.at(iptl).xProd());
      vtx_y.push_back(pythia.event.at(iptl).yProd());
      vtx_z.push_back(pythia.event.at(iptl).zProd());
      vtx_t.push_back(pythia.event.at(iptl).tProd());
      dvtx_x.push_back(pythia.event.at(iptl).xDec());
      dvtx_y.push_back(pythia.event.at(iptl).yDec());
      dvtx_z.push_back(pythia.event.at(iptl).zDec());
      dvtx_t.push_back(pythia.event.at(iptl).tDec());





    }
    
    met = sqrt(metx*metx + mety*mety);

    tree->Fill(); // Fill all events
}}}
  //pythia.stat();

  file->cd();
  tree->Write();
  file->Close();

  return 0;
}
