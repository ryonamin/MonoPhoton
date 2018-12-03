#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include <iostream>

using namespace std;

#include "eventselection.C"

void Event::process(string fname)
{
  TFile fin(fname.c_str());
  TTree* data = static_cast<TTree*>(fin.Get("evtdata"));

  int nevents = data->GetEntries();

  int nPhoton_MC = 0;
  int nPhoton_Rec = 0;

  // event loop
  //nevents = 100; // for testing
  //Event event(data);
  SetTree(data);
  cerr << "# of events = " << nevents << endl;
  for (int ev = 0; ev < nevents; ev++) {
    callGetEntry(ev);
    
    //cerr << "Accepted event. ev = " << ev << endl;
    //data->GetEntry(ev); 

    if (!isAcceptableEvent(ev)) continue; 
    // check if corresponding pfo exists.
    // pfo loop
    for (int i = 0; i < npfos; i++) {
      // e+/e- pt distribution 
      if (abs(pfo_pdg[i])==11) { // select electron/positron 
        TLorentzVector p4(pfo_px[i],pfo_py[i],pfo_pz[i],pfo_e[i]);
        if (mcr_isoverlay[i])                outputs.hPt_ep_ol->Fill(p4.Pt());
        //else if (mcr_isOriginatedFromISR[i]) outputs.hPt_ep_isr->Fill(p4.Pt());
        else if ((mcr_pdg[mcr_parentIndex[i][0]]==22&&(mcr_parentIndex[i][0]>=8&&mcr_parentIndex[i][0]<=10))||
                 (mcr_pdg[mcr_parentIndex[i][1]]==22&&(mcr_parentIndex[i][1]>=8&&mcr_parentIndex[i][1]<=10)) ) outputs.hPt_ep_isr->Fill(p4.Pt());
        else                                       outputs.hPt_ep_other->Fill(p4.Pt());
      }

      // Charged PFO pt distribution 
      if (fabs(pfo_chrg[i])>0) { // select charged pfo 
        TLorentzVector p4(pfo_px[i],pfo_py[i],pfo_pz[i],pfo_e[i]);
        if (mcr_isoverlay[i])                outputs.hPt_pfo_ol->Fill(p4.Pt());
        //else if (mcr_isOriginatedFromISR[i]) outputs.hPt_pfo_isr->Fill(p4.Pt());
        else if ((mcr_pdg[mcr_parentIndex[i][0]]==22&&(mcr_parentIndex[i][0]>=8&&mcr_parentIndex[i][0]<=10))||
                 (mcr_pdg[mcr_parentIndex[i][1]]==22&&(mcr_parentIndex[i][1]>=8&&mcr_parentIndex[i][1]<=10)) ) outputs.hPt_pfo_isr->Fill(p4.Pt());
        else                                       outputs.hPt_pfo_other->Fill(p4.Pt());
      }

    }

//    if (!isPassedPtCut(ev)) continue; 

    // pfo loop
    float esum = 0.;
    float esum_wo_pi_n = 0.;
    for (int i = 0; i < npfos; i++) {

      // photon E distribution 
      if (abs(pfo_pdg[i])==22&&!mcr_isOriginatedFromISR[i]) { // select photon 
        if (mcr_isoverlay[i])                outputs.hE_photon_ol->Fill(pfo_e[i]);
        else                                       outputs.hE_photon_rest->Fill(pfo_e[i]);
      }

      // V0 E distribution 
      if (abs(pfo_pdg[i])==310||abs(pfo_pdg[i])==3122) { // select V0 
        if (mcr_isoverlay[i])                outputs.hE_V0_ol->Fill(pfo_e[i]);
        else                                       outputs.hE_V0_rest->Fill(pfo_e[i]);
      }

      // neutron E distribution 
      if (abs(pfo_pdg[i])==2112) { // select neutron 
        if (mcr_isoverlay[i])                outputs.hE_neutron_ol->Fill(pfo_e[i]);
        else                                       outputs.hE_neutron_rest->Fill(pfo_e[i]);
      }

      // electron E distribution 
      if (abs(pfo_pdg[i])==11) { // select electron 
        if (mcr_isoverlay[i])                outputs.hE_electron_ol->Fill(pfo_e[i]);
        else                                       outputs.hE_electron_rest->Fill(pfo_e[i]);
      }

      // muon E distribution 
      if (abs(pfo_pdg[i])==13) { // select muon 
        if (mcr_isoverlay[i])                outputs.hE_muon_ol->Fill(pfo_e[i]);
        else                                       outputs.hE_muon_rest->Fill(pfo_e[i]);
      }

      // pion E distribution 
      if (abs(pfo_pdg[i])==211) { // select muon 
        if (mcr_isoverlay[i])                outputs.hE_pion_ol->Fill(pfo_e[i]);
        else                                       outputs.hE_pion_rest->Fill(pfo_e[i]);
      }

      if (pfo_e[i]>5.) { // individual energy cut
        esum += pfo_e[i];
        if (!(abs(pfo_pdg[i])==211||abs(pfo_pdg[i])==111||abs(pfo_pdg[i])==2112)) {
          esum_wo_pi_n += pfo_e[i];
        }
      }

    }
    outputs.hEvis_pfo->Fill(esum);
    outputs.hEvis_pfo_wo_pi_n->Fill(esum_wo_pi_n);

//    if (!isPassedECut(ev)) continue; 

    //cerr << "passed event. ev = " << ev << endl;
    outputs.hNBcalClusters->Fill(nbcalhits);
    if      (getNISRPhotons()==1) outputs.hNBcalClusters1ISR->Fill(nbcalhits);
    else if (getNISRPhotons()>1)  outputs.hNBcalClustersMultiISR->Fill(nbcalhits);
  }
}

void mcp_pfo_test()
{
  Event nung;
  nung.outputs.hPt_ep_isr             = new TH1F("hPt_ep_nung_isr",";Pt [GeV/c];",100,0,250);
  nung.outputs.hPt_ep_ol              = new TH1F("hPt_ep_nung_ol",";Pt [GeV/c];",100,0,250);
  nung.outputs.hPt_ep_other           = new TH1F("hPt_ep_nung_other",";Pt [GeV/c];",100,0,250);
  nung.outputs.hPt_pfo_isr            = new TH1F("hPt_pfo_nung_isr",";Pt [GeV/c];",100,0,250);
  nung.outputs.hPt_pfo_ol             = new TH1F("hPt_pfo_nung_ol",";Pt [GeV/c];",100,0,250);
  nung.outputs.hPt_pfo_other          = new TH1F("hPt_pfo_nung_other",";Pt [GeV/c];",100,0,250);
  nung.outputs.hE_photon_rest         = new TH1F("hE_photon_nung_rest",";E [GeV];",100,0,350);
  nung.outputs.hE_photon_ol           = new TH1F("hE_photon_nung_ol",";E [GeV];",100,0,350);
  nung.outputs.hE_photon_electron     = new TH1F("hE_photon_nung_electron",";E [GeV];",100,0,350);
  nung.outputs.hE_V0_rest             = new TH1F("hE_V0_nung_rest",";E [GeV];",100,0,350);
  nung.outputs.hE_V0_ol               = new TH1F("hE_V0_nung_ol",";E [GeV];",100,0,350);
  nung.outputs.hE_neutron_rest        = new TH1F("hE_neutron_nung_rest",";E [GeV];",100,0,350);
  nung.outputs.hE_neutron_ol          = new TH1F("hE_neutron_nung_ol",";E [GeV];",100,0,350);
  nung.outputs.hE_electron_rest       = new TH1F("hE_electron_nung_rest",";E [GeV];",100,0,350);
  nung.outputs.hE_electron_ol         = new TH1F("hE_electron_nung_ol",";E [GeV];",100,0,350);
  nung.outputs.hE_muon_rest           = new TH1F("hE_muon_nung_rest",";E [GeV];",100,0,350);
  nung.outputs.hE_muon_ol             = new TH1F("hE_muon_nung_ol",";E [GeV];",100,0,350);
  nung.outputs.hE_pion_rest           = new TH1F("hE_pion_nung_rest",";E [GeV];",100,0,350);
  nung.outputs.hE_pion_ol             = new TH1F("hE_pion_nung_ol",";E [GeV];",100,0,350);
  nung.outputs.hEvis_pfo              = new TH1F("hEvis_pfo_nung",";E [GeV];",100,0,750);
  nung.outputs.hEvis_pfo_wo_pi_n      = new TH1F("hEvis_pfo_wo_pi_n_nung",";E [GeV];",100,0,750);
  nung.outputs.hNBcalClusters         = new TH1F("hNBcalClusters_nung",";# of BCal clusters;",6,0,6);
  nung.outputs.hNBcalClusters1ISR     = new TH1F("hNBcalClusters1ISR_nung",";# of BCal clusters;",6,0,6);
  nung.outputs.hNBcalClustersMultiISR = new TH1F("hNBcalClustersMultiISR_nung",";# of BCal clusters;",6,0,6);
  nung.outputs.hPt_ep_isr->SetFillColor(kAzure+6);
  nung.outputs.hPt_ep_ol->SetFillColor(kYellow+1);
  nung.outputs.hPt_ep_other->SetFillColor(kGreen+3);
  nung.outputs.hPt_pfo_isr->SetFillColor(kAzure+6);
  nung.outputs.hPt_pfo_ol->SetFillColor(kYellow+1);
  nung.outputs.hPt_pfo_other->SetFillColor(kGreen+3);
  nung.outputs.hE_photon_rest->SetFillColor(kGreen+3);
  nung.outputs.hE_photon_ol->SetFillColor(kYellow+1);
  nung.outputs.hE_V0_rest->SetFillColor(kGreen+3);
  nung.outputs.hE_V0_ol->SetFillColor(kYellow+1);
  nung.outputs.hE_neutron_rest->SetFillColor(kGreen+3);
  nung.outputs.hE_neutron_ol->SetFillColor(kYellow+1);
  nung.outputs.hE_electron_rest->SetFillColor(kGreen+3);
  nung.outputs.hE_electron_ol->SetFillColor(kYellow+1);
  nung.outputs.hE_muon_rest->SetFillColor(kGreen+3);
  nung.outputs.hE_muon_ol->SetFillColor(kYellow+1);
  nung.outputs.hE_pion_rest->SetFillColor(kGreen+3);
  nung.outputs.hE_pion_ol->SetFillColor(kYellow+1);
  nung.outputs.hEvis_pfo->SetFillColor(kGreen+3);
  nung.outputs.hEvis_pfo_wo_pi_n->SetFillColor(kGreen+3);
  nung.outputs.hNBcalClusters1ISR->SetFillColor(kYellow+1);
  nung.outputs.hNBcalClustersMultiISR->SetFillColor(kGreen+3);

  //for (int i = 1; i < 100; i++) {
  for (int i = 1; i < 2; i++) {
    stringstream fname;
    fname << "/hsm/ilc/users/yonamine/physics/monophoton/run_l5/nung/root/l5_500GeV.nung_" << i << ".root" << ends;
    nung.process(fname.str().data());
  }
  //nung.process("nunug/root/l5_500GeV.nung_1.root");

  Event bhabhang;
  bhabhang.outputs.hPt_ep_isr   = new TH1F("hPt_ep_bhabhang_isr",";Pt [GeV/c];",100,0,250);
  bhabhang.outputs.hPt_ep_ol    = new TH1F("hPt_ep_bhabhang_ol",";Pt [GeV/c];",100,0,250);
  bhabhang.outputs.hPt_ep_other = new TH1F("hPt_ep_bhabhang_other",";Pt [GeV/c];",100,0,250);
  bhabhang.outputs.hPt_pfo_isr   = new TH1F("hPt_pfo_bhabhang_isr",";Pt [GeV/c];",100,0,250);
  bhabhang.outputs.hPt_pfo_ol    = new TH1F("hPt_pfo_bhabhang_ol",";Pt [GeV/c];",100,0,250);
  bhabhang.outputs.hPt_pfo_other = new TH1F("hPt_pfo_bhabhang_other",";Pt [GeV/c];",100,0,250);
  bhabhang.outputs.hE_photon_rest     = new TH1F("hE_photon_bhabhang_rest",";E [GeV];",100,0,350);
  bhabhang.outputs.hE_photon_ol       = new TH1F("hE_photon_bhabhang_ol",";E [GeV];",100,0,350);
  bhabhang.outputs.hE_photon_electron = new TH1F("hE_photon_bhabhang_electron",";E [GeV];",100,0,350);
  bhabhang.outputs.hE_V0_rest         = new TH1F("hE_V0_bhabhang_rest",";E [GeV];",100,0,350);
  bhabhang.outputs.hE_V0_ol           = new TH1F("hE_V0_bhabhang_ol",";E [GeV];",100,0,350);
  bhabhang.outputs.hE_neutron_rest    = new TH1F("hE_neutron_bhabhang_rest",";E [GeV];",100,0,350);
  bhabhang.outputs.hE_neutron_ol      = new TH1F("hE_neutron_bhabhang_ol",";E [GeV];",100,0,350);
  bhabhang.outputs.hE_electron_rest   = new TH1F("hE_electron_bhabhang_rest",";E [GeV];",100,0,350);
  bhabhang.outputs.hE_electron_ol     = new TH1F("hE_electron_bhabhang_ol",";E [GeV];",100,0,350);
  bhabhang.outputs.hE_muon_rest       = new TH1F("hE_muon_bhabhang_rest",";E [GeV];",100,0,350);
  bhabhang.outputs.hE_muon_ol         = new TH1F("hE_muon_bhabhang_ol",";E [GeV];",100,0,350);
  bhabhang.outputs.hE_pion_rest       = new TH1F("hE_pion_bhabhang_rest",";E [GeV];",100,0,350);
  bhabhang.outputs.hE_pion_ol         = new TH1F("hE_pion_bhabhang_ol",";E [GeV];",100,0,350);
  bhabhang.outputs.hEvis_pfo          = new TH1F("hEvis_pfo_bhabhang",";E [GeV];",100,0,750);
  bhabhang.outputs.hEvis_pfo_wo_pi_n  = new TH1F("hEvis_pfo_wo_pi_n_bhabhang",";E [GeV];",100,0,750);
  bhabhang.outputs.hNBcalClusters = new TH1F("hNBcalClusters_bhabhang",";# of BCal clusters;",6,0,6);
  bhabhang.outputs.hNBcalClusters1ISR = new TH1F("hNBcalClusters1ISR_bhabhang",";# of BCal clusters;",6,0,6);
  bhabhang.outputs.hNBcalClustersMultiISR = new TH1F("hNBcalClustersMultiISR_bhabhang",";# of BCal clusters;",6,0,6);
  bhabhang.outputs.hPt_ep_isr->SetFillColor(kOrange+4);
  bhabhang.outputs.hPt_ep_ol->SetFillColor(kMagenta+1);
  bhabhang.outputs.hPt_ep_other->SetFillColor(kOrange+6);
  bhabhang.outputs.hPt_pfo_isr->SetFillColor(kOrange+4);
  bhabhang.outputs.hPt_pfo_ol->SetFillColor(kMagenta+1);
  bhabhang.outputs.hPt_pfo_other->SetFillColor(kOrange+6);
  bhabhang.outputs.hE_photon_rest->SetFillColor(kOrange+6);
  bhabhang.outputs.hE_photon_electron->SetFillColor(kOrange+4);
  bhabhang.outputs.hE_photon_ol->SetFillColor(kMagenta+1);
  bhabhang.outputs.hE_V0_rest->SetFillColor(kOrange+6);
  bhabhang.outputs.hE_V0_ol->SetFillColor(kMagenta+1);
  bhabhang.outputs.hE_neutron_rest->SetFillColor(kOrange+6);
  bhabhang.outputs.hE_neutron_ol->SetFillColor(kMagenta+1);
  bhabhang.outputs.hE_electron_rest->SetFillColor(kOrange+6);
  bhabhang.outputs.hE_electron_ol->SetFillColor(kMagenta+1);
  bhabhang.outputs.hE_muon_rest->SetFillColor(kOrange+6);
  bhabhang.outputs.hE_muon_ol->SetFillColor(kMagenta+1);
  bhabhang.outputs.hE_pion_rest->SetFillColor(kOrange+6);
  bhabhang.outputs.hE_pion_ol->SetFillColor(kMagenta+1);
  bhabhang.outputs.hEvis_pfo->SetFillColor(kOrange+6);
  bhabhang.outputs.hEvis_pfo_wo_pi_n->SetFillColor(kOrange+6);
  bhabhang.outputs.hNBcalClusters->SetFillColor(kOrange+6);

  //for (int i = 1; i < 131; i++) {
  for (int i = 1; i < 3; i++) {
    stringstream fname;
    fname << "/hsm/ilc/users/yonamine/physics/monophoton/run_l5/bhabhang/root/l5_500GeV.bhabhang_" << i << ".root" << ends;
    bhabhang.process(fname.str().data());
  }
  //bhabhang.process("bhabhang/root/l5_500GeV.bhabhang_52.root"); 

  TCanvas* c1 = new TCanvas("c1","",600,400);
  gPad->SetLogy();
  THStack* hPt_e_all = new THStack("hPt_e_all","Electron/Positron;P_{T} [GeV];");
  hPt_e_all->Add(nung.outputs.hPt_ep_other);
  hPt_e_all->Add(nung.outputs.hPt_ep_ol);
  hPt_e_all->Add(nung.outputs.hPt_ep_isr);
  hPt_e_all->Add(bhabhang.outputs.hPt_ep_ol);
  hPt_e_all->Add(bhabhang.outputs.hPt_ep_isr);
  hPt_e_all->Add(bhabhang.outputs.hPt_ep_other);
  hPt_e_all->Draw();

  TCanvas* c2 = new TCanvas("c2","",600,400);
  gPad->SetLogy();
  THStack* hPt_pfo_all = new THStack("hPt_pfo_all","PFO;P_{T} [GeV];");
  hPt_pfo_all->Add(nung.outputs.hPt_pfo_other);
  hPt_pfo_all->Add(nung.outputs.hPt_pfo_ol);
  hPt_pfo_all->Add(nung.outputs.hPt_pfo_isr);
  hPt_pfo_all->Add(bhabhang.outputs.hPt_pfo_ol);
  hPt_pfo_all->Add(bhabhang.outputs.hPt_pfo_isr);
  hPt_pfo_all->Add(bhabhang.outputs.hPt_pfo_other);
  hPt_pfo_all->Draw();

  TCanvas* c3 = new TCanvas("c3","",600,400);
  gPad->SetLogy();
  THStack* hE_photon_all = new THStack("hE_photon_all","Photon;E [GeV];");
  hE_photon_all->Add(nung.outputs.hE_photon_rest);
  hE_photon_all->Add(nung.outputs.hE_photon_ol);
  hE_photon_all->Add(bhabhang.outputs.hE_photon_rest);
  hE_photon_all->Add(bhabhang.outputs.hE_photon_electron);
  hE_photon_all->Add(bhabhang.outputs.hE_photon_ol);
  hE_photon_all->Draw();

  TCanvas* c4 = new TCanvas("c4","",600,400);
  gPad->SetLogy();
  THStack* hE_V0_all = new THStack("hE_V0_all","V0;E [GeV];");
  hE_V0_all->Add(nung.outputs.hE_V0_rest);
  hE_V0_all->Add(nung.outputs.hE_V0_ol);
  hE_V0_all->Add(bhabhang.outputs.hE_V0_rest);
  hE_V0_all->Add(bhabhang.outputs.hE_V0_ol);
  hE_V0_all->Draw();

  TCanvas* c5 = new TCanvas("c5","",600,400);
  gPad->SetLogy();
  THStack* hE_neutron_all = new THStack("hE_neutron_all","Neutron;E [GeV];");
  hE_neutron_all->Add(nung.outputs.hE_neutron_rest);
  hE_neutron_all->Add(nung.outputs.hE_neutron_ol);
  hE_neutron_all->Add(bhabhang.outputs.hE_neutron_rest);
  hE_neutron_all->Add(bhabhang.outputs.hE_neutron_ol);
  hE_neutron_all->Draw();

  TCanvas* c6 = new TCanvas("c6","",600,400);
  gPad->SetLogy();
  THStack* hE_electron_all = new THStack("hE_electron_all","Electron;E [GeV];");
  hE_electron_all->Add(nung.outputs.hE_electron_rest);
  hE_electron_all->Add(nung.outputs.hE_electron_ol);
  hE_electron_all->Add(bhabhang.outputs.hE_electron_rest);
  hE_electron_all->Add(bhabhang.outputs.hE_electron_ol);
  hE_electron_all->Draw();

  TCanvas* c7 = new TCanvas("c7","",600,400);
  gPad->SetLogy();
  THStack* hE_muon_all = new THStack("hE_muon_all","Muon;E [GeV];");
  hE_muon_all->Add(nung.outputs.hE_muon_rest);
  hE_muon_all->Add(nung.outputs.hE_muon_ol);
  hE_muon_all->Add(bhabhang.outputs.hE_muon_rest);
  hE_muon_all->Add(bhabhang.outputs.hE_muon_ol);
  hE_muon_all->Draw();

  TCanvas* c8 = new TCanvas("c8","",600,400);
  gPad->SetLogy();
  THStack* hE_pion_all = new THStack("hE_pion_all","Pion;E [GeV];");
  hE_pion_all->Add(nung.outputs.hE_pion_rest);
  hE_pion_all->Add(nung.outputs.hE_pion_ol);
  hE_pion_all->Add(bhabhang.outputs.hE_pion_rest);
  hE_pion_all->Add(bhabhang.outputs.hE_pion_ol);
  hE_pion_all->Draw();

  TCanvas* c9 = new TCanvas("c9","",600,400);
  gPad->SetLogy();
  THStack* hEvis_all = new THStack("hEvis_all","PFO;Evis [GeV];");
  hEvis_all->Add(nung.outputs.hEvis_pfo);
  hEvis_all->Add(bhabhang.outputs.hEvis_pfo);
  hEvis_all->Draw();

  TCanvas* c10 = new TCanvas("c10","",600,400);
  gPad->SetLogy();
  THStack* hEvis_wo_pi_n_all = new THStack("hEvis_wo_pi_n_all","PFO w/o Pion, Nutron;Evis [GeV];");
  hEvis_wo_pi_n_all->Add(nung.outputs.hEvis_pfo_wo_pi_n);
  hEvis_wo_pi_n_all->Add(bhabhang.outputs.hEvis_pfo_wo_pi_n);
  hEvis_wo_pi_n_all->Draw();

  TCanvas* c11 = new TCanvas("c11","",600,400);
  gPad->SetLogy();
  THStack* hNBcalClusters_all = new THStack("hNBcalClusters_all",";# of BCal clusters;");
  hNBcalClusters_all->Add(nung.outputs.hNBcalClusters1ISR);
  hNBcalClusters_all->Add(nung.outputs.hNBcalClustersMultiISR);
  hNBcalClusters_all->Add(bhabhang.outputs.hNBcalClusters);
  hNBcalClusters_all->Draw(); 
}

