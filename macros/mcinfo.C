
int nmcps;
const int NMAXMCPS = 500;
float mcp_e[NMAXMCPS], mcp_phi[NMAXMCPS], mcp_theta[NMAXMCPS];
bool mcp_iscreatedinsim[NMAXMCPS];
int mcp_pdg[NMAXMCPS];

void mcinfo ()
{
  TFile* fin = new TFile("output.root");
  TTree* T = static_cast<TTree*>(fin->Get("evtdata"));

  TH1D* he_gen = new TH1D("he_gen","",100,0,15);
  TH1D* he_sim = new TH1D("he_sim","",100,0,15);
  TH1D* hphi_gen = new TH1D("hphi_gen","",100,-3.2,3.2);
  TH1D* hphi_sim = new TH1D("hphi_sim","",100,-3.2,3.2);
  TH1D* htheta_gen = new TH1D("htheta_gen","",100,-0.3,0.3);
  TH1D* htheta_sim = new TH1D("htheta_sim","",100,-0.3,0.3);
  TH1D* hpdg_gen = new TH1D("hpdg_gen","",51,-15,35);
  TH1D* hpdg_sim = new TH1D("hpdg_sim","",51,-15,35);
  TCanvas* ce = new TCanvas("ce","",600,400);
  TCanvas* cphi = new TCanvas("cphi","",600,400);
  TCanvas* ctheta = new TCanvas("ctheta","",600,400);
  TCanvas* cpdg = new TCanvas("cpdg","",600,400);

  T->SetBranchAddress("nmcps",&nmcps);
  T->SetBranchAddress("mcp_e",mcp_e);
  T->SetBranchAddress("mcp_pdg",mcp_pdg);
  T->SetBranchAddress("mcp_phi",mcp_phi);
  T->SetBranchAddress("mcp_theta",mcp_theta);
  T->SetBranchAddress("mcp_iscreatedinsim",mcp_iscreatedinsim);

  int nevents = T->GetEntries();

  for (int ev = 0; ev < nevents; ev++) {
    
    T->GetEntry(ev);
    for ( int ip = 0; ip < nmcps; ip++ ) {
      if (!mcp_iscreatedinsim[ip]) {
        he_gen->Fill(mcp_e[ip]);
        hphi_gen->Fill(mcp_phi[ip]);
        htheta_gen->Fill(mcp_theta[ip]);
        hpdg_gen->Fill(mcp_pdg[ip]);
      } else {
        he_sim->Fill(mcp_e[ip]);
        hphi_sim->Fill(mcp_phi[ip]);
        htheta_sim->Fill(mcp_theta[ip]);
        hpdg_sim->Fill(mcp_pdg[ip]);
      }
    }
  }
  ce->cd();
  he_gen->SetLineColor(4);
  he_gen->Draw();
  he_sim->SetLineColor(2);
  he_sim->Draw("same");
  cphi->cd();
  hphi_gen->SetLineColor(4);
  hphi_gen->SetMaximum(1100);
  hphi_gen->Draw();
  hphi_sim->SetLineColor(2);
  hphi_sim->Draw("same");
  ctheta->cd();
  htheta_gen->SetLineColor(4);
  htheta_gen->SetMaximum(8000);
  htheta_gen->Draw();
  htheta_sim->SetLineColor(2);
  htheta_sim->Draw("same");
  cpdg->cd();
  hpdg_gen->SetLineColor(4);
  hpdg_gen->Draw();
  hpdg_sim->SetLineColor(2);
  hpdg_sim->Draw("same");
  stringstream fnamee,fnamephi,fnametheta, fnamepdg;
  fnamee << "mcp_e.pdf" << ends;
  fnamephi << "mcp_phi.pdf" << ends;
  fnametheta << "mcp_theta.pdf" << ends;
  fnamepdg << "mcp_pdg.pdf" << ends;
  ce->Print(fnamee.str().data());
  cphi->Print(fnamephi.str().data());
  ctheta->Print(fnametheta.str().data());
  cpdg->Print(fnamepdg.str().data());

}
