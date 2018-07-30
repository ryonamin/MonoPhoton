
int npfos;
const int NMAXPFOS = 500;
float pfo_e[NMAXPFOS], pfo_theta[NMAXPFOS];
float mcr_e[NMAXPFOS], mcr_theta[NMAXPFOS];
bool mcr_iscreatedinsim[NMAXPFOS];

void pfoe_theta ()
{
  TFile* fin = new TFile("output.root");
  TTree* T = static_cast<TTree*>(fin->Get("evtdata"));

  TH2D* hpfoetheta_gen_1pfo = new TH2D("hpfoetheta_gen_1pfo","",100,-0.3,0.3,100,0,15);
  TH2D* hpfoetheta_gen_other = new TH2D("hpfoetheta_gen_other","",100,-0.3,0.3,100,0,15);
  TH2D* hpfoetheta_sim = new TH2D("hpfoetheta_sim","",100,-0.3,0.3,100,0,15);
  TH2D* hmcretheta_gen = new TH2D("hmcretheta_gen","",100,-0.3,0.3,100,0,15);
  TH2D* hmcretheta_sim = new TH2D("hmcretheta_sim","",100,-0.3,0.3,100,0,15);
  TCanvas* c = new TCanvas("c","",600,400);

  T->SetBranchAddress("npfos",&npfos);
  T->SetBranchAddress("pfo_e",pfo_e);
  T->SetBranchAddress("pfo_theta",pfo_theta);
  T->SetBranchAddress("mcr_e",mcr_e);
  T->SetBranchAddress("mcr_theta",mcr_theta);
  T->SetBranchAddress("mcr_iscreatedinsim",mcr_iscreatedinsim);

  int nevents = T->GetEntries();

  for (int ev = 0; ev < nevents; ev++) {
    
    T->GetEntry(ev);
    for ( int ip = 0; ip < npfos; ip++ ) {
      if (!mcr_iscreatedinsim[ip]) {
        if (npfos==1) hpfoetheta_gen_1pfo->Fill(pfo_theta[ip],pfo_e[ip]);
        else hpfoetheta_gen_other->Fill(pfo_theta[ip],pfo_e[ip]);
        hmcretheta_gen->Fill(mcr_theta[ip],mcr_e[ip]);
      } else {
        hpfoetheta_sim->Fill(pfo_theta[ip],pfo_e[ip]);
        hmcretheta_sim->Fill(mcr_theta[ip],mcr_e[ip]);
      }
    }
  }
  c->cd();
  hpfoetheta_gen_1pfo->SetMarkerColor(4);
  hpfoetheta_gen_other->SetMarkerColor(2);
  hmcretheta_gen->SetMarkerColor(3);
  hmcretheta_sim->SetMarkerColor(5);
  hmcretheta_sim->SetStats(0);
  hmcretheta_sim->SetTitle(";theta [rad];E [GeV]");
  hmcretheta_sim->Draw();
  hpfoetheta_gen_1pfo->Draw("same");
  hpfoetheta_gen_other->Draw("same");
  hmcretheta_gen->Draw("same");
  stringstream fname;
  fname << "pfo_etheta.pdf" << ends;
  c->Print(fname.str().data());
}
