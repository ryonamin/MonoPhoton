
int npfos;
const int NMAXPFOS = 500;
float pfo_e[NMAXPFOS], pfo_phi[NMAXPFOS], mcr_e[NMAXPFOS];
bool mcr_iscreatedinsim[NMAXPFOS];

void DrawSlice(float phistart, float phiend, TTree* T = 0, TH1D* h = 0, TCanvas* c=0)
{
  h->Reset();
  int nevents = T->GetEntries();

  for (int ev = 0; ev < nevents; ev++) {
    
    T->GetEntry(ev);
    for ( int ip = 0; ip < npfos; ip++ ) {
      float phi = pfo_phi[ip];
      if (phi>=phistart&&phi<phiend) {
        if (!mcr_iscreatedinsim[ip]) h->Fill(pfo_e[ip]-mcr_e[ip]);
      }
    }
  }
  c->cd();
  h->Draw();
  stringstream fname;
  fname << "resol_" << phistart << "_" << phiend << ".pdf" << ends;
  c->Print(fname.str().data());
}

void resol_phi ()
{
  TFile* fin = new TFile("output.root");
  TTree* data = static_cast<TTree*>(fin->Get("evtdata"));

  TH1D* h = new TH1D("h","",1000,-20,20);
  TCanvas* c = new TCanvas("c","",600,400);

  data->SetBranchAddress("npfos",&npfos);
  data->SetBranchAddress("pfo_e",pfo_e);
  data->SetBranchAddress("pfo_phi",pfo_phi);
  data->SetBranchAddress("mcr_e",mcr_e);
  data->SetBranchAddress("mcr_iscreatedinsim",mcr_iscreatedinsim);

  const int nbin = 12;
  const float dphi = 2*TMath::Pi()/nbin; 
  for (int i = 0; i < nbin; i++) {
    float phistart = (i-nbin/2)*dphi; // -2pi to 2pi
    float phiend = (i-nbin/2+1)*dphi;
    DrawSlice(phistart,phiend,data,h,c);
  }

}
