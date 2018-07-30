
int npfos;
const int NMAXPFOS = 500;
float pfo_e[NMAXPFOS], pfo_phi[NMAXPFOS], mcr_e[NMAXPFOS];

void DrawSlice(float phistart, float phiend, TTree* T = 0, TH1D* h = 0, TCanvas* c=0)
{
  h->Reset();
  int nevents = T->GetEntries();

  for (int ev = 0; ev < nevents; ev++) {
    
    T->GetEntry(ev);
    for ( int ip = 0; ip < npfos; ip++ ) {
      float phi = pfo_phi[ip];
      if (phi>=phistart&&phi<phiend) {
        //h->Fill(pfo_e[ip]-mcr_e[ip]);
        h->Fill(pfo_e[ip]);
      }
    }
  }
  c->cd();
  h->Draw();
  stringstream fname;
  fname << "pfoe_phi." << phistart << "_" << phiend << ".pdf" << ends;
  c->Print(fname.str().data());
}

void pfoe_phi ()
{
  TFile* fin = new TFile("output.root");
  TTree* data = static_cast<TTree*>(fin->Get("evtdata"));

  TH1D* h = new TH1D("h","",100,0,12);
  TCanvas* c = new TCanvas("c","",600,400);

  data->SetBranchAddress("npfos",&npfos);
  data->SetBranchAddress("pfo_e",pfo_e);
  data->SetBranchAddress("pfo_phi",pfo_phi);
  data->SetBranchAddress("mcr_e",mcr_e);

  const int nbin = 12;
  const float dphi = 2*TMath::Pi()/nbin; 
  for (int i = 0; i < nbin; i++) {
    float phistart = (i-nbin/2)*dphi; // -2pi to 2pi
    float phiend = (i-nbin/2+1)*dphi;
    DrawSlice(phistart,phiend,data,h,c);
  }

}
