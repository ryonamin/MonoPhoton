
int npfos;
const int NMAXPFOS = 500;
float pfo_phi[NMAXPFOS];

void DrawSlice(TTree* T = 0, TH1D* h = 0, TCanvas* c=0)
{
  h->Reset();
  int nevents = T->GetEntries();

  for (int ev = 0; ev < nevents; ev++) {
    
    T->GetEntry(ev);
    for ( int ip = 0; ip < npfos; ip++ ) {
      h->Fill(pfo_phi[ip]);
    }
  }
  c->cd();
  h->Draw();
  stringstream fname;
  fname << "phi.pdf" << ends;
  c->Print(fname.str().data());
}

void phi ()
{
  TFile* fin = new TFile("output.root");
  TTree* data = static_cast<TTree*>(fin->Get("evtdata"));

  TH1D* h = new TH1D("h","",100,-4,4);
  TCanvas* c = new TCanvas("c","",600,400);

  data->SetBranchAddress("npfos",&npfos);
  data->SetBranchAddress("pfo_phi",pfo_phi);

  DrawSlice(data,h,c);

}
