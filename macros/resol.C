
int npfos;
const int NMAXPFOS = 500;
float pfo_e[NMAXPFOS], mcr_e[NMAXPFOS];
bool mcr_iscreatedinsim[NMAXPFOS];

void DrawSlice(TTree* T = 0, TH1D* h1 = 0, TH1D* h2 = 0, TH1D* h3 = 0, TCanvas* c=0)
{
  h1->Reset();
  h2->Reset();
  int nevents = T->GetEntries();

  for (int ev = 0; ev < nevents; ev++) {
    
    T->GetEntry(ev);
    for ( int ip = 0; ip < npfos; ip++ ) {
      h1->Fill(pfo_e[ip]-mcr_e[ip]);
      if ( !mcr_iscreatedinsim[ip]) h2->Fill(pfo_e[ip]-mcr_e[ip]);
      else h3->Fill(pfo_e[ip]-mcr_e[ip]);
    }
  }
  c->cd();
  c->SetLogy();
  h1->SetTitle(";E_{PFO}-E_{MC} [GeV];");
  h1->Draw();
  h1->SetLineStyle(7);
  h2->SetLineColor(4);
  h2->Draw("same");
  h3->SetLineColor(2);
  h3->Draw("same");
  TLegend *leg = new TLegend(0.15,0.85,0.45,0.7);
  leg->SetTextFont(132);
  leg->SetTextSize(0.03); 
  leg->SetTextAlign(12);
  leg->AddEntry(h1,"All PFOs","f"); 
  leg->AddEntry(h2,"PFOs created in generator","f"); 
  leg->AddEntry(h3,"PFOs created in simulation","f"); 
  leg->Draw();
  stringstream fname;
  fname << "resol.pdf" << ends;
  c->Print(fname.str().data());
}

void resol ()
{
  TFile* fin = new TFile("output.root");
  TTree* data = static_cast<TTree*>(fin->Get("evtdata"));

  TH1D* hall = new TH1D("hall","",300,-15,15);
  TH1D* hgen = new TH1D("hgen","",300,-15,15);
  TH1D* hsim = new TH1D("hsim","",300,-15,15);
  TCanvas* c = new TCanvas("c","",600,400);

  data->SetBranchAddress("npfos",&npfos);
  data->SetBranchAddress("pfo_e",pfo_e);
  data->SetBranchAddress("mcr_e",mcr_e);
  data->SetBranchAddress("mcr_iscreatedinsim",mcr_iscreatedinsim);

  DrawSlice(data,hall,hgen,hsim,c);

}
