
int npfos;
const int NMAXPFOS = 500;
float pfo_x[NMAXPFOS], pfo_y[NMAXPFOS], pfo_z[NMAXPFOS];

void DrawSlice(float zstart = -3000., float zend = 3000., TTree* T = 0, TH2D* h = 0, TCanvas* c=0)
{
  h->Reset();
  int nevents = T->GetEntries();

  for (int ev = 0; ev < nevents; ev++) {
    
    T->GetEntry(ev);
    for ( int ip = 0; ip < npfos; ip++ ) {
      float x = pfo_x[ip];
      float y = pfo_y[ip];
      float z = pfo_z[ip];
      if (z >= zstart && z<zend) {
        h->Fill(x,y);
      }
    }
  }
  c->cd();
  h->Draw("COL");
  stringstream fname;
  fname << "slice_" << zstart << "_" << zend << ".pdf" << ends;
  c->Print(fname.str().data());
}

void hit ()
{
  TFile* fin = new TFile("output.root");
  TTree* data = static_cast<TTree*>(fin->Get("evtdata"));

  TH2D* h = new TH2D("h","",100,-1000,1000,100,-1000,1000);
  TCanvas* c = new TCanvas("c","",600,400);

  data->SetBranchAddress("npfos",&npfos);
  data->SetBranchAddress("pfo_cal_x",pfo_x);
  data->SetBranchAddress("pfo_cal_y",pfo_y);
  data->SetBranchAddress("pfo_cal_z",pfo_z);

#if 0
  float zstart = 2350.; 
  float dz = 50.;
#else
  float zstart = -2850.; 
  float dz = 50.;
#endif
  for (int i = 0; i < 10; i++) {
    DrawSlice(zstart+i*dz,zstart+(i+1)*dz,data,h,c);
  }

}
