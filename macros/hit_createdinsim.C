
int npfos;
const int NMAXPFOS = 500;
float x[NMAXPFOS], y[NMAXPFOS], z[NMAXPFOS];
float mcr_e[NMAXPFOS], mcr_theta[NMAXPFOS];
bool mcr_iscreatedinsim[NMAXPFOS];

void hit_createdinsim ()
{
  TFile* fin = new TFile("output.root");
  TTree* T = static_cast<TTree*>(fin->Get("evtdata"));

  TH2D* hxy = new TH2D("hxy","",1000,-1000,1000,1000,-1000,1000);
  TH2D* hrz = new TH2D("hrz","",1000,-3000,3000,1000,0,1000);
  TCanvas* cxy = new TCanvas("cxy","",600,400);
  TCanvas* crz = new TCanvas("crz","",600,400);

  T->SetBranchAddress("npfos",&npfos);
  T->SetBranchAddress("pfo_cal_x",x);
  T->SetBranchAddress("pfo_cal_y",y);
  T->SetBranchAddress("pfo_cal_z",z);
  T->SetBranchAddress("mcr_iscreatedinsim",mcr_iscreatedinsim);

  int nevents = T->GetEntries();

  for (int ev = 0; ev < nevents; ev++) {
    
    T->GetEntry(ev);
    for ( int ip = 0; ip < npfos; ip++ ) {
      if (mcr_iscreatedinsim[ip]) {
        hxy->Fill(x[ip],y[ip]);
        hrz->Fill(z[ip],TMath::Sqrt(x[ip]*x[ip]+y[ip]*y[ip]));
      }
    }
  }
  cxy->cd();
  hxy->Draw();
  crz->cd();
  hrz->Draw();
  stringstream fnamexy,fnamerz;
  fnamexy << "hit_xy_createdinsim.pdf" << ends;
  fnamerz << "hit_rz_createdinsim.pdf" << ends;
  cxy->Print(fnamexy.str().data());
  crz->Print(fnamerz.str().data());
}
