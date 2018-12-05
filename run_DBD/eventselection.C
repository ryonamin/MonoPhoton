class Event
{
   public:
	//Event(TTree* T) { 
	void SetTree(TTree* T) { 
           _T = T;
           // MCParticle
           _T->SetBranchAddress("nmcps",&nmcps);
           _T->SetBranchAddress("mcp_pdg",mcp_pdg);
           _T->SetBranchAddress("mcp_index",mcp_index);
           _T->SetBranchAddress("mcp_nparents",mcp_nparents);
           _T->SetBranchAddress("mcp_parentIndex",mcp_parentIndex);
           _T->SetBranchAddress("mcp_ndaughters",mcp_ndaughters);
           _T->SetBranchAddress("mcp_daughterIndex",mcp_daughterIndex);
           _T->SetBranchAddress("mcp_e",mcp_e);
           _T->SetBranchAddress("mcp_phi",mcp_phi);
           _T->SetBranchAddress("mcp_theta",mcp_theta);
           _T->SetBranchAddress("mcp_endx",mcp_endx);
           _T->SetBranchAddress("mcp_endy",mcp_endy);
           _T->SetBranchAddress("mcp_endz",mcp_endz);
           // PFO
           _T->SetBranchAddress("npfos",&npfos);
           _T->SetBranchAddress("pfo_pdg",pfo_pdg);
           _T->SetBranchAddress("pfo_chrg",pfo_chrg);
           _T->SetBranchAddress("pfo_e",pfo_e);
           _T->SetBranchAddress("pfo_px",pfo_px);
           _T->SetBranchAddress("pfo_py",pfo_py);
           _T->SetBranchAddress("pfo_pz",pfo_pz);
           _T->SetBranchAddress("pfo_phi",pfo_phi);
           _T->SetBranchAddress("pfo_theta",pfo_theta);
           _T->SetBranchAddress("pfo_cal_x",pfo_cal_x);
           _T->SetBranchAddress("pfo_cal_y",pfo_cal_y);
           _T->SetBranchAddress("pfo_cal_z",pfo_cal_z);
           _T->SetBranchAddress("pfo_ecal_e",pfo_ecal_e);
           _T->SetBranchAddress("pfo_hcal_e",pfo_hcal_e);
           _T->SetBranchAddress("emaxphoton_e",&emaxphoton_e);
           _T->SetBranchAddress("ptmaxphoton_e",&ptmaxphoton_e);
           _T->SetBranchAddress("emaxphoton_pt_bcalcoord",&emaxphoton_pt_bcalcoord);
           _T->SetBranchAddress("ptmaxphoton_pt_bcalcoord",&ptmaxphoton_pt_bcalcoord);
           _T->SetBranchAddress("emaxphoton_phi_bcalcoord",&emaxphoton_phi_bcalcoord);
           _T->SetBranchAddress("ptmaxphoton_phi_bcalcoord",&ptmaxphoton_phi_bcalcoord);
           _T->SetBranchAddress("emaxphoton_theta_bcalcoord",&emaxphoton_theta_bcalcoord);
           _T->SetBranchAddress("ptmaxphoton_theta_bcalcoord",&ptmaxphoton_theta_bcalcoord);
           // RecoMCTruthLink 
           _T->SetBranchAddress("mcr_pdg",mcr_pdg);
           _T->SetBranchAddress("mcr_weight",mcr_weight);
           _T->SetBranchAddress("mcr_index",mcr_index);
           _T->SetBranchAddress("mcr_nparents",mcr_nparents);
           _T->SetBranchAddress("mcr_parentIndex",mcr_parentIndex);
           _T->SetBranchAddress("mcr_e",mcr_e);
           _T->SetBranchAddress("mcr_phi",mcr_phi);
           _T->SetBranchAddress("mcr_theta",mcr_theta);
           _T->SetBranchAddress("mcr_isoverlay",mcr_isoverlay);
           _T->SetBranchAddress("mcr_isOriginatedFromISR",mcr_isOriginatedFromISR);
           // BCal cluster
           _T->SetBranchAddress("nbcalclrs",&nbcalclrs);
           _T->SetBranchAddress("bcal_e",bcal_e);
           _T->SetBranchAddress("bcal_phi",bcal_phi);
           _T->SetBranchAddress("bcal_theta",bcal_theta);
           _T->SetBranchAddress("bcal_x",bcal_x);
           _T->SetBranchAddress("bcal_y",bcal_y);
           _T->SetBranchAddress("bcal_z",bcal_z);
        }

        void callGetEntry(int ev) {
          if (currentEv==ev) return;
          _T->GetEntry(ev);
          currentEv = ev; 

          int nAllPhoton_MC_per_evt = 0;
          int nAcceptablePhoton_MC_per_evt = 0;
          signal_index = -1.;
          signal_e = -1.;
          nISRPhotons = 0;
          for (int i = 0; i < nmcps; i++) {
            if (mcp_pdg[i]==22 && mcp_parentIndex[i][0]==2&&mcp_parentIndex[i][1]==3) { // select only initial ones
              nAllPhoton_MC_per_evt++;
              if (mcp_ndaughters[i] < 2 && mcp_e[i]>2. && mcp_theta[i]>7/180.*TMath::Pi() && mcp_theta[i]<173/180.*TMath::Pi()) { 
                // select photon (non converted photon only)
                nAcceptablePhoton_MC_per_evt++;
                signal_index = i;
                signal_e = mcp_e[i];
                signal_phi = mcp_phi[i];
                signal_theta = mcp_theta[i];
                nISRPhotons++;
              }
            }
          }
          _isAcceptableEvent = (nAllPhoton_MC_per_evt==1&&nAcceptablePhoton_MC_per_evt==1); 

        }

        bool isAcceptableEvent(int ev) {
          callGetEntry(ev);
          return _isAcceptableEvent;
        }

        bool isPassedPtCut(int ev) {
          callGetEntry(ev);
          for (int i = 0; i < npfos; i++) {
            TLorentzVector p4(pfo_px[i],pfo_py[i],pfo_pz[i],pfo_e[i]);
            if (abs(pfo_pdg[i])==11) { // electron/positron cases
              if (p4.Pt()>0.5) return false;
            } else { // other pfos
              if (p4.Pt()>3.) return false;
            } 
          }
          return true;
        }

        bool isPassedECut(int ev) {
          callGetEntry(ev);
          float esum = 0.;
          float esum_wo_pi_n = 0.;
          for (int i = 0; i < npfos; i++) {
            if (pfo_e[i]<5.) return false; // individual energy cut
            esum += pfo_e[i];
            if (!(abs(pfo_pdg[i])==211||abs(pfo_pdg[i])==111||abs(pfo_pdg[i])==2112)) {
              esum_wo_pi_n += pfo_e[i];
            }
          }
          if (esum>30.) return false;
          if (esum_wo_pi_n>10.) return false;
          return true;
        }

        int getSignalIndex() const { 
          if (signal_index>0) return signal_index; 
          else return -1.;
        }

        float getSignalE() const { 
          if (signal_e>0) return signal_e; 
          else return -1.;
        }

        float getSignalTheta() const { 
          if (signal_e>0) return signal_theta; 
          else return -1.;
        }

        int getNISRPhotons() const { return nISRPhotons; }

        void process(string fname);

        struct Outputs
        {
          TH1F* hE_photon;
          TH1F* hNrecNgen_photon;
          TH2F* hNrecNgenEmc_photon;
          TH2F* hNrecNgenCostheta_photon;
          TGraphErrors* gNrecNgenEmc_photon;
          TGraphErrors* gNrecNgenCostheta_photon;

          TH1F* hPt_ep_isr;
          TH1F* hPt_ep_ol;
          TH1F* hPt_ep_other;
          TH1F* hPt_pfo_isr;
          TH1F* hPt_pfo_ol;
          TH1F* hPt_pfo_other;
          TH1F* hE_photon_ol;
          TH1F* hE_photon_electron;
          TH1F* hE_photon_rest;
          TH1F* hE_V0_ol;
          TH1F* hE_V0_rest;
          TH1F* hE_neutron_ol;
          TH1F* hE_neutron_rest;
          TH1F* hE_electron_ol;
          TH1F* hE_electron_rest;
          TH1F* hE_muon_ol;
          TH1F* hE_muon_rest;
          TH1F* hE_pion_ol;
          TH1F* hE_pion_rest;
          TH1F* hEvis_pfo;
          TH1F* hEvis_pfo_wo_pi_n;
          TH1F* hNBcalClusters;
          TH1F* hNBcalClusters1ISR;
          TH1F* hNBcalClustersMultiISR;
        };

        Outputs outputs;
   private:
	TTree* _T;
        int currentEv = -1;
        static const int NMAX = 10000;
        // MCParticle
	int signal_index;
        float signal_e, signal_phi, signal_theta;
        bool _isAcceptableEvent;
        bool _isSignalLikeEvent;
        int nISRPhotons;
   public:
        int nmcps, mcp_pdg[NMAX], mcp_index[NMAX];
        int mcp_nparents[NMAX], mcp_parentIndex[NMAX][10];
        int mcp_ndaughters[NMAX], mcp_daughterIndex[NMAX][10];
        float mcp_e[NMAX], mcp_phi[NMAX], mcp_theta[NMAX];
        float mcp_endx[NMAX], mcp_endy[NMAX], mcp_endz[NMAX];
        int npfos, pfo_pdg[NMAX],mcr_pdg[NMAX], mcr_index[NMAX], mcr_nparents[NMAX], mcr_parentIndex[NMAX][10];
        float pfo_e[NMAX], pfo_px[NMAX], pfo_py[NMAX], pfo_pz[NMAX], pfo_phi[NMAX], pfo_theta[NMAX];
        float pfo_cal_x[NMAX],  pfo_cal_y[NMAX],  pfo_cal_z[NMAX];
        float pfo_ecal_e[NMAX], pfo_hcal_e[NMAX];
        float pfo_chrg[NMAX];
        float mcr_e[NMAX], mcr_phi[NMAX], mcr_theta[NMAX], mcr_weight[NMAX];
        bool mcr_isOriginatedFromISR[NMAX], mcr_isoverlay[NMAX];
        int nbcalclrs;
        float bcal_e[NMAX], bcal_phi[NMAX], bcal_theta[NMAX];
        float bcal_x[NMAX], bcal_y[NMAX], bcal_z[NMAX];
        float emaxphoton_e, ptmaxphoton_pt_bcalcoord, emaxphoton_pt_bcalcoord, ptmaxphoton_e;
        float ptmaxphoton_phi_bcalcoord, emaxphoton_phi_bcalcoord;
        float ptmaxphoton_theta_bcalcoord, emaxphoton_theta_bcalcoord;
};

