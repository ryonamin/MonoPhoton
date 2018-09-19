#ifndef MonoPhotonProcessor_h
#define MonoPhotonProcessor_h 1

#include "marlin/Processor.h"
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/MCParticle.h>
#include <EVENT/Vertex.h>
#include "lcio.h"
#include <string>

#include "UTIL/LCRelationNavigator.h"

class TFile;
class TTree;
class TVector3;
const int NMAX_PFOS   = 500;
const int NMAX_MCPS   = 10000;
const int NMAX_CLRS   = 500;
const int NMAX_PARENTS   = 10;
const int NMAX_DAUGHTERS = 10;

using namespace lcio ;
using namespace marlin ;


/**  Example processor for marlin.
 * 
 *  If compiled with MARLIN_USE_AIDA 
 *  it creates a histogram (cloud) of the MCParticle energies.
 * 
 *  <h4>Input - Prerequisites</h4>
 *  Needs the collection of MCParticles.
 *
 *  <h4>Output</h4> 
 *  A histogram.
 * 
 * @param CollectionName Name of the MCParticle collection
 * 
 * @author F. Gaede, DESY
 * @version $Id: MonoPhotonProcessor.h,v 1.4 2005-10-11 12:57:39 gaede Exp $ 
 */

class MonoPhotonProcessor : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new MonoPhotonProcessor ; }
  
  
  MonoPhotonProcessor() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  
 protected:

 /** Prepare NTuple 
  */ 
  void makeNTuple() ;

 /** utility functions 
  */ 
  TVector3 getIP(EVENT::MCParticle* p); 

  /** MCParticle pointer and index 
   */
  std::map<MCParticle*,int> _mcpmap;

  /** Get most probable MCParticle pointer form PFOs 
   */
  MCParticle* getBestMCParticleOf(ReconstructedParticle* p, LCRelationNavigator* nav);

  /** Input collection name.
   */
  std::string _colMCP ;
  std::string _colPFO ;
  std::string _colMCPFORelation ;
  LCRelationNavigator* _navpfo;

  int _nRun ;
  int _nEvt ;

  TFile* _otfile;
  TTree* _evtdata;

  struct EVTFILLDATA {
    int    evt ;       // event number
    int    npfos ;     // # of PFOs 
    float  pfo_e[NMAX_PFOS];
    float  pfo_px[NMAX_PFOS];
    float  pfo_py[NMAX_PFOS];
    float  pfo_pz[NMAX_PFOS];
    float  pfo_phi[NMAX_PFOS];
    float  pfo_theta[NMAX_PFOS];
    //float  pfo_startx[NMAX_PFOS];
    //float  pfo_starty[NMAX_PFOS];
    //float  pfo_startz[NMAX_PFOS];
    //float  pfo_endx[NMAX_PFOS];
    //float  pfo_endy[NMAX_PFOS];
    //float  pfo_endz[NMAX_PFOS];
    float  pfo_chrg[NMAX_PFOS];
    int    pfo_pdg[NMAX_PFOS];
    int    pfo_goodnessOfPid[NMAX_PFOS];
    int    pfo_ntrk[NMAX_PFOS];
    float  pfo_d0[NMAX_PFOS];
    float  pfo_d0sig[NMAX_PFOS];
    float  pfo_z0[NMAX_PFOS];
    float  pfo_z0sig[NMAX_PFOS];
    float  pfo_trkphi[NMAX_PFOS];
    float  pfo_omega[NMAX_PFOS];
    float  pfo_tanlambda[NMAX_PFOS];
    int    pfo_nclus[NMAX_PFOS]; 
    float  pfo_cal_x[NMAX_PFOS];
    float  pfo_cal_y[NMAX_PFOS];
    float  pfo_cal_z[NMAX_PFOS];
    float  pfo_ecal_e[NMAX_PFOS];
    float  pfo_hcal_e[NMAX_PFOS];
    float  pfo_yoke_e[NMAX_PFOS];
    float  pfo_lcal_e[NMAX_PFOS];
    float  pfo_lhcal_e[NMAX_PFOS];
    float  pfo_bcal_e[NMAX_PFOS];

    // MC Relation 
    int    nmcr[NMAX_PFOS];
    float  mcr_weight[NMAX_PFOS];
    int    mcr_index[NMAX_PFOS];
    float  mcr_e[NMAX_PFOS];
    float  mcr_px[NMAX_PFOS];
    float  mcr_py[NMAX_PFOS];
    float  mcr_pz[NMAX_PFOS];
    float  mcr_phi[NMAX_PFOS];
    float  mcr_theta[NMAX_PFOS];
    float  mcr_chrg[NMAX_PFOS];
    float  mcr_startx[NMAX_PFOS];
    float  mcr_starty[NMAX_PFOS];
    float  mcr_startz[NMAX_PFOS];
    float  mcr_endx[NMAX_PFOS];
    float  mcr_endy[NMAX_PFOS];
    float  mcr_endz[NMAX_PFOS];
    int    mcr_pdg[NMAX_PFOS];
    int    mcr_nparents[NMAX_PFOS];
    int    mcr_parentIndex[NMAX_PFOS][NMAX_PARENTS];
    int    mcr_ndaughters[NMAX_PFOS];
    int    mcr_daughterIndex[NMAX_PFOS][NMAX_DAUGHTERS];
    int    mcr_genstatus[NMAX_PFOS];
    int    mcr_simstatus[NMAX_PFOS];
    bool   mcr_iscreatedinsim[NMAX_PFOS];

    // MC Info
    int    nmcps ;     // # of MC Particles 
    float  ipx;
    float  ipy;
    float  ipz;
    int    mcp_index[NMAX_MCPS];
    float  mcp_e[NMAX_MCPS];
    float  mcp_px[NMAX_MCPS];
    float  mcp_py[NMAX_MCPS];
    float  mcp_pz[NMAX_MCPS];
    float  mcp_phi[NMAX_MCPS];
    float  mcp_theta[NMAX_MCPS];
    float  mcp_chrg[NMAX_MCPS];
    float  mcp_startx[NMAX_MCPS];
    float  mcp_starty[NMAX_MCPS];
    float  mcp_startz[NMAX_MCPS];
    float  mcp_endx[NMAX_MCPS];
    float  mcp_endy[NMAX_MCPS];
    float  mcp_endz[NMAX_MCPS];
    int    mcp_pdg[NMAX_MCPS];
    int    mcp_nparents[NMAX_MCPS];
    int    mcp_parentIndex[NMAX_MCPS][NMAX_PARENTS];
    int    mcp_ndaughters[NMAX_MCPS];
    int    mcp_daughterIndex[NMAX_MCPS][NMAX_DAUGHTERS];
    int    mcp_genstatus[NMAX_MCPS];
    int    mcp_simstatus[NMAX_MCPS];
    bool   mcp_iscreatedinsim[NMAX_MCPS];

    int    nclrhits ;    // # of Cal hit clusters;
    float  clr_x[NMAX_CLRS];
    float  clr_y[NMAX_CLRS];
    float  clr_z[NMAX_CLRS];
  };

  EVTFILLDATA _data;

  // output root file name
  std::string _rootfilename;

} ;

#endif



