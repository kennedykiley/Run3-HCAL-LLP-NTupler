//--------------------------------------------------------------------------------------------------
// $Id $
//
// ElectronMVAEstimatorRun2NonTrig
//
// Helper Class for applying MVA electron ID selection
//
// Authors: D.Benedetti, E.DiMaro, S.Xie
//--------------------------------------------------------------------------------------------------




#ifndef ElectronMVAEstimatorRun2NonTrig_H
#define ElectronMVAEstimatorRun2NonTrig_H

#include <vector>
#include <TROOT.h>
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

using namespace std;

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "DataFormats/PatCandidates/interface/Electron.h"

using namespace reco;

class ElectronMVAEstimatorRun2NonTrig{
  public:
    ElectronMVAEstimatorRun2NonTrig();
    ~ElectronMVAEstimatorRun2NonTrig(); 
  
    enum MVAType {
      kPHYS14  = 0,                      // MVA for PHYS14
    };
  
    const uint nCategories = 6;
    enum mvaCategories {
      UNDEFINED = -1,
      CAT_EB1_PT5to10  = 0,
      CAT_EB2_PT5to10  = 1,
      CAT_EE_PT5to10   = 2,
      CAT_EB1_PT10plus = 3,
      CAT_EB2_PT10plus = 4,
      CAT_EE_PT10plus  = 5
    };
    
    void     initialize( std::string methodName,
                         std::string weightsfile,
                         ElectronMVAEstimatorRun2NonTrig::MVAType type);
    void     initialize( std::string methodName,
                         ElectronMVAEstimatorRun2NonTrig::MVAType type,
                         Bool_t useBinnedVersion,
                         std::vector<std::string> weightsfiles );
    
    Bool_t   isInitialized() const { return fisInitialized; }
    UInt_t   GetMVABin(double eta,double pt ) const;
    
    void bindVariables();
    

    Double_t mvaValue(const reco::GsfElectron& ele, 
                      const reco::Vertex& vertex, 
                      const TransientTrackBuilder& transientTrackBuilder,
                      noZS::EcalClusterLazyTools myEcalCluster,
		      edm::Handle<vector<reco::Conversion> > conversions,
		      const math::XYZPoint beamspotPosition,
                      bool printDebug = kFALSE);


    Double_t mvaValue(const pat::Electron& ele,
		      edm::Handle<vector<reco::Conversion> > conversions,
		      const math::XYZPoint beamspotPosition,
                      bool printDebug);


  private:

    std::vector<TMVA::Reader*> fTMVAReader;
    std::string                fMethodname;
    Bool_t                     fisInitialized;
    MVAType                    fMVAType;

    Float_t                    fMVAVar_see;
    Float_t                    fMVAVar_spp;
    Float_t                    fMVAVar_OneMinusE1x5E5x5;
    Float_t                    fMVAVar_R9;
    Float_t                    fMVAVar_etawidth;
    Float_t                    fMVAVar_phiwidth;
    Float_t                    fMVAVar_HoE;
    Float_t                    fMVAVar_PreShowerOverRaw;
    Float_t                    fMVAVar_kfhits;    
    Float_t                    fMVAVar_kfchi2;
    Float_t                    fMVAVar_gsfchi2;
    Float_t                    fMVAVar_fbrem;
    Float_t                    fMVAVar_gsfhits;
    Float_t                    fMVAVar_expectedMissingInnerHits;
    Float_t                    fMVAVar_convVtxFitProbability;
    Float_t                    fMVAVar_EoP;
    Float_t                    fMVAVar_eleEoPout;
    Float_t                    fMVAVar_IoEmIoP;
    Float_t                    fMVAVar_deta;
    Float_t                    fMVAVar_dphi;
    Float_t                    fMVAVar_detacalo;
    Float_t                    fMVAVar_pt;    
    Float_t                    fMVAVar_isBarrel; 
    Float_t                    fMVAVar_isEndcap; 
    Float_t                    fMVAVar_SCeta; 
    Float_t                    fMVAVar_eClass; 
    Float_t                    fMVAVar_pfRelIso; 
    Float_t                    fMVAVar_expectedInnerHits; 
    Float_t                    fMVAVar_vtxconv; 
    Float_t                    fMVAVar_mcEventWeight; 
    Float_t                    fMVAVar_mcCBmatchingCategory; 
 
};

#endif
