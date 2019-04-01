//--------------------------------------------------------------------------------------------------
// $Id $
//
// EGammaMvaPhotonEstimator
//
// Helper Class for evaluating MVA photon ID
//--------------------------------------------------------------------------------------------------




#ifndef EGammaMvaPhotonEstimator_H
#define EGammaMvaPhotonEstimator_H

#include <vector>
#include <TROOT.h>
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

using namespace std;

#include "DataFormats/PatCandidates/interface/Photon.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

/* #include "DataFormats/TrackReco/interface/Track.h" */
/* #include "DataFormats/GsfTrackReco/interface/GsfTrack.h" */
/* #include "DataFormats/TrackReco/interface/TrackFwd.h" */
/* #include "DataFormats/EgammaCandidates/interface/GsfElectron.h" */
/* #include "DataFormats/EgammaReco/interface/SuperCluster.h" */
/* #include "DataFormats/MuonReco/interface/Muon.h" */
/* #include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h" */
/* #include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h" */
/* #include "DataFormats/Common/interface/RefToPtr.h" */
/* #include "DataFormats/VertexReco/interface/Vertex.h" */
/* #include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h" */
/* #include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h" */
/* #include "TrackingTools/IPTools/interface/IPTools.h" */
/* #include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h" */
/* #include "DataFormats/Common/interface/RefToPtr.h" */

/* #include "DataFormats/PatCandidates/interface/Electron.h" */

using namespace reco;

class EGammaMvaPhotonEstimator{
  public:
    EGammaMvaPhotonEstimator();
    ~EGammaMvaPhotonEstimator(); 
  
    enum MVAType {
      kPhotonMVATypeDefault = 0                     // MVA for triggering electrons
    };
  
    void     initialize( std::vector<std::string> methodName ,
			 std::vector<std::string> weightsfiles ,
			 EGammaMvaPhotonEstimator::MVAType type);
    
    Bool_t   isInitialized() const { return fisInitialized; }
    UInt_t   GetMVABin(double sceta ) const;
    
    Double_t mvaValue(const pat::Photon& pho, 				
		      double varRho,
		      double varPhotonIsoRaw,
		      double varChargedIsoRaw,
		      double varWorstVertexChargedIsoRaw,
		      noZS::EcalClusterLazyTools *myEcalClusterLazyTool,
                      bool printDebug = kFALSE);

  private:

    std::vector<TMVA::Reader*> fTMVAReader;
    std::vector<string>        fMethodnames;
    Bool_t                     fisInitialized;
    MVAType                    fMVAType;
    UInt_t                     fNMVABins;

    Float_t                    fMVAVar_Phi;
    Float_t                    fMVAVar_R9;
    Float_t                    fMVAVar_Sieie;    
    Float_t                    fMVAVar_Sieip; 
    Float_t                    fMVAVar_E1x3OverE5x5;
    Float_t                    fMVAVar_E2x2OverE5x5;
    Float_t                    fMVAVar_E2x5MaxOverE5x5;
    Float_t                    fMVAVar_Sceta;
    Float_t                    fMVAVar_RawE;
    Float_t                    fMVAVar_ScEtaWidth;
    Float_t                    fMVAVar_ScPhiWidth;
    Float_t                    fMVAVar_Rho;
    Float_t                    fMVAVar_PhotonIsoRaw;
    Float_t                    fMVAVar_ChargedIsoRaw;
    Float_t                    fMVAVar_WorstVertexChargedIsoRaw;
    Float_t                    fMVAVar_ESEnOverRawE;
    Float_t                    fMVAVar_ESEffSigmaRR;

    Float_t                    fMVAVar_pt;
    Float_t                    fMVAVar_eta;
 
};

#endif
