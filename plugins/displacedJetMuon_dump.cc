// -*- C++ -*-
// Class:      displacedJetMuon_dump
/*
  Description: Base class for miniAOD analysis with CRAB
*/
//         Author:  Cristián Peña and Si Xie.
//         Created:  Thu, 26 March 2019 15:00:06 GMT

#include "displacedJetMuon_dump.h"
#include <fastjet/internal/base.hh>
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/PseudoJet.hh"


//------ Constructors and destructor ------//
displacedJetMuon_dump::displacedJetMuon_dump(const edm::ParameterSet& iConfig):
  //get inputs from config file
  isData_(iConfig.getParameter<bool> ("isData")),
  //isFourJet_(iConfig.getParameter<bool> ("isFourJet")),
  useGen_(iConfig.getParameter<bool> ("useGen")),
  isFastsim_(iConfig.getParameter<bool> ("isFastsim")),
  //isQCD_(iConfig.getParameter<bool> ("isQCD")),
  enableTriggerInfo_(iConfig.getParameter<bool> ("enableTriggerInfo")),
  enableCaloJet_(iConfig.getParameter<bool> ("enableCaloJet")),
  enableGenLLPInfo_(iConfig.getParameter<bool> ("enableGenLLPInfo")),
  enableEcalRechits_(iConfig.getParameter<bool> ("enableEcalRechits")),
  readGenVertexTime_(iConfig.getParameter<bool> ("readGenVertexTime")),
  verticesToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  tracksTag_(consumes<edm::View<reco::Track> >(iConfig.getParameter<edm::InputTag>("tracks"))),
  cscSegmentInputToken_(consumes<CSCSegmentCollection>(edm::InputTag("cscSegments"))),
  dtSegmentInputToken_(consumes<DTRecSegment4DCollection>(edm::InputTag("dt4DCosmicSegments"))),
  rpcRecHitInputToken_(consumes<RPCRecHitCollection>(edm::InputTag("rpcRecHits"))),
  muonsToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  electronsToken_(consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  tausToken_(consumes<reco::PFTauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
  photonsToken_(consumes<reco::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
  jetsCaloToken_(consumes<reco::CaloJetCollection>(iConfig.getParameter<edm::InputTag>("jetsCalo"))),
  // jetsPFToken_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jetsPF"))),
  jetsToken_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
  jetsPuppiToken_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jetsPuppi"))),
  jetsAK8Token_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jetsAK8"))),
  PFCandsToken_(consumes<reco::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
  PFClustersToken_(consumes<reco::PFClusterCollection>(iConfig.getParameter<edm::InputTag>("pfClusters"))),
  //genParticlesToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"))),
  //genParticlesToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"))),
  genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))),
  genJetsToken_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJets"))),
  triggerBitsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBits"))),
  hepMCToken_(consumes<edm::HepMCProduct>(iConfig.getParameter<edm::InputTag>("hepMC"))),
  //triggerObjectsToken_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
  //triggerPrescalesToken_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("triggerPrescales"))),
  genMetCaloToken_(consumes<reco::GenMETCollection>(iConfig.getParameter<edm::InputTag>("genMetsCalo"))),
  genMetTrueToken_(consumes<reco::GenMETCollection>(iConfig.getParameter<edm::InputTag>("genMetsTrue"))),
  metToken_(consumes<reco::PFMETCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
//  metNoHFToken_(consumes<reco::PFMETCollection>(iConfig.getParameter<edm::InputTag>("metsNoHF"))),
  metPuppiToken_(consumes<reco::PFMETCollection>(iConfig.getParameter<edm::InputTag>("metsPuppi"))),
  metFilterBitsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("metFilterBits"))),
  //hbheNoiseFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("hbheNoiseFilter"))),
  //hbheTightNoiseFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("hbheTightNoiseFilter"))),
  //hbheIsoNoiseFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("hbheIsoNoiseFilter"))),
  //badChargedCandidateFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("BadChargedCandidateFilter"))),
  //badMuonFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("BadMuonFilter"))),
//  lheRunInfoTag_(iConfig.getParameter<edm::InputTag>("lheInfo")),
//  lheRunInfoToken_(consumes<LHERunInfoProduct,edm::InRun>(lheRunInfoTag_)),
//  lheInfoToken_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheInfo"))),
  genInfoToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genInfo"))),
  genLumiHeaderToken_(consumes<GenLumiInfoHeader,edm::InLumi>(edm::InputTag("generator",""))),
  puInfoToken_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("puInfo"))),
  //hcalNoiseInfoToken_(consumes<HcalNoiseSummary>(iConfig.getParameter<edm::InputTag>("hcalNoiseInfo"))),
  secondaryVerticesToken_(consumes<vector<reco::VertexCompositePtrCandidate> >(iConfig.getParameter<edm::InputTag>("secondaryVertices"))),
  rhoAllToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoAll"))),
  rhoFastjetAllToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetAll"))),
  rhoFastjetAllCaloToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetAllCalo"))),
  rhoFastjetCentralCaloToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetCentralCalo"))),
  rhoFastjetCentralChargedPileUpToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetCentralChargedPileUp"))),
  rhoFastjetCentralNeutralToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetCentralNeutral"))),
  beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
  ebRecHitsToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > >(iConfig.getParameter<edm::InputTag>("ebRecHits"))),
  eeRecHitsToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > >(iConfig.getParameter<edm::InputTag>("eeRecHits"))),
  esRecHitsToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > >(iConfig.getParameter<edm::InputTag>("esRecHits"))),
  ebeeClustersToken_(consumes<vector<reco::CaloCluster> >(iConfig.getParameter<edm::InputTag>("ebeeClusters"))),
  esClustersToken_(consumes<vector<reco::CaloCluster> >(iConfig.getParameter<edm::InputTag>("esClusters"))),
  conversionsToken_(consumes<vector<reco::Conversion> >(iConfig.getParameter<edm::InputTag>("conversions"))),
  singleLegConversionsToken_(consumes<vector<reco::Conversion> >(iConfig.getParameter<edm::InputTag>("singleLegConversions"))),
  gedGsfElectronCoresToken_(consumes<vector<reco::GsfElectronCore> >(iConfig.getParameter<edm::InputTag>("gedGsfElectronCores"))),
  gedPhotonCoresToken_(consumes<vector<reco::PhotonCore> >(iConfig.getParameter<edm::InputTag>("gedPhotonCores"))),
  generalTrackToken_(consumes<std::vector<reco::Track>>(edm::InputTag("generalTracks")))
  //superClustersToken_(consumes<vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("superClusters"))),
  //  lostTracksToken_(consumes<vector<reco::PFCandidate> >(iConfig.getParameter<edm::InputTag>("lostTracks")))
  // mvaGeneralPurposeValuesMapToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaGeneralPurposeValuesMap"))),
  // mvaGeneralPurposeCategoriesMapToken_(consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("mvaGeneralPurposeCategoriesMap"))),
  // mvaHZZValuesMapToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaHZZValuesMap"))),
  // mvaHZZCategoriesMapToken_(consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("mvaHZZCategoriesMap")))
{
 
  if(readGenVertexTime_) genParticles_t0_Token_ = consumes<float>(iConfig.getParameter<edm::InputTag>("genParticles_t0"));
 
}

displacedJetMuon_dump::~displacedJetMuon_dump()
{
};

//***********************************************
//Enable output ntuple branches
//***********************************************

//------ Load the miniAOD objects and reset tree variables for each event ------//
void displacedJetMuon_dump::loadEvent(const edm::Event& iEvent)//load all miniAOD objects for the current event
{
  iEvent.getByToken(triggerBitsToken_, triggerBits);
  iEvent.getByToken(hepMCToken_, hepMC);
  iEvent.getByToken(triggerBitsToken_, triggerBits);
  iEvent.getByToken(metFilterBitsToken_, metFilterBits);
  iEvent.getByToken(verticesToken_, vertices);
  iEvent.getByToken(cscSegmentInputToken_,cscSegments);
  iEvent.getByToken(dtSegmentInputToken_,dtSegments);
  iEvent.getByToken(rpcRecHitInputToken_,rpcRecHits);
  iEvent.getByToken(tracksTag_,tracks);
  iEvent.getByToken(PFCandsToken_, pfCands);
  iEvent.getByToken(PFClustersToken_, pfClusters);
  iEvent.getByToken(muonsToken_, muons);
  iEvent.getByToken(electronsToken_, electrons);
  iEvent.getByToken(photonsToken_, photons);
  iEvent.getByToken(tausToken_, taus);
  iEvent.getByToken(jetsCaloToken_, jetsCalo);
  // iEvent.getByToken(jetsPFToken_, jetsPF);
  iEvent.getByToken(jetsToken_, jets);
  iEvent.getByToken(jetsPuppiToken_, jetsPuppi);
  iEvent.getByToken(jetsAK8Token_, jetsAK8);
  iEvent.getByToken(genMetCaloToken_, genMetsCalo);
  iEvent.getByToken(genMetTrueToken_, genMetsTrue);
  iEvent.getByToken(metToken_, mets);
  //iEvent.getByToken(metNoHFToken_, metsNoHF);
  iEvent.getByToken(metPuppiToken_, metsPuppi);
//  iEvent.getByToken(hcalNoiseInfoToken_,hcalNoiseInfo);
  iEvent.getByToken(secondaryVerticesToken_,secondaryVertices);
  iEvent.getByToken(rhoAllToken_,rhoAll);
  iEvent.getByToken(rhoFastjetAllToken_,rhoFastjetAll);
  iEvent.getByToken(rhoFastjetAllCaloToken_,rhoFastjetAllCalo);
  iEvent.getByToken(rhoFastjetCentralCaloToken_,rhoFastjetCentralCalo);
  iEvent.getByToken(rhoFastjetCentralChargedPileUpToken_,rhoFastjetCentralChargedPileUp);
  iEvent.getByToken(rhoFastjetCentralNeutralToken_,rhoFastjetCentralNeutral);
  iEvent.getByToken(beamSpotToken_,beamSpot);
  iEvent.getByToken(ebRecHitsToken_,ebRecHits);
  iEvent.getByToken(eeRecHitsToken_,eeRecHits);
  iEvent.getByToken(esRecHitsToken_,esRecHits);
  iEvent.getByToken(ebeeClustersToken_,ebeeClusters);
  iEvent.getByToken(esClustersToken_,esClusters);
  iEvent.getByToken(conversionsToken_,conversions);
  iEvent.getByToken(singleLegConversionsToken_,singleLegConversions);
  iEvent.getByToken(gedGsfElectronCoresToken_,gedGsfElectronCores);
  iEvent.getByToken(gedPhotonCoresToken_, gedPhotonCores);
  iEvent.getByToken(generalTrackToken_,generalTracks);
//  iEvent.getByToken(superClustersToken_,superClusters);
//  iEvent.getByToken(lostTracksToken_,lostTracks);
//  iEvent.getByToken(hbheNoiseFilterToken_, hbheNoiseFilter);
//  iEvent.getByToken(hbheTightNoiseFilterToken_, hbheTightNoiseFilter);
//  iEvent.getByToken(hbheIsoNoiseFilterToken_, hbheIsoNoiseFilter);
  //iEvent.getByToken(badChargedCandidateFilterToken_, badChargedCandidateFilter);
  //iEvent.getByToken(badMuonFilterToken_, badMuonFilter);
  if(readGenVertexTime_) iEvent.getByToken(genParticles_t0_Token_,genParticles_t0);
  if (useGen_) {
//    iEvent.getByToken(genParticlesToken_,genParticles);
    iEvent.getByToken(genParticlesToken_,genParticles);
    iEvent.getByToken(genJetsToken_,genJets);

    //for Spring16 fastsim, this has been changed and removed
//    if (!isFastsim_) iEvent.getByToken(lheInfoToken_, lheInfo);

    iEvent.getByToken(genInfoToken_,genInfo);
    iEvent.getByToken(puInfoToken_,puInfo);
  }


}



void displacedJetMuon_dump::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {


}


//------ Method called for each lumi block ------//
void displacedJetMuon_dump::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) {

} 


//**********************************************************
// // Dump Version to be used later for clustered objects
//**********************************************************
void displacedJetMuon_dump::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;
  using namespace fastjet;

  //initialize
  loadEvent(iEvent); //loads objects and resets tree branches



  edm::ESHandle<CSCGeometry> cscG;
  edm::ESHandle<DTGeometry> dtG;
  edm::ESHandle<RPCGeometry> rpcG;

  iSetup.get<MuonGeometryRecord>().get(cscG);
  iSetup.get<MuonGeometryRecord>().get(dtG);
  iSetup.get<MuonGeometryRecord>().get(rpcG);

  
  vector< const CSCSegment*> myCSCSegments;
  vector<PseudoJet> mySegments;
  vector<vector<float> > trackSegments;
  int nCscSegments = 0;
  for (const CSCSegment cscSegment : *cscSegments) {

    myCSCSegments.push_back( &cscSegment );

    float globX = 0.;
    float globY = 0.;
    float globZ = 0.;
        

    CSCDetId id  = (CSCDetId)(cscSegment).cscDetId();
    LocalPoint segPos = (cscSegment).localPosition();
    LocalVector segDirection = (cscSegment).localDirection();
    const CSCChamber* cscchamber = cscG->chamber(id);
    if (cscchamber) {
      GlobalPoint globalPosition = cscchamber->toGlobal(segPos);
      GlobalVector globalDirection = cscchamber->toGlobal(segDirection);
      globX = globalPosition.x();
      globY = globalPosition.y();
      globZ = globalPosition.z();  

      cout << "CSC Segment: " << globX << " " << globY << " " << globZ 
	   << " ( " << globalPosition.eta() << " , " << globalPosition.phi() << ") : "
	   << cscSegment.nRecHits() << " " << cscSegment.chi2() << " " 
	   << " : " 
	   << globalDirection.x() << " " << globalDirection.y() << " " << globalDirection.z() 
	   << " ( " << globalDirection.eta() << " , " << globalDirection.phi() << " ) " 
	   << " : "
	   <<  cscSegment.time() << " " 
	//<< cscSegment.parameters()[0] << " " << cscSegment.parameters()[1] << " " << cscSegment.parameters()[2] << " " << cscSegment.parameters()[3] << " | " 
	//<< cscSegment.parametersError()
	   << "\n";
      
      vector<float> trackSegment;
      trackSegment.push_back( globX );
      trackSegment.push_back( globY );
      trackSegment.push_back( globZ );
      trackSegment.push_back( globX + globalDirection.x());
      trackSegment.push_back( globY + globalDirection.y());
      trackSegment.push_back( globZ + globalDirection.z());
      trackSegments.push_back(trackSegment);
      nCscSegments++;

      // look at flags of the rechits
      const std::vector<CSCRecHit2D> cscrechits2d = cscSegment.specificRecHits();
      int cscNRecHits_flagged = 0;
      for (const CSCRecHit2D recHit2d : cscrechits2d) {
	if (!(recHit2d.quality()==1)) continue;
	if(recHit2d.badStrip()) continue;
	if (recHit2d.badWireGroup()) continue;
	
	// std::cout<<cscSegment.nRecHits()<<", " << recHit2d.quality()<<", "<<recHit2d.badStrip()<<", "<<recHit2d.badWireGroup()<<", "<<recHit2d.errorWithinStrip()<<", "<<recHit2d.energyDepositedInLayer()<<std::endl;
	
	//cout << "---> Rechit : " << recHit2d	
      }     
    }
    mySegments.push_back( PseudoJet ( globX, globY, globZ, 0 ) );
  }



  //Check for distance of closest approach
  cout << "nCscSegments = " << nCscSegments << "\n";
  cout << "\n\n\n Check for closest approach\n";
  cout << "Number of Segments: " << trackSegments.size() << "\n";
  // for (uint i=0; i<trackSegments.size() ;i++) {
  //   cout << "Track Segment " << i << " : " 
  // 	 << trackSegments[i][0] << "," << trackSegments[i][1] << "," << trackSegments[i][2] << " | " 
  // 	 << trackSegments[i][3] - trackSegments[i][0]  << "," << trackSegments[i][4] - trackSegments[i][1]  << "," << trackSegments[i][5] -trackSegments[i][2]  
  // 	 << "\n";
  //   for (uint j=0; j<trackSegments.size() ;j++) {
  //     if (i==j) continue;
  //     cout << "   Closest approach to Track Segment " << j << " : " 
  // 	   << trackSegments[j][0] << "," << trackSegments[j][1] << "," << trackSegments[j][2] << " | " 
  // 	   << trackSegments[j][3] - trackSegments[j][0]  << "," << trackSegments[j][4] - trackSegments[j][1]  << "," << trackSegments[j][5] -trackSegments[j][2]  
  // 	   << " | ";
	
  //     //u is the direction of line i
  //     double u_x = trackSegments[i][3] - trackSegments[i][0];
  //     double u_y = trackSegments[i][4] - trackSegments[i][1];
  //     double u_z = trackSegments[i][5] - trackSegments[i][2];
  //     //v is the direction of line j
  //     double v_x = trackSegments[j][3] - trackSegments[j][0];
  //     double v_y = trackSegments[j][4] - trackSegments[j][1];
  //     double v_z = trackSegments[j][5] - trackSegments[j][2];
  //     //w is the vector between the anchor points on the two lines
  //     double w_x = trackSegments[i][0] - trackSegments[j][0];
  //     double w_y = trackSegments[i][1] - trackSegments[j][1];
  //     double w_z = trackSegments[i][2] - trackSegments[j][2];

  //     double a = u_x*u_x + u_y*u_y + u_z*u_z;         // u.u , always >= 0	
  //     double b = u_x*v_x + u_y*v_y + u_z*v_z;         // u.v
  //     double c = v_x*v_x + v_y*v_y + v_z*u_z;         // v.v , always >= 0
  //     double d = u_x*w_x + u_y*w_y + u_z*w_z;         // u.w
  //     double e = v_x*w_x + v_y*w_y + v_z*w_z;         // v.w
  //     double D = a*c - b*b;        // always >= 0
  //     double sc, tc;

  //     // compute the line parameters of the two closest points
  //     if (D < 0.000001) {          // the lines are almost parallel
  // 	sc = 0.0;
  // 	tc = (b>c ? d/b : e/c);    // use the largest denominator
  //     }
  //     else {
  // 	sc = (b*e - c*d) / D;
  // 	tc = (a*e - b*d) / D;
  //     }
	
  //     // get the difference of the two closest points : vector DP = L1(sc) - L2(tc)
  //     double dP_x = w_x + (sc * u_x) - (tc * v_x);
  //     double dP_y = w_y + (sc * u_y) - (tc * v_y);
  //     double dP_z = w_z + (sc * u_z) - (tc * v_z);
  //     double closestApproach = sqrt(dP_x*dP_x + dP_y*dP_y + dP_z*dP_z); //norm(dP)

  //     double midPoint_x = 0.5* (trackSegments[i][0] + (sc * u_x) + trackSegments[j][0] + (tc * v_x));
  //     double midPoint_y = 0.5* (trackSegments[i][1] + (sc * u_y) + trackSegments[j][1] + (tc * v_y));
  //     double midPoint_z = 0.5* (trackSegments[i][2] + (sc * u_z) + trackSegments[j][2] + (tc * v_z));

  //     cout << closestApproach << " | ";
  //     cout << midPoint_x << " , " << midPoint_y << " , " << midPoint_z << " ";
  //     cout << "\n";	
  //   }
  //   cout << " done \n";     
  // }

    cout << "\n******\n\n\n\n";



   
  // // choose a jet definition
  // double R = 0.1;
  // JetDefinition jet_def(cambridge_algorithm, R);
  
  // // run the clustering, extract the jets
  // ClusterSequence cs(mySegments, jet_def);
  // vector<PseudoJet> mySegmentClusters = sorted_by_pt(cs.inclusive_jets());
  
  // // print out some infos
  // cout << "Clustering with " << jet_def.description() << "\n";
  
  // cout <<   "        pt y phi" << "\n";
  // for (unsigned i = 0; i < mySegmentClusters.size(); i++) {
  //   cout << "mySegmentCluster " << i << ": "<< mySegmentClusters[i].px() << " " 
  // 	 << mySegmentClusters[i].py() << " " << mySegmentClusters[i].pz() << "\n";
  //   vector<PseudoJet> mySegmentClusterConstituents = mySegmentClusters[i].constituents();
  //   for (unsigned j = 0; j < mySegmentClusterConstituents.size(); j++) {
  //     cout << "    constituent segment" << j << ": " << mySegmentClusterConstituents[j].px() << " " << mySegmentClusterConstituents[j].py() << " " << mySegmentClusterConstituents[j].pz()
  // 	   << "\n";
  //   }
  // }

};
//**********************************************************
//**********************************************************


//------ Method called once each job just before starting event loop ------//
void displacedJetMuon_dump::beginJob()
{

}

//------ Method called once each job just after ending the event loop ------//
void displacedJetMuon_dump::endJob(){};





//define this as a plug-in
DEFINE_FWK_MODULE(displacedJetMuon_dump);
