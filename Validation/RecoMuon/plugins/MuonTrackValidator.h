#ifndef MuonTrackValidator_h
#define MuonTrackValidator_h

/** \class MuonTrackValidator
* Class that produces histograms to validate Muon Track Reconstruction performances
*
*/
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "Validation/RecoMuon/plugins/MuonTrackValidatorBase.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"

#include "DataFormats/HeavyIonEvent/interface/Centrality.h"

class MuonTrackValidator : public DQMEDAnalyzer, protected MuonTrackValidatorBase {
 public:
  /// Constructor
  MuonTrackValidator(const edm::ParameterSet& pset):MuonTrackValidatorBase(pset){
    dirName_ = pset.getParameter<std::string>("dirName");
    associatormap = pset.getParameter< edm::InputTag >("associatormap");
    UseAssociators = pset.getParameter< bool >("UseAssociators");
    tpSelector = TrackingParticleSelector(pset.getParameter<double>("ptMinTP"),
					  pset.getParameter<double>("minRapidityTP"),
					  pset.getParameter<double>("maxRapidityTP"),
					  pset.getParameter<double>("tipTP"),
					  pset.getParameter<double>("lipTP"),
					  pset.getParameter<int>("minHitTP"),
					  pset.getParameter<bool>("signalOnlyTP"),
					  pset.getParameter<bool>("chargedOnlyTP"),
					  pset.getParameter<bool>("stableOnlyTP"),
					  pset.getParameter<std::vector<int> >("pdgIdTP"));
    cosmictpSelector = CosmicTrackingParticleSelector(pset.getParameter<double>("ptMinTP"),
						      pset.getParameter<double>("minRapidityTP"),
						      pset.getParameter<double>("maxRapidityTP"),
						      pset.getParameter<double>("tipTP"),
						      pset.getParameter<double>("lipTP"),
						      pset.getParameter<int>("minHitTP"),
						      pset.getParameter<bool>("chargedOnlyTP"),
						      pset.getParameter<std::vector<int> >("pdgIdTP"));
    
    minPhi = pset.getParameter<double>("minPhi");
    maxPhi = pset.getParameter<double>("maxPhi");
    nintPhi = pset.getParameter<int>("nintPhi");
    useGsf = pset.getParameter<bool>("useGsf");
    BiDirectional_RecoToSim_association = pset.getParameter<bool>("BiDirectional_RecoToSim_association");

    //// Weighting numbers for different pT bin MC samples in PbPb
    pTWeight = pset.getParameter<double>("pTWeight");
    _centralityBinTag=consumes<int>(pset.getParameter<edm::InputTag>("CentralityBinSrc"));

    _centralityRanges = pset.getParameter< std::vector<double> >("centralityRanges");
    nintcent = _centralityRanges.size() - 1;
    centralityRanges = new float[nintcent+1];
    for (int idx=0; idx<=nintcent; idx++) centralityRanges[idx] = _centralityRanges[idx]/2.5;

    _pTRanges = pset.getParameter< std::vector<double> >("pTRanges");
    if (!_pTRanges.empty()) {
      nintpT = _pTRanges.size() - 1;
      std::cout << "nintpT is changed to: " << nintpT << std::endl;
      pTRanges = new float[nintpT+1];
      for (int idx=0; idx<=nintpT; idx++) pTRanges[idx] = _pTRanges[idx];
    } else {
      pTRanges = new float[1];
      pTRanges[0] = -1;
    }

    _rapArr = pset.getParameter< std::vector<double> >("rapArr");
    _rapArrRes = pset.getParameter< std::vector<double> >("rapArrRes");
    _ptArrRes = pset.getParameter< std::vector<double> >("ptArrRes");
    nintRapArr = _rapArr.size() -1;
    nintRapArrRes = _rapArrRes.size() -1;
    nintPtArrRes = _ptArrRes.size() -1;
    rapArr = new float[nintRapArr+1];
    rapArrRes = new float[nintRapArrRes+1];
    ptArrRes =  new float[nintPtArrRes+1];
    for (int idx=0; idx<=nintRapArr; idx++) rapArr[idx] = _rapArr[idx];
    for (int idx=0; idx<=nintRapArrRes; idx++) rapArrRes[idx] = _rapArrRes[idx];
    for (int idx=0; idx<=nintPtArrRes; idx++) ptArrRes[idx] = _ptArrRes[idx];

    // dump cfg parameters
    edm::LogVerbatim("MuonTrackValidator") << "constructing MuonTrackValidator: " << pset.dump();
    
    // Declare consumes (also for the base class)
    bsSrc_Token = consumes<reco::BeamSpot>(bsSrc);
    tp_effic_Token = consumes<TrackingParticleCollection>(label_tp_effic);
    tp_fake_Token = consumes<TrackingParticleCollection>(label_tp_fake);
    for (unsigned int www=0;www<label.size();www++){
      track_Collection_Token.push_back(consumes<edm::View<reco::Track> >(label[www]));
    }
    simToRecoCollection_Token = consumes<reco::SimToRecoCollection>(associatormap);
    recoToSimCollection_Token = consumes<reco::RecoToSimCollection>(associatormap);

    _simHitTpMapTag = mayConsume<SimHitTPAssociationProducer::SimHitTPAssociationList>(pset.getParameter<edm::InputTag>("simHitTpMapTag"));

    MABH = false;
    if (!UseAssociators) {
      // flag MuonAssociatorByHits
      if (associators[0] == "MuonAssociationByHits") MABH = true;
      // reset string associators to the map label
      associators.clear();
      associators.push_back(associatormap.label());
      edm::LogVerbatim("MuonTrackValidator") << "--> associators reset to: " <<associators[0];
    } else {
      for (auto const& associator :associators) {
        consumes<reco::TrackToTrackingParticleAssociator>(edm::InputTag(associator));
      }
    }
    
    // inform on which SimHits will be counted
    if (usetracker) edm::LogVerbatim("MuonTrackValidator")
      <<"\n usetracker = TRUE : Tracker SimHits WILL be counted";
    else edm::LogVerbatim("MuonTrackValidator")
      <<"\n usetracker = FALSE : Tracker SimHits WILL NOT be counted";
    if (usemuon) edm::LogVerbatim("MuonTrackValidator")
      <<" usemuon = TRUE : Muon SimHits WILL be counted";
    else edm::LogVerbatim("MuonTrackValidator")
      <<" usemuon = FALSE : Muon SimHits WILL NOT be counted"<<std::endl;
    
    // loop over the reco::Track collections to validate: check for inconsistent input settings
    for (unsigned int www=0;www<label.size();www++) {
      std::string recoTracksLabel = label[www].label();
      std::string recoTracksInstance = label[www].instance();
      
      // tracks with hits only on tracker
      if (recoTracksLabel=="generalTracks" ||
	  (recoTracksLabel.find("cutsRecoTracks") != std::string::npos) ||
	  recoTracksLabel=="ctfWithMaterialTracksP5LHCNavigation" ||
	  recoTracksLabel=="hltL3TkTracksFromL2" ||
	  (recoTracksLabel=="hltL3Muons" && recoTracksInstance=="L2Seeded"))
	{
	  if (usemuon) {
	    edm::LogWarning("MuonTrackValidator")
	      <<"\n*** WARNING : inconsistent input tracksTag = "<<label[www]
	      <<"\n with usemuon == true"<<"\n ---> please change to usemuon == false ";
	  }
	  if (!usetracker) {
	    edm::LogWarning("MuonTrackValidator")
	      <<"\n*** WARNING : inconsistent input tracksTag = "<<label[www]
	      <<"\n with usetracker == false"<<"\n ---> please change to usetracker == true ";
	  }	
	}
      
      // tracks with hits only on muon detectors
      else if (recoTracksLabel=="standAloneMuons" ||
	       recoTracksLabel=="standAloneSETMuons" ||
	       recoTracksLabel=="cosmicMuons" ||
	       recoTracksLabel=="hltL2Muons")
	{
	  if (usetracker) {
	    edm::LogWarning("MuonTrackValidator")
	      <<"\n*** WARNING : inconsistent input tracksTag = "<<label[www]
	      <<"\n with usetracker == true"<<"\n ---> please change to usetracker == false ";
	  }
	  if (!usemuon) {
	    edm::LogWarning("MuonTrackValidator")
	      <<"\n*** WARNING : inconsistent input tracksTag = "<<label[www]
	      <<"\n with usemuon == false"<<"\n ---> please change to usemuon == true ";
	  }
	}
      
    } // for (unsigned int www=0;www<label.size();www++)
  } // end of constructor
  
  /// Destructor
  virtual ~MuonTrackValidator(){
    delete[] centralityRanges;
    delete[] pTRanges;
    delete[] rapArr;
    delete[] rapArrRes;
    delete[] ptArrRes;
  }

  /// Method called before the event loop
  //  void beginRun(edm::Run const&, edm::EventSetup const&);
  /// Method called once per event
  void analyze(const edm::Event&, const edm::EventSetup& );
  /// Method called at the end of the event loop
  void endRun(edm::Run const&, edm::EventSetup const&);
  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&);

private:
  /// retrieval of reconstructed momentum components from reco::Track (== mean values for GSF)
  void getRecoMomentum (const reco::Track& track, double& pt, double& ptError,
			double& qoverp, double& qoverpError, double& lambda, double& lambdaError,
			double& phi, double& phiError ) const;
  /// retrieval of reconstructed momentum components based on the mode of a reco::GsfTrack
  void getRecoMomentum (const reco::GsfTrack& gsfTrack, double& pt, double& ptError,
			double& qoverp, double& qoverpError, double& lambda, double& lambdaError,
			double& phi, double& phiError) const;
  Double_t findNcoll(int hiBin);
  void setUpAsymmVectors();

 private:
  std::string dirName_;
  edm::InputTag associatormap;
  edm::EDGetTokenT<reco::SimToRecoCollection> simToRecoCollection_Token;
  edm::EDGetTokenT<reco::RecoToSimCollection> recoToSimCollection_Token;
  edm::EDGetTokenT<SimHitTPAssociationProducer::SimHitTPAssociationList> _simHitTpMapTag;

  bool UseAssociators;
  double minPhi, maxPhi;
  int nintPhi;
  bool useGsf;
  // select tracking particles
  //(i.e. "denominator" of the efficiency ratio)
  TrackingParticleSelector tpSelector;	
  CosmicTrackingParticleSelector cosmictpSelector;

  // flag new validation logic (bidirectional RecoToSim association)
  bool BiDirectional_RecoToSim_association;
  // flag MuonAssociatorByHits
  bool MABH;
  
  //1D
  std::vector<MonitorElement*> h_nchi2, h_nchi2_prob, h_losthits;

  //2D
  std::vector<MonitorElement*> chi2_vs_nhits, etares_vs_eta;
  std::vector<MonitorElement*> h_ptshifteta;
  std::vector<MonitorElement*> ptres_vs_phi, chi2_vs_phi, nhits_vs_phi, phires_vs_phi;

  //Profile2D
  std::vector<MonitorElement*> ptmean_vs_eta_phi, phimean_vs_eta_phi;

  //assoc chi2
  std::vector<MonitorElement*> h_assochi2, h_assochi2_prob;

  //chi2 and # lost hits vs eta: to be used with doProfileX
  std::vector<MonitorElement*> chi2_vs_eta, nlosthits_vs_eta;
  std::vector<MonitorElement*> h_chi2meanh, h_losthits_eta;
  std::vector<MonitorElement*> h_hits_phi;
  std::vector<MonitorElement*> h_chi2meanhitsh, h_chi2mean_vs_phi;

  //resolution of track params: to be used with fitslicesytool
  std::vector<MonitorElement*> dxyres_vs_eta, ptres_vs_eta, dzres_vs_eta, phires_vs_eta, cotThetares_vs_eta;
  std::vector<MonitorElement*> dxyres_vs_pt, ptres_vs_pt, dzres_vs_pt, phires_vs_pt, cotThetares_vs_pt;

  //pulls of track params vs eta: to be used with fitslicesytool
  std::vector<MonitorElement*> dxypull_vs_eta, ptpull_vs_eta, dzpull_vs_eta, phipull_vs_eta, thetapull_vs_eta;
  std::vector<MonitorElement*> ptpull_vs_phi, phipull_vs_phi, thetapull_vs_phi;
  std::vector<MonitorElement*> h_dxypulleta, h_ptpulleta, h_dzpulleta, h_phipulleta, h_thetapulleta;
  std::vector<MonitorElement*> h_ptpullphi, h_phipullphi, h_thetapullphi;

  // Weighting for different pT bins of MC samples (in pbpb)
  double pTWeight;

  // Centrality related
  edm::EDGetTokenT<int> _centralityBinTag;
  int nintcent, centBin;

  std::vector<double> _centralityRanges;
  std::vector<double> _pTRanges;
  float *centralityRanges;
  float *pTRanges;

  std::vector<double> _rapArr;
  std::vector<double> _rapArrRes;
  std::vector<double> _ptArrRes;
  float *rapArr;
  float *rapArrRes;
  float *ptArrRes;
  int nintRapArr;
  int nintRapArrRes;
  int nintPtArrRes;

};


#endif
