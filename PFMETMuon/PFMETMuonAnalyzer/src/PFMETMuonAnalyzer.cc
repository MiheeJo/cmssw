// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// ana
#include "PFMETMuon/PFMETMuonAnalyzer/interface/PFMETMuonAnalyzer.h"


using namespace std;
using namespace edm;
using namespace reco;
using namespace pat;

//
// constructors and destructor
//
PFMETMuonAnalyzer::PFMETMuonAnalyzer(const edm::ParameterSet& iConfig):
  hltPrescaleProvider(iConfig, consumesCollector(), *this)
{
  // Event source
  // Event Info
  centralityTagToken_ = consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("CentralitySrc"));
  centralityBinTagToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("CentralityBinSrc"));

  pfCandidateLabel_ = iConfig.getParameter<edm::InputTag>("pfCandidateLabel");
  pfCandidatePF_ = consumes<reco::PFCandidateCollection>(pfCandidateLabel_);
  pfCandidateView_ = consumes<reco::CandidateView>(pfCandidateLabel_);
  pfPtMin_ = iConfig.getParameter<double>("pfPtMin");
  genPtMin_ = iConfig.getParameter<double>("genPtMin");
  jetPtMin_ = iConfig.getParameter<double>("jetPtMin");

  etaBins_ = iConfig.getParameter<int>("etaBins");
  fourierOrder_ = iConfig.getParameter<int>("fourierOrder");

  doVS_ = iConfig.getUntrackedParameter<bool>("doVS",false);
  if(doVS_){
    edm::InputTag vsTag = iConfig.getParameter<edm::InputTag>("bkg");
    srcVorFloat_ = consumes<std::vector<float> >(vsTag);
    srcVorMap_ = consumes<reco::VoronoiMap>(vsTag);
  }

  // debug
  verbosity_ = iConfig.getUntrackedParameter<int>("verbosity", false);

  doJets_ = iConfig.getUntrackedParameter<bool>("doJets",false);
  if(doJets_){
    jetLabel_ = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetLabel"));
  }
  doUEraw_ = iConfig.getUntrackedParameter<bool>("doUEraw",false);

  doMC_ = iConfig.getUntrackedParameter<bool>("doMC",false);
  isHI_ = iConfig.getUntrackedParameter<bool>("isHI",false);
  isPA_ = iConfig.getUntrackedParameter<bool>("isPA",true);

  if(doMC_){
    genLabel_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genLabel"));
  }
  skipCharged_ = iConfig.getUntrackedParameter<bool>("skipCharged",false);
  qualityString_ = iConfig.getParameter<std::string>("trackQuality");

  pfMETLabel_ = consumes<reco::PFMETCollection>(iConfig.getParameter<edm::InputTag>("pfMETLabel"));
  trackLabel_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("trackLabel"));
  muonLabel_ = consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonLabel"));
  vtxLabel_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vtxLabel"));

  // Trigger information
  triggerResultsLabel_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResultsLabel"));
  theTriggerNames = iConfig.getParameter< std::vector<string> >("triggerPathNames");
  HLTLastFilters = iConfig.getParameter< std::vector<string> >("triggerFilterNames");

}


PFMETMuonAnalyzer::~PFMETMuonAnalyzer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
PFMETMuonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  pfEvt_.Clear();

  // Initialize trigger and prescale info for each event
  for(std::map< std::string, int >::iterator clearIt= mapTriggerNameToIntFired_.begin(); clearIt != mapTriggerNameToIntFired_.end(); clearIt++){
    clearIt->second=0;
  }
  for(std::map< std::string, int >::iterator clearIt= mapTriggerNameToPrescaleFac_.begin(); clearIt != mapTriggerNameToPrescaleFac_.end(); clearIt++){
    clearIt->second=-1;
  }

  // Trigger information filled up
  hltReport(iEvent, iSetup);

  pfEvt_.runNb = iEvent.id().run();
  pfEvt_.eventNb = iEvent.id().event();
  pfEvt_.lumiSection = iEvent.luminosityBlock();

  // Centrality
  edm::Handle<reco::Centrality> centrality;
  edm::Handle<int> cbin_;
  if (isHI_ || isPA_)  {
    iEvent.getByToken(centralityTagToken_, centrality); 
    iEvent.getByToken(centralityBinTagToken_, cbin_);
  }
  if (centrality.isValid() && cbin_.isValid()) {
    pfEvt_.CentBin           = *cbin_; // Bin number of centrality (0-200 for HI)

    pfEvt_.Npix              = centrality->multiplicityPixel();
    pfEvt_.NpixelTracks      = centrality->NpixelTracks();
    pfEvt_.Ntracks           = centrality->Ntracks();
    pfEvt_.NtracksPtCut      = centrality->NtracksPtCut();
    pfEvt_.NtracksEtaCut     = centrality->NtracksEtaCut();
    pfEvt_.NtracksEtaPtCut   = centrality->NtracksEtaPtCut();

    pfEvt_.SumET_HF          = centrality->EtHFtowerSum();
    pfEvt_.SumET_HFplus      = centrality->EtHFtowerSumPlus();
    pfEvt_.SumET_HFminus     = centrality->EtHFtowerSumMinus();
    pfEvt_.SumET_HFplusEta4  = centrality->EtHFtruncatedPlus();
    pfEvt_.SumET_HFminusEta4 = centrality->EtHFtruncatedMinus();

    pfEvt_.SumET_HFhit       = centrality->EtHFhitSum(); 
    pfEvt_.SumET_HFhitPlus   = centrality->EtHFhitSumPlus();
    pfEvt_.SumET_HFhitMinus  = centrality->EtHFhitSumMinus();

    pfEvt_.SumET_ZDC         = centrality->zdcSum();
    pfEvt_.SumET_ZDCplus     = centrality->zdcSumPlus();
    pfEvt_.SumET_ZDCminus    = centrality->zdcSumMinus();

    pfEvt_.SumET_EEplus      = centrality->EtEESumPlus();
    pfEvt_.SumET_EEminus     = centrality->EtEESumMinus();
    pfEvt_.SumET_EE          = centrality->EtEESum();
    pfEvt_.SumET_EB          = centrality->EtEBSum();
    pfEvt_.SumET_ET          = centrality->EtMidRapiditySum();
  }

  edm::Handle<reco::VertexCollection> privtxs;
  iEvent.getByToken(vtxLabel_, privtxs);

  reco::VertexCollection::const_iterator privtx;

  if ( privtxs->begin() != privtxs->end() ) {
    privtx = privtxs->begin();
    pfEvt_.RefVtx = privtx->position();
    pfEvt_.RefVtx_x = pfEvt_.RefVtx.X();
    pfEvt_.RefVtx_y = pfEvt_.RefVtx.Y();
    pfEvt_.RefVtx_z = pfEvt_.RefVtx.Z();
    pfEvt_.RefVtx_xError = privtx->xError();
    pfEvt_.RefVtx_yError = privtx->yError();
    pfEvt_.RefVtx_zError = privtx->zError();
  } else {
    pfEvt_.RefVtx.SetXYZ(0.,0.,0.);
    pfEvt_.RefVtx_x = 0;
    pfEvt_.RefVtx_y = 0;
    pfEvt_.RefVtx_z = 0;
    pfEvt_.RefVtx_xError = 0.0;
    pfEvt_.RefVtx_yError = 0.0;
    pfEvt_.RefVtx_zError = 0.0;
  }
  pfEvt_.nPV = privtxs->size();

  // Fill PF info
  edm::Handle<reco::PFCandidateCollection> pfCandidates;
  iEvent.getByToken(pfCandidatePF_,pfCandidates);
  iEvent.getByToken(pfCandidateView_,candidates_);
  const reco::PFCandidateCollection *pfCandidateColl = pfCandidates.product();
  if (doVS_) {
   iEvent.getByToken(srcVorMap_,backgrounds_);
   iEvent.getByToken(srcVorFloat_,vn_);
   UEParameters vnUE(vn_.product(),fourierOrder_,etaBins_);
   const std::vector<float>& vue = vnUE.get_raw();

   for(int ieta = 0; ieta < etaBins_; ++ieta){
     pfEvt_.sumpt[ieta] = vnUE.get_sum_pt(ieta);
     for(int ifour = 0; ifour < fourierOrder_; ++ifour){
       pfEvt_.vn[ifour * etaBins_ + ieta] = vnUE.get_vn(ifour,ieta);
       pfEvt_.psin[ifour * etaBins_ + ieta] = vnUE.get_psin(ifour,ieta);
     }
   }

   for(int iue = 0; iue < etaBins_*fourierOrder_*2*3; ++iue){
     pfEvt_.ueraw[iue] = vue[iue];
   }
  }

  for(unsigned icand=0;icand<pfCandidateColl->size(); icand++) {
    const reco::PFCandidate pfCandidate = pfCandidateColl->at(icand);
    reco::CandidateViewRef ref(candidates_,icand);

    double vsPtInitial=-999, vsPt=-999, vsArea = -999;

    if (doVS_) {
      const reco::VoronoiBackground& voronoi = (*backgrounds_)[ref];
      vsPt = voronoi.pt();
      vsPtInitial = voronoi.pt_subtracted();
      vsArea = voronoi.area();
    }

    double pt =  pfCandidate.pt();
    double energy = pfCandidate.energy();
    if (pt<=pfPtMin_) continue;

    int id = pfCandidate.particleId();
    if (skipCharged_ && (abs(id) == 1 || abs(id) == 3)) continue;

    bool matched2Jet = false;
    if (doJets_) {
      pfEvt_.jetMatchIndex_.push_back( -1 ); // default index is -1
      
      edm::Handle<pat::JetCollection> jets;
      iEvent.getByToken(jetLabel_,jets);
      const pat::JetCollection *jetColl = &(*jets);

      for (unsigned ijet=0;ijet<jetColl->size(); ijet++) {
        const pat::Jet jet = jetColl->at(ijet);
                
        if (jet.pt()>jetPtMin_) {
          std::vector<reco::PFCandidatePtr> pfConstituents = jet.getPFConstituents();
          for (std::vector<reco::PFCandidatePtr>::const_iterator ibegin=pfConstituents.begin(), iend=pfConstituents.end(), iconstituent=ibegin; iconstituent!=iend; ++iconstituent) {
            
            reco::PFCandidatePtr candptr(pfCandidates, icand);
            edm::Ptr<reco::PFCandidate> pfBackRef ( *iconstituent );
              
            // couldn't figure out the matching by ref, so just do it like this:
            if (pfBackRef->pt() == pfCandidate.pt() && pfBackRef->eta()== pfCandidate.eta() && pfBackRef->particleId()== pfCandidate.particleId()) {
              
              pfEvt_.jetMatchIndex_[pfEvt_.nPFpart_] = ijet; // when a match is found, change -1 to matched jet index
              matched2Jet=true;
              break;
            }
          }
        }
        if(matched2Jet==true) break;
      }
    } // end of if(doJets)

    pfEvt_.pfId_.push_back( id );
    pfEvt_.pfPt_.push_back( rndSF(pt,4) );
    pfEvt_.pfEnergy_.push_back( rndSF(energy,4) );
    pfEvt_.pfVsPt_.push_back( rndSF(vsPt,4) );
    pfEvt_.pfVsPtInitial_.push_back( rndSF(vsPtInitial,4) );
    pfEvt_.pfArea_.push_back( rndSF(vsArea,4) );
    pfEvt_.pfEta_.push_back( rndDP(pfCandidate.eta(),3) );
    pfEvt_.pfPhi_.push_back( rndDP(pfCandidate.phi(),3) );
    pfEvt_.pfTheta_.push_back( rndDP(pfCandidate.theta(),3) ); 
    pfEvt_.pfEt_.push_back( rndSF(pfCandidate.et(),4) );
    pfEvt_.pfCharge_.push_back( pfCandidate.charge() ); 

    // More information on charged particle(1) and muon(3)
    Float_t TrackerMuon=0, TrackerMuonPt=0, GlobalMuonPt=0;
    Int_t TrackHits=0;
    Float_t Dxy=0, Dz=0, Chi2=0;
    Float_t MuonPx=0, MuonPy=0, MuonPz=0;
    Float_t ChargedPx=0, ChargedPy=0, ChargedPz=0, ChargedTrackRefPt=0;

    if(abs(id) == 3){ // in case of muons
      MuonPx = pfCandidate.px();
      MuonPy = pfCandidate.py();
      MuonPz = pfCandidate.pz();

      const reco::MuonRef muonRef = pfCandidate.muonRef();  
      if ( muonRef->isTrackerMuon() ){
        reco::TrackRef trackRef = muonRef->track();
        
        if (trackRef.isNonnull()) {
          TrackerMuon = 1;
          TrackerMuonPt = rndSF(trackRef->pt(),4) ;
          TrackHits = trackRef->numberOfValidHits();
          Dxy = rndDP(trackRef->dxy(pfEvt_.RefVtx),4); 
          Dz  = rndDP(trackRef->dz(pfEvt_.RefVtx),4); 
          Chi2 = rndDP(trackRef->normalizedChi2(),4);

          if ( muonRef->isGlobalMuon() ){
            reco::TrackRef globalmuon = muonRef->globalTrack();
            GlobalMuonPt = rndSF(globalmuon->pt(),4);
          }
        }
      }
    }
    
    if(abs(id) == 1){
      ChargedPx = pfCandidate.px();
      ChargedPy = pfCandidate.py();
      ChargedPz = pfCandidate.pz();
      const reco::TrackRef trackRef = pfCandidate.trackRef(); 
      if (trackRef.isNonnull()) {
        ChargedTrackRefPt = rndSF(trackRef->pt(),4);
      }
    }

    pfEvt_.pfMuonPx_.push_back( MuonPx );
    pfEvt_.pfMuonPy_.push_back( MuonPy );
    pfEvt_.pfMuonPz_.push_back( MuonPz );
    pfEvt_.pfTrackerMuon_.push_back( TrackerMuon );
    pfEvt_.pfTrackerMuonPt_.push_back( TrackerMuonPt );
    pfEvt_.pfTrackHits_.push_back( TrackHits );
    pfEvt_.pfDxy_.push_back( Dxy );
    pfEvt_.pfDz_.push_back( Dz );
    pfEvt_.pfChi2_.push_back( Chi2 );
    pfEvt_.pfGlobalMuonPt_.push_back( GlobalMuonPt );
    pfEvt_.pfChargedPx_.push_back( ChargedPx );
    pfEvt_.pfChargedPy_.push_back( ChargedPy );
    pfEvt_.pfChargedPz_.push_back( ChargedPz );
    pfEvt_.pfChargedTrackRefPt_.push_back( ChargedTrackRefPt );

    pfEvt_.nPFpart_++;
  }

  // Fill GEN info
  if(doMC_){
    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByToken(genLabel_,genParticles);
    const reco::GenParticleCollection* genColl= &(*genParticles);

    for(unsigned igen=0;igen<genColl->size(); igen++) {

      const reco::GenParticle gen = genColl->at(igen);
      double eta = gen.eta();
      double pt = gen.pt();

      if(gen.status()==1 && fabs(eta)<3.0 && pt> genPtMin_){
        pfEvt_.genPDGId_.push_back( gen.pdgId() );
        pfEvt_.genPt_.push_back( rndSF(pt,4) );
        pfEvt_.genEta_.push_back( rndDP(eta,3) );
        pfEvt_.genPhi_.push_back( rndDP(gen.phi(),3) );
        pfEvt_.nGENpart_++;
      }
    }
  }

  // Fill Jet info
  if(doJets_){
    edm::Handle<pat::JetCollection> jets;
    iEvent.getByToken(jetLabel_,jets);
    const pat::JetCollection *jetColl = &(*jets);

    for(unsigned ijet=0;ijet<jetColl->size(); ijet++) {
      const pat::Jet jet = jetColl->at(ijet);

      double pt =  jet.pt();
      double energy =  jet.energy();
      if(pt>jetPtMin_){
        pfEvt_.jetPt_.push_back( pt );
        pfEvt_.jetEnergy_.push_back( energy );
        pfEvt_.jetEta_.push_back( jet.eta() );
        pfEvt_.jetPhi_.push_back( jet.phi() );
	pfEvt_.jetMass_.push_back( jet.mass() );
	pfEvt_.jetPU_.push_back( jet.pileup() );
        pfEvt_.njets_++;
      }
    }
  }

  // Fill generalTracks
  edm::Handle<reco::TrackCollection> trackCollection;
  iEvent.getByToken(trackLabel_, trackCollection);
  const reco::TrackCollection* trackColl= &(*trackCollection);
 
  for(unsigned itrack=0;itrack<trackColl->size(); itrack++) { 
    const reco::Track tra = trackColl->at(itrack); 

    if (tra.quality(reco::TrackBase::qualityByName(qualityString_))) pfEvt_.traQual_.push_back(1);
    else pfEvt_.traQual_.push_back(0);
    
    pfEvt_.traCharge_.push_back(tra.charge());
    pfEvt_.traPt_.push_back(rndSF(tra.pt(),4));
    pfEvt_.traEta_.push_back(rndDP(tra.eta(),3));
    pfEvt_.traPhi_.push_back(rndDP(tra.phi(),3));
    pfEvt_.traAlgo_.push_back(tra.algo());
    pfEvt_.traHits_.push_back(tra.numberOfValidHits());
    pfEvt_.nTRACKpart_++;
  }

  // Fill MET Object
  Handle<reco::PFMETCollection>   recoPfMETHandle;
  iEvent.getByToken(pfMETLabel_, recoPfMETHandle);
  const reco::PFMET& pfmet = recoPfMETHandle->at(0);

  if (recoPfMETHandle.isValid()) {
    pfEvt_.recoPfMET_ = rndSF(pfmet.et(),4);
    pfEvt_.recoPfMETPhi_ = rndDP(pfmet.phi(),3);
    pfEvt_.recoPfMETsumEt_  = pfmet.sumEt();
    pfEvt_.recoPfMETmEtSig_ = pfmet.mEtSig();
    pfEvt_.recoPfMETSig_    = pfmet.significance();
  } 


  // Fill single muon information
  Handle<pat::MuonCollection> muonCollection; 
  iEvent.getByToken(muonLabel_, muonCollection);
  for (unsigned imuon=0;imuon<muonCollection->size(); imuon++) { 
    const pat::Muon& muon = muonCollection->at(imuon);

    // Inner track information of a muon
    reco::TrackRef iTrack = muon.innerTrack();
    if (iTrack.isNonnull()) {
      pfEvt_.muSelectionType_.push_back( muonIDmask(muon) );
      pfEvt_.muHighPurity_.push_back( iTrack->quality(reco::TrackBase::highPurity) );
      pfEvt_.muIsGoodMuon_.push_back( muon::isGoodMuon(muon, muon::TMOneStationTight) );
      pfEvt_.muTrkMuArb_.push_back( muon.muonID("TrackerMuonArbitrated") );
      pfEvt_.muTMOneStaTight_.push_back( muon.muonID("TMOneStationTight") );
      pfEvt_.muNTrkHits_.push_back( iTrack->found() );
      pfEvt_.muNormChi2Inner_.push_back( iTrack->normalizedChi2() );
      pfEvt_.muNPixValHits_.push_back( iTrack->hitPattern().numberOfValidPixelHits() );
      pfEvt_.muNPixWMea_.push_back( iTrack->hitPattern().pixelLayersWithMeasurement() );
      pfEvt_.muNTrkWMea_.push_back( iTrack->hitPattern().trackerLayersWithMeasurement() );
      pfEvt_.muStationsMatched_.push_back( muon.numberOfMatchedStations() );
      pfEvt_.muDxy_.push_back( iTrack->dxy(pfEvt_.RefVtx) );
      pfEvt_.muDxyErr_.push_back( iTrack->dxyError() );
      pfEvt_.muDz_.push_back( iTrack->dz(pfEvt_.RefVtx) );
      pfEvt_.muDzErr_.push_back( iTrack->dzError() );
      pfEvt_.muPtInner_.push_back( iTrack->pt() );
      pfEvt_.muPtErrInner_.push_back( iTrack->ptError() );

      Int_t muNMuValHits=0;
      Float_t muNormChi2Global=999, muPtGlobal=-1, muPtErrGlobal=-1;
      
      // Outer track information of a muon
      if (muon.isGlobalMuon()) {
        reco::TrackRef gTrack = muon.globalTrack();
        if (gTrack.isNonnull()) {
          muNMuValHits = gTrack->hitPattern().numberOfValidMuonHits();
          muNormChi2Global = gTrack->normalizedChi2();
          muPtGlobal = gTrack->pt();
          muPtErrGlobal = gTrack->ptError();
        }
      }

      pfEvt_.muNMuValHits_.push_back( muNMuValHits );
      pfEvt_.muNormChi2Global_.push_back( muNormChi2Global );
      pfEvt_.muPtGlobal_.push_back( muPtGlobal );
      pfEvt_.muPtErrGlobal_.push_back( muPtErrGlobal );

      pfEvt_.muPt_.push_back( muon.pt() ); 
      pfEvt_.muPx_.push_back( muon.px() );
      pfEvt_.muPy_.push_back( muon.py() );
      pfEvt_.muPz_.push_back( muon.pz() );     
      pfEvt_.muEta_.push_back( muon.eta() );      
      pfEvt_.muPhi_.push_back( muon.phi() );
      pfEvt_.muCharge_.push_back( muon.charge() );
      pfEvt_.muTrackIso_.push_back( muon.trackIso() );
      pfEvt_.muCaloIso_.push_back( muon.caloIso() );
      pfEvt_.muEcalIso_.push_back( muon.ecalIso() );
      pfEvt_.muHcalIso_.push_back( muon.hcalIso() );
      
      // Muon triggers are matched to muons?
      ULong64_t muTrig=0;
      for (unsigned int iTr = 0; iTr<HLTLastFilters.size(); ++iTr) {
        const pat::TriggerObjectStandAloneCollection mu1HLTMatchesFilter = muon.triggerObjectMatchesByFilter( HLTLastFilters[iTr] );
        const pat::TriggerObjectStandAloneCollection mu1HLTMatchesPath = muon.triggerObjectMatchesByPath( theTriggerNames.at(iTr), true, false );
        
        if (mu1HLTMatchesFilter.size() > 0) {
          muTrig += pow(2,iTr);
        }
      }
      pfEvt_.muTrig_.push_back( muTrig );
      
      pfEvt_.nMUpart_++;
    } // end of iTrack.isNonnull

  }

  // Trigger information
  for (unsigned int iTr = 0 ; iTr < theTriggerNames.size() ; iTr++) {
    if (mapTriggerNameToIntFired_[theTriggerNames.at(iTr)] == 3) {
      pfEvt_.HLTriggers += pow(2,iTr);
    }
    pfEvt_.trigPrescale.push_back(mapTriggerNameToPrescaleFac_[theTriggerNames.at(iTr)]);
  }

  // All done
  pfTree_->Fill();
}

void PFMETMuonAnalyzer::beginJob()
{
  // -- trees --
  pfTree_ = fs->make<TTree>("pfTree","W analysis tree");
  pfEvt_.SetTree(pfTree_);
  pfEvt_.doMC = doMC_;
  pfEvt_.doJets = doJets_;

  pfEvt_.SetBranches(etaBins_, fourierOrder_, doUEraw_);

  // -- init trigger info map --
  for(std::vector<std::string>::iterator it = theTriggerNames.begin(); it != theTriggerNames.end(); ++it){
    mapTriggerNameToIntFired_[*it] = -9999;
    mapTriggerNameToPrescaleFac_[*it] = -1;
  }
}

// ------------ method called once each job just after ending the event loop  ------------
void
PFMETMuonAnalyzer::endJob() {
}

void PFMETMuonAnalyzer::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {
  //init HLTConfigProvider
  EDConsumerBase::Labels labelTriggerResults;
  EDConsumerBase::labelsForToken(triggerResultsLabel_, labelTriggerResults); 
  const std::string pro = labelTriggerResults.process;

  //bool init(const edm::Run& iRun, const edm::EventSetup& iSetup, const std::string& processName, bool& changed);
  bool changed = true;
  hltConfigInit = false;
  if( hltConfig.init(iRun, iSetup, pro, changed) ) hltConfigInit = true;

  changed = true;
  hltPrescaleInit = false;
  if( hltPrescaleProvider.init(iRun, iSetup, pro, changed) ) hltPrescaleInit = true;
}

// constructors
TreePFCandEventData::TreePFCandEventData(){
}


// set branches
void TreePFCandEventData::SetBranches(int etaBins, int fourierOrder, bool doUEraw)
{
  // --event level--
  tree_->Branch("runNb",  &runNb,       "runNb/i");
  tree_->Branch("eventNb",&eventNb,     "eventNb/i");
  tree_->Branch("LS",     &lumiSection, "LS/i"); 

  tree_->Branch("CentBin",&(this->CentBin),"CentBin/I");
  tree_->Branch("Npix",&(this->Npix),"Npix/I");
  tree_->Branch("NpixelTracks",&(this->NpixelTracks),"NpixelTracks/I");
  tree_->Branch("Ntracks",&(this->Ntracks),"Ntracks/I");
  tree_->Branch("NtracksPtCut",&(this->NtracksPtCut),"NtracksPtCut/I");
  tree_->Branch("NtracksEtaCut",&(this->NtracksEtaCut),"NtracksEtaCut/I");
  tree_->Branch("NtracksEtaPtCut",&(this->NtracksEtaPtCut),"NtracksEtaPtCut/I");
  tree_->Branch("SumET_HF",&(this->SumET_HF),"SumET_HF/F");
  tree_->Branch("SumET_HFplus",&(this->SumET_HFplus),"SumET_HFplus/F");
  tree_->Branch("SumET_HFminus",&(this->SumET_HFminus),"SumET_HFminus/F");
  tree_->Branch("SumET_HFplusEta4",&(this->SumET_HFplusEta4),"SumET_HFplusEta4/F");
  tree_->Branch("SumET_HFminusEta4",&(this->SumET_HFminusEta4),"SumET_HFminusEta4/F");
  tree_->Branch("SumET_HFhit",&(this->SumET_HFhit),"SumET_HFhit/F");
  tree_->Branch("SumET_HFhitPlus",&(this->SumET_HFhitPlus),"SumET_HFhitPlus/F");
  tree_->Branch("SumET_HFhitMinus",&(this->SumET_HFhitMinus),"SumET_HFhitMinus/F");
  tree_->Branch("SumET_ZDC",&(this->SumET_ZDC),"SumET_ZDC/F");
  tree_->Branch("SumET_ZDCplus",&(this->SumET_ZDCplus),"SumET_ZDCplus/F");
  tree_->Branch("SumET_ZDCminus",&(this->SumET_ZDCminus),"SumET_ZDCminus/F");
  tree_->Branch("SumET_EEplus",&(this->SumET_EEplus),"SumET_EEplus/F"); 
  tree_->Branch("SumET_EEminus",&(this->SumET_EEminus),"SumET_EEminus/F");
  tree_->Branch("SumET_EE",&(this->SumET_EE),"SumET_EE/F");
  tree_->Branch("SumET_EB",&(this->SumET_EB),"SumET_EB/F");
  tree_->Branch("SumET_ET",&(this->SumET_ET),"SumET_ET/F");

  tree_->Branch("nPV",&(this->nPV),"nPV/F");
  tree_->Branch("RefVtx_x",&(this->RefVtx_x),"RefVtx_x/F");
  tree_->Branch("RefVtx_y",&(this->RefVtx_y),"RefVtx_y/F");
  tree_->Branch("RefVtx_z",&(this->RefVtx_z),"RefVtx_z/F");
  tree_->Branch("RefVtx_xError",&(this->RefVtx_xError),"RefVtx_xError/F");
  tree_->Branch("RefVtx_yError",&(this->RefVtx_yError),"RefVtx_yError/F");
  tree_->Branch("RefVtx_zError",&(this->RefVtx_zError),"RefVtx_zError/F");

  // -- particle info --
  tree_->Branch("nPFpart",&(this->nPFpart_),"nPFpart/I");
  tree_->Branch("pfId",&(this->pfId_));
  tree_->Branch("pfPt",&(this->pfPt_));
  tree_->Branch("pfEnergy",&(this->pfEnergy_));
  tree_->Branch("pfVsPt",&(this->pfVsPt_));
  tree_->Branch("pfVsPtInitial",&(this->pfVsPtInitial_));
  tree_->Branch("pfArea",&(this->pfArea_));

  tree_->Branch("pfEta",&(this->pfEta_));
  tree_->Branch("pfPhi",&(this->pfPhi_));
  tree_->Branch("pfCharge",&(this->pfCharge_));
  tree_->Branch("pfTheta",&(this->pfTheta_));
  tree_->Branch("pfEt",&(this->pfEt_));

  // -- jet info --
  if(doJets){
    tree_->Branch("njets",&(this->njets_),"njets/I");
    tree_->Branch("jetEnergy",&(this->jetEnergy_));
    tree_->Branch("jetPt",&(this->jetPt_));
    tree_->Branch("jetEta",&(this->jetEta_));
    tree_->Branch("jetPhi",&(this->jetPhi_));
    tree_->Branch("jetMass",&(this->jetMass_));
    tree_->Branch("jetMatchIndex",&(this->jetMatchIndex_));
    tree_->Branch("jetPU",&(this->jetPU_));
  }

  tree_->Branch("vn",this->vn,Form("vn[%d][%d]/F",fourierOrder,etaBins));
  tree_->Branch("psin",this->psin,Form("vpsi[%d][%d]/F",fourierOrder,etaBins));
  tree_->Branch("sumpt",this->sumpt,Form("sumpt[%d]/F",etaBins));
  if(doUEraw){
    tree_->Branch("ueraw",this->ueraw,Form("ueraw[%d]/F",(fourierOrder*etaBins*2*3)));
  }
  
  // -- particle info --
  tree_->Branch("pfMuonPx",&(this->pfMuonPx_));
  tree_->Branch("pfMuonPy",&(this->pfMuonPy_));
  tree_->Branch("pfMuonPz",&(this->pfMuonPz_));
  tree_->Branch("pfTrackerMuon",&(this->pfTrackerMuon_));
  tree_->Branch("pfTrackerMuonPt",&(this->pfTrackerMuonPt_));
  tree_->Branch("pfTrackHits",&(this->pfTrackHits_));
  tree_->Branch("pfDxy",&(this->pfDxy_));
  tree_->Branch("pfDz",&(this->pfDz_));
  tree_->Branch("pfChi2",&(this->pfChi2_));
  tree_->Branch("pfGlobalMuonPt",&(this->pfGlobalMuonPt_));
  tree_->Branch("pfChargedPx",&(this->pfChargedPx_));
  tree_->Branch("pfChargedPy",&(this->pfChargedPy_));
  tree_->Branch("pfChargedPz",&(this->pfChargedPz_));
  tree_->Branch("pfChargedTrackRefPt",&(this->pfChargedTrackRefPt_));
  
  // -- gen info --
  if(doMC){
    tree_->Branch("nGENpart",&(this->nGENpart_),"nGENpart/I");
    tree_->Branch("genPDGId",&(this->genPDGId_));
    tree_->Branch("genPt",&(this->genPt_));
    tree_->Branch("genEta",&(this->genEta_));
    tree_->Branch("genPhi",&(this->genPhi_));
  }

  // -- generalTracks info --
  tree_->Branch("nTRACKpart",&(this->nTRACKpart_),"nTRACKpart/I");
  tree_->Branch("traQual",&(this->traQual_));
  tree_->Branch("traCharge",&(this->traCharge_));
  tree_->Branch("traPt",&(this->traPt_));
  tree_->Branch("traEta",&(this->traEta_));
  tree_->Branch("traPhi",&(this->traPhi_));
  tree_->Branch("traAlgo",&(this->traAlgo_));
  tree_->Branch("traHits",&(this->traHits_));

  // -- MET info --
  tree_->Branch("recoPfMET",&(this->recoPfMET_),"recoPfMET/F");
  tree_->Branch("recoPfMETPhi",&(this->recoPfMETPhi_),"recoPfMETPhi/F");
  tree_->Branch("recoPfMETsumEt",&(this->recoPfMETsumEt_),"recoPfMETsumEt/F");
  tree_->Branch("recoPfMETmEtSig",&(this->recoPfMETmEtSig_),"recoPfMETmEtSig/F");
  tree_->Branch("recoPfMETSig",&(this->recoPfMETSig_),"recoPfMETSig/F");
  
  // -- Muon info (pat::muons) --
  tree_->Branch("nMUpart",&(this->nMUpart_),"nMUpart/I");
  tree_->Branch("muPx",&(this->muPx_));
  tree_->Branch("muPy",&(this->muPy_));
  tree_->Branch("muPz",&(this->muPz_));
  tree_->Branch("muCharge",&(this->muCharge_));
  tree_->Branch("muSelectionType",&(this->muSelectionType_));
  tree_->Branch("muTrackIso",&(this->muTrackIso_));
  tree_->Branch("muCaloIso",&(this->muCaloIso_));
  tree_->Branch("muEcalIso",&(this->muEcalIso_));
  tree_->Branch("muHcalIso",&(this->muHcalIso_));
  tree_->Branch("muHighPurity",&(this->muHighPurity_));
  tree_->Branch("muIsGoodMuon",&(this->muIsGoodMuon_));
  tree_->Branch("muTrkMuArb",&(this->muTrkMuArb_));
  tree_->Branch("muTMOneStaTight",&(this->muTMOneStaTight_));
  tree_->Branch("muNTrkHits",&(this->muNTrkHits_));
  tree_->Branch("muNPixValHits",&(this->muNPixValHits_));
  tree_->Branch("muNPixWMea",&(this->muNPixWMea_));
  tree_->Branch("muNTrkWMea",&(this->muNTrkWMea_));
  tree_->Branch("muStationsMatched",&(this->muStationsMatched_));
  tree_->Branch("muNMuValHits",&(this->muNMuValHits_));
  tree_->Branch("muDxy",&(this->muDxy_));
  tree_->Branch("muDxyErr",&(this->muDxyErr_));
  tree_->Branch("muDz",&(this->muDz_));
  tree_->Branch("muDzErr",&(this->muDzErr_));
  tree_->Branch("muPtInner",&(this->muPtInner_));
  tree_->Branch("muPtErrInner",&(this->muPtErrInner_));
  tree_->Branch("muPtGlobal",&(this->muPtGlobal_));
  tree_->Branch("muPtErrGlobal",&(this->muPtErrGlobal_));
  tree_->Branch("muNormChi2Inner",&(this->muNormChi2Inner_));
  tree_->Branch("muNormChi2Global",&(this->muNormChi2Global_));
  
  // -- Trigger info --
  tree_->Branch("muTrig",&(this->muTrig_));
  tree_->Branch("HLTriggers",&(this->HLTriggers));
  tree_->Branch("trigPrescale",&(this->trigPrescale));
  


}
void TreePFCandEventData::Clear()
{
  // for every event, below variables will be cleaned
  nPFpart_      = 0;
  nGENpart_     = 0;
  njets_        = 0;
  nTRACKpart_   = 0;
  nMUpart_      = 0;

  pfId_.clear();
  genPDGId_.clear();
  pfEnergy_.clear();
  jetEnergy_.clear();
  pfPt_.clear();
  genPt_.clear();
  jetPt_.clear();
  pfEta_.clear();
  genEta_.clear();
  jetEta_.clear();
  pfPhi_.clear();
  genPhi_.clear();
  jetPhi_.clear();
  jetMass_.clear();
  jetPU_.clear();
  jetMatchIndex_.clear();
  pfTheta_.clear();
  pfEt_.clear();
  pfCharge_.clear();
  pfVsPt_.clear();
  pfVsPtInitial_.clear();
  pfArea_.clear();

  pfMuonPx_.clear();
  pfMuonPy_.clear();
  pfMuonPz_.clear();
  pfTrackerMuon_.clear();
  pfTrackerMuonPt_.clear();
  pfTrackHits_.clear();
  pfDxy_.clear();
  pfDz_.clear();
  pfChi2_.clear();
  pfGlobalMuonPt_.clear();
  pfChargedPx_.clear();
  pfChargedPy_.clear();
  pfChargedPz_.clear();
  pfChargedTrackRefPt_.clear();

  traCharge_.clear();
  traPt_.clear();
  traEta_.clear();
  traPhi_.clear();
  traAlgo_.clear();
  traHits_.clear();

  recoPfMET_ = 0;
  recoPfMETPhi_ = 0;
  recoPfMETsumEt_  = 0;
  recoPfMETmEtSig_ = 0;
  recoPfMETSig_    = 0;

  muPt_.clear();
  muPx_.clear();
  muPy_.clear();
  muPz_.clear();
  muEta_.clear();
  muPhi_.clear();
  muCharge_.clear();
  muTrackIso_.clear();
  muCaloIso_.clear();
  muEcalIso_.clear();
  muHcalIso_.clear();
  muTrig_.clear();
  
  muSelectionType_.clear();
  muHighPurity_.clear();
  muIsGoodMuon_.clear();
  muTrkMuArb_.clear();
  muTMOneStaTight_.clear();
  muNTrkHits_.clear();
  muNormChi2Inner_.clear();
  muNPixValHits_.clear();
  muNPixWMea_.clear();
  muNTrkWMea_.clear();
  muStationsMatched_.clear();
  muDxy_.clear();
  muDxyErr_.clear();
  muDz_.clear();
  muDzErr_.clear();
  muPtInner_.clear();
  muPtErrInner_.clear();

  muNMuValHits_.clear();
  muNormChi2Global_.clear();
  muPtGlobal_.clear();
  muPtErrGlobal_.clear();

}

void PFMETMuonAnalyzer::hltReport(const edm::Event &iEvent ,const edm::EventSetup& iSetup)
{
  std::map<std::string, bool> mapTriggernameToTriggerFired;
  std::map<std::string, unsigned int> mapTriggernameToHLTbit;

  for(std::vector<std::string>::const_iterator it=theTriggerNames.begin(); it !=theTriggerNames.end(); ++it){
    mapTriggernameToTriggerFired[*it]=false;
    mapTriggernameToHLTbit[*it]=1000;
  }

  // HLTConfigProvider
  if ( hltConfigInit ) {
    //! Use HLTConfigProvider
    const unsigned int n= hltConfig.size();
    for (std::map<std::string, unsigned int>::iterator it = mapTriggernameToHLTbit.begin(); it != mapTriggernameToHLTbit.end(); it++) {
      unsigned int triggerIndex= hltConfig.triggerIndex( it->first );
      if (it->first == "NoTrigger") continue;
      if (triggerIndex >= n) {
              std::cout << "[PFMETMuonAnalyzer::hltReport] --- TriggerName " << it->first << " not available in config!" << std::endl;
      }
      else {
        it->second= triggerIndex;
        //      std::cout << "[PFMETMuonAnalyzer::hltReport] --- TriggerName " << it->first << " available in config!" << std::endl;
      }
    }
  }
    
  // Get Trigger Results
  try {
    iEvent.getByToken( triggerResultsLabel_, collTriggerResults );
    //    std::cout << "[PFMETMuonAnalyzer::hltReport] --- J/psi TriggerResult is present in current event" << std::endl;
  }
  catch(...) {
    //    std::cout << "[PFMETMuonAnalyzer::hltReport] --- J/psi TriggerResults NOT present in current event" << std::endl;
  }
  if ( collTriggerResults.isValid() && (collTriggerResults->size()==hltConfig.size()) ){
    //    std::cout << "[PFMETMuonAnalyzer::hltReport] --- J/psi TriggerResults IS valid in current event" << std::endl;
      
    // loop over Trigger Results to check if paths was fired
    for(std::vector< std::string >::iterator itHLTNames= theTriggerNames.begin(); itHLTNames != theTriggerNames.end(); itHLTNames++){
      const std::string triggerPathName =  *itHLTNames;
      if ( mapTriggernameToHLTbit[triggerPathName] < 1000 ) {
        if (collTriggerResults->accept( mapTriggernameToHLTbit[triggerPathName] ) ){
          mapTriggerNameToIntFired_[triggerPathName] = 3;
        }
        if (doMC_) {
          mapTriggerNameToPrescaleFac_[triggerPathName] = 1;
        } else {
          //-------prescale factor------------
          if ( hltPrescaleInit && hltPrescaleProvider.prescaleSet(iEvent,iSetup)>=0 ) {
            std::pair<std::vector<std::pair<std::string,int> >,int> detailedPrescaleInfo = hltPrescaleProvider.prescaleValuesInDetail(iEvent, iSetup, triggerPathName);
            //get HLT prescale info from hltPrescaleProvider     
            const int hltPrescale = detailedPrescaleInfo.second;
            //get L1 prescale info from hltPrescaleProvider
            int l1Prescale = -1;     
            if (detailedPrescaleInfo.first.size()==1) {
              l1Prescale = detailedPrescaleInfo.first.at(0).second;
            }
            else if (detailedPrescaleInfo.first.size()>1) {
              std::cout << "[PFMETMuonAnalyzer::hltReport] --- TriggerName " << triggerPathName << " has complex L1 seed " << hltConfig.hltL1GTSeeds(triggerPathName).at(0).second << std::endl;
              std::cout << "[PFMETMuonAnalyzer::hltReport] --- Need to define a proper way to compute the total L1 prescale, default L1 prescale value set to 1 "  << std::endl;
            }
            else {
              std::cout << "[PFMETMuonAnalyzer::hltReport] --- L1 prescale was NOT found for TriggerName " << triggerPathName  << " , default L1 prescale value set to 1 " <<  std::endl;
            }
            //compute the total prescale = HLT prescale * L1 prescale
            mapTriggerNameToPrescaleFac_[triggerPathName] = hltPrescale * l1Prescale;
          }
        }
      }
    }
  } else std::cout << "[PFMETMuonAnalyzer::hltReport] --- TriggerResults NOT valid in current event" << std::endl;

  return;
}

DEFINE_FWK_MODULE(PFMETMuonAnalyzer);
