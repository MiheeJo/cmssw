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
  usePfMuonsOnly_ = iConfig.getUntrackedParameter<bool> ("usePfMuonsOnly", true);

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


  // Fill MET Object
  Handle<reco::PFMETCollection>   recoPFMETHandle;
  iEvent.getByToken(pfMETLabel_, recoPFMETHandle);
  const reco::PFMET& pfmet = recoPFMETHandle->at(0);

  if (recoPFMETHandle.isValid()) {
    pfEvt_.recoPFMET = pfmet.et();
    pfEvt_.recoPFMETPhi = pfmet.phi();
    pfEvt_.recoPFMETsumEt  = pfmet.sumEt();
    pfEvt_.recoPFMETmEtSig = pfmet.mEtSig();
    pfEvt_.recoPFMETSig    = pfmet.significance();
  } 


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
    if (pt<pfPtMin_) continue;

    int id = pfCandidate.particleId();
    if (skipCharged_ && (abs(id) == 1 || abs(id) == 3)) continue;

    bool matched2Jet = false;
    if (doJets_) {
      pfEvt_.jetMatchIndex.push_back( -1 ); // default index is -1
      
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
              
              pfEvt_.jetMatchIndex[pfEvt_.nPFpart] = ijet; // when a match is found, change -1 to matched jet index
              matched2Jet=true;
              break;
            }
          }
        }
        if(matched2Jet==true) break;
      }
    } // end of if(doJets)

    pfEvt_.pfId.push_back( id );
    pfEvt_.pfPt.push_back( rndSF(pt,4) );
    pfEvt_.pfEnergy.push_back( rndSF(energy,4) );
    pfEvt_.pfVsPt.push_back( rndSF(vsPt,4) );
    pfEvt_.pfVsPtInitial.push_back( rndSF(vsPtInitial,4) );
    pfEvt_.pfArea.push_back( rndSF(vsArea,4) );
    pfEvt_.pfEta.push_back( rndDP(pfCandidate.eta(),3) );
    pfEvt_.pfPhi.push_back( rndDP(pfCandidate.phi(),3) );
    pfEvt_.pfTheta.push_back( rndDP(pfCandidate.theta(),3) ); 
    pfEvt_.pfEt.push_back( rndSF(pfCandidate.et(),4) );
    pfEvt_.pfCharge.push_back( pfCandidate.charge() ); 

    // More information on charged particle(1) and muon(3)
    Float_t TrackerMuon=0, TrackerMuonPt=0, GlobalMuonPt=0;
    Int_t TrackHits=0;
    Float_t Dxy=0, Dz=0, Chi2=0;
    Float_t MuonPt=0, MuonPx=0, MuonPy=0, MuonPz=0;
    Float_t MuonPhi=0, MuonEta=0, MuonTM=0;
    Int_t   MuonCharge=0;
    Float_t ChargedPt=0, ChargedPx=0, ChargedPy=0,ChargedPz=0, ChargedPhi=0, ChargedEta=0, ChargedTrackRefPt=0;

    if ( id==reco::PFCandidate::mu ) {
      MuonPt = pfCandidate.pt();
      MuonPx = pfCandidate.px();
      MuonPy = pfCandidate.py();
      MuonPz = pfCandidate.pz();
      MuonPhi = pfCandidate.phi();
      MuonEta = pfCandidate.eta();
      MuonCharge = pfCandidate.charge();

      MuonTM = TMath::Sqrt( 2*pfCandidate.pt()*pfEvt_.recoPFMET*(1-TMath::Cos(pfCandidate.phi()-pfEvt_.recoPFMETPhi)) );
      
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
    
    if ( id==reco::PFCandidate::h ) {
      ChargedPx = pfCandidate.px();
      ChargedPy = pfCandidate.py();
      ChargedPz = pfCandidate.pz();
      ChargedPhi = pfCandidate.phi();
      ChargedEta = pfCandidate.eta();
      const reco::TrackRef trackRef = pfCandidate.trackRef(); 
      if (trackRef.isNonnull()) {
        ChargedTrackRefPt = rndSF(trackRef->pt(),4);
      }
    }

    pfEvt_.pfMuonMt.push_back( MuonTM ); 
    pfEvt_.pfMuonPt.push_back( MuonPt );
    pfEvt_.pfMuonPx.push_back( MuonPx );
    pfEvt_.pfMuonPy.push_back( MuonPy );
    pfEvt_.pfMuonPz.push_back( MuonPz );
    pfEvt_.pfMuonPhi.push_back( MuonPhi );
    pfEvt_.pfMuonEta.push_back( MuonEta );
    pfEvt_.pfMuonCharge.push_back( MuonCharge );
    pfEvt_.pfTrackerMuon.push_back( TrackerMuon );
    pfEvt_.pfTrackHits.push_back( TrackHits );
    pfEvt_.pfTrackerMuonPt.push_back( TrackerMuonPt );
    pfEvt_.pfGlobalMuonPt.push_back( GlobalMuonPt );
    pfEvt_.pfDxy.push_back( Dxy );
    pfEvt_.pfDz.push_back( Dz );
    pfEvt_.pfChi2.push_back( Chi2 );
    pfEvt_.pfChargedPt.push_back( ChargedPt );
    pfEvt_.pfChargedPx.push_back( ChargedPx );
    pfEvt_.pfChargedPy.push_back( ChargedPy );
    pfEvt_.pfChargedPz.push_back( ChargedPz );
    pfEvt_.pfChargedPhi.push_back( ChargedPhi );
    pfEvt_.pfChargedEta.push_back( ChargedEta );
    pfEvt_.pfChargedTrackRefPt.push_back( ChargedTrackRefPt );

    pfEvt_.nPFpart++;
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
        pfEvt_.genPDGId.push_back( gen.pdgId() );
        pfEvt_.genPt.push_back( rndSF(pt,4) );
        pfEvt_.genEta.push_back( rndDP(eta,3) );
        pfEvt_.genPhi.push_back( rndDP(gen.phi(),3) );
        pfEvt_.nGENpart++;
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
        pfEvt_.jetPt.push_back( pt );
        pfEvt_.jetEnergy.push_back( energy );
        pfEvt_.jetEta.push_back( jet.eta() );
        pfEvt_.jetPhi.push_back( jet.phi() );
	pfEvt_.jetMass.push_back( jet.mass() );
	pfEvt_.jetPU.push_back( jet.pileup() );
        pfEvt_.njets++;
      }
    }
  }

  // Fill generalTracks
  edm::Handle<reco::TrackCollection> trackCollection;
  iEvent.getByToken(trackLabel_, trackCollection);
  const reco::TrackCollection* trackColl= &(*trackCollection);
 
  for(unsigned itrack=0;itrack<trackColl->size(); itrack++) { 
    const reco::Track tra = trackColl->at(itrack); 

    if (tra.quality(reco::TrackBase::qualityByName(qualityString_))) pfEvt_.traQual.push_back(1);
    else pfEvt_.traQual.push_back(0);
    
    pfEvt_.traCharge.push_back(tra.charge());
    pfEvt_.traPt.push_back(rndSF(tra.pt(),4));
    pfEvt_.traEta.push_back(rndDP(tra.eta(),3));
    pfEvt_.traPhi.push_back(rndDP(tra.phi(),3));
    pfEvt_.traAlgo.push_back(tra.algo());
    pfEvt_.traHits.push_back(tra.numberOfValidHits());
    pfEvt_.nTRACKpart++;
  }


  // Fill single muon information
  Handle<pat::MuonCollection> muonCollection; 
  iEvent.getByToken(muonLabel_, muonCollection);
  for (unsigned imuon=0;imuon<muonCollection->size(); imuon++) { 
    const pat::Muon& muon = muonCollection->at(imuon);

    // Inner track information of a muon
    reco::TrackRef iTrack = muon.innerTrack();
    if (iTrack.isNonnull()) {
      // When this muon is NOT PF muon, take values from reco::muon collection
      float muon_pt = muon.pt();
      float muon_px = muon.px();
      float muon_py = muon.py();
      float muon_pz = muon.pz();
      float muon_eta = muon.eta();
      float muon_phi = muon.phi();
      float muon_charge = muon.charge();

      muonIDmask(muon);
      pfEvt_.muSelectionType.push_back( muonIDmask(muon) );
      pfEvt_.muHighPurity.push_back( iTrack->quality(reco::TrackBase::highPurity) );
      pfEvt_.muIsTightMuon.push_back( muon::isTightMuon(muon, *privtx) );
      pfEvt_.muIsGoodMuon.push_back( muon::isGoodMuon(muon, muon::TMOneStationTight) );
      pfEvt_.muTrkMuArb.push_back( muon.muonID("TrackerMuonArbitrated") );
      pfEvt_.muTMOneStaTight.push_back( muon.muonID("TMOneStationTight") );
      pfEvt_.muNTrkHits.push_back( iTrack->found() );
      pfEvt_.muNormChi2Inner.push_back( iTrack->normalizedChi2() );
      pfEvt_.muNPixValHits.push_back( iTrack->hitPattern().numberOfValidPixelHits() );
      pfEvt_.muNPixWMea.push_back( iTrack->hitPattern().pixelLayersWithMeasurement() );
      pfEvt_.muNTrkWMea.push_back( iTrack->hitPattern().trackerLayersWithMeasurement() );
      pfEvt_.muStationsMatched.push_back( muon.numberOfMatchedStations() );
      pfEvt_.muDxy.push_back( iTrack->dxy(pfEvt_.RefVtx) );
      pfEvt_.muDxyErr.push_back( iTrack->dxyError() );
      pfEvt_.muDz.push_back( iTrack->dz(pfEvt_.RefVtx) );
      pfEvt_.muDzErr.push_back( iTrack->dzError() );
      pfEvt_.muPtInner.push_back( iTrack->pt() );
      pfEvt_.muPtErrInner.push_back( iTrack->ptError() );

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

      pfEvt_.muNMuValHits.push_back( muNMuValHits );
      pfEvt_.muNormChi2Global.push_back( muNormChi2Global );
      pfEvt_.muPtGlobal.push_back( muPtGlobal );
      pfEvt_.muPtErrGlobal.push_back( muPtErrGlobal );

      // Muon triggers are matched to muons?
      ULong64_t muTrig=0;
      for (unsigned int iTr = 0; iTr<HLTLastFilters.size(); ++iTr) {
        const pat::TriggerObjectStandAloneCollection mu1HLTMatchesFilter = muon.triggerObjectMatchesByFilter( HLTLastFilters[iTr] );
        const pat::TriggerObjectStandAloneCollection mu1HLTMatchesPath = muon.triggerObjectMatchesByPath( theTriggerNames.at(iTr), true, false );
        
        if (mu1HLTMatchesFilter.size() > 0) {
          muTrig += pow(2,iTr);
//          cout << "muTrig: " << iTr << " " << muTrig << endl;
        }
      }
      pfEvt_.muTrig.push_back( muTrig );
      
      // Loop over PFCollection to get sum of energy within a cone around muon
      // Interpret "track" as charged particles (e,mu, chraged hadrons)
      // Interpret "em" as photons and also as electromagnetic energy in HF.
      // Interpret "had" as neutral hadrons and also as hadronic energy in HF.
      
      // Ask for PfMuon consistency if requested
      bool muonFound = false;

      // Set isolations
      float iso03_vetoPt=0;
      float iso03_sumPt=0, iso04_sumPt=0, iso05_sumPt=0;
      float iso03_emEt=0, iso04_emEt=0, iso05_emEt=0;
      float iso03_hadEt=0, iso04_hadEt=0, iso05_hadEt=0;
      int iso03_nTracks=0, iso04_nTracks=0, iso05_nTracks=0;
      bool muNotPFMuon = false;
      
      // Loop on all PF candidates
      for (unsigned icand=0; icand<pfCandidateColl->size(); icand++) {
        const reco::PFCandidate pfCandidate = pfCandidateColl->at(icand);
        if (pfCandidate.pt()<pfPtMin_) continue; // Don't use low pT tracks for a sum of pT in a cone 

        // Check the muon is in the PF collection when required
        bool thisIsTheMuon = false;
        if ( iTrack.isNonnull() && pfCandidate.particleId()==reco::PFCandidate::mu &&
             (pfCandidate.pt()-muon_pt)<0.0001 && (pfCandidate.eta()-muon_eta)<0.0001 &&
             (pfCandidate.phi()-muon_phi)<0.0001 && pfCandidate.charge()==muon_charge
           ) {
          thisIsTheMuon = true;
          muonFound = true;
          // Take muon information when it matches with PFCandidate collection
          if (usePfMuonsOnly_) {
            muon_pt = pfCandidate.pt();
            muon_px = pfCandidate.px();
            muon_py = pfCandidate.py();
            muon_pz = pfCandidate.pz();
            muon_eta = pfCandidate.eta();
            muon_phi = pfCandidate.phi();
            muon_charge = pfCandidate.charge();
          }
        }

        // Get dR. Nothing to add if dR>0.5
        double deltaR = reco::deltaR(muon.eta(),muon.phi(),pfCandidate.eta(),pfCandidate.phi());
        //if (deltaR>0.5) continue;

        // Fill "tracker" components
        if (   pfCandidate.particleId()==reco::PFCandidate::h
            || pfCandidate.particleId()==reco::PFCandidate::e
            || pfCandidate.particleId()==reco::PFCandidate::mu ) {
              if (!thisIsTheMuon) {
                iso05_sumPt += pfCandidate.pt();
                iso05_nTracks++;
                if (deltaR<0.4) {
                  iso04_sumPt += pfCandidate.pt();
                  iso04_nTracks++;
                }
                if (deltaR<0.3) {
                  iso03_sumPt += pfCandidate.pt();
                  iso03_nTracks++;
                } else { // Veto sumPt for cross-check
                  iso03_vetoPt += pfCandidate.pt();
                }
              }
        // Fill "em" components
        } else if (   pfCandidate.particleId()==reco::PFCandidate::gamma
                   || pfCandidate.particleId()==reco::PFCandidate::egamma_HF) {
                      iso05_emEt += pfCandidate.pt();
                      if (deltaR<0.4) iso04_emEt += pfCandidate.pt();
                      if (deltaR<0.3) iso03_emEt += pfCandidate.pt();
        // Fill "had" components
        } else if (   pfCandidate.particleId()==reco::PFCandidate::h0
                   || pfCandidate.particleId()==reco::PFCandidate::h_HF) {
                      iso05_hadEt += pfCandidate.pt();
                      if (deltaR<0.4) iso04_hadEt += pfCandidate.pt();
                      if (deltaR<0.3) iso03_hadEt += pfCandidate.pt();
        }
      }

      // Do not take this muon (under explicit request) if it is not a PfMuon
      if (usePfMuonsOnly_ && (!muonFound)) muNotPFMuon = true;

      pfEvt_.muIso03_vetoPt.push_back( iso03_vetoPt );
      pfEvt_.muIso03_sumPt.push_back( iso03_sumPt );
      pfEvt_.muIso04_sumPt.push_back( iso04_sumPt );
      pfEvt_.muIso05_sumPt.push_back( iso05_sumPt );
      pfEvt_.muIso03_emEt.push_back( iso03_emEt );
      pfEvt_.muIso04_emEt.push_back( iso04_emEt );
      pfEvt_.muIso05_emEt.push_back( iso05_emEt );
      pfEvt_.muIso03_hadEt.push_back( iso03_hadEt );
      pfEvt_.muIso04_hadEt.push_back( iso04_hadEt );
      pfEvt_.muIso05_hadEt.push_back( iso05_hadEt );
      pfEvt_.muIso03_nTracks.push_back( iso03_nTracks );
      pfEvt_.muIso04_nTracks.push_back( iso04_nTracks );
      pfEvt_.muIso05_nTracks.push_back( iso05_nTracks );
      pfEvt_.muNotPFMuon.push_back( muNotPFMuon );

      // Muon POG standard isolation cuts and variables
      pfEvt_.muSumChargedHadronPt.push_back( muon.pfIsolationR03().sumChargedHadronPt );
      pfEvt_.muSumNeutralHadronEt.push_back( muon.pfIsolationR03().sumNeutralHadronEt );
      pfEvt_.muSumPhotonEt.push_back( muon.pfIsolationR03().sumPhotonEt );
      pfEvt_.muSumPUPt.push_back( muon.pfIsolationR03().sumPUPt );
      float PFBasedIso = (muon.pfIsolationR04().sumChargedHadronPt + max(0., static_cast<double>(muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt)))/muon_pt;
      pfEvt_.muPFBasedDBetaIso.push_back( PFBasedIso );
      pfEvt_.muTrackIso.push_back( muon.trackIso() );
      pfEvt_.muCaloIso.push_back( muon.caloIso() );
      pfEvt_.muEcalIso.push_back( muon.ecalIso() );
      pfEvt_.muHcalIso.push_back( muon.hcalIso() );


      float transverseMass = TMath::Sqrt( 2*muon_pt*pfEvt_.recoPFMET*(1-TMath::Cos(muon_phi-pfEvt_.recoPFMETPhi)) );
      pfEvt_.muMt.push_back( transverseMass ); 
      pfEvt_.muPt.push_back( muon_pt ); 
      pfEvt_.muPx.push_back( muon_px );
      pfEvt_.muPy.push_back( muon_py );
      pfEvt_.muPz.push_back( muon_pz );     
      pfEvt_.muEta.push_back( muon_eta );      
      pfEvt_.muPhi.push_back( muon_phi );
      pfEvt_.muCharge.push_back( muon_charge );

      pfEvt_.nMUpart++;
    } // end of iTrack.isNonnull

  } // end of muonCollection loop

  // Trigger information for each event
  for (unsigned int iTr = 0 ; iTr < theTriggerNames.size() ; iTr++) {
    if (mapTriggerNameToIntFired_[theTriggerNames.at(iTr)] == 3) {
      pfEvt_.HLTriggers += pow(2,iTr);
//      cout << "HLTriggers: " << iTr << " " << pfEvt_.HLTriggers << endl;
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
  tree_->Branch("nPFpart",&(this->nPFpart),"nPFpart/I");
  tree_->Branch("pfId",&(this->pfId));
  tree_->Branch("pfPt",&(this->pfPt));
  tree_->Branch("pfEnergy",&(this->pfEnergy));
  tree_->Branch("pfVsPt",&(this->pfVsPt));
  tree_->Branch("pfVsPtInitial",&(this->pfVsPtInitial));
  tree_->Branch("pfArea",&(this->pfArea));

  tree_->Branch("pfEta",&(this->pfEta));
  tree_->Branch("pfPhi",&(this->pfPhi));
  tree_->Branch("pfCharge",&(this->pfCharge));
  tree_->Branch("pfTheta",&(this->pfTheta));
  tree_->Branch("pfEt",&(this->pfEt));

  // -- jet info --
  if(doJets){
    tree_->Branch("njets",&(this->njets),"njets/I");
    tree_->Branch("jetEnergy",&(this->jetEnergy));
    tree_->Branch("jetPt",&(this->jetPt));
    tree_->Branch("jetEta",&(this->jetEta));
    tree_->Branch("jetPhi",&(this->jetPhi));
    tree_->Branch("jetMass",&(this->jetMass));
    tree_->Branch("jetMatchIndex",&(this->jetMatchIndex));
    tree_->Branch("jetPU",&(this->jetPU));
  }

  tree_->Branch("vn",this->vn,Form("vn[%d][%d]/F",fourierOrder,etaBins));
  tree_->Branch("psin",this->psin,Form("vpsi[%d][%d]/F",fourierOrder,etaBins));
  tree_->Branch("sumpt",this->sumpt,Form("sumpt[%d]/F",etaBins));
  if(doUEraw){
    tree_->Branch("ueraw",this->ueraw,Form("ueraw[%d]/F",(fourierOrder*etaBins*2*3)));
  }
  
  // -- particle info --
  tree_->Branch("pfMuonMt",&(this->pfMuonMt));
  tree_->Branch("pfMuonPt",&(this->pfMuonPt));
  tree_->Branch("pfMuonPx",&(this->pfMuonPx));
  tree_->Branch("pfMuonPy",&(this->pfMuonPy));
  tree_->Branch("pfMuonPz",&(this->pfMuonPz));
  tree_->Branch("pfMuonPhi",&(this->pfMuonPhi));
  tree_->Branch("pfMuonEta",&(this->pfMuonEta));
  tree_->Branch("pfMuonCharge",&(this->pfMuonCharge));
  tree_->Branch("pfTrackerMuon",&(this->pfTrackerMuon));
  tree_->Branch("pfTrackerMuonPt",&(this->pfTrackerMuonPt));
  tree_->Branch("pfTrackHits",&(this->pfTrackHits));
  tree_->Branch("pfDxy",&(this->pfDxy));
  tree_->Branch("pfDz",&(this->pfDz));
  tree_->Branch("pfChi2",&(this->pfChi2));
  tree_->Branch("pfGlobalMuonPt",&(this->pfGlobalMuonPt));
  tree_->Branch("pfChargedPt",&(this->pfChargedPt));
  tree_->Branch("pfChargedPx",&(this->pfChargedPx));
  tree_->Branch("pfChargedPy",&(this->pfChargedPy));
  tree_->Branch("pfChargedPz",&(this->pfChargedPz));
  tree_->Branch("pfChargedPhi",&(this->pfChargedPhi));
  tree_->Branch("pfChargedEta",&(this->pfChargedEta));
  tree_->Branch("pfChargedTrackRefPt",&(this->pfChargedTrackRefPt));
  
  // -- gen info --
  if(doMC){
    tree_->Branch("nGENpart",&(this->nGENpart),"nGENpart/I");
    tree_->Branch("genPDGId",&(this->genPDGId));
    tree_->Branch("genPt",&(this->genPt));
    tree_->Branch("genEta",&(this->genEta));
    tree_->Branch("genPhi",&(this->genPhi));
  }

  // -- generalTracks info --
  tree_->Branch("nTRACKpart",&(this->nTRACKpart),"nTRACKpart/I");
  tree_->Branch("traQual",&(this->traQual));
  tree_->Branch("traCharge",&(this->traCharge));
  tree_->Branch("traPt",&(this->traPt));
  tree_->Branch("traEta",&(this->traEta));
  tree_->Branch("traPhi",&(this->traPhi));
  tree_->Branch("traAlgo",&(this->traAlgo));
  tree_->Branch("traHits",&(this->traHits));

  // -- MET info --
  tree_->Branch("recoPFMET",&(this->recoPFMET),"recoPFMET/F");
  tree_->Branch("recoPFMETPhi",&(this->recoPFMETPhi),"recoPFMETPhi/F");
  tree_->Branch("recoPFMETsumEt",&(this->recoPFMETsumEt),"recoPFMETsumEt/F");
  tree_->Branch("recoPFMETmEtSig",&(this->recoPFMETmEtSig),"recoPFMETmEtSig/F");
  tree_->Branch("recoPFMETSig",&(this->recoPFMETSig),"recoPFMETSig/F");
  
  // -- Muon info (pat::muons) --
  tree_->Branch("nMUpart",&(this->nMUpart),"nMUpart/I");
  tree_->Branch("muPx",&(this->muPx));
  tree_->Branch("muPy",&(this->muPy));
  tree_->Branch("muPz",&(this->muPz));
  tree_->Branch("muPt",&(this->muPt));
  tree_->Branch("muMt",&(this->muMt));
  tree_->Branch("muEta",&(this->muEta));
  tree_->Branch("muPhi",&(this->muPhi));
  tree_->Branch("muCharge",&(this->muCharge));
  tree_->Branch("muSelectionType",&(this->muSelectionType));
  tree_->Branch("muTrackIso",&(this->muTrackIso));
  tree_->Branch("muCaloIso",&(this->muCaloIso));
  tree_->Branch("muEcalIso",&(this->muEcalIso));
  tree_->Branch("muHcalIso",&(this->muHcalIso));
  tree_->Branch("muSumChargedHadronPt",&(this->muSumChargedHadronPt));
  tree_->Branch("muSumNeutralHadronEt",&(this->muSumNeutralHadronEt));
  tree_->Branch("muSumPhotonEt",&(this->muSumPhotonEt));
  tree_->Branch("muSumPUPt",&(this->muSumPUPt));
  tree_->Branch("muPFBasedDBetaIso",&(this->muPFBasedDBetaIso));
  tree_->Branch("muHighPurity",&(this->muHighPurity));
  tree_->Branch("muIsTightMuon",&(this->muIsTightMuon));
  tree_->Branch("muIsGoodMuon",&(this->muIsGoodMuon));
  tree_->Branch("muTrkMuArb",&(this->muTrkMuArb));
  tree_->Branch("muTMOneStaTight",&(this->muTMOneStaTight));
  tree_->Branch("muNTrkHits",&(this->muNTrkHits));
  tree_->Branch("muNPixValHits",&(this->muNPixValHits));
  tree_->Branch("muNPixWMea",&(this->muNPixWMea));
  tree_->Branch("muNTrkWMea",&(this->muNTrkWMea));
  tree_->Branch("muStationsMatched",&(this->muStationsMatched));
  tree_->Branch("muNMuValHits",&(this->muNMuValHits));
  tree_->Branch("muDxy",&(this->muDxy));
  tree_->Branch("muDxyErr",&(this->muDxyErr));
  tree_->Branch("muDz",&(this->muDz));
  tree_->Branch("muDzErr",&(this->muDzErr));
  tree_->Branch("muPtInner",&(this->muPtInner));
  tree_->Branch("muPtErrInner",&(this->muPtErrInner));
  tree_->Branch("muPtGlobal",&(this->muPtGlobal));
  tree_->Branch("muPtErrGlobal",&(this->muPtErrGlobal));
  tree_->Branch("muNormChi2Inner",&(this->muNormChi2Inner));
  tree_->Branch("muNormChi2Global",&(this->muNormChi2Global));
  tree_->Branch("muIso03_vetoPt",&(this->muIso03_vetoPt));
  tree_->Branch("muIso03_sumPt",&(this->muIso03_sumPt));
  tree_->Branch("muIso04_sumPt",&(this->muIso04_sumPt));
  tree_->Branch("muIso05_sumPt",&(this->muIso05_sumPt));
  tree_->Branch("muIso03_emEt",&(this->muIso03_emEt));
  tree_->Branch("muIso04_emEt",&(this->muIso04_emEt));
  tree_->Branch("muIso05_emEt",&(this->muIso05_emEt));
  tree_->Branch("muIso03_hadEt",&(this->muIso03_hadEt));
  tree_->Branch("muIso04_hadEt",&(this->muIso04_hadEt));
  tree_->Branch("muIso05_hadEt",&(this->muIso05_hadEt));
  tree_->Branch("muIso03_nTracks",&(this->muIso03_nTracks));
  tree_->Branch("muIso04_nTracks",&(this->muIso04_nTracks));
  tree_->Branch("muIso05_nTracks",&(this->muIso05_nTracks));
  tree_->Branch("muNotPFMuon",&(this->muNotPFMuon));
  
  // -- Trigger info --
  tree_->Branch("muTrig",&(this->muTrig));
  tree_->Branch("HLTriggers",&(this->HLTriggers),"HLTriggers/l");
  tree_->Branch("trigPrescale",&(this->trigPrescale));
  
}


void TreePFCandEventData::Clear()
{
  // for every event, below variables will be cleaned
  runNb = 0;
  eventNb = 0;
  lumiSection = 0;

  CentBin = 0;
  Npix = 0;
  NpixelTracks = 0;
  Ntracks = 0;
  NtracksPtCut = 0;
  NtracksEtaCut = 0;
  NtracksEtaPtCut = 0;
  SumET_HF = 0;
  SumET_HFplus = 0;
  SumET_HFminus = 0;
  SumET_HFplusEta4 = 0;
  SumET_HFminusEta4 = 0;
  SumET_HFhit = 0;
  SumET_HFhitPlus = 0;
  SumET_HFhitMinus = 0;
  SumET_ZDC = 0;
  SumET_ZDCplus = 0;
  SumET_ZDCminus = 0;
  SumET_EEplus = 0;
  SumET_EEminus = 0;
  SumET_EE = 0;
  SumET_EB = 0;
  SumET_ET = 0;
 
  nPV = 0;
  RefVtx_x = 0;
  RefVtx_y = 0;
  RefVtx_z = 0;
  RefVtx_xError = 0;
  RefVtx_yError = 0;
  RefVtx_zError = 0;

  // -- particle info --
  nPFpart = 0;
  pfId.clear();
  pfPt.clear();
  pfEnergy.clear();
  pfVsPt.clear();
  pfVsPtInitial.clear();
  pfArea.clear();
 
  pfEta.clear();
  pfPhi.clear();
  pfCharge.clear();
  pfTheta.clear();
  pfEt.clear();

  // -- jet info --
  if(doJets){
    njets = 0;
    jetEnergy.clear();
    jetPt.clear();
    jetEta.clear();
    jetPhi.clear();
    jetMass.clear();
    jetMatchIndex.clear();
    jetPU.clear();
  }

  // -- particle info --
  pfMuonPt.clear();
  pfMuonMt.clear();
  pfMuonPx.clear();
  pfMuonPy.clear();
  pfMuonPz.clear();
  pfMuonEta.clear();
  pfMuonPhi.clear();
  pfMuonCharge.clear();
  pfTrackerMuon.clear();
  pfTrackerMuonPt.clear();
  pfTrackHits.clear();
  pfDxy.clear();
  pfDz.clear();
  pfChi2.clear();
  pfGlobalMuonPt.clear();
  pfChargedPt.clear();
  pfChargedPx.clear();
  pfChargedPy.clear();
  pfChargedPz.clear();
  pfChargedPhi.clear();
  pfChargedEta.clear();
  pfChargedTrackRefPt.clear();
  
  // -- gen info --
  if(doMC){
    nGENpart = 0;
    genPDGId.clear();
    genPt.clear();
    genEta.clear();
    genPhi.clear();
  }

  // -- generalTracks info --
  nTRACKpart = 0;
  traQual.clear();
  traCharge.clear();
  traPt.clear();
  traEta.clear();
  traPhi.clear();
  traAlgo.clear();
  traHits.clear();

  // -- MET info --
  recoPFMET = 0;
  recoPFMETPhi = 0;
  recoPFMETsumEt = 0;
  recoPFMETmEtSig = 0;
  recoPFMETSig = 0;
  
  // -- Muon info (pat::muons) --
  nMUpart = 0;
  muPx.clear();
  muPy.clear();
  muPz.clear();
  muPt.clear();
  muMt.clear();
  muEta.clear();
  muPhi.clear();
  muCharge.clear();
  muSelectionType.clear();
  muTrackIso.clear();
  muCaloIso.clear();
  muEcalIso.clear();
  muHcalIso.clear();
  muSumChargedHadronPt.clear();
  muSumNeutralHadronEt.clear();
  muSumPhotonEt.clear();
  muSumPUPt.clear();
  muPFBasedDBetaIso.clear();
  muHighPurity.clear();
  muIsTightMuon.clear();
  muIsGoodMuon.clear();
  muTrkMuArb.clear();
  muTMOneStaTight.clear();
  muNTrkHits.clear();
  muNPixValHits.clear();
  muNPixWMea.clear();
  muNTrkWMea.clear();
  muStationsMatched.clear();
  muNMuValHits.clear();
  muDxy.clear();
  muDxyErr.clear();
  muDz.clear();
  muDzErr.clear();
  muPtInner.clear();
  muPtErrInner.clear();
  muPtGlobal.clear();
  muPtErrGlobal.clear();
  muNormChi2Inner.clear();
  muNormChi2Global.clear();
  muIso03_sumPt.clear();
  muIso04_sumPt.clear();
  muIso05_sumPt.clear();
  muIso03_emEt.clear();
  muIso04_emEt.clear();
  muIso05_emEt.clear();
  muIso03_hadEt.clear();
  muIso04_hadEt.clear();
  muIso05_hadEt.clear();
  muIso03_nTracks.clear();
  muIso04_nTracks.clear();
  muIso05_nTracks.clear();
  muNotPFMuon.clear();

  // -- Trigger info --
  muTrig.clear();
  HLTriggers = 0;
  trigPrescale.clear();
 
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
