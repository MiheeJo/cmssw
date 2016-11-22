// system include files
#include <memory>
#include <iostream>
#include <vector>
#include <cmath>

// stl
#include <algorithm>
#include <utility>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HeavyIonEvent/interface/VoronoiBackground.h"

#include "RecoJets/JetAlgorithms/interface/JetAlgoHelper.h"
#include "RecoHI/HiJetAlgos/interface/UEParameters.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"
#include "SimDataFormats/HiGenData/interface/GenHIEvent.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TMath.h"
#include "TStopwatch.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"

//rounds to first few nonzero sig figs
float rndSF(float value, int nSignificantDigits) 
{
  if(value==0) return 0; 
  
  float dSign = (value > 0.0) ? 1 : -1; 
  value *= dSign; 

  int nOffset = static_cast<int>(log10(value)); 
  if(nOffset>=0) ++nOffset;  

  float dScale = pow(10.0,nSignificantDigits-nOffset); 

  return dSign * static_cast<float>(round(value*dScale))/dScale;    
}

//keeps first n digits after the decimal place
inline float rndDP(float value, int nPlaces)
{
  return float(int(value*pow(10,nPlaces))/pow(10,nPlaces));
}

int muonIDmask(const pat::Muon& muon)
{
   int mask = 0;
   int type;
   for (type=muon::All; type<=muon::RPCMuLoose; type++)
      if (muon::isGoodMuon(muon,(muon::SelectionType) type))
         mask = mask | (int) pow(2, type);

   return mask;
}
//----------------------------------------------------------------
//
// DiJet ana Event Data Tree definition
//
class TreePFCandEventData
{
public:
  // ===== Class Methods =====
  void SetDefaults();
  TreePFCandEventData();
  void SetTree(TTree * t) { tree_=t; }
  void SetBranches(int etaBins, int fourierOrder, bool doUEraw = 0);
  void Clear();
  bool doJets;
  bool doMC;
  int muonIDmask(const pat::Muon* muon);

  Float_t         jdphi_;

  unsigned int runNb;
  unsigned int eventNb;
  unsigned int lumiSection;

  // -- centrality variables --
  Int_t           CentBin;
  Int_t           Npix, NpixelTracks, Ntracks, NtracksPtCut, NtracksEtaCut, NtracksEtaPtCut;
  Float_t         SumET_HF, SumET_HFplus, SumET_HFminus, SumET_HFplusEta4, SumET_HFminusEta4;
  Float_t         SumET_HFhit, SumET_HFhitPlus, SumET_HFhitMinus;
  Float_t         SumET_ZDC, SumET_ZDCplus, SumET_ZDCminus;
  Float_t         SumET_EEplus, SumET_EEminus;
  Float_t         SumET_EE, SumET_EB, SumET_ET;

  // -- Primary Vetex --
  Float_t         nPV;
  math::XYZPoint  RefVtx;
  Float_t         RefVtx_x, RefVtx_y, RefVtx_z;
  Float_t         RefVtx_xError, RefVtx_yError, RefVtx_zError;

  // -- particle info & GEN info --
  Int_t                        nPFpart_, nGENpart_, njets_;
  std::vector<Int_t>           pfId_, genPDGId_;
  std::vector<Float_t>         pfEnergy_, jetEnergy_;
  std::vector<Float_t>         pfPt_, genPt_,  jetPt_;
  std::vector<Float_t>         pfEta_, genEta_,  jetEta_;
  std::vector<Float_t>         pfPhi_, genPhi_,  jetPhi_;
  std::vector<Float_t>         jetMass_, jetPU_;
  std::vector<Int_t>           jetMatchIndex_, pfCharge_;
  std::vector<Float_t>         pfTheta_, pfEt_;
  std::vector<Float_t>         pfVsPt_, pfVsPtInitial_, pfArea_;
  Float_t                      vn[200];
  Float_t                      psin[200];
  Float_t                      sumpt[20];
  Float_t                      ueraw[1200];
  // (particle flow charged hadrons and muons)
  std::vector<Float_t>         pfMuonPx_, pfMuonPy_, pfMuonPz_;
  std::vector<bool>            pfTrackerMuon_;
  std::vector<Float_t>         pfTrackerMuonPt_;
  std::vector<Int_t>           pfTrackHits_;
  std::vector<Float_t>         pfDxy_, pfDz_, pfChi2_;
  std::vector<Float_t>         pfGlobalMuonPt_;
  std::vector<Float_t>         pfChargedPx_, pfChargedPy_, pfChargedPz_;
  std::vector<Float_t>         pfChargedTrackRefPt_;

  // -- generalTracks info --
  Int_t                        nTRACKpart_;
  std::vector<Int_t>           traQual_, traCharge_;
  std::vector<Float_t>         traPt_,   traEta_,  traPhi_;
  std::vector<Int_t>           traAlgo_,  traHits_;
  // track algorithm enum:
  // https://cmssdt.cern.ch/SDT/doxygen/CMSSW_8_0_24/doc/html/da/d0c/TrackBase_8h_source.html#l00099

  // -- MET info --
  Float_t                      recoPfMET_, recoPfMETPhi_, recoPfMETsumEt_;
  Float_t                      recoPfMETmEtSig_, recoPfMETSig_;

  // -- Muon info (pat::muons) --
  Int_t                        nMUpart_;
  std::vector<Float_t>         muPx_, muPy_, muPz_;
  std::vector<Float_t>         muPt_, muEta_, muPhi_;
  std::vector<Int_t>           muCharge_, mu_SelectionType_;
  std::vector<Float_t>         muTrackIso_, muCaloIso_, muEcalIso_, muHcalIso_;
  std::vector<bool>            muHighPurity_, muIsTightMuon_, muIsGoodMuon_, muTrkMuArb_, muTMOneStaTight_;
  std::vector<Int_t>           muSelectionType_, muNTrkHits_, muNPixValHits_, muNPixWMea_, muNTrkWMea_, muStationsMatched_, muNMuValHits_;
  std::vector<Float_t>         muDxy_, muDxyErr_, muDz_, muDzErr_;
  std::vector<Float_t>         muPtInner_, muPtErrInner_, muPtGlobal_, muPtErrGlobal_;
  std::vector<Float_t>         muNormChi2Inner_, muNormChi2Global_;

  // -- Trigger info --
  std::vector<ULong64_t>       muTrig_;
  ULong64_t                    HLTriggers;
  std::vector<Int_t>           trigPrescale;

private:
  TTree*                       tree_;
};

class PFMETMuonAnalyzer : public edm::EDAnalyzer {
public:
  explicit PFMETMuonAnalyzer(const edm::ParameterSet&);
  ~PFMETMuonAnalyzer();

  // class methods


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  void beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup);
  void hltReport(const edm::Event &iEvent, const edm::EventSetup& iSetup);

  // ----------member data ---------------------------
  edm::Service<TFileService> fs;
  edm::Handle<reco::VoronoiMap> backgrounds_;
  edm::Handle<std::vector<float> > vn_;
  edm::Handle<reco::CandidateView> candidates_;
  edm::Handle<edm::TriggerResults> collTriggerResults;

  // === Ana setup ===

  // Event Info
  edm::InputTag pfCandidateLabel_;
  edm::EDGetTokenT<reco::PFCandidateCollection> pfCandidatePF_;
  edm::EDGetTokenT<reco::CandidateView> pfCandidateView_;
  edm::EDGetTokenT<reco::GenParticleCollection> genLabel_;
  edm::EDGetTokenT<pat::JetCollection> jetLabel_;
  edm::EDGetTokenT<std::vector<float> > srcVorFloat_;
  edm::EDGetTokenT<reco::VoronoiMap> srcVorMap_;
  edm::EDGetTokenT<reco::VertexCollection> vtxLabel_;
  edm::EDGetTokenT<reco::TrackCollection> trackLabel_;
  edm::EDGetTokenT<reco::PFMETCollection> pfMETLabel_;
  edm::EDGetTokenT<pat::MuonCollection> muonLabel_;
  edm::EDGetTokenT<reco::Centrality> centralityTagToken_;
  edm::EDGetTokenT<int> centralityBinTagToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerResultsLabel_;

  TTree	  *pfTree_;
  TreePFCandEventData pfEvt_;

  // cuts
  Double_t        pfPtMin_;
  Double_t        jetPtMin_;
  Double_t        genPtMin_;

  int           fourierOrder_;
  int           etaBins_;

  // debug
  Int_t	  verbosity_;

  bool        doJets_;
  bool        doMC_;
  bool        isHI_;
  bool        isPA_;
  bool        doVS_;
  bool        doUEraw_;
  bool        skipCharged_;
  std::string qualityString_;
  
  // Trigger prescales, names, information
  HLTConfigProvider       hltConfig;
  bool                    hltConfigInit;
  HLTPrescaleProvider     hltPrescaleProvider;
  bool                    hltPrescaleInit;

  std::vector<std::string>     theTriggerNames, HLTLastFilters;
  std::map<std::string, int>   mapTriggerNameToIntFired_;
  std::map<std::string, int>   mapTriggerNameToPrescaleFac_;
};

//    int id = pfCandidate.particleId();
//    enum ParticleType :
//      X=0,     // undefined
//      h,       // charged hadron
//      e,       // electron 
//      mu,      // muon 
//      gamma,   // photon
//      h0,      // neutral hadron
//      h_HF,        // HF tower identified as a hadron
//      egamma_HF    // HF tower identified as an EM particle
