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
      if (muon::isGoodMuon(muon,(muon::SelectionType) type)) {
         mask = mask | (int) pow(2, type);
      }
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
  Int_t                        nPFpart, nGENpart, njets;
  std::vector<Int_t>           pfId, genPDGId;
  std::vector<Float_t>         pfEnergy, jetEnergy;
  std::vector<Float_t>         pfPt, genPt,  jetPt;
  std::vector<Float_t>         pfEta, genEta,  jetEta;
  std::vector<Float_t>         pfPhi, genPhi,  jetPhi;
  std::vector<Float_t>         jetMass, jetPU;
  std::vector<Int_t>           jetMatchIndex, pfCharge;
  std::vector<Float_t>         pfTheta, pfEt;
  std::vector<Float_t>         pfVsPt, pfVsPtInitial, pfArea;
  Float_t                      vn[200];
  Float_t                      psin[200];
  Float_t                      sumpt[20];
  Float_t                      ueraw[1200];
  // (particle flow charged hadrons and muons)
  std::vector<Float_t>         pfMuonPt, pfMuonPx, pfMuonPy, pfMuonPz, pfMuonPhi, pfMuonEta, pfMuonMt;
  std::vector<Int_t>           pfMuonCharge;
  std::vector<bool>            pfTrackerMuon;
  std::vector<Float_t>         pfTrackerMuonPt;
  std::vector<Int_t>           pfTrackHits;
  std::vector<Float_t>         pfDxy, pfDz, pfChi2;
  std::vector<Float_t>         pfGlobalMuonPt;
  std::vector<Float_t>         pfChargedPt, pfChargedPx, pfChargedPy, pfChargedPz, pfChargedPhi, pfChargedEta;
  std::vector<Float_t>         pfChargedTrackRefPt;

  // -- generalTracks info --
  Int_t                        nTRACKpart;
  std::vector<Int_t>           traQual, traCharge;
  std::vector<Float_t>         traPt,   traEta,  traPhi;
  std::vector<Int_t>           traAlgo,  traHits;
  // track algorithm enum:
  // https://cmssdt.cern.ch/SDT/doxygen/CMSSW_8_0_24/doc/html/da/d0c/TrackBase_8h_source.html#l00099

  // -- MET info --
  Float_t                      recoPFMET, recoPFMETPhi, recoPFMETsumEt;
  Float_t                      recoPFMETmEtSig, recoPFMETSig;

  // -- Muon info (pat::muons) --
  Int_t                        nMUpart;
  std::vector<Float_t>         muPx, muPy, muPz;
  std::vector<Float_t>         muMt, muPt, muEta, muPhi;
  std::vector<Int_t>           muCharge, mu_SelectionType;
  std::vector<Float_t>         muTrackIso, muCaloIso, muEcalIso, muHcalIso; //R0.3 default pp isolation setup
  std::vector<Float_t>         muSumChargedHadronPt, muSumNeutralHadronEt, muSumPhotonEt, muSumPUPt, muPFBasedDBetaIso; //R0.4 (MU POG default PF-based isolation)
  std::vector<bool>            muHighPurity, muIsTightMuon, muIsGoodMuon, muTrkMuArb, muTMOneStaTight;
  std::vector<Int_t>           muSelectionType, muNTrkHits, muNPixValHits, muNPixWMea, muNTrkWMea, muStationsMatched, muNMuValHits;
  std::vector<Float_t>         muDxy, muDxyErr, muDz, muDzErr;
  std::vector<Float_t>         muPtInner, muPtErrInner, muPtGlobal, muPtErrGlobal;
  std::vector<Float_t>         muNormChi2Inner, muNormChi2Global;
  // PF types : 
  // http://cmslxr.fnal.gov/source/DataFormats/ParticleFlowCandidate/interface/PFCandidate.h?v=CMSSW_8_0_24#0044
  // Isolation code reference: https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/ElectroWeakAnalysis/Utilities/src/MuonWithPFIsoProducer.cc
  std::vector<Float_t>         muIso03_vetoPt, muIso03_sumPt, muIso04_sumPt, muIso05_sumPt, muIso03_emEt, muIso04_emEt, muIso05_emEt, muIso03_hadEt, muIso04_hadEt, muIso05_hadEt;
  std::vector<Int_t>           muIso03_nTracks, muIso04_nTracks, muIso05_nTracks;
  std::vector<bool>            muNotPFMuon;

  // -- Trigger info --
  ULong64_t                    HLTriggers;
  std::vector<Int_t>           trigPrescale; // size of this array will be same as number of trigger names
  std::vector<ULong64_t>       muTrig;

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
  Int_t	      verbosity_;

  bool        usePfMuonsOnly_;
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
