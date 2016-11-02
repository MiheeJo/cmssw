// -*- C++ -*-
//
// Package:    simpleMuonAnalyzer/MuonAnalyzer
// Class:      MuonAnalyzer
// 
/**\class MuonAnalyzer MuonAnalyzer.cc simpleMuonAnalyzer/MuonAnalyzer/plugins/MuonAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mihee Jo
//         Created:  Mon, 19 Sep 2016 16:05:31 GMT
//
//


// system include files
#include <memory>
#include <iostream>
#include <sstream>
#include <utility>

#include <TTree.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TClonesArray.h>
#include <TMath.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class MuonAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit MuonAnalyzer(const edm::ParameterSet&);
    ~MuonAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
    
    // initialization
    void InitEvent();
    void InitTree();
    
    // cuts (reco::Muon cuts are only available)
    bool isSoftMuon(const reco::Muon *aMuon) ;
    bool isMuonInAccept(const reco::Muon *aMuon, std::string muonType);
    
    // genparticles
    void fillGenInfo();
    bool isAbHadron(int pdgID);
    bool isAMixedbHadron(int pdgID, int momPdgID);
    reco::GenParticleRef findMotherRef(reco::GenParticleRef GenParticleMother, int GenParticlePDG);
    std::pair<int, std::pair<float, float> >  findGenMCInfo(const reco::GenParticle *genJpsi);
    reco::GenParticleRef findDaughterRef(reco::GenParticleRef GenParticleDaughter, int GenParticlePDG);
    
    // miscellaneous
    TLorentzVector lorentzMomentum(const reco::Candidate::LorentzVector& p);
    
    
    // ----------member data ---------------------------
    edm::Service<TFileService> fs;
    TTree* myTree;
 
    static const int Max_QQ_size = 100;
    static const int Max_mu_size = 10000;

    // centrality
    int centBin;
    
    TClonesArray* Reco_mu_4mom;
    int Reco_mu_size; // Number of reconstructed muons
    int Reco_mu_charge[Max_mu_size]; // Vector of charge of muons
    int Reco_mu_type[Max_mu_size]; // Vector of type of muon (glb&trk=0, !glb&trk=1)
    
    bool Reco_mu_highPurity[Max_mu_size];    // Vector of high purity flag  
    
    int Reco_mu_nPixValHits[Max_mu_size];  // Number of valid pixel hits in sta muons
    int Reco_mu_nMuValHits[Max_mu_size];  // Number of valid muon hits in sta muons
    int Reco_mu_nTrkHits[Max_mu_size];  // track hits global muons
    int Reco_mu_nPixWMea[Max_mu_size];  // pixel layers with measurement for inner track muons
    int Reco_mu_nTrkWMea[Max_mu_size];  // track layers with measurement for inner track muons
    int Reco_mu_StationsMatched[Max_mu_size];  // number of stations matched for inner track muons
    float Reco_mu_normChi2_inner[Max_mu_size];  // chi2/ndof for inner track muons
    float Reco_mu_normChi2_global[Max_mu_size];  // chi2/ndof for global muons
    float Reco_mu_dxy[Max_mu_size];  // dxy for inner track muons
    float Reco_mu_dxyErr[Max_mu_size];  // dxy error for inner track muons
    float Reco_mu_dz[Max_mu_size];  // dz for inner track muons
    float Reco_mu_dzErr[Max_mu_size];  // dz error for inner track muons
    float Reco_mu_pt_inner[Max_mu_size];  // pT for inner track muons
    float Reco_mu_pt_global[Max_mu_size];  // pT for global muons
    float Reco_mu_ptErr_inner[Max_mu_size];  // pT error for inner track muons
    float Reco_mu_ptErr_global[Max_mu_size];  // pT error for global muons
    
    std::string _muType; // type of muon (Glb, Trk) 

    TClonesArray* Gen_mu_4mom;
    TClonesArray* Gen_QQ_4mom;
    TClonesArray* Gen_QQ_mupl_4mom;
    TClonesArray* Gen_QQ_mumi_4mom;
    int Gen_QQ_size; // number of generated Onia
    float Gen_QQ_ctau[Max_QQ_size];    // ctau: flight time
    float Gen_QQ_ctau3D[Max_QQ_size]; // ctau3D: 3D flight time    
    int Gen_QQ_momId[Max_QQ_size]; // pdgId of mother particle of 2 muons
    int Gen_mu_size; // number of generated muons
    int Gen_mu_charge[Max_mu_size]; // muon charge

    int _oniaPDG;

    unsigned int runNb;
    unsigned int eventNb;
    unsigned int lumiSection;
    math::XYZPoint RefVtx;

    // switches
    bool _isMC;
    bool _isHI;
    bool _checkIDCuts;
    bool _checkAcceptance;

    // handles
    edm::Handle<reco::MuonCollection> collMuon;
    edm::Handle<reco::GenParticleCollection> collGenParticles;
    edm::Handle<reco::VertexCollection> collPriVtx;

    edm::EDGetTokenT<reco::MuonCollection> _muonToken;
    edm::EDGetTokenT<reco::GenParticleCollection> _genParticleToken;
    edm::EDGetTokenT<reco::VertexCollection> _thePVsToken;
    edm::EDGetTokenT<reco::Centrality> _centralityTagToken;
    edm::EDGetTokenT<int> _centralityBinTagToken;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MuonAnalyzer::MuonAnalyzer(const edm::ParameterSet& iConfig):
  _muType(iConfig.getParameter<std::string>("muonType")),
  _oniaPDG(iConfig.getParameter<int>("oniaPDG")),
  _isMC(iConfig.getParameter<bool>("isMC")),
  _isHI(iConfig.getParameter<bool>("isHI")),
  _checkIDCuts(iConfig.getParameter<bool>("checkIDCuts")),
  _checkAcceptance(iConfig.getParameter<bool>("checkAcceptance")),
  _muonToken(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("srcMuon"))),
  _genParticleToken(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))),
  _thePVsToken(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertexTag"))),
  _centralityTagToken(consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag> ("CentralitySrc"))),
  _centralityBinTagToken(consumes<int>(iConfig.getParameter<edm::InputTag> ("CentralityBinSrc")))
{
  //now do what ever initialization is needed
  usesResource("TFileService");
}


MuonAnalyzer::~MuonAnalyzer()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MuonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;

  InitEvent();

  iEvent.getByToken(_muonToken,collMuon);
  iEvent.getByToken(_thePVsToken,collPriVtx);
  reco::VertexCollection::const_iterator privtx;
  privtx = collPriVtx->begin();
  RefVtx = privtx->position();

  if (_isMC) {
    iEvent.getByToken(_genParticleToken,collGenParticles);
  }

  if (_isHI) {
    // Obtain various detailed centrality information from 'centrality'
    edm::Handle<reco::Centrality> centrality;
    iEvent.getByToken(_centralityTagToken, centrality); 

    // Obtain centrality bin number. (200 bins available)
    edm::Handle<int> cbin_;
    iEvent.getByToken(_centralityBinTagToken, cbin_);
    centBin = *cbin_;
  }

  if (Gen_QQ_size >= Max_QQ_size) {
    std::cout << "Too many GEN dimuons: " << Gen_QQ_size << std::endl;
    std::cout << "Maximum allowed: " << Max_QQ_size << std::endl;
    return;
  }

  if (Gen_mu_size >= Max_mu_size) {
    std::cout << "Too many GEN muons: " << Gen_mu_size << std::endl;
    std::cout << "Maximum allowed: " << Max_mu_size << std::endl;
    return;
  }
  
  if (Reco_mu_size >= Max_mu_size) {
    std::cout << "Too many RECO muons: " << Reco_mu_size << std::endl;
    std::cout << "Maximum allowed: " << Max_mu_size << std::endl;
    return;
  }

  // fillGenInfo()
  if (collGenParticles.isValid()) {
    for(std::vector<reco::GenParticle>::const_iterator it=collGenParticles->begin();
        it!=collGenParticles->end();++it) {
      const reco::GenParticle* gen = &(*it);

      if (abs(gen->pdgId()) == _oniaPDG  && gen->status() == 2 &&
          gen->numberOfDaughters() >= 2) {

        reco::GenParticleRef genMuon1 = findDaughterRef(gen->daughterRef(0), gen->pdgId());
        reco::GenParticleRef genMuon2 = findDaughterRef(gen->daughterRef(1), gen->pdgId());

        if ( abs(genMuon1->pdgId()) == 13 &&
             abs(genMuon2->pdgId()) == 13 &&
             genMuon1->status() == 1 &&
             genMuon2->status() == 1 ) {
          
          std::pair<int, std::pair<float, float> > MCinfo = findGenMCInfo(gen);
          Gen_QQ_momId[Gen_QQ_size] = MCinfo.first;
          Gen_QQ_ctau[Gen_QQ_size] = 10.0*MCinfo.second.first;
          Gen_QQ_ctau3D[Gen_QQ_size] = 10.0*MCinfo.second.second;
          
          TLorentzVector vJpsi = lorentzMomentum(gen->p4());
          new((*Gen_QQ_4mom)[Gen_QQ_size])TLorentzVector(vJpsi);

          TLorentzVector vMuon1 = lorentzMomentum(genMuon1->p4());
          TLorentzVector vMuon2 = lorentzMomentum(genMuon2->p4());
            
          if (genMuon1->charge() > genMuon2->charge()) {
            new((*Gen_QQ_mupl_4mom)[Gen_QQ_size])TLorentzVector(vMuon1);
            new((*Gen_QQ_mumi_4mom)[Gen_QQ_size])TLorentzVector(vMuon2);
          } else {
            new((*Gen_QQ_mupl_4mom)[Gen_QQ_size])TLorentzVector(vMuon2);
            new((*Gen_QQ_mumi_4mom)[Gen_QQ_size])TLorentzVector(vMuon1);
          }

          Gen_QQ_size++;
        }
      }

      if (abs(gen->pdgId()) == 13  && gen->status() == 1) {
        Gen_mu_charge[Gen_mu_size] = gen->charge();

        TLorentzVector vMuon = lorentzMomentum(gen->p4());
        new((*Gen_mu_4mom)[Gen_mu_size])TLorentzVector(vMuon);

        Gen_mu_size++;
      }
    } // End of "collGenParticles" loop
  } // if (collGenParticles.isValid()) 
  // END OF fillGenInfo()
  
  if (collMuon.isValid()) {
    for (std::vector<reco::Muon>::const_iterator it=collMuon->begin();
        it!=collMuon->end();++it) {
      
      const reco::Muon *muon = &(*it);
        
      // Skip reco single muon if this is not in acceptance region and doesn't pass the cut
      if ( (_checkIDCuts && isSoftMuon(muon)) || (_checkAcceptance && isMuonInAccept(muon, _muType)) ) {
        // fillTreeMuon()
        Reco_mu_charge[Reco_mu_size] = muon->charge();
        int iType=0;
        if (!(muon->innerTrack().isNonnull())) continue;
        // if isGLB && isTRK -> iType=0, isTRK -> iType=1,
        if (muon->isTrackerMuon()) iType=1;
        if (muon->isGlobalMuon() && muon->isTrackerMuon()) iType=0;
        Reco_mu_type[Reco_mu_size] = iType;
        
        TLorentzVector vMuon = lorentzMomentum(muon->p4());
        new((*Reco_mu_4mom)[Reco_mu_size])TLorentzVector(vMuon);

        reco::TrackRef iTrack = muon->innerTrack();
        
        Reco_mu_highPurity[Reco_mu_size] = iTrack->quality(reco::TrackBase::highPurity);
        Reco_mu_nTrkHits[Reco_mu_size] = iTrack->found();
        Reco_mu_normChi2_inner[Reco_mu_size] = iTrack->normalizedChi2();
        Reco_mu_nPixValHits[Reco_mu_size] = iTrack->hitPattern().numberOfValidPixelHits();
        Reco_mu_nPixWMea[Reco_mu_size] = iTrack->hitPattern().pixelLayersWithMeasurement();
        Reco_mu_nTrkWMea[Reco_mu_size] = iTrack->hitPattern().trackerLayersWithMeasurement();
        Reco_mu_StationsMatched[Reco_mu_size] = muon->numberOfMatchedStations();
        Reco_mu_dxy[Reco_mu_size] = iTrack->dxy(RefVtx);
        Reco_mu_dxyErr[Reco_mu_size] = iTrack->dxyError();
        Reco_mu_dz[Reco_mu_size] = iTrack->dz(RefVtx);
        Reco_mu_dzErr[Reco_mu_size] = iTrack->dzError();
        Reco_mu_pt_inner[Reco_mu_size] = iTrack->pt();
        Reco_mu_ptErr_inner[Reco_mu_size] = iTrack->ptError();
        
        if (muon->isGlobalMuon()) {
          reco::TrackRef gTrack = muon->globalTrack();
          Reco_mu_nMuValHits[Reco_mu_size] = gTrack->hitPattern().numberOfValidMuonHits();
          Reco_mu_normChi2_global[Reco_mu_size] = gTrack->normalizedChi2();
          Reco_mu_pt_global[Reco_mu_size] = gTrack->pt();
          Reco_mu_ptErr_global[Reco_mu_size] = gTrack->ptError();
        } else {
          Reco_mu_nMuValHits[Reco_mu_size] = -1;
          Reco_mu_normChi2_global[Reco_mu_size] = 999;
          Reco_mu_pt_global[Reco_mu_size] = -1;
          Reco_mu_ptErr_global[Reco_mu_size] = -1;
        }
        Reco_mu_size++;
      }
    } // END OF for (std::vector<reco::Muon>::const_iterator it=collMuon->begin(); it!=collMuon->end();++it) 
  }
  
  std::cout << "RECO muons: " << Reco_mu_size << std::endl;
  myTree->Fill();
  // END OF fillTreeMuon()
}

void MuonAnalyzer::InitEvent() {
  std::cout << "InitEvent" << std::endl;
  Reco_mu_size = 0;
  Gen_mu_size = 0;
  Gen_QQ_size = 0;

  Reco_mu_4mom->Clear();
  if (_isMC) {
    Gen_mu_4mom->Clear();
    Gen_QQ_4mom->Clear();
    Gen_QQ_mupl_4mom->Clear();
    Gen_QQ_mumi_4mom->Clear();
  }
}

void MuonAnalyzer::InitTree() {
   Reco_mu_4mom = new TClonesArray("TLorentzVector", 200);
   if (_isMC) {
     Gen_mu_4mom = new TClonesArray("TLorentzVector", 100);
     Gen_QQ_4mom = new TClonesArray("TLorentzVector", 10);
     Gen_QQ_mupl_4mom = new TClonesArray("TLorentzVector", 10);
     Gen_QQ_mumi_4mom = new TClonesArray("TLorentzVector", 10);
   }

   myTree = fs->make<TTree>("myTree","My TTree of dimuons");

   myTree->Branch("eventNb", &eventNb,   "eventNb/i");
   myTree->Branch("runNb",   &runNb,     "runNb/i");
   myTree->Branch("LS", &lumiSection, "LS/i"); 
   if (_isHI) myTree->Branch("Centrality", &centBin, "Centrality/I");

   myTree->Branch("Reco_mu_size", &Reco_mu_size,  "Reco_mu_size/I");
   myTree->Branch("Reco_mu_type", Reco_mu_type,   "Reco_mu_type[Reco_mu_size]/I");
   myTree->Branch("Reco_mu_charge", Reco_mu_charge,   "Reco_mu_charge[Reco_mu_size]/I");
   myTree->Branch("Reco_mu_4mom", "TClonesArray", &Reco_mu_4mom, 32000, 0);
   myTree->Branch("Reco_mu_highPurity", Reco_mu_highPurity,   "Reco_mu_highPurity[Reco_mu_size]/O");
   myTree->Branch("Reco_mu_nPixValHits", Reco_mu_nPixValHits,   "Reco_mu_nPixValHits[Reco_mu_size]/I");
   myTree->Branch("Reco_mu_nMuValHits", Reco_mu_nMuValHits,   "Reco_mu_nMuValHits[Reco_mu_size]/I");
   myTree->Branch("Reco_mu_nTrkHits",Reco_mu_nTrkHits, "Reco_mu_nTrkHits[Reco_mu_size]/I");
   myTree->Branch("Reco_mu_normChi2_inner",Reco_mu_normChi2_inner, "Reco_mu_normChi2_inner[Reco_mu_size]/F");
   myTree->Branch("Reco_mu_normChi2_global",Reco_mu_normChi2_global, "Reco_mu_normChi2_global[Reco_mu_size]/F");
   myTree->Branch("Reco_mu_nPixWMea",Reco_mu_nPixWMea, "Reco_mu_nPixWMea[Reco_mu_size]/I");
   myTree->Branch("Reco_mu_nTrkWMea",Reco_mu_nTrkWMea, "Reco_mu_nTrkWMea[Reco_mu_size]/I");
   myTree->Branch("Reco_mu_StationsMatched",Reco_mu_StationsMatched, "Reco_mu_StationsMatched[Reco_mu_size]/I");
   myTree->Branch("Reco_mu_dxy",Reco_mu_dxy, "Reco_mu_dxy[Reco_mu_size]/F");
   myTree->Branch("Reco_mu_dxyErr",Reco_mu_dxyErr, "Reco_mu_dxyErr[Reco_mu_size]/F");
   myTree->Branch("Reco_mu_dz",Reco_mu_dz, "Reco_mu_dz[Reco_mu_size]/F");
   myTree->Branch("Reco_mu_dzErr",Reco_mu_dzErr, "Reco_mu_dzErr[Reco_mu_size]/F");
   myTree->Branch("Reco_mu_pt_inner",Reco_mu_pt_inner, "Reco_mu_pt_inner[Reco_mu_size]/F");
   myTree->Branch("Reco_mu_pt_global",Reco_mu_pt_global, "Reco_mu_pt_global[Reco_mu_size]/F");
   myTree->Branch("Reco_mu_ptErr_inner",Reco_mu_ptErr_inner, "Reco_mu_ptErr_inner[Reco_mu_size]/F");
   myTree->Branch("Reco_mu_ptErr_global",Reco_mu_ptErr_global, "Reco_mu_ptErr_global[Reco_mu_size]/F");

  if (_isMC) {
    myTree->Branch("Gen_QQ_size",      &Gen_QQ_size,    "Gen_QQ_size/I");
    myTree->Branch("Gen_QQ_4mom",      "TClonesArray", &Gen_QQ_4mom, 32000, 0);
    myTree->Branch("Gen_QQ_momId",      Gen_QQ_momId,    "Gen_QQ_momId[Gen_QQ_size]/I");
    myTree->Branch("Gen_QQ_ctau",      Gen_QQ_ctau,    "Gen_QQ_ctau[Gen_QQ_size]/F");
    myTree->Branch("Gen_QQ_ctau3D",      Gen_QQ_ctau3D,    "Gen_QQ_ctau3D[Gen_QQ_size]/F");  
    myTree->Branch("Gen_QQ_mupl_4mom", "TClonesArray", &Gen_QQ_mupl_4mom, 32000, 0); 
    myTree->Branch("Gen_QQ_mumi_4mom", "TClonesArray", &Gen_QQ_mumi_4mom, 32000, 0);

    myTree->Branch("Gen_mu_size",   &Gen_mu_size,  "Gen_mu_size/I");
    myTree->Branch("Gen_mu_charge", Gen_mu_charge, "Gen_mu_charge[Gen_mu_size]/I");
    myTree->Branch("Gen_mu_4mom",   "TClonesArray", &Gen_mu_4mom, 32000, 0);
  }

}

TLorentzVector MuonAnalyzer::lorentzMomentum(const reco::Candidate::LorentzVector& p) {
  TLorentzVector res;
  res.SetPtEtaPhiM(p.pt(), p.eta(), p.phi(), p.mass());

  return res;
}

bool MuonAnalyzer::isAbHadron(int pdgID) {
  if (abs(pdgID) == 511 || abs(pdgID) == 521 || abs(pdgID) == 531 || abs(pdgID) == 5122) return true;
  return false;
}

bool MuonAnalyzer::isAMixedbHadron(int pdgID, int momPdgID) {
  if ((abs(pdgID) == 511 && abs(momPdgID) == 511 && pdgID*momPdgID < 0) || 
      (abs(pdgID) == 531 && abs(momPdgID) == 531 && pdgID*momPdgID < 0)) 
      return true;
  return false;
}

reco::GenParticleRef   
MuonAnalyzer::findMotherRef(reco::GenParticleRef GenParticleMother, int GenParticlePDG) {
  for(int i=0; i<1000; ++i) {
    if (GenParticleMother.isNonnull() && (GenParticleMother->pdgId()==GenParticlePDG) && GenParticleMother->numberOfMothers()>0) {        
      GenParticleMother = GenParticleMother->motherRef();
    } else break;
  }
  return GenParticleMother;
}

std::pair<int, std::pair<float, float> >  
MuonAnalyzer::findGenMCInfo(const reco::GenParticle* genJpsi) {

  int momJpsiID = 0;
  float trueLife = -99.;
  float trueLife3D = -99.;

  if (genJpsi->numberOfMothers()>0) {
    TVector3 trueVtx(0.0,0.0,0.0);
    TVector3 trueP(0.0,0.0,0.0);
    TVector3 trueVtxMom(0.0,0.0,0.0);

    trueVtx.SetXYZ(genJpsi->vertex().x(),genJpsi->vertex().y(),genJpsi->vertex().z());
    trueP.SetXYZ(genJpsi->momentum().x(),genJpsi->momentum().y(),genJpsi->momentum().z());

    bool aBhadron = false;
    reco::GenParticleRef Jpsimom = findMotherRef(genJpsi->motherRef(), genJpsi->pdgId());      
    if (Jpsimom.isNull()) {
      std::pair<float, float> trueLifePair = std::make_pair(trueLife, trueLife3D);
      std::pair<int, std::pair<float, float>> result = std::make_pair(momJpsiID, trueLifePair);
      return result;
    } 
    else if (Jpsimom->numberOfMothers()<=0) {
      if (isAbHadron(Jpsimom->pdgId())) {  
        momJpsiID = Jpsimom->pdgId();
        trueVtxMom.SetXYZ(Jpsimom->vertex().x(),Jpsimom->vertex().y(),Jpsimom->vertex().z());
        aBhadron = true;
      }
    } 
    else {
      reco::GenParticleRef Jpsigrandmom = findMotherRef(Jpsimom->motherRef(), Jpsimom->pdgId());   
      if (isAbHadron(Jpsimom->pdgId())) {       
        if (Jpsigrandmom.isNonnull() && isAMixedbHadron(Jpsimom->pdgId(),Jpsigrandmom->pdgId())) {       
          momJpsiID = Jpsigrandmom->pdgId();
          trueVtxMom.SetXYZ(Jpsigrandmom->vertex().x(),Jpsigrandmom->vertex().y(),Jpsigrandmom->vertex().z());
        } 
        else {                  
          momJpsiID = Jpsimom->pdgId();
          trueVtxMom.SetXYZ(Jpsimom->vertex().x(),Jpsimom->vertex().y(),Jpsimom->vertex().z());
        }
        aBhadron = true;
      } 
      else if (Jpsigrandmom.isNonnull() && isAbHadron(Jpsigrandmom->pdgId()))  {  
        if (Jpsigrandmom->numberOfMothers()<=0) {
          momJpsiID = Jpsigrandmom->pdgId();
          trueVtxMom.SetXYZ(Jpsigrandmom->vertex().x(),Jpsigrandmom->vertex().y(),Jpsigrandmom->vertex().z());
        } 
        else { 
          reco::GenParticleRef JpsiGrandgrandmom = findMotherRef(Jpsigrandmom->motherRef(), Jpsigrandmom->pdgId());
          if (JpsiGrandgrandmom.isNonnull() && isAMixedbHadron(Jpsigrandmom->pdgId(),JpsiGrandgrandmom->pdgId())) {
            momJpsiID = JpsiGrandgrandmom->pdgId();
            trueVtxMom.SetXYZ(JpsiGrandgrandmom->vertex().x(),JpsiGrandgrandmom->vertex().y(),JpsiGrandgrandmom->vertex().z());
          } 
          else {
            momJpsiID = Jpsigrandmom->pdgId();
            trueVtxMom.SetXYZ(Jpsigrandmom->vertex().x(),Jpsigrandmom->vertex().y(),Jpsigrandmom->vertex().z());
          }
        }
        aBhadron = true;
      }
    }
    if (!aBhadron) {
      momJpsiID = Jpsimom->pdgId();
      trueVtxMom.SetXYZ(Jpsimom->vertex().x(),Jpsimom->vertex().y(),Jpsimom->vertex().z()); 
    }
    
    TVector3 vdiff = trueVtx - trueVtxMom;
    trueLife = vdiff.Perp()*3.096916/trueP.Perp();
    trueLife3D = vdiff.Mag()*3.096916/trueP.Mag();
  }

  std::pair<float, float> trueLifePair = std::make_pair(trueLife, trueLife3D);
  std::pair<int, std::pair<float, float> > result = std::make_pair(momJpsiID, trueLifePair);
  return result;
}

reco::GenParticleRef  
MuonAnalyzer::findDaughterRef(reco::GenParticleRef GenParticleDaughter, int GenParticlePDG) {
  reco::GenParticleRef GenParticleTmp = GenParticleDaughter;   
  for(int j=0; j<1000; ++j) { 
    if (GenParticleTmp.isNonnull() && ((GenParticleTmp->pdgId()==GenParticlePDG) || (GenParticleTmp->pdgId()==GenParticleDaughter->pdgId())) && GenParticleTmp->numberOfDaughters()>0) { 
      GenParticleDaughter = GenParticleTmp;
      GenParticleTmp      = GenParticleDaughter->daughterRef(0);
    } else break;
  }
  if (GenParticleTmp.isNonnull() && (GenParticleTmp->pdgId()==GenParticleDaughter->pdgId())) {
    GenParticleDaughter = GenParticleTmp;
  }   
  return GenParticleDaughter;
}

bool MuonAnalyzer::isSoftMuon(const reco::Muon *aMuon) {
  if (!(aMuon->innerTrack().isNonnull())) return false;
  else
    return (
            aMuon->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5   &&
            aMuon->innerTrack()->hitPattern().pixelLayersWithMeasurement()   > 0   &&
            fabs(aMuon->innerTrack()->dxy(RefVtx)) < 0.3 &&
            fabs(aMuon->innerTrack()->dz(RefVtx)) < 20. 
           );
}

bool MuonAnalyzer::isMuonInAccept(const reco::Muon* aMuon, const std::string muonType) {
  if (muonType == (std::string)("GLB")) {
    return (fabs(aMuon->eta()) < 2.4 &&
            ((fabs(aMuon->eta()) < 1.2 && aMuon->pt() >= 3.5) ||
             (1.2 <= fabs(aMuon->eta()) && fabs(aMuon->eta()) < 2.1 && aMuon->pt() >= 5.77-1.89*fabs(aMuon->eta())) ||
             (2.1 <= fabs(aMuon->eta()) && aMuon->pt() >= 1.8)));
  }
  else if (muonType == (std::string)("TRK")) {
    return (fabs(aMuon->eta()) < 2.4 &&
            ((fabs(aMuon->eta()) < 1.3 && aMuon->pt() >= 3.3) ||
             (1.3 <= fabs(aMuon->eta()) && fabs(aMuon->eta()) < 2.2 && aMuon->p() >= 2.9) ||
             (2.2 <= fabs(aMuon->eta()) && aMuon->pt() >= 0.8)));
  }
  else {
    std::cout << "ERROR: Incorrect Muon Type ===>>> Regards it as GLB" << std::endl;
    return (fabs(aMuon->eta()) < 2.4 &&
            ((fabs(aMuon->eta()) < 1.2 && aMuon->pt() >= 3.5) ||
             (1.2 <= fabs(aMuon->eta()) && fabs(aMuon->eta()) < 2.1 && aMuon->pt() >= 5.77-1.89*fabs(aMuon->eta())) ||
             (2.1 <= fabs(aMuon->eta()) && aMuon->pt() >= 1.8)));
  }
  
  return false;
}


// ------------ method called once each job just before starting event loop  ------------
void MuonAnalyzer::beginJob() {
  InitTree();
}

// ------------ method called once each job just after ending the event loop  ------------
void MuonAnalyzer::endJob() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonAnalyzer);
