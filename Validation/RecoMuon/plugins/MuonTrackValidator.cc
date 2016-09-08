#include "Validation/RecoMuon/plugins/MuonTrackValidator.h"
#include "DQMServices/ClientConfig/interface/FitSlicesYTool.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "SimTracker/TrackAssociation/plugins/ParametersDefinerForTPESProducer.h"
#include "SimTracker/TrackAssociation/plugins/CosmicParametersDefinerForTPESProducer.h"

#include "TMath.h"
#include <TF1.h>

using namespace std;
using namespace edm;

Double_t MuonTrackValidator::findNcoll(int hiBin) {
    const int nbins = 200;
    const Double_t Ncoll[nbins] = {1976.95, 1944.02, 1927.29, 1891.9, 1845.3, 1807.2, 1760.45, 1729.18, 1674.8, 1630.3, 1590.52, 1561.72, 1516.1, 1486.5, 1444.68, 1410.88, 1376.4, 1347.32, 1309.71, 1279.98, 1255.31, 1219.89, 1195.13, 1165.96, 1138.92, 1113.37, 1082.26, 1062.42, 1030.6, 1009.96, 980.229, 955.443, 936.501, 915.97, 892.063, 871.289, 847.364, 825.127, 806.584, 789.163, 765.42, 751.187, 733.001, 708.31, 690.972, 677.711, 660.682, 640.431, 623.839, 607.456, 593.307, 576.364, 560.967, 548.909, 530.475, 519.575, 505.105, 490.027, 478.133, 462.372, 451.115, 442.642, 425.76, 416.364, 405.154, 392.688, 380.565, 371.167, 360.28, 348.239, 340.587, 328.746, 320.268, 311.752, 300.742, 292.172, 281.361, 274.249, 267.025, 258.625, 249.931, 240.497, 235.423, 228.63, 219.854, 214.004, 205.425, 199.114, 193.618, 185.644, 180.923, 174.289, 169.641, 161.016, 157.398, 152.151, 147.425, 140.933, 135.924, 132.365, 127.017, 122.127, 117.817, 113.076, 109.055, 105.16, 101.323, 98.098, 95.0548, 90.729, 87.6495, 84.0899, 80.2237, 77.2201, 74.8848, 71.3554, 68.7745, 65.9911, 63.4136, 61.3859, 58.1903, 56.4155, 53.8486, 52.0196, 49.2921, 47.0735, 45.4345, 43.8434, 41.7181, 39.8988, 38.2262, 36.4435, 34.8984, 33.4664, 31.8056, 30.351, 29.2074, 27.6924, 26.7754, 25.4965, 24.2802, 22.9651, 22.0059, 21.0915, 19.9129, 19.1041, 18.1487, 17.3218, 16.5957, 15.5323, 14.8035, 14.2514, 13.3782, 12.8667, 12.2891, 11.61, 11.0026, 10.3747, 9.90294, 9.42648, 8.85324, 8.50121, 7.89834, 7.65197, 7.22768, 6.7755, 6.34855, 5.98336, 5.76555, 5.38056, 5.11024, 4.7748, 4.59117, 4.23247, 4.00814, 3.79607, 3.68702, 3.3767, 3.16309, 2.98282, 2.8095, 2.65875, 2.50561, 2.32516, 2.16357, 2.03235, 1.84061, 1.72628, 1.62305, 1.48916, 1.38784, 1.28366, 1.24693, 1.18552, 1.16085, 1.12596, 1.09298, 1.07402, 1.06105, 1.02954};
    return Ncoll[hiBin];
};

void MuonTrackValidator::bookHistograms(DQMStore::IBooker& ibooker, edm::Run const&, edm::EventSetup const& setup) {

  int j=0;
  for (unsigned int ww=0;ww<associators.size();ww++){
    for (unsigned int www=0;www<label.size();www++){

      std::cout << "##bookHistograms## ww: " << ww << "/" << associators.size()
                << " www: " << www << "/" << label.size() << std::endl;
      ibooker.cd();
      InputTag algo = label[www];
      string dirName=dirName_;
      if (algo.process()!="")
        dirName+=algo.process()+"_";
      if(algo.label()!="")
        dirName+=algo.label();
      if(algo.instance()!="")
        dirName+=("_"+algo.instance());
      if (dirName.find("Tracks")<dirName.length()){
        dirName.replace(dirName.find("Tracks"),6,"Trks");
      }
      if (dirName.find("UpdatedAtVtx")<dirName.length()){
        dirName.replace(dirName.find("UpdatedAtVtx"),12,"UpdAtVtx");
      }
      string assoc= associators[ww];
      if (assoc.find("tpToTkmuTrackAssociation")<assoc.length()){
        dirName+="_TkAsso";
      }
      std::replace(dirName.begin(), dirName.end(), ':', '_');
      ibooker.setCurrentFolder(dirName.c_str());

      setUpVectors();
      if (!_pTRanges.empty()) setUpAsymmVectors();
      centBin = -1; // initialize centrality bin
      std::cout << "nintpT : " << nintpT << " " << pTintervals[j].size() << std::endl;
      std::cout << "pTRanges: ";
      for (int idx=0; idx<=nintpT ; idx++) std::cout << pTRanges[idx] << " " << pTintervals[j][idx] << " ";
      std::cout << std::endl;

      ibooker.goUp();
      string subDirName = dirName + "/simulation";
      ibooker.setCurrentFolder(subDirName.c_str());
      if (!_pTRanges.empty()) h_ptSIM.push_back( ibooker.book1D("ptSIM", "generated p_{t}", nintpT, pTRanges) );
      else h_ptSIM.push_back( ibooker.book1D("ptSIM", "generated p_{t}", nintpT, minpT, maxpT ) );
      h_etaSIM.push_back( ibooker.book1D("etaSIM", "generated pseudorapidity", nint, min, max ) );
      h_tracksSIM.push_back( ibooker.book1D("tracksSIM","number of simulated tracks",200,-0.5,99.5) );
      h_vertposSIM.push_back( ibooker.book1D("vertposSIM","Transverse position of sim vertices",100,0.,120.) );
      
      ibooker.cd();
      ibooker.setCurrentFolder(dirName.c_str());
      h_tracks.push_back( ibooker.book1D("tracks","number of reconstructed tracks",200,-0.5,19.5) );
      h_fakes.push_back( ibooker.book1D("fakes","number of fake reco tracks",20,-0.5,19.5) );
      h_charge.push_back( ibooker.book1D("charge","charge",3,-1.5,1.5) );
      h_hits.push_back( ibooker.book1D("hits", "number of hits per track", nintHit,minHit,maxHit ) );
      h_losthits.push_back( ibooker.book1D("losthits", "number of lost hits per track", nintHit,minHit,maxHit) );
      h_nchi2.push_back( ibooker.book1D("chi2", "normalized #chi^{2}", 200, 0, 20 ) );
      h_nchi2_prob.push_back( ibooker.book1D("chi2_prob", "normalized #chi^{2} probability",100,0,1));

      /// this are needed to calculate efficiency during tha harvesting for the automated validation
      h_recocent.push_back( ibooker.book1D("num_reco_cent","N of reco track vs cent",nintcent,centralityRanges) );
      h_assoccent.push_back( ibooker.book1D("num_assoc(simToReco)_cent","N of associated tracks (simToReco) vs cent",nintcent,centralityRanges) );
      h_assoc2cent.push_back( ibooker.book1D("num_assoc(recoToSim)_cent","N of associated (recoToSim) tracks vs cent",nintcent,centralityRanges) );
      h_simulcent.push_back( ibooker.book1D("num_simul_cent","N of simulated tracks vs cent",nintcent,centralityRanges) );

      h_recoeta.push_back( ibooker.book1D("num_reco_eta","N of reco track vs eta",nint,min,max) );
      h_assoceta.push_back( ibooker.book1D("num_assoc(simToReco)_eta","N of associated tracks (simToReco) vs eta",nint,min,max) );
      h_assoc2eta.push_back( ibooker.book1D("num_assoc(recoToSim)_eta","N of associated (recoToSim) tracks vs eta",nint,min,max) );
      h_simuleta.push_back( ibooker.book1D("num_simul_eta","N of simulated tracks vs eta",nint,min,max) );

      if (!_pTRanges.empty()) {
        h_recopT.push_back( ibooker.book1D("num_reco_pT","N of reco track vs pT",nintpT,pTRanges) );
        h_assocpT.push_back( ibooker.book1D("num_assoc(simToReco)_pT","N of associated tracks (simToReco) vs pT",nintpT,pTRanges) );
        h_assoc2pT.push_back( ibooker.book1D("num_assoc(recoToSim)_pT","N of associated (recoToSim) tracks vs pT",nintpT,pTRanges) );
        h_simulpT.push_back( ibooker.book1D("num_simul_pT","N of simulated tracks vs pT",nintpT,pTRanges) );
      } else {
        h_recopT.push_back( ibooker.book1D("num_reco_pT","N of reco track vs pT",nintpT,minpT,maxpT) );
        h_assocpT.push_back( ibooker.book1D("num_assoc(simToReco)_pT","N of associated tracks (simToReco) vs pT",nintpT,minpT,maxpT) );
        h_assoc2pT.push_back( ibooker.book1D("num_assoc(recoToSim)_pT","N of associated (recoToSim) tracks vs pT",nintpT,minpT,maxpT) );
        h_simulpT.push_back( ibooker.book1D("num_simul_pT","N of simulated tracks vs pT",nintpT,minpT,maxpT) );
      }

      for (int idx1=0; idx1<nintRapArr; idx1++) {
        for (int idx2=0; idx2<nintPtArrRes; idx2++) {
          h_recocent.push_back( ibooker.book1D(Form("num_reco_cent_Rap%.0f-%.0f_Pt%.0f-%.0f",rapArr[idx1]*10,rapArr[idx1+1]*10,ptArrRes[idx2]*10,ptArrRes[idx2+1]*10),"N of reco track vs cent",nintcent,centralityRanges) );
          h_assoccent.push_back( ibooker.book1D(Form("num_assoc(simToReco)_cent_Rap%.0f-%.0f_Pt%.0f-%.0f",rapArr[idx1]*10,rapArr[idx1+1]*10,ptArrRes[idx2]*10,ptArrRes[idx2+1]*10),"N of associated tracks (simToReco) vs cent",nintcent,centralityRanges) );
          h_assoc2cent.push_back( ibooker.book1D(Form("num_assoc(recoToSim)_cent_Rap%.0f-%.0f_Pt%.0f-%.0f",rapArr[idx1]*10,rapArr[idx1+1]*10,ptArrRes[idx2]*10,ptArrRes[idx2+1]*10),"N of associated (recoToSim) tracks vs cent",nintcent,centralityRanges) );
          h_simulcent.push_back( ibooker.book1D(Form("num_simul_cent_Rap%.0f-%.0f_Pt%.0f-%.0f",rapArr[idx1]*10,rapArr[idx1+1]*10,ptArrRes[idx2]*10,ptArrRes[idx2+1]*10),"N of simulated tracks vs cent",nintcent,centralityRanges) );
          if (!_pTRanges.empty()) {
            h_recopT.push_back( ibooker.book1D(Form("num_reco_pT_Rap%.0f-%.0f_Pt%.0f-%.0f",rapArr[idx1]*10,rapArr[idx1+1]*10,ptArrRes[idx2]*10,ptArrRes[idx2+1]*10),"N of reco track vs pT",nintpT,pTRanges) );
            h_assocpT.push_back( ibooker.book1D(Form("num_assoc(simToReco)_pT_Rap%.0f-%.0f_Pt%.0f-%.0f",rapArr[idx1]*10,rapArr[idx1+1]*10,ptArrRes[idx2]*10,ptArrRes[idx2+1]*10),"N of associated tracks (simToReco) vs pT",nintpT,pTRanges) );
            h_assoc2pT.push_back( ibooker.book1D(Form("num_assoc(recoToSim)_pT_Rap%.0f-%.0f_Pt%.0f-%.0f",rapArr[idx1]*10,rapArr[idx1+1]*10,ptArrRes[idx2]*10,ptArrRes[idx2+1]*10),"N of associated (recoToSim) tracks vs pT",nintpT,pTRanges) );
            h_simulpT.push_back( ibooker.book1D(Form("num_simul_pT_Rap%.0f-%.0f_Pt%.0f-%.0f",rapArr[idx1]*10,rapArr[idx1+1]*10,ptArrRes[idx2]*10,ptArrRes[idx2+1]*10),"N of simulated tracks vs pT",nintpT,pTRanges) );
          } else {
            h_recopT.push_back( ibooker.book1D(Form("num_reco_pT_Rap%.0f-%.0f_Pt%.0f-%.0f",rapArr[idx1]*10,rapArr[idx1+1]*10,ptArrRes[idx2]*10,ptArrRes[idx2+1]*10),"N of reco track vs pT",nintpT,ptArrRes[idx2]*10,ptArrRes[idx2+1]*10) );
            h_assocpT.push_back( ibooker.book1D(Form("num_assoc(simToReco)_pT_Rap%.0f-%.0f_Pt%.0f-%.0f",rapArr[idx1]*10,rapArr[idx1+1]*10,ptArrRes[idx2]*10,ptArrRes[idx2+1]*10),"N of associated tracks (simToReco) vs pT",nintpT,ptArrRes[idx2]*10,ptArrRes[idx2+1]*10) );
            h_assoc2pT.push_back( ibooker.book1D(Form("num_assoc(recoToSim)_pT_Rap%.0f-%.0f_Pt%.0f-%.0f",rapArr[idx1]*10,rapArr[idx1+1]*10,ptArrRes[idx2]*10,ptArrRes[idx2+1]*10),"N of associated (recoToSim) tracks vs pT",nintpT,ptArrRes[idx2]*10,ptArrRes[idx2+1]*10) );
            h_simulpT.push_back( ibooker.book1D(Form("num_simul_pT_Rap%.0f-%.0f_Pt%.0f-%.0f",rapArr[idx1]*10,rapArr[idx1+1]*10,ptArrRes[idx2]*10,ptArrRes[idx2+1]*10),"N of simulated tracks vs pT",nintpT,ptArrRes[idx2]*10,ptArrRes[idx2+1]*10) );
          }
        }
      }
      //
      h_recohit.push_back( ibooker.book1D("num_reco_hit","N of reco track vs hit",nintHit,minHit,maxHit) );
      h_assochit.push_back( ibooker.book1D("num_assoc(simToReco)_hit","N of associated tracks (simToReco) vs hit",nintHit,minHit,maxHit) );
      h_assoc2hit.push_back( ibooker.book1D("num_assoc(recoToSim)_hit","N of associated (recoToSim) tracks vs hit",nintHit,minHit,maxHit) );
      h_simulhit.push_back( ibooker.book1D("num_simul_hit","N of simulated tracks vs hit",nintHit,minHit,maxHit) );
      //
      h_recophi.push_back( ibooker.book1D("num_reco_phi","N of reco track vs phi",nintPhi,minPhi,maxPhi) );
      h_assocphi.push_back( ibooker.book1D("num_assoc(simToReco)_phi","N of associated tracks (simToReco) vs phi",nintPhi,minPhi,maxPhi) );
      h_assoc2phi.push_back( ibooker.book1D("num_assoc(recoToSim)_phi","N of associated (recoToSim) tracks vs phi",nintPhi,minPhi,maxPhi) );
      h_simulphi.push_back( ibooker.book1D("num_simul_phi","N of simulated tracks vs phi",nintPhi,minPhi,maxPhi) );

      h_recodxy.push_back( ibooker.book1D("num_reco_dxy","N of reco track vs dxy",nintDxy,minDxy,maxDxy) );
      h_assocdxy.push_back( ibooker.book1D("num_assoc(simToReco)_dxy","N of associated tracks (simToReco) vs dxy",nintDxy,minDxy,maxDxy) );
      h_assoc2dxy.push_back( ibooker.book1D("num_assoc(recoToSim)_dxy","N of associated (recoToSim) tracks vs dxy",nintDxy,minDxy,maxDxy) );
      h_simuldxy.push_back( ibooker.book1D("num_simul_dxy","N of simulated tracks vs dxy",nintDxy,minDxy,maxDxy) );
      
      h_recodz.push_back( ibooker.book1D("num_reco_dz","N of reco track vs dz",nintDz,minDz,maxDz) );
      h_assocdz.push_back( ibooker.book1D("num_assoc(simToReco)_dz","N of associated tracks (simToReco) vs dz",nintDz,minDz,maxDz) );
      h_assoc2dz.push_back( ibooker.book1D("num_assoc(recoToSim)_dz","N of associated (recoToSim) tracks vs dz",nintDz,minDz,maxDz) );
      h_simuldz.push_back( ibooker.book1D("num_simul_dz","N of simulated tracks vs dz",nintDz,minDz,maxDz) );

      h_assocvertpos.push_back( ibooker.book1D("num_assoc(simToReco)_vertpos","N of associated tracks (simToReco) vs transverse vert position",nintVertpos,minVertpos,maxVertpos) );
      h_simulvertpos.push_back( ibooker.book1D("num_simul_vertpos","N of simulated tracks vs transverse vert position",nintVertpos,minVertpos,maxVertpos) );

      h_assoczpos.push_back( ibooker.book1D("num_assoc(simToReco)_zpos","N of associated tracks (simToReco) vs z vert position",nintZpos,minZpos,maxZpos) );
      h_simulzpos.push_back( ibooker.book1D("num_simul_zpos","N of simulated tracks vs z vert position",nintZpos,minZpos,maxZpos) );

      /////////////////////////////////

      h_eta.push_back( ibooker.book1D("eta", "pseudorapidity residue", 1000, -0.1, 0.1 ) );
      h_pt.push_back( ibooker.book1D("pullPt", "pull of p_{t}", 100, -10, 10 ) );
      ///// Resolution plots in eta and pT
      h_eta2.push_back( dbe_->book1D("etaSimRecRes", "(Rec - Sim)/Sim #eta", 200, -2, 2 ) );
      h_pt2.push_back( dbe_->book1D("ptSimRecRes", "(Rec - Sim)/Sim p_{T}", 200, -2, 2 ) );
      for (int idx1=0; idx1<nintRapArrRes; idx1++) {
        for (int idx2=0; idx2<nintPtArrRes; idx2++) {
          h_eta2.push_back( dbe_->book1D(Form("etaSimRecRes_Rap%.0f-%.0f_Pt%.0f-%.0f",rapArrRes[idx1]*10,rapArrRes[idx1+1]*10,ptArrRes[idx2]*10,ptArrRes[idx2+1]*10), "(Rec - Sim)/Sim #eta", 200, -2, 2 ) );
          h_pt2.push_back( dbe_->book1D(Form("ptSimRecRes_Rap%.0f-%.0f_Pt%.0f-%.0f",rapArrRes[idx1]*10,rapArrRes[idx1+1]*10,ptArrRes[idx2]*10,ptArrRes[idx2+1]*10), "(Rec - Sim)/Sim p_{T}", 200, -2, 2 ) );
        }
      }
      ///// End of resolution plots in eta and pT
      h_pullTheta.push_back( ibooker.book1D("pullTheta","pull of #theta parameter",250,-25,25) );
      h_pullPhi.push_back( ibooker.book1D("pullPhi","pull of #phi parameter",250,-25,25) );
      h_pullDxy.push_back( ibooker.book1D("pullDxy","pull of dxy parameter",250,-25,25) );
      h_pullDz.push_back( ibooker.book1D("pullDz","pull of dz parameter",250,-25,25) );
      h_pullQoverp.push_back( ibooker.book1D("pullQoverp","pull of qoverp parameter",250,-25,25) );
      
      if (associators[ww]=="trackAssociatorByChi2"){
        h_assochi2.push_back( ibooker.book1D("assocChi2","track association #chi^{2}",1000000,0,100000) );
        h_assochi2_prob.push_back(ibooker.book1D("assocChi2_prob","probability of association #chi^{2}",100,0,1));
      } else if (associators[ww]=="trackAssociatorByHits"){
        h_assocFraction.push_back( ibooker.book1D("assocFraction","fraction of shared hits",200,0,2) );
        h_assocSharedHit.push_back(ibooker.book1D("assocSharedHit","number of shared hits",20,0,20));
      }
      
      chi2_vs_nhits.push_back( ibooker.book2D("chi2_vs_nhits","#chi^{2} vs nhits",25,0,25,100,0,10) );
      h_chi2meanhitsh.push_back( ibooker.bookProfile("chi2mean_vs_nhits","mean #chi^{2} vs nhits",25,0,25,100,0,10) );

      etares_vs_eta.push_back( ibooker.book2D("etares_vs_eta","etaresidue vs eta",nint,min,max,200,-0.1,0.1) );
      nrec_vs_nsim.push_back( ibooker.book2D("nrec_vs_nsim","nrec vs nsim",20,-0.5,19.5,20,-0.5,19.5) );

      chi2_vs_eta.push_back( ibooker.book2D("chi2_vs_eta","chi2_vs_eta",nint,min,max, 200, 0, 20 ));
      h_chi2meanh.push_back( ibooker.bookProfile("chi2mean","mean #chi^{2} vs #eta",nint,min,max, 200, 0, 20) );
      chi2_vs_phi.push_back( ibooker.book2D("chi2_vs_phi","#chi^{2} vs #phi",nintPhi,minPhi,maxPhi, 200, 0, 20 ) );
      h_chi2mean_vs_phi.push_back( ibooker.bookProfile("chi2mean_vs_phi","mean of #chi^{2} vs #phi",nintPhi,minPhi,maxPhi, 200, 0, 20) );

      nhits_vs_eta.push_back( ibooker.book2D("nhits_vs_eta","nhits vs eta",nint,min,max,nintHit,minHit,maxHit) );
      nDThits_vs_eta.push_back( ibooker.book2D("nDThits_vs_eta","# DT hits vs eta",nint,min,max,nintHit,minHit,maxHit) );
      nCSChits_vs_eta.push_back( ibooker.book2D("nCSChits_vs_eta","# CSC hits vs eta",nint,min,max,nintHit,minHit,maxHit) );
      nRPChits_vs_eta.push_back( ibooker.book2D("nRPChits_vs_eta","# RPC hits vs eta",nint,min,max,nintHit,minHit,maxHit) );

      h_DThits_eta.push_back( ibooker.bookProfile("DThits_eta","mean # DT hits vs eta",nint,min,max,nintHit,minHit,maxHit) );
      h_CSChits_eta.push_back( ibooker.bookProfile("CSChits_eta","mean # CSC hits vs eta",nint,min,max,nintHit,minHit,maxHit) );
      h_RPChits_eta.push_back( ibooker.bookProfile("RPChits_eta","mean # RPC hits vs eta",nint,min,max,nintHit,minHit,maxHit) );
      h_hits_eta.push_back( ibooker.bookProfile("hits_eta","mean #hits vs eta",nint,min,max,nintHit,minHit,maxHit) );
      nhits_vs_phi.push_back( ibooker.book2D("nhits_vs_phi","#hits vs #phi",nintPhi,minPhi,maxPhi,nintHit,minHit,maxHit) );
      h_hits_phi.push_back( ibooker.bookProfile("hits_phi","mean #hits vs #phi",nintPhi,minPhi,maxPhi, nintHit,minHit,maxHit) );

      nlosthits_vs_eta.push_back( ibooker.book2D("nlosthits_vs_eta","nlosthits vs eta",nint,min,max,nintHit,minHit,maxHit) );
      h_losthits_eta.push_back( ibooker.bookProfile("losthits_eta","losthits_eta",nint,min,max,nintHit,minHit,maxHit) );

      ptres_vs_eta.push_back(ibooker.book2D("ptres_vs_eta","ptres_vs_eta",nint,min,max, ptRes_nbin, ptRes_rangeMin, ptRes_rangeMax));
      ptres_vs_phi.push_back( ibooker.book2D("ptres_vs_phi","p_{t} res vs #phi",nintPhi,minPhi,maxPhi, ptRes_nbin, ptRes_rangeMin, ptRes_rangeMax));
      if (!_pTRanges.empty()) ptres_vs_pt.push_back(ibooker.book2D("ptres_vs_pt","ptres_vs_pt",nintpT,pTRanges[0],pTRanges[nintpT], ptRes_nbin, ptRes_rangeMin, ptRes_rangeMax));
      else ptres_vs_pt.push_back(ibooker.book2D("ptres_vs_pt","ptres_vs_pt",nintpT,minpT,maxpT, ptRes_nbin, ptRes_rangeMin, ptRes_rangeMax));

      cotThetares_vs_eta.push_back(ibooker.book2D("cotThetares_vs_eta","cotThetares_vs_eta",nint,min,max,cotThetaRes_nbin, cotThetaRes_rangeMin, cotThetaRes_rangeMax));
      if (!_pTRanges.empty()) cotThetares_vs_pt.push_back(ibooker.book2D("cotThetares_vs_pt","cotThetares_vs_pt",nintpT,pTRanges[0],pTRanges[nintpT], cotThetaRes_nbin, cotThetaRes_rangeMin, cotThetaRes_rangeMax));
      else cotThetares_vs_pt.push_back(ibooker.book2D("cotThetares_vs_pt","cotThetares_vs_pt",nintpT,minpT,maxpT, cotThetaRes_nbin, cotThetaRes_rangeMin, cotThetaRes_rangeMax));

      phires_vs_eta.push_back(ibooker.book2D("phires_vs_eta","phires_vs_eta",nint,min,max, phiRes_nbin, phiRes_rangeMin, phiRes_rangeMax));
      if (!_pTRanges.empty()) phires_vs_pt.push_back(ibooker.book2D("phires_vs_pt","phires_vs_pt",nintpT,pTRanges[0],pTRanges[nintpT], phiRes_nbin, phiRes_rangeMin, phiRes_rangeMax));
      else phires_vs_pt.push_back(ibooker.book2D("phires_vs_pt","phires_vs_pt",nintpT,minpT,maxpT, phiRes_nbin, phiRes_rangeMin, phiRes_rangeMax));
      phires_vs_phi.push_back(ibooker.book2D("phires_vs_phi","#phi res vs #phi",nintPhi,minPhi,maxPhi,phiRes_nbin, phiRes_rangeMin, phiRes_rangeMax));

      dxyres_vs_eta.push_back(ibooker.book2D("dxyres_vs_eta","dxyres_vs_eta",nint,min,max,dxyRes_nbin, dxyRes_rangeMin, dxyRes_rangeMax));
      if (!_pTRanges.empty()) dxyres_vs_pt.push_back( ibooker.book2D("dxyres_vs_pt","dxyres_vs_pt",nintpT,pTRanges[0],pTRanges[nintpT], dxyRes_nbin, dxyRes_rangeMin, dxyRes_rangeMax));
      else dxyres_vs_pt.push_back( ibooker.book2D("dxyres_vs_pt","dxyres_vs_pt",nintpT,minpT,maxpT,dxyRes_nbin, dxyRes_rangeMin, dxyRes_rangeMax));

      dzres_vs_eta.push_back(ibooker.book2D("dzres_vs_eta","dzres_vs_eta",nint,min,max,dzRes_nbin, dzRes_rangeMin, dzRes_rangeMax));
      if (!_pTRanges.empty()) dzres_vs_pt.push_back(ibooker.book2D("dzres_vs_pt","dzres_vs_pt",nintpT,pTRanges[0],pTRanges[nintpT], dzRes_nbin, dzRes_rangeMin, dzRes_rangeMax));
      else dzres_vs_pt.push_back(ibooker.book2D("dzres_vs_pt","dzres_vs_pt",nintpT,minpT,maxpT,dzRes_nbin, dzRes_rangeMin, dzRes_rangeMax));

      ptmean_vs_eta_phi.push_back(ibooker.bookProfile2D("ptmean_vs_eta_phi","mean p_{t} vs #eta and #phi",nintPhi,minPhi,maxPhi,nint,min,max,1000,0,1000));
      phimean_vs_eta_phi.push_back(ibooker.bookProfile2D("phimean_vs_eta_phi","mean #phi vs #eta and #phi",nintPhi,minPhi,maxPhi,nint,min,max,nintPhi,minPhi,maxPhi));

      //pulls of track params vs eta: to be used with fitslicesytool
      dxypull_vs_eta.push_back(ibooker.book2D("dxypull_vs_eta","dxypull_vs_eta",nint,min,max,100,-10,10));
      ptpull_vs_eta.push_back(ibooker.book2D("ptpull_vs_eta","ptpull_vs_eta",nint,min,max,100,-10,10));
      dzpull_vs_eta.push_back(ibooker.book2D("dzpull_vs_eta","dzpull_vs_eta",nint,min,max,100,-10,10));
      phipull_vs_eta.push_back(ibooker.book2D("phipull_vs_eta","phipull_vs_eta",nint,min,max,100,-10,10));
      thetapull_vs_eta.push_back(ibooker.book2D("thetapull_vs_eta","thetapull_vs_eta",nint,min,max,100,-10,10));

      //pulls of track params vs phi
      ptpull_vs_phi.push_back(ibooker.book2D("ptpull_vs_phi","p_{t} pull vs #phi",nintPhi,minPhi,maxPhi,100,-10,10));
      phipull_vs_phi.push_back(ibooker.book2D("phipull_vs_phi","#phi pull vs #phi",nintPhi,minPhi,maxPhi,100,-10,10));
      thetapull_vs_phi.push_back(ibooker.book2D("thetapull_vs_phi","#theta pull vs #phi",nintPhi,minPhi,maxPhi,100,-10,10));

      nrecHit_vs_nsimHit_sim2rec.push_back( ibooker.book2D("nrecHit_vs_nsimHit_sim2rec","nrecHit vs nsimHit (Sim2RecAssoc)",nintHit,minHit,maxHit, nintHit,minHit,maxHit ));
      nrecHit_vs_nsimHit_rec2sim.push_back( ibooker.book2D("nrecHit_vs_nsimHit_rec2sim","nrecHit vs nsimHit (Rec2simAssoc)",nintHit,minHit,maxHit, nintHit,minHit,maxHit ));
      
      if (MABH) {
        h_PurityVsQuality.push_back( ibooker.book2D("PurityVsQuality","Purity vs Quality (MABH)",20,0.01,1.01,20,0.01,1.01) );
        h_assoceta_Quality05.push_back( ibooker.book1D("num_assoc(simToReco)_eta_Q05","N of associated tracks (simToReco) vs eta (Quality>0.5)",nint,min,max) );
        h_assoceta_Quality075.push_back( ibooker.book1D("num_assoc(simToReco)_eta_Q075","N of associated tracks (simToReco) vs eta (Quality>0.75)",nint,min,max) );
        if (!_pTRanges.empty()) {
          h_assocpT_Quality05.push_back( ibooker.book1D("num_assoc(simToReco)_pT_Q05","N of associated tracks (simToReco) vs pT (Quality>0.5)",nintpT,pTRanges) );
          h_assocpT_Quality075.push_back( ibooker.book1D("num_assoc(simToReco)_pT_Q075","N of associated tracks (simToReco) vs pT (Quality>0.75)",nintpT,pTRanges) );
        } else {
          h_assocpT_Quality05.push_back( ibooker.book1D("num_assoc(simToReco)_pT_Q05","N of associated tracks (simToReco) vs pT (Quality>0.5)",nintpT,minpT,maxpT) );
          h_assocpT_Quality075.push_back( ibooker.book1D("num_assoc(simToReco)_pT_Q075","N of associated tracks (simToReco) vs pT (Quality>0.75)",nintpT,minpT,maxpT) );
        }
        h_assocphi_Quality05.push_back( ibooker.book1D("num_assoc(simToReco)_phi_Q05","N of associated tracks (simToReco) vs phi (Quality>0.5)",nintPhi,minPhi,maxPhi) );
        h_assocphi_Quality075.push_back( ibooker.book1D("num_assoc(simToReco)_phi_Q075","N of associated tracks (simToReco) vs phi (Quality>0.75)",nintPhi,minPhi,maxPhi) );
      }

      if(useLogPt){
        BinLogX(dzres_vs_pt[j]->getTH2F());
        BinLogX(dxyres_vs_pt[j]->getTH2F());
        BinLogX(phires_vs_pt[j]->getTH2F());
        BinLogX(cotThetares_vs_pt[j]->getTH2F());
        BinLogX(ptres_vs_pt[j]->getTH2F());
        BinLogX(h_recopT[j]->getTH1F());
        BinLogX(h_assocpT[j]->getTH1F());
        BinLogX(h_assoc2pT[j]->getTH1F());
        BinLogX(h_simulpT[j]->getTH1F());
        if (MABH) {
          BinLogX(h_assocpT_Quality05[j]->getTH1F());
          BinLogX(h_assocpT_Quality075[j]->getTH1F());
        }
        j++;
      }

    }
  }
}

void MuonTrackValidator::analyze(const edm::Event& event, const edm::EventSetup& setup){
  using namespace reco;

  edm::LogInfo("MuonTrackValidator") << "\n====================================================" << "\n"
                                     << "Analyzing new event" << "\n"
                                     << "====================================================\n" << "\n";

  edm::ESHandle<ParametersDefinerForTP> Lhc_parametersDefinerTP;
  std::unique_ptr<ParametersDefinerForTP> Cosmic_parametersDefinerTP;

  if(parametersDefiner=="LhcParametersDefinerForTP") {
    setup.get<TrackAssociatorRecord>().get(parametersDefiner, Lhc_parametersDefinerTP);
  }
  else if(parametersDefiner=="CosmicParametersDefinerForTP") {
    edm::ESHandle<CosmicParametersDefinerForTP>  _Cosmic_parametersDefinerTP;
    setup.get<TrackAssociatorRecord>().get(parametersDefiner, _Cosmic_parametersDefinerTP);

     //Since we modify the object, we must clone it
    Cosmic_parametersDefinerTP = _Cosmic_parametersDefinerTP->clone();

    edm::Handle<SimHitTPAssociationProducer::SimHitTPAssociationList> simHitsTPAssoc;
    //warning: make sure the TP collection used in the map is the same used here
    event.getByToken(_simHitTpMapTag,simHitsTPAssoc);
    Cosmic_parametersDefinerTP->initEvent(simHitsTPAssoc);
    cosmictpSelector.initEvent(simHitsTPAssoc);
  }
  else {
    edm::LogError("MuonTrackValidator")
      << "Unexpected label: parametersDefiner = "<< parametersDefiner.c_str() << "\n";
  }

  edm::Handle<TrackingParticleCollection> TPCollectionHeff ;
  event.getByToken(tp_effic_Token,TPCollectionHeff);
  TrackingParticleCollection const & tPCeff = *(TPCollectionHeff.product());

  edm::Handle<TrackingParticleCollection> TPCollectionHfake ;
  event.getByToken(tp_fake_Token,TPCollectionHfake);
  
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  event.getByToken(bsSrc_Token,recoBeamSpotHandle);
  reco::BeamSpot bs = *recoBeamSpotHandle;

  std::vector<const reco::TrackToTrackingParticleAssociator*> associator;
  if (UseAssociators) {
    edm::Handle<reco::TrackToTrackingParticleAssociator> theAssociator;
    for (unsigned int w=0;w<associators.size();w++) {
      event.getByLabel(associators[w],theAssociator);
      associator.push_back( theAssociator.product() );
    }
  }
  
  // centrality settings loading
  edm::Handle<int> cbin_;
  event.getByToken(_centralityBinTag,cbin_);
  centBin = *cbin_;

  double weighting = findNcoll(centBin) * pTWeight;
  edm::LogVerbatim("MuonTrackValidator") << "centrality: " << centBin << "\t"
                                         << "NcollWeight: " << findNcoll(centBin) << "\n";

  int w=0;
  for (unsigned int ww=0;ww<associators.size();ww++){
    for (unsigned int www=0;www<label.size();www++){
      //
      //get collections from the event
      //
      edm::Handle<edm::View<Track> > trackCollection;

      reco::RecoToSimCollection recSimColl;
      reco::SimToRecoCollection simRecColl;
      unsigned int trackCollectionSize = 0;

      if(!event.getByToken(track_Collection_Token[www], trackCollection)&&ignoremissingtkcollection_) {

        recSimColl.post_insert();
        simRecColl.post_insert();
        
      } else {
        trackCollectionSize = trackCollection->size();

        //associate tracks
        if(UseAssociators){
          edm::LogVerbatim("MuonTrackValidator") << "Analyzing "
                                                 << label[www].process()<<":"
                                                 << label[www].label()<<":"
                                                 << label[www].instance()<<" with "
                                                 << associators[ww].c_str() <<"\n";
          
          LogTrace("MuonTrackValidator") << "Calling associateRecoToSim method" << "\n";
          recSimColl=associator[ww]->associateRecoToSim(trackCollection,
                                                        TPCollectionHfake);
          LogTrace("MuonTrackValidator") << "Calling associateSimToReco method" << "\n";
          simRecColl=associator[ww]->associateSimToReco(trackCollection,
                                                        TPCollectionHeff);
        } else {
          edm::LogVerbatim("MuonTrackValidator") << "Analyzing "
                                                 << label[www].process()<<":"
                                                 << label[www].label()<<":"
                                                 << label[www].instance()<<" with "
                                                 << associatormap.process()<<":"
                                                 << associatormap.label()<<":"
                                                 << associatormap.instance()<<"\n";
          
          Handle<reco::SimToRecoCollection > simtorecoCollectionH;
          event.getByToken(simToRecoCollection_Token,simtorecoCollectionH);
          simRecColl= *(simtorecoCollectionH.product());
          
          Handle<reco::RecoToSimCollection > recotosimCollectionH;
          event.getByToken(recoToSimCollection_Token,recotosimCollectionH);
          recSimColl= *(recotosimCollectionH.product());
        }
        
      }
      
      //
      //fill simulation histograms
      //compute number of tracks per eta interval
      //
      edm::LogVerbatim("MuonTrackValidator") << "\n# of TrackingParticles: " << tPCeff.size() << "\n";
      int ats = 0;
      int st = 0;
      for (TrackingParticleCollection::size_type i=0; i<tPCeff.size(); i++){
        bool TP_is_matched = false;
        double quality = 0.;
        bool Quality05 = false;
        bool Quality075 = false;
        
        TrackingParticleRef tpr(TPCollectionHeff, i);
        TrackingParticle* tp = const_cast<TrackingParticle*>(tpr.get());
        
        TrackingParticle::Vector momentumTP;
        TrackingParticle::Point vertexTP;
        double dxySim = 0;
        double dzSim = 0;
        
        //If the TrackingParticle is collision-like, get the momentum and vertex at production state
        //and the impact parameters w.r.t. PCA
        if(parametersDefiner=="LhcParametersDefinerForTP")
        {
          LogTrace("MuonTrackValidator") <<"TrackingParticle "<< i;
          if(! tpSelector(*tp)) continue;
          momentumTP = tp->momentum();
          vertexTP = tp->vertex();
          TrackingParticle::Vector momentum = Lhc_parametersDefinerTP->momentum(event,setup,tpr);
          TrackingParticle::Point vertex = Lhc_parametersDefinerTP->vertex(event,setup,tpr);
          dxySim = (-vertex.x()*sin(momentum.phi())+vertex.y()*cos(momentum.phi()));
          dzSim = vertex.z() - (vertex.x()*momentum.x()+vertex.y()*momentum.y()) /
                                sqrt(momentum.perp2()) * momentum.z()/sqrt(momentum.perp2());
        }
        //for cosmics get the momentum and vertex at PCA
        else if(parametersDefiner=="CosmicParametersDefinerForTP")
        {
          edm::LogVerbatim("MuonTrackValidator") <<"TrackingParticle "<< i;
          if(! cosmictpSelector(tpr,&bs,event,setup)) continue;       
          momentumTP = Cosmic_parametersDefinerTP->momentum(event,setup,tpr);
          vertexTP = Cosmic_parametersDefinerTP->vertex(event,setup,tpr);
          dxySim = (-vertexTP.x()*sin(momentumTP.phi())+vertexTP.y()*cos(momentumTP.phi()));
          dzSim = vertexTP.z() - (vertexTP.x()*momentumTP.x()+vertexTP.y()*momentumTP.y()) /
                                  sqrt(momentumTP.perp2()) * momentumTP.z()/sqrt(momentumTP.perp2());
        }
        edm::LogVerbatim("MuonTrackValidator") <<"--------------------Selected TrackingParticle #"<<tpr.key();
        edm::LogVerbatim("MuonTrackValidator") <<"momentumTP: pt = "<<sqrt(momentumTP.perp2())<<", pz = "<<momentumTP.z() 
                                               <<", \t vertexTP: radius = "<<sqrt(vertexTP.perp2())<<  ",  z = "<<vertexTP.z() <<"\n";
        st++;

        h_ptSIM[w]->Fill(sqrt(momentumTP.perp2()),weighting);
        h_etaSIM[w]->Fill(momentumTP.eta(),weighting);
        h_vertposSIM[w]->Fill(sqrt(vertexTP.perp2()),weighting);
        
        std::vector<std::pair<RefToBase<Track>, double> > rt;
        if(simRecColl.find(tpr) != simRecColl.end()) {
          rt = (std::vector<std::pair<RefToBase<Track>, double> >) simRecColl[tpr];
          if (rt.size()!=0) {
            RefToBase<Track> assoc_recoTrack = rt.begin()->first;
            edm::LogVerbatim("MuonTrackValidator")<<"-----------------------------associated Track #"<<assoc_recoTrack.key();
            TP_is_matched = true;
            ats++;
            quality = rt.begin()->second;
            edm::LogVerbatim("MuonTrackValidator") << "TrackingParticle #" <<tpr.key()
                                                   << " with pt=" << sqrt(momentumTP.perp2())
                                                   << " associated with quality:" << quality <<"\n";
            if (MABH) {
              if (quality > 0.75) {
                Quality075 = true;
                Quality05 = true;
              }
              else if (quality > 0.5) {
                Quality05 = true;
              }
            }
          }
        }else{
          edm::LogVerbatim("MuonTrackValidator")
            << "TrackingParticle #" << tpr.key()
            << " with pt,eta,phi: "
            << sqrt(momentumTP.perp2()) << " , "
            << momentumTP.eta() << " , "
            << momentumTP.phi() << " , "
            << " NOT associated to any reco::Track" << "\n";
        }
        
        // merged over all eta
        h_simulcent[w]->Fill(centBin,weighting);
        if (TP_is_matched) h_assoccent[w]->Fill(centBin,weighting);
        for (int idx1=0; idx1<nintRapArr; idx1++) {
          if (getEta(momentumTP.eta())>rapArr[idx1] && getEta(momentumTP.eta())<rapArr[idx1+1]) {
            for (int idx2=0; idx2<nintPtArrRes; idx2++) {
              if (getPt(sqrt(momentumTP.perp2()))>ptArrRes[idx2] && getPt(sqrt(momentumTP.perp2()))<ptArrRes[idx2+1]) {
                int idx = nintPtArrRes*idx1 + idx2+1;
                // h_simulcent[0] is for merged case: differential ones should start from idx=1
                h_simulcent[idx]->Fill(centBin,weighting);
                if (TP_is_matched) h_assoccent[idx]->Fill(centBin,weighting);
              }
            }
          }
        }

        for (unsigned int f=0; f<etaintervals[w].size()-1; f++){
          if (getEta(momentumTP.eta())>etaintervals[w][f]&&
              getEta(momentumTP.eta())<etaintervals[w][f+1]) {
            totSIMeta[w][f] = totSIMeta[w][f] + weighting;
            if (TP_is_matched) {
              totASSeta[w][f] = totASSeta[w][f] + weighting;

              if (MABH) {
                if (Quality075) {
                  totASSeta_Quality075[w][f] = totASSeta_Quality075[w][f] + weighting;
                  totASSeta_Quality05[w][f] = totASSeta_Quality05[w][f] + weighting;
                }
                else if (Quality05) {
                  totASSeta_Quality05[w][f] = totASSeta_Quality05[w][f] + weighting;
                }
              }
            }
          }
        } // END for (unsigned int f=0; f<etaintervals[w].size()-1; f++)
        
        for (unsigned int f=0; f<phiintervals[w].size()-1; f++){
          if (momentumTP.phi() > phiintervals[w][f]&&
              momentumTP.phi() <phiintervals[w][f+1]) {
            totSIM_phi[w][f] = totSIM_phi[w][f] + weighting;
            if (TP_is_matched) {
              totASS_phi[w][f] = totASS_phi[w][f] + weighting;
              
              if (MABH) {
                if (Quality075) {
                  totASS_phi_Quality075[w][f] = totASS_phi_Quality075[w][f] + weighting;
                  totASS_phi_Quality05[w][f] = totASS_phi_Quality05[w][f] + weighting;
                }
                else if (Quality05) {
                  totASS_phi_Quality05[w][f] = totASS_phi_Quality05[w][f] + weighting;
                }
              }
            }
          }
        } // END for (unsigned int f=0; f<phiintervals[w].size()-1; f++)
        
        // merged over all eta
        for (unsigned int f=0; f<pTintervals[w].size()-1; f++){
          if (getPt(sqrt(momentumTP.perp2()))>pTintervals[w][f]&&
              getPt(sqrt(momentumTP.perp2()))<pTintervals[w][f+1]) {
            totSIMpT[w][f] = totSIMpT[w][f] + weighting;
            if (TP_is_matched) {
              totASSpT[w][f] = totASSpT[w][f] + weighting;
              
              if (MABH) {
                if (Quality075) {
                  totASSpT_Quality075[w][f] = totASSpT_Quality075[w][f] + weighting;
                  totASSpT_Quality05[w][f] = totASSpT_Quality05[w][f] + weighting;
                }
                else if (Quality05) { 
                  totASSpT_Quality05[w][f] = totASSpT_Quality05[w][f] + weighting;
                }
              }
            }
          }
        } // END for (unsigned int f=0; f<pTintervals[w].size()-1; f++)


        for (int idx1=0; idx1<nintRapArr; idx1++) {
          if (getEta(momentumTP.eta())>rapArr[idx1] && getEta(momentumTP.eta())<rapArr[idx1+1]) {
            for (int idx2=0; idx2<nintPtArrRes; idx2++) {
              if (getPt(sqrt(momentumTP.perp2()))>ptArrRes[idx2] && getPt(sqrt(momentumTP.perp2()))<ptArrRes[idx2+1]) {
                int idx = nintPtArrRes*idx1 + idx2+1;
                // h_simulpT[0] is for merged case: differential ones should start from idx=1
                h_simulpT[idx]->Fill(getPt(sqrt(momentumTP.perp2())),weighting);
                if (TP_is_matched) h_assocpT[idx]->Fill(getPt(sqrt(momentumTP.perp2())),weighting);
              }
            }
          }
        }

        for (unsigned int f=0; f<dxyintervals[w].size()-1; f++){
          if (dxySim>dxyintervals[w][f]&&
              dxySim<dxyintervals[w][f+1]) {
            totSIM_dxy[w][f] = totSIM_dxy[w][f] + weighting;
            if (TP_is_matched) {
              totASS_dxy[w][f] = totASS_dxy[w][f] + weighting;
            }
          }
        } // END for (unsigned int f=0; f<dxyintervals[w].size()-1; f++)

        for (unsigned int f=0; f<dzintervals[w].size()-1; f++){
          if (dzSim>dzintervals[w][f]&&
              dzSim<dzintervals[w][f+1]) {
            totSIM_dz[w][f] = totSIM_dz[w][f] + weighting;
            if (TP_is_matched) {
              totASS_dz[w][f] = totASS_dz[w][f] + weighting;
            }
          }
        } // END for (unsigned int f=0; f<dzintervals[w].size()-1; f++)

        for (unsigned int f=0; f<vertposintervals[w].size()-1; f++){
          if (sqrt(vertexTP.perp2())>vertposintervals[w][f]&&
              sqrt(vertexTP.perp2())<vertposintervals[w][f+1]) {
            totSIM_vertpos[w][f] = totSIM_vertpos[w][f] + weighting;
            if (TP_is_matched) {
              totASS_vertpos[w][f] = totASS_vertpos[w][f] + weighting;
            }
          }
        } // END for (unsigned int f=0; f<vertposintervals[w].size()-1; f++)

        for (unsigned int f=0; f<zposintervals[w].size()-1; f++){
          if (vertexTP.z()>zposintervals[w][f]&&
              vertexTP.z()<zposintervals[w][f+1]) {
            totSIM_zpos[w][f] = totSIM_zpos[w][f] + weighting;
            if (TP_is_matched) {
              totASS_zpos[w][f] = totASS_zpos[w][f] + weighting;
            }
          }
        } // END for (unsigned int f=0; f<zposintervals[w].size()-1; f++)
        
        int nSimHits = 0;
        if (usetracker && usemuon) {
          nSimHits= tpr.get()->numberOfHits();
        }
        else if (!usetracker && usemuon) {
          nSimHits= tpr.get()->numberOfHits() - tpr.get()->numberOfTrackerHits();
        }
        else if (usetracker && !usemuon) {
          nSimHits=tpr.get()->numberOfTrackerHits();
        }
        
        int tmp = std::min(nSimHits,int(maxHit-1));
        edm::LogVerbatim("MuonTrackValidator") << "\t N simhits = "<< nSimHits<<"\n";

        totSIM_hit[w][tmp] = totSIM_hit[w][tmp] + weighting;
        if (TP_is_matched) totASS_hit[w][tmp] = totASS_hit[w][tmp] + weighting;

        if (TP_is_matched)
        {
          RefToBase<Track> assoctrack = rt.begin()->first;
          nrecHit_vs_nsimHit_sim2rec[w]->Fill( assoctrack->numberOfValidHits(),nSimHits);
        }
      } // End for (TrackingParticleCollection::size_type i=0; i<tPCeff.size(); i++)
      if (st!=0) h_tracksSIM[w]->Fill(st,weighting);
      // END for simToReco matching
      
      //
      //fill reconstructed track histograms
      //
      edm::LogVerbatim("MuonTrackValidator") << "\n# of reco::Tracks with "
                                             << label[www].process()<<":"
                                             << label[www].label()<<":"
                                             << label[www].instance()
                                             << ": " << trackCollectionSize << "\n";
      int at = 0;
      int rT = 0;
      for(edm::View<Track>::size_type i=0; i<trackCollectionSize; ++i){
        bool Track_is_matched = false;
        RefToBase<Track> track(trackCollection, i);
        rT++;

        std::vector<std::pair<TrackingParticleRef, double> > tp;
        TrackingParticleRef tpr;

        // new logic (bidirectional)
        if (BiDirectional_RecoToSim_association) {      
          edm::LogVerbatim("MuonTrackValidator")<<"----------------------------------------Track #"<< track.key();

          if(recSimColl.find(track) != recSimColl.end()) {
            tp = recSimColl[track];     
            if (tp.size() != 0) {
              tpr = tp.begin()->first;  
              // RtS and StR must associate the same pair !
              if(simRecColl.find(tpr) != simRecColl.end()) {
                std::vector<std::pair<RefToBase<Track>, double> > track_checkback = simRecColl[tpr];
                RefToBase<Track> assoc_track_checkback;
                assoc_track_checkback = track_checkback.begin()->first;

                if ( assoc_track_checkback.key() == track.key() ) {
                  edm::LogVerbatim("MuonTrackValidator")<<"------------------associated TrackingParticle #"<<tpr.key();
                  Track_is_matched = true;
                  at++;
                  double Purity = tp.begin()->second;
                  double Quality = track_checkback.begin()->second;
                  edm::LogVerbatim("MuonTrackValidator") << "reco::Track #" << track.key() << " with pt=" << track->pt()
                                                         << " associated with quality:" << Purity <<"\n";
                  if (MABH) h_PurityVsQuality[w]->Fill(Quality,Purity,weighting);
                }
              }
            }
          }
          
          if (!Track_is_matched)
            edm::LogVerbatim("MuonTrackValidator")
              << "reco::Track #" << track.key() << " with pt=" << track->pt() << " NOT associated to any TrackingParticle" << "\n";
        }
        // old logic (bugged for collision scenario, still valid for cosmics 2 legs reco)
        else {
          if(recSimColl.find(track) != recSimColl.end()){
            tp = recSimColl[track];
            if (tp.size()!=0) {
              Track_is_matched = true;
              tpr = tp.begin()->first;
              at++;
              edm::LogVerbatim("MuonTrackValidator") << "reco::Track #" << track.key() << " with pt=" << track->pt()
                                                     << " associated with quality:" << tp.begin()->second <<"\n";
            }
          } else {
            edm::LogVerbatim("MuonTrackValidator") << "reco::Track #" << track.key() << " with pt=" << track->pt()
                                                   << " NOT associated to any TrackingParticle" << "\n";        
          }
        }
        
        // merged over all eta
        h_recocent[w]->Fill(centBin,weighting);
        if (Track_is_matched) h_assoc2cent[w]->Fill(centBin,weighting);
        for (int idx1=0; idx1<nintRapArr; idx1++) {
          if (getEta(track->momentum().eta())>rapArr[idx1] && getEta(track->momentum().eta())<rapArr[idx1+1]) {
            for (int idx2=0; idx2<nintPtArrRes; idx2++) {
              if (track->pt()>ptArrRes[idx2] && track->pt()<ptArrRes[idx2+1]) {
                int idx = nintPtArrRes*idx1 + idx2+1;
                // h_recocent[0] is for merged case: differential ones should start from idx=1
                h_recocent[idx]->Fill(centBin,weighting);
                if (Track_is_matched) h_assoc2cent[idx]->Fill(centBin,weighting);
              }
            }
          }
        }
        
        //Compute fake rate vs eta
        for (unsigned int f=0; f<etaintervals[w].size()-1; f++){
          if (getEta(track->momentum().eta())>etaintervals[w][f]&&
              getEta(track->momentum().eta())<etaintervals[w][f+1]) {
            totRECeta[w][f] = totRECeta[w][f] + weighting; 
            if (Track_is_matched) {
              totASS2eta[w][f] = totASS2eta[w][f] + weighting;
            }
          }
        } // End for (unsigned int f=0; f<etaintervals[w].size()-1; f++)

        for (unsigned int f=0; f<phiintervals[w].size()-1; f++){
          if (track->momentum().phi()>phiintervals[w][f]&&
              track->momentum().phi()<phiintervals[w][f+1]) {
            totREC_phi[w][f] = totREC_phi[w][f] + weighting; 
            if (Track_is_matched) {
              totASS2_phi[w][f] = totASS2_phi[w][f] + weighting;
            }
          }
        } // End for (unsigned int f=0; f<phiintervals[w].size()-1; f++)

        for (unsigned int f=0; f<pTintervals[w].size()-1; f++){
          if (getPt(sqrt(track->momentum().perp2()))>pTintervals[w][f]&&
              getPt(sqrt(track->momentum().perp2()))<pTintervals[w][f+1]) {
            totRECpT[w][f] = totRECpT[w][f] + weighting; 
            if (Track_is_matched) {
              totASS2pT[w][f] = totASS2pT[w][f] + weighting;
            }         
          }
        } // End for (unsigned int f=0; f<pTintervals[w].size()-1; f++)

        for (int idx1=0; idx1<nintRapArr; idx1++) {
          if (getEta(track->momentum().eta())>rapArr[idx1] && getEta(track->momentum().eta())<rapArr[idx1+1]) {
            for (int idx2=0; idx2<nintPtArrRes; idx2++) {
              if (getPt(sqrt(track->momentum().perp2()))>ptArrRes[idx2] && getPt(sqrt(track->momentum().perp2()))<ptArrRes[idx2+1]) {
                int idx = nintPtArrRes*idx1 + idx2+1;
                // h_simulpT[0] is for merged case: differential ones should start from idx=1
                h_recopT[idx]->Fill(getPt(sqrt(track->momentum().perp2())),weighting);
                if (Track_is_matched) h_assoc2pT[idx]->Fill(getPt(sqrt(track->momentum().perp2())),weighting);
              }
            }
          }
        } // End for (unsigned int f=0; f<pTintervals[w].size()-1; f++)

        for (unsigned int f=0; f<dxyintervals[w].size()-1; f++){
          if (track->dxy(bs.position())>dxyintervals[w][f]&&
              track->dxy(bs.position())<dxyintervals[w][f+1]) {
            totREC_dxy[w][f] = totREC_dxy[w][f] + weighting; 
            if (Track_is_matched) {
              totASS2_dxy[w][f] = totASS2_dxy[w][f] + weighting;
            }
          }
        } // End for (unsigned int f=0; f<dxyintervals[w].size()-1; f++)

        for (unsigned int f=0; f<dzintervals[w].size()-1; f++){
          if (track->dz(bs.position())>dzintervals[w][f]&&
              track->dz(bs.position())<dzintervals[w][f+1]) {
            totREC_dz[w][f] = totREC_dz[w][f] + weighting; 
            if (Track_is_matched) {
              totASS2_dz[w][f] = totASS2_dz[w][f] + weighting;
            }
          }
        } // End for (unsigned int f=0; f<dzintervals[w].size()-1; f++)

        int tmp = std::min((int)track->found(),int(maxHit-1));
        totREC_hit[w][tmp] = totREC_hit[w][tmp] + weighting;
        if (Track_is_matched) totASS2_hit[w][tmp] = totASS2_hit[w][tmp] + weighting;

        edm::LogVerbatim("MuonTrackValidator") << "\t N valid rechits = "<< (int)track->found() <<"\n";

        // Fill other histos
        TrackingParticle* tpp = const_cast<TrackingParticle*>(tpr.get());
        // TrackingParticle parameters at point of closest approach to the beamline
        TrackingParticle::Vector momentumTP;
        TrackingParticle::Point vertexTP;
        
        if (parametersDefiner=="LhcParametersDefinerForTP") {
          // following reco plots are made only from tracks associated to selected signal TPs 
          if (! (Track_is_matched && tpSelector(*tpp)) ) continue;
          else {
            momentumTP = Lhc_parametersDefinerTP->momentum(event,setup,tpr) ;
            vertexTP = Lhc_parametersDefinerTP->vertex(event,setup,tpr);
          }
        }
        else if (parametersDefiner=="CosmicParametersDefinerForTP") {
          // following reco plots are made only from tracks associated to selected signal TPs 
          if (! (Track_is_matched && cosmictpSelector(tpr,&bs,event,setup)) ) continue;
          else {
            momentumTP = Cosmic_parametersDefinerTP->momentum(event,setup,tpr) ;
            vertexTP = Cosmic_parametersDefinerTP->vertex(event,setup,tpr);         
          }
        }
        
        if (associators[ww]=="trackAssociatorByChi2"){
          //association chi2
          double assocChi2 = -tp.begin()->second;//in association map is stored -chi2
          h_assochi2[www]->Fill(assocChi2,weighting);
          h_assochi2_prob[www]->Fill(TMath::Prob((assocChi2)*5,5),weighting);
        }
        else if (associators[ww]=="trackAssociatorByHits"){
          double fraction = tp.begin()->second;
          h_assocFraction[www]->Fill(fraction,weighting);
          h_assocSharedHit[www]->Fill(fraction*track->numberOfValidHits(),weighting);
        }
        
        //nchi2 and hits global distributions
        h_nchi2[w]->Fill(track->normalizedChi2(),weighting);
        h_nchi2_prob[w]->Fill(TMath::Prob(track->chi2(),(int)track->ndof()),weighting);
        h_hits[w]->Fill(track->numberOfValidHits(),weighting);
        h_losthits[w]->Fill(track->numberOfLostHits(),weighting);
        chi2_vs_nhits[w]->Fill(track->numberOfValidHits(),track->normalizedChi2(),weighting);
        h_charge[w]->Fill( track->charge() ,weighting);
        
        double ptSim = sqrt(momentumTP.perp2());
        double qoverpSim = tpr->charge()/sqrt(momentumTP.x()*momentumTP.x()+momentumTP.y()*momentumTP.y()+momentumTP.z()*momentumTP.z());
        double thetaSim = momentumTP.theta();
        double lambdaSim = M_PI/2-momentumTP.theta();
        double phiSim = momentumTP.phi();
        double dxySim = (-vertexTP.x()*sin(momentumTP.phi())+vertexTP.y()*cos(momentumTP.phi()));
        double dzSim = vertexTP.z() - (vertexTP.x()*momentumTP.x()+vertexTP.y()*momentumTP.y()) /
                                       sqrt(momentumTP.perp2()) * momentumTP.z()/sqrt(momentumTP.perp2());
        
        // removed unused variable, left this in case it has side effects
        track->parameters();
        
        double qoverpRec(0);
        double qoverpErrorRec(0);
        double ptRec(0);
        double ptErrorRec(0);
        double lambdaRec(0);
        double lambdaErrorRec(0);
        double phiRec(0);
        double phiErrorRec(0);
        
        //loop to decide whether to take gsfTrack (utilisation of mode-function) or common track
        const GsfTrack* gsfTrack(0);
        if(useGsf){
          gsfTrack = dynamic_cast<const GsfTrack*>(&(*track));
          if (gsfTrack==0) edm::LogInfo("MuonTrackValidator") << "Trying to access mode for a non-GsfTrack";
        }
        
        if (gsfTrack) {
          // get values from mode
          getRecoMomentum(*gsfTrack, ptRec, ptErrorRec, qoverpRec, qoverpErrorRec,
                          lambdaRec,lambdaErrorRec, phiRec, phiErrorRec);
        }
        
        else {
          // get values from track (without mode)
          getRecoMomentum(*track, ptRec, ptErrorRec, qoverpRec, qoverpErrorRec,
                          lambdaRec,lambdaErrorRec, phiRec, phiErrorRec);
        }
        
        double thetaRec = track->theta();
        double ptError = ptErrorRec;
        double ptres = ptRec - ptSim;
        double etares = track->eta()-momentumTP.Eta();
        double dxyRec = track->dxy(bs.position());
        double dzRec = track->dz(bs.position());
        // eta residue; pt, k, theta, phi, dxy, dz pulls
        double qoverpPull=(qoverpRec-qoverpSim)/qoverpErrorRec;
        double thetaPull=(lambdaRec-lambdaSim)/lambdaErrorRec;
        double phiDiff = phiRec - phiSim;
        if (abs(phiDiff) > M_PI) {
          if (phiDiff >0.) phiDiff = phiDiff - 2.*M_PI;
          else phiDiff = phiDiff + 2.*M_PI;
        }
        double phiPull=phiDiff/phiErrorRec;
        double dxyPull=(dxyRec-dxySim)/track->dxyError();
        double dzPull=(dzRec-dzSim)/track->dzError();
        
        double contrib_Qoverp = ((qoverpRec-qoverpSim)/qoverpErrorRec)*((qoverpRec-qoverpSim)/qoverpErrorRec)/5;
        double contrib_dxy = ((dxyRec-dxySim)/track->dxyError())*((dxyRec-dxySim)/track->dxyError())/5;
        double contrib_dz = ((dzRec-dzSim)/track->dzError())*((dzRec-dzSim)/track->dzError())/5;
        double contrib_theta = ((lambdaRec-lambdaSim)/lambdaErrorRec)*((lambdaRec-lambdaSim)/lambdaErrorRec)/5;
        double contrib_phi = (phiDiff/phiErrorRec)*(phiDiff/phiErrorRec)/5;
        
        edm::LogVerbatim("MuonTrackValidator") << "assocChi2=" << tp.begin()->second << "\n"
                                               << "" << "\n"
                                               << "ptREC=" << ptRec << "\n"
                                               << "etaREC=" << track->eta() << "\n"
                                               << "qoverpREC=" << qoverpRec << "\n"
                                               << "dxyREC=" << dxyRec << "\n"
                                               << "dzREC=" << dzRec << "\n"
                                               << "thetaREC=" << track->theta() << "\n"
                                               << "phiREC=" << phiRec << "\n"
                                               << "" << "\n"
                                               << "qoverpError()=" << qoverpErrorRec << "\n"
                                               << "dxyError()=" << track->dxyError() << "\n"
                                               << "dzError()=" << track->dzError() << "\n"
                                               << "thetaError()=" << lambdaErrorRec << "\n"
                                               << "phiError()=" << phiErrorRec << "\n"
                                               << "" << "\n"
                                               << "ptSIM=" << ptSim << "\n"
                                               << "etaSIM=" << momentumTP.Eta() << "\n"
                                               << "qoverpSIM=" << qoverpSim << "\n"
                                               << "dxySIM=" << dxySim << "\n"
                                               << "dzSIM=" << dzSim << "\n"
                                               << "thetaSIM=" << M_PI/2-lambdaSim << "\n"
                                               << "phiSIM=" << phiSim << "\n"
                                               << "" << "\n"
                                               << "contrib_Qoverp=" << contrib_Qoverp << "\n"
                                               << "contrib_dxy=" << contrib_dxy << "\n"
                                               << "contrib_dz=" << contrib_dz << "\n"
                                               << "contrib_theta=" << contrib_theta << "\n"
                                               << "contrib_phi=" << contrib_phi << "\n"
                                               << "" << "\n"
                                               <<"chi2PULL="<<contrib_Qoverp+contrib_dxy+contrib_dz+contrib_theta+contrib_phi<<"\n";
        
        h_pullQoverp[w]->Fill(qoverpPull,weighting);
        h_pullTheta[w]->Fill(thetaPull,weighting);
        h_pullPhi[w]->Fill(phiPull,weighting);
        h_pullDxy[w]->Fill(dxyPull,weighting);
        h_pullDz[w]->Fill(dzPull,weighting);
        
        h_pt[w]->Fill(ptres/ptError,weighting);
        h_eta[w]->Fill(etares,weighting);
        etares_vs_eta[w]->Fill(getEta(track->eta()),etares,weighting);
        
        //chi2 and #hit vs eta: fill 2D histos
        chi2_vs_eta[w]->Fill(getEta(track->eta()),track->normalizedChi2(),weighting);
        nhits_vs_eta[w]->Fill(getEta(track->eta()),track->numberOfValidHits(),weighting);
        nDThits_vs_eta[w]->Fill(getEta(track->eta()),track->hitPattern().numberOfValidMuonDTHits(),weighting);
        nCSChits_vs_eta[w]->Fill(getEta(track->eta()),track->hitPattern().numberOfValidMuonCSCHits(),weighting);
        nRPChits_vs_eta[w]->Fill(getEta(track->eta()),track->hitPattern().numberOfValidMuonRPCHits(),weighting);
        
        nlosthits_vs_eta[w]->Fill(getEta(track->eta()),track->numberOfLostHits(),weighting);
        
        //resolution of track params: fill 2D histos
        dxyres_vs_eta[w]->Fill(getEta(track->eta()),dxyRec-dxySim,weighting);
        ptres_vs_eta[w]->Fill(getEta(track->eta()),(ptRec-ptSim)/ptRec,weighting);
        dzres_vs_eta[w]->Fill(getEta(track->eta()),dzRec-dzSim,weighting);
        phires_vs_eta[w]->Fill(getEta(track->eta()),phiDiff,weighting);
        cotThetares_vs_eta[w]->Fill(getEta(track->eta()), cos(thetaRec)/sin(thetaRec) - cos(thetaSim)/sin(thetaSim),weighting);
        
        //same as before but vs pT
        dxyres_vs_pt[w]->Fill(getPt(ptRec),dxyRec-dxySim,weighting);
        ptres_vs_pt[w]->Fill(getPt(ptRec),(ptRec-ptSim)/ptRec,weighting);
        dzres_vs_pt[w]->Fill(getPt(ptRec),dzRec-dzSim,weighting);
        phires_vs_pt[w]->Fill(getPt(ptRec),phiDiff,weighting);
        cotThetares_vs_pt[w]->Fill(getPt(ptRec), cos(thetaRec)/sin(thetaRec) - cos(thetaSim)/sin(thetaSim),weighting);
        
        //pulls of track params vs eta: fill 2D histos
        dxypull_vs_eta[w]->Fill(getEta(track->eta()),dxyPull,weighting);
        ptpull_vs_eta[w]->Fill(getEta(track->eta()),ptres/ptError,weighting);
        dzpull_vs_eta[w]->Fill(getEta(track->eta()),dzPull,weighting);
        phipull_vs_eta[w]->Fill(getEta(track->eta()),phiPull,weighting);
        thetapull_vs_eta[w]->Fill(getEta(track->eta()),thetaPull,weighting);
        
        //plots vs phi
        nhits_vs_phi[w]->Fill(phiRec,track->numberOfValidHits(),weighting);
        chi2_vs_phi[w]->Fill(phiRec,track->normalizedChi2(),weighting);
        ptmean_vs_eta_phi[w]->Fill(phiRec,getEta(track->eta()),ptRec,weighting);
        phimean_vs_eta_phi[w]->Fill(phiRec,getEta(track->eta()),phiRec,weighting);
        ptres_vs_phi[w]->Fill(phiRec,(ptRec-ptSim)/ptRec,weighting);
        phires_vs_phi[w]->Fill(phiRec,phiDiff,weighting);
        ptpull_vs_phi[w]->Fill(phiRec,ptres/ptError,weighting);
        phipull_vs_phi[w]->Fill(phiRec,phiPull,weighting);
        thetapull_vs_phi[w]->Fill(phiRec,thetaPull,weighting);
        
        int nSimHits = 0;
        if (usetracker && usemuon) {
          nSimHits= tpr.get()->numberOfHits();
        }
        else if (!usetracker && usemuon) {
          nSimHits= tpr.get()->numberOfHits() - tpr.get()->numberOfTrackerHits();
        }
        else if (usetracker && !usemuon) {
          nSimHits=tpr.get()->numberOfTrackerHits();
        }
        
        nrecHit_vs_nsimHit_rec2sim[w]->Fill(track->numberOfValidHits(), nSimHits,weighting);
        
      } // End of for(View<Track>::size_type i=0; i<trackCollectionSize; ++i)
      
      if (at!=0) h_tracks[w]->Fill(at,weighting);
      h_fakes[w]->Fill(rT-at,weighting);
      edm::LogVerbatim("MuonTrackValidator") << "Total Simulated: " << st << "\n"
                                             << "Total Associated (simToReco): " << ats << "\n"
                                             << "Total Reconstructed: " << rT << "\n"
                                             << "Total Associated (recoToSim): " << at << "\n"
                                             << "Total Fakes: " << rT-at << "\n";
      nrec_vs_nsim[w]->Fill(rT,st,weighting);
      w++;
      std::cout << "w: " << w << std::endl;
    } // End of for (unsigned int www=0;www<label.size();www++)
  } //END of for (unsigned int ww=0;ww<associators.size();ww++)
  std::cout << "analyze::::: LAST w: " << w << std::endl;
}

void MuonTrackValidator::endRun(Run const&, EventSetup const&) {
  
  int w=0;
  for (unsigned int ww=0;ww<associators.size();ww++){
    for (unsigned int www=0;www<label.size();www++){
      
      //chi2 and #hit vs eta: get mean from 2D histos
      std::cout << "endRun::::: Before doProfileX(): " << chi2_vs_eta[w] << "\t" << h_chi2meanh[w] << "\n";
      doProfileX(chi2_vs_eta[w],h_chi2meanh[w]);
      doProfileX(nhits_vs_eta[w],h_hits_eta[w]);
      doProfileX(nDThits_vs_eta[w],h_DThits_eta[w]);
      doProfileX(nCSChits_vs_eta[w],h_CSChits_eta[w]);
      doProfileX(nRPChits_vs_eta[w],h_RPChits_eta[w]);
      
      doProfileX(nlosthits_vs_eta[w],h_losthits_eta[w]);
      //vs phi
      doProfileX(chi2_vs_nhits[w],h_chi2meanhitsh[w]);
      doProfileX(chi2_vs_phi[w],h_chi2mean_vs_phi[w]);
      doProfileX(nhits_vs_phi[w],h_hits_phi[w]);

      std::cout << "fillPlotFromVector(0): " << h_recoeta[w] << " " << totRECeta[w][0] << "\n";
      fillPlotFromVector(h_recoeta[w],totRECeta[w]);
      std::cout << "fillPlotFromVector(1): " << h_simuleta[w] << " " << totSIMeta[w][0] << "\n";
      fillPlotFromVector(h_simuleta[w],totSIMeta[w]);
      std::cout << "fillPlotFromVector(2): " << h_assoceta[w] << " " << totASSeta[w][0] << "\n";
      fillPlotFromVector(h_assoceta[w],totASSeta[w]);
      std::cout << "fillPlotFromVector(3): " << h_assoc2eta[w] << " " << totASS2eta[w][0] << "\n";
      fillPlotFromVector(h_assoc2eta[w],totASS2eta[w]);

      std::cout << "fillPlotFromVector(4): " << h_recopT[w] << " " << totRECpT[w][0] << "\n";
      fillPlotFromVector(h_recopT[w],totRECpT[w]);
      fillPlotFromVector(h_simulpT[w],totSIMpT[w]);
      fillPlotFromVector(h_assocpT[w],totASSpT[w]);
      fillPlotFromVector(h_assoc2pT[w],totASS2pT[w]);

      std::cout << "fillPlotFromVector(5): " << h_recohit[w] << " " << totREC_hit[w][0] << "\n";
      fillPlotFromVector(h_recohit[w],totREC_hit[w]);
      fillPlotFromVector(h_simulhit[w],totSIM_hit[w]);
      fillPlotFromVector(h_assochit[w],totASS_hit[w]);
      fillPlotFromVector(h_assoc2hit[w],totASS2_hit[w]);

      fillPlotFromVector(h_recophi[w],totREC_phi[w]);
      fillPlotFromVector(h_simulphi[w],totSIM_phi[w]);
      fillPlotFromVector(h_assocphi[w],totASS_phi[w]);
      fillPlotFromVector(h_assoc2phi[w],totASS2_phi[w]);

      fillPlotFromVector(h_recodxy[w],totREC_dxy[w]);
      fillPlotFromVector(h_simuldxy[w],totSIM_dxy[w]);
      fillPlotFromVector(h_assocdxy[w],totASS_dxy[w]);
      fillPlotFromVector(h_assoc2dxy[w],totASS2_dxy[w]);

      fillPlotFromVector(h_recodz[w],totREC_dz[w]);
      fillPlotFromVector(h_simuldz[w],totSIM_dz[w]);
      fillPlotFromVector(h_assocdz[w],totASS_dz[w]);
      fillPlotFromVector(h_assoc2dz[w],totASS2_dz[w]);

      fillPlotFromVector(h_simulvertpos[w],totSIM_vertpos[w]);
      fillPlotFromVector(h_assocvertpos[w],totASS_vertpos[w]);

      fillPlotFromVector(h_simulzpos[w],totSIM_zpos[w]);
      fillPlotFromVector(h_assoczpos[w],totASS_zpos[w]);
      
      if (MABH) {
        fillPlotFromVector(h_assoceta_Quality05[w] ,totASSeta_Quality05[w]);
        fillPlotFromVector(h_assoceta_Quality075[w],totASSeta_Quality075[w]);
        fillPlotFromVector(h_assocpT_Quality05[w] ,totASSpT_Quality05[w]);
        fillPlotFromVector(h_assocpT_Quality075[w],totASSpT_Quality075[w]);
        fillPlotFromVector(h_assocphi_Quality05[w] ,totASS_phi_Quality05[w]);
        fillPlotFromVector(h_assocphi_Quality075[w],totASS_phi_Quality075[w]);
      }
      
      std::cout << "Before w++ " << w << "\n";
      w++;
    }
  }
  
  std::cout << "Will be saved " << out.size() << "\t" << dbe_ << std::endl;
  if ( out.size() != 0 && dbe_ ) {
    dbe_->save(out);
    std::cout << "endRun::::: Saved " << w << std::endl;
  }
}

void MuonTrackValidator::getRecoMomentum (const reco::Track& track,
                                          double& pt, double& ptError, double& qoverp, double& qoverpError,
                                          double& lambda,double& lambdaError, double& phi, double& phiError ) const {
  pt = track.pt();
  ptError = track.ptError();
  qoverp = track.qoverp();
  qoverpError = track.qoverpError();
  lambda = track.lambda();
  lambdaError = track.lambdaError();
  phi = track.phi();
  phiError = track.phiError();
}

void MuonTrackValidator::getRecoMomentum (const reco::GsfTrack& gsfTrack,
                                          double& pt, double& ptError, double& qoverp, double& qoverpError,
                                          double& lambda,double& lambdaError, double& phi, double& phiError ) const {
  pt = gsfTrack.ptMode();
  ptError = gsfTrack.ptModeError();
  qoverp = gsfTrack.qoverpMode();
  qoverpError = gsfTrack.qoverpModeError();
  lambda = gsfTrack.lambdaMode();
  lambdaError = gsfTrack.lambdaModeError();
  phi = gsfTrack.phiMode();
  phiError = gsfTrack.phiModeError(); 
}

void MuonTrackValidator::setUpAsymmVectors () {
  if (!_pTRanges.empty()) {
    std::vector<double> pTintervalsv;
    for (int idx=0; idx<=nintpT; idx++) pTintervalsv.push_back(pTRanges[idx]);

    pTintervals.clear();
    pTintervals.push_back(pTintervalsv);
  }
}
