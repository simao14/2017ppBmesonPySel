#include "TROOT.h"
#include "TH1.h"
#include "TTree.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TSystem.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom.h"
#include <iostream>
#include <fstream>

using namespace std;

using std::cout;
using std::endl;

enum Tracking{
  loose = 0,
  standard = 1,
  tight = 2,
};

std::map<Tracking, double> ptErr{
  {Tracking::standard, 0.1},
  {Tracking::tight, 0.05},
  {Tracking::loose, 0.15}
};

std::map<Tracking, double> chi2Nlayer{
  {Tracking::standard, 0.18},
  {Tracking::tight, 0.15},
  {Tracking::loose, 0.18}
};

void smallfile(int meson_n){

    cout<<"Running"<<endl;
    TString var;
    TString min;
    TString max;

    if (meson_n==0){var="BP";min="3";max="60";}
    else {var="Bs";min="3";max="50";}

    TFile* fileS = new TFile(Form("/data3/smcosta/2017ppBmesonTMVA/%s/sample/%sMC_%s_%s.root",var.Data(),var.Data(),min.Data(),max.Data()));
    TFile* fileB = new TFile(Form("/data3/smcosta/2017ppBmesonTMVA/%s/sample/%sData_%s_%s.root",var.Data(),var.Data(),min.Data(),max.Data()));

    auto outS = Form("/data3/smcosta/2017ppBmesonTMVA/%s/sample/%sMC_%s_%s_small.root",var.Data(),var.Data(),min.Data(),max.Data());
    auto outB = Form("/data3/smcosta/2017ppBmesonTMVA/%s/sample/%sData_%s_%s_small.root",var.Data(),var.Data(),min.Data(),max.Data());
    auto outG = Form("/data3/smcosta/2017ppBmesonTMVA/%s/sample/%sGen.root",var.Data(),var.Data());

    TFile* outfS = new TFile( outS,"RECREATE");
    TFile* outfB = new TFile( outB,"RECREATE");
    TFile* outfG = new TFile( outG,"RECREATE");


    TTree* background;
    TTree* backgroundroot;
    TTree* backgroundhl;
    TTree* backgroundskim;
    TTree* backgroundnew = new TTree("background2", "background full");

    if (meson_n==0){ background = (TTree*) fileB->Get("Bfinder/ntKp");}
    else { background = (TTree*) fileB->Get("Bfinder/ntphi");}
    backgroundroot = (TTree*) fileB->Get("Bfinder/root");
    backgroundhl = (TTree*) fileB->Get("hltanalysis/HltTree");
    backgroundskim = (TTree*) fileB->Get("skimanalysis/HltTree");
    

    TTree* signal;
    TTree* gen;
    TTree* signalroot;
    TTree* signalhl;
    TTree* signalskim;
    TTree* signalhi;
    TTree* signalnew = new TTree("signal2", "signal full");
    if (meson_n==0){ signal = (TTree*) fileS->Get("Bfinder/ntKp");}
    else { signal = (TTree*) fileS->Get("Bfinder/ntphi");}
    signalroot = (TTree*) fileS->Get("Bfinder/root");
    signalhl = (TTree*) fileS->Get("hltanalysis/HltTree");
    signalskim = (TTree*) fileS->Get("skimanalysis/HltTree");
    signalhi = (TTree*) fileS->Get("hiEvtAnalyzer/HiTree");
    gen = (TTree*) fileS->Get("Bfinder/ntGen");


    TTree* gennew = gen->CloneTree(0);

    const int ncand=1000;
    std::array<Float_t,ncand> Btrk1Dz1B;
    std::array<Float_t,ncand> Btrk1DzError1B;
    std::array<Float_t,ncand> Btrk2Dz1B;
    std::array<Float_t,ncand> Btrk2DzError1B;
    std::array<Float_t,ncand> Btrk1Dxy1B;
    std::array<Float_t,ncand> Btrk1DxyError1B;
    std::array<Float_t,ncand> Btrk2Dxy1B;
    std::array<Float_t,ncand> Btrk2DxyError1B;
    std::array<Float_t,ncand> BtktkmassB;
    std::array<Float_t,ncand> BsvpvDistanceB;
    std::array<Float_t,ncand> BsvpvDisErrB;
    std::array<Float_t,ncand> Bd0B;
    std::array<Float_t,ncand> BmassB;
    std::array<Float_t,ncand> BptB;
    std::array<Bool_t,ncand> Bmu1isTriggeredB;
    std::array<Bool_t,ncand> Bmu2isTriggeredB;
    std::array<Float_t,ncand> Btrk1PtB;
    std::array<Float_t,ncand> Btrk1PtErrB;
    std::array<Float_t,ncand> Btrk2PtErrB;
    std::array<Float_t,ncand> Bchi2clB;
    std::array<Float_t,ncand> Btrk1EtaB;
    std::array<Float_t,ncand> ByB;
    std::array<Float_t,ncand> BmumumassB;
    std::array<Float_t,ncand> Bmu1etaB;
    std::array<Float_t,ncand> Bmu1ptB;
    std::array<Float_t,ncand> Bmu2etaB;
    std::array<Float_t,ncand> Bmu2ptB;
    std::array<Bool_t,ncand> Bmu1TMOneStationTightB;
    std::array<Bool_t,ncand> Bmu2TMOneStationTightB;
    std::array<Int_t,ncand> Bmu1InPixelLayerB;
    std::array<Int_t,ncand> Bmu1InStripLayerB;
    std::array<Int_t,ncand> Bmu2InPixelLayerB;
    std::array<Int_t,ncand> Bmu2InStripLayerB;
    std::array<Float_t,ncand> Bmu1dxyPVB;
    std::array<Float_t,ncand> Bmu2dxyPVB;
    std::array<Float_t,ncand> Bmu1dzPVB;
    std::array<Float_t,ncand> Bmu2dzPVB;
    std::array<Bool_t,ncand> Bmu1isTrackerMuonB;
    std::array<Bool_t,ncand> Bmu2isTrackerMuonB;
    std::array<Bool_t,ncand> Bmu1isGlobalMuonB;
    std::array<Bool_t,ncand> Bmu2isGlobalMuonB;
    std::array<Bool_t,ncand> Btrk1highPurityB;
    std::array<Bool_t,ncand> Btrk2highPurityB;
    std::array<Float_t,ncand> Btrk1PixelHitB;
    std::array<Float_t,ncand> Btrk2PixelHitB;
    std::array<Float_t,ncand> Btrk1StripHitB;
    std::array<Float_t,ncand> Btrk2StripHitB;
    std::array<Float_t,ncand> Btrk1Chi2ndfB;
    std::array<Float_t,ncand> Btrk2Chi2ndfB;
    std::array<Float_t,ncand> Btrk2PtB;
    std::array<Float_t,ncand> BalphaB;
    std::array<Float_t,ncand> Btrk2EtaB;
    std::array<Float_t,ncand> Btrk1nStripLayerB;
    std::array<Float_t,ncand> Btrk2nStripLayerB;
    std::array<Float_t,ncand> Btrk1nPixelLayerB;
    std::array<Float_t,ncand> Btrk2nPixelLayerB;

   Float_t Btrk1Dz1flatB;
   Float_t Btrk1DzError1flatB;
   Float_t Btrk2Dz1flatB;
   Float_t Btrk2DzError1flatB;
   Float_t Btrk1Dxy1flatB;
   Float_t Btrk1DxyError1flatB;
   Float_t Btrk2Dxy1flatB;
   Float_t Btrk2DxyError1flatB;
   Float_t BtktkmassflatB;
   Float_t BsvpvDistanceflatB;
   Float_t BsvpvDisErrflatB;
   Float_t Bd0flatB;
   Float_t BmassflatB;
   Float_t BptflatB;
   Bool_t Bmu1isTriggeredflatB;
   Bool_t Bmu2isTriggeredflatB;
   Float_t Btrk1PtflatB;
   Float_t Btrk1PtErrflatB;
   Float_t Btrk2PtErrflatB;
   Float_t Bchi2clflatB;
   Float_t Btrk1EtaflatB;
   Float_t ByflatB;
   Float_t BmumumassflatB;
   Float_t Bmu1etaflatB;
   Float_t Bmu1ptflatB;
   Float_t Bmu2etaflatB;
   Float_t Bmu2ptflatB;
   Bool_t Bmu1TMOneStationTightflatB;
   Bool_t Bmu2TMOneStationTightflatB;
   Int_t Bmu1InPixelLayerflatB;
   Int_t Bmu1InStripLayerflatB;
   Int_t Bmu2InPixelLayerflatB;
   Int_t Bmu2InStripLayerflatB;
   Float_t Bmu1dxyPVflatB;
   Float_t Bmu2dxyPVflatB;
   Float_t Bmu1dzPVflatB;
   Float_t Bmu2dzPVflatB;
   Bool_t Bmu1isTrackerMuonflatB;
   Bool_t Bmu2isTrackerMuonflatB;
   Bool_t Bmu1isGlobalMuonflatB;
   Bool_t Bmu2isGlobalMuonflatB;
   Bool_t Btrk1highPurityflatB;
   Bool_t Btrk2highPurityflatB;
   Float_t Btrk1PixelHitflatB;
   Float_t Btrk2PixelHitflatB;
   Float_t Btrk1StripHitflatB;
   Float_t Btrk2StripHitflatB;
   Float_t Btrk1Chi2ndfflatB;
   Float_t Btrk2Chi2ndfflatB;
   Float_t Btrk2PtflatB;
   Float_t BalphaflatB;
   Float_t Btrk2EtaflatB;
   Float_t Btrk1nStripLayerflatB;
   Float_t Btrk2nStripLayerflatB;
   Float_t Btrk1nPixelLayerflatB;
   Float_t Btrk2nPixelLayerflatB;

    Float_t PVzB;
    Int_t BsizeB;
    Int_t nMultB;
    Int_t trackB;
    

    Int_t pPAprimaryVertexFilterB;
    Int_t pBeamScrapingFilterB;
    Int_t HLT_HIL1DoubleMu0_v1B;

    std::array<Float_t,ncand> Btrk1Dz1S;
    std::array<Float_t,ncand> Btrk1DzError1S;
    std::array<Float_t,ncand> Btrk2Dz1S;
    std::array<Float_t,ncand> Btrk2DzError1S;
    std::array<Float_t,ncand> Btrk1Dxy1S;
    std::array<Float_t,ncand> Btrk1DxyError1S;
    std::array<Float_t,ncand> Btrk2Dxy1S;
    std::array<Float_t,ncand> Btrk2DxyError1S;
    std::array<Float_t,ncand> BtktkmassS;
    std::array<Float_t,ncand> BsvpvDistanceS;
    std::array<Float_t,ncand> BsvpvDisErrS;
    std::array<Float_t,ncand> Bd0S;
    std::array<Float_t,ncand> BmassS;
    std::array<Float_t,ncand> BptS;
    std::array<Bool_t,ncand> Bmu1isTriggeredS;
    std::array<Bool_t,ncand> Bmu2isTriggeredS;
    std::array<Float_t,ncand> Btrk1PtS;
    std::array<Float_t,ncand> Btrk1PtErrS;
    std::array<Float_t,ncand> Btrk2PtErrS;
    std::array<Float_t,ncand> Bchi2clS;
    std::array<Float_t,ncand> Btrk1EtaS;
    std::array<Float_t,ncand> ByS;
    std::array<Float_t,ncand> BmumumassS;
    std::array<Float_t,ncand> Bmu1etaS;
    std::array<Float_t,ncand> Bmu1ptS;
    std::array<Float_t,ncand> Bmu2etaS;
    std::array<Float_t,ncand> Bmu2ptS;
    std::array<Bool_t,ncand> Bmu1TMOneStationTightS;
    std::array<Bool_t,ncand> Bmu2TMOneStationTightS;
    std::array<Int_t,ncand> Bmu1InPixelLayerS;
    std::array<Int_t,ncand> Bmu1InStripLayerS;
    std::array<Int_t,ncand> Bmu2InPixelLayerS;
    std::array<Int_t,ncand> Bmu2InStripLayerS;
    std::array<Float_t,ncand> Bmu1dxyPVS;
    std::array<Float_t,ncand> Bmu2dxyPVS;
    std::array<Float_t,ncand> Bmu1dzPVS;
    std::array<Float_t,ncand> Bmu2dzPVS;
    std::array<Bool_t,ncand> Bmu1isTrackerMuonS;
    std::array<Bool_t,ncand> Bmu2isTrackerMuonS;
    std::array<Bool_t,ncand> Bmu1isGlobalMuonS;
    std::array<Bool_t,ncand> Bmu2isGlobalMuonS;
    std::array<Bool_t,ncand> Btrk1highPurityS;
    std::array<Bool_t,ncand> Btrk2highPurityS;
    std::array<Float_t,ncand> Btrk1PixelHitS;
    std::array<Float_t,ncand> Btrk2PixelHitS;
    std::array<Float_t,ncand> Btrk1StripHitS;
    std::array<Float_t,ncand> Btrk2StripHitS;
    std::array<Float_t,ncand> Btrk1Chi2ndfS;
    std::array<Float_t,ncand> Btrk2Chi2ndfS;
    std::array<Float_t,ncand> Btrk2PtS;
    std::array<Float_t,ncand> BalphaS;
    std::array<Float_t,ncand> Btrk2EtaS;
    std::array<Float_t,ncand> BgenS;
    std::array<Float_t,ncand> BgenptS;
    std::array<Float_t,ncand> BgenyS;
    std::array<Float_t,ncand> Btrk1nStripLayerS;
    std::array<Float_t,ncand> Btrk2nStripLayerS;
    std::array<Float_t,ncand> Btrk1nPixelLayerS;
    std::array<Float_t,ncand> Btrk2nPixelLayerS;

    Float_t Btrk1Dz1flatS;
    Float_t Btrk1DzError1flatS;
    Float_t Btrk2Dz1flatS;
    Float_t Btrk2DzError1flatS;
    Float_t Btrk1Dxy1flatS;
    Float_t Btrk1DxyError1flatS;
    Float_t Btrk2Dxy1flatS;
    Float_t Btrk2DxyError1flatS;
    Float_t BtktkmassflatS;
    Float_t BsvpvDistanceflatS;
    Float_t BsvpvDisErrflatS;
    Float_t Bd0flatS;
    Float_t BmassflatS;
    Float_t BptflatS;
    Bool_t Bmu1isTriggeredflatS;
    Bool_t Bmu2isTriggeredflatS;
    Float_t Btrk1PtflatS;
    Float_t Btrk1PtErrflatS;
    Float_t Btrk2PtErrflatS;
    Float_t Bchi2clflatS;
    Float_t Btrk1EtaflatS;
    Float_t ByflatS;
    Float_t BmumumassflatS;
    Float_t Bmu1etaflatS;
    Float_t Bmu1ptflatS;
    Float_t Bmu2etaflatS;
    Float_t Bmu2ptflatS;
    Bool_t Bmu1TMOneStationTightflatS;
    Bool_t Bmu2TMOneStationTightflatS;
    Int_t Bmu1InPixelLayerflatS;
    Int_t Bmu1InStripLayerflatS;
    Int_t Bmu2InPixelLayerflatS;
    Int_t Bmu2InStripLayerflatS;
    Float_t Bmu1dxyPVflatS;
    Float_t Bmu2dxyPVflatS;
    Float_t Bmu1dzPVflatS;
    Float_t Bmu2dzPVflatS;
    Bool_t Bmu1isTrackerMuonflatS;
    Bool_t Bmu2isTrackerMuonflatS;
    Bool_t Bmu1isGlobalMuonflatS;
    Bool_t Bmu2isGlobalMuonflatS;
    Bool_t Btrk1highPurityflatS;
    Bool_t Btrk2highPurityflatS;
    Float_t Btrk1PixelHitflatS;
    Float_t Btrk2PixelHitflatS;
    Float_t Btrk1StripHitflatS;
    Float_t Btrk2StripHitflatS;
    Float_t Btrk1Chi2ndfflatS;
    Float_t Btrk2Chi2ndfflatS;
    Float_t Btrk2PtflatS;
    Float_t BalphaflatS;
    Float_t Btrk2EtaflatS;
    Float_t Btrk1nStripLayerflatS;
    Float_t Btrk2nStripLayerflatS;
    Float_t Btrk1nPixelLayerflatS;
    Float_t Btrk2nPixelLayerflatS;
    Float_t BgenflatS;
    Float_t BgenptflatS;
    Float_t BgenyflatS;

    Float_t PVzS;
    Int_t BsizeS;
    Int_t nMultS;
    Int_t trackS;

    Int_t pPAprimaryVertexFilterS;
    Int_t pBeamScrapingFilterS;
    Int_t HLT_HIL1DoubleMu0_v1S;

    Float_t weightS;
   
    signalhi->SetBranchAddress("weight",&weightS);
    signalnew->Branch("weight",&weightS);
    gennew->Branch("weight",&weightS);

    backgroundnew->Branch("track",&trackB,"track/I");
    signalnew->Branch("track",&trackS,"track/I");

    background->SetBranchAddress("PVz",&PVzB);
    backgroundnew->Branch("PVz",&PVzB);

    signal->SetBranchAddress("PVz",&PVzS);
    signalnew->Branch("PVz",&PVzS);
    gennew->Branch("PVz",&PVzS);

    background->SetBranchAddress("Bsize",&BsizeB);
    backgroundnew->Branch("Bsize",&BsizeB,"Bsize/I");

    signal->SetBranchAddress("Bsize",&BsizeS);
    signalnew->Branch("Bsize",&BsizeS,"Bsize/I");

    backgroundroot->SetBranchAddress("EvtInfo.nMult",&nMultB);
    backgroundnew->Branch("nMult",&nMultB,"nMult/I");

    signalroot->SetBranchAddress("EvtInfo.nMult",&nMultS);
    signalnew->Branch("nMult",&nMultS,"nMult/I");

    signal->SetBranchAddress("Bgen",&BgenS);
    signalnew->Branch("Bgen",&BgenflatS,"Bgen/F");

    signal->SetBranchAddress("Bgenpt",&BgenptS);
    signalnew->Branch("Bgenpt",&BgenptflatS,"Bgenpt/F");

    signal->SetBranchAddress("Bgeny",&BgenyS);
    signalnew->Branch("Bgeny",&BgenyflatS,"Bgeny/F");

    background->SetBranchAddress("Btrk1Dz1",&Btrk1Dz1B);
    backgroundnew->Branch("Btrk1Dz1",&Btrk1Dz1flatB,"Btrk1Dz1/F");

    signal->SetBranchAddress("Btrk1Dz1",&Btrk1Dz1S);
    signalnew->Branch("Btrk1Dz1",&Btrk1Dz1flatS,"Btrk1Dz1/F");

    background->SetBranchAddress("Btrk1DzError1",&Btrk1DzError1B);
    backgroundnew->Branch("Btrk1DzError1",&Btrk1DzError1flatB,"Btrk1DzError1/F");

    signal->SetBranchAddress("Btrk1DzError1",&Btrk1DzError1S);
    signalnew->Branch("Btrk1DzError1",&Btrk1DzError1flatS,"Btrk1DzError1/F");

    background->SetBranchAddress("Btrk2Dz1",&Btrk2Dz1B);
    backgroundnew->Branch("Btrk2Dz1",&Btrk2Dz1flatB,"Btrk2Dz1/F");

    signal->SetBranchAddress("Btrk2Dz1",&Btrk2Dz1S);
    signalnew->Branch("Btrk2Dz1",&Btrk2Dz1flatS,"Btrk2Dz1/F");

    background->SetBranchAddress("Btrk2DzError1",&Btrk2DzError1B);
    backgroundnew->Branch("Btrk2DzError1",&Btrk2DzError1flatB,"Btrk2DzError1/F");

    signal->SetBranchAddress("Btrk2DzError1",&Btrk2DzError1S);
    signalnew->Branch("Btrk2DzError1",&Btrk2DzError1flatS,"Btrk2DzError1/F");

    background->SetBranchAddress("Btrk1Dxy1",&Btrk1Dxy1B);
    backgroundnew->Branch("Btrk1Dxy1",&Btrk1Dxy1flatB,"Btrk1Dxy1/F");

    signal->SetBranchAddress("Btrk1Dxy1",&Btrk1Dxy1S);
    signalnew->Branch("Btrk1Dxy1",&Btrk1Dxy1flatS,"Btrk1Dxy1/F");

    background->SetBranchAddress("Btrk1DxyError1",&Btrk1DxyError1B);
    backgroundnew->Branch("Btrk1DxyError1",&Btrk1DxyError1flatB,"Btrk1DxyError1/F");

    signal->SetBranchAddress("Btrk1DxyError1",&Btrk1DxyError1S);
    signalnew->Branch("Btrk1DxyError1",&Btrk1DxyError1flatS,"Btrk1DxyError1/F");

    background->SetBranchAddress("Btrk2Dxy1",&Btrk2Dxy1B);
    backgroundnew->Branch("Btrk2Dxy1",&Btrk2Dxy1flatB,"Btrk2Dxy1/F");

    signal->SetBranchAddress("Btrk2Dxy1",&Btrk2Dxy1S);
    signalnew->Branch("Btrk2Dxy1",&Btrk2Dxy1flatS,"Btrk2Dxy1/F");

    background->SetBranchAddress("Btrk2DxyError1",&Btrk2DxyError1B);
    backgroundnew->Branch("Btrk2DxyError1",&Btrk2DxyError1flatB,"Btrk2DxyError1/F");

    signal->SetBranchAddress("Btrk2DxyError1",&Btrk2DxyError1S);
    signalnew->Branch("Btrk2DxyError1",&Btrk2DxyError1flatS,"Btrk2DxyError1/F");

    background->SetBranchAddress("Btktkmass",&BtktkmassB);
    backgroundnew->Branch("Btktkmass",&BtktkmassflatB,"Btktkmass/F");

    signal->SetBranchAddress("Btktkmass",&BtktkmassS);
    signalnew->Branch("Btktkmass",&BtktkmassflatS,"Btktkmass/F");

    background->SetBranchAddress("BsvpvDistance",&BsvpvDistanceB);
    backgroundnew->Branch("BsvpvDistance",&BsvpvDistanceflatB,"BsvpvDistance/F");

    signal->SetBranchAddress("BsvpvDistance",&BsvpvDistanceS);
    signalnew->Branch("BsvpvDistance",&BsvpvDistanceflatS,"BsvpvDistance/F");

    background->SetBranchAddress("BsvpvDisErr",&BsvpvDisErrB);
    backgroundnew->Branch("BsvpvDisErr",&BsvpvDisErrflatB,"BsvpvDisErr/F");

    signal->SetBranchAddress("BsvpvDisErr",&BsvpvDisErrS);
    signalnew->Branch("BsvpvDisErr",&BsvpvDisErrflatS,"BsvpvDisErr/F");

    background->SetBranchAddress("Bd0",&Bd0B);
    backgroundnew->Branch("Bd0",&Bd0flatB,"Bd0/F");

    signal->SetBranchAddress("Bd0",&Bd0S);
    signalnew->Branch("Bd0",&Bd0flatS,"Bd0/F");

    background->SetBranchAddress("Bmass",&BmassB);
    backgroundnew->Branch("Bmass",&BmassflatB,"Bmass/F");

    signal->SetBranchAddress("Bmass",&BmassS);
    signalnew->Branch("Bmass",&BmassflatS,"Bmass/F");

    background->SetBranchAddress("Bpt",&BptB);
    backgroundnew->Branch("Bpt",&BptflatB,"Bpt/F");

    signal->SetBranchAddress("Bpt",&BptS);
    signalnew->Branch("Bpt",&BptflatS,"Bpt/F");

    background->SetBranchAddress("Btrk1Pt",&Btrk1PtB);
    backgroundnew->Branch("Btrk1Pt",&Btrk1PtflatB,"Btrk1Pt/F");

    signal->SetBranchAddress("Btrk1Pt",&Btrk1PtS);
    signalnew->Branch("Btrk1Pt",&Btrk1PtflatS,"Btrk1Pt/F");

    background->SetBranchAddress("Btrk1PtErr",&Btrk1PtErrB);
    backgroundnew->Branch("Btrk1PtErr",&Btrk1PtErrflatB,"Btrk1PtErr/F");

    signal->SetBranchAddress("Btrk1PtErr",&Btrk1PtErrS);
    signalnew->Branch("Btrk1PtErr",&Btrk1PtErrflatS,"Btrk1PtErr/F");

    background->SetBranchAddress("Btrk2PtErr",&Btrk2PtErrB);
    backgroundnew->Branch("Btrk2PtErr",&Btrk2PtErrflatB,"Btrk2PtErr/F");

    signal->SetBranchAddress("Btrk2PtErr",&Btrk2PtErrS);
    signalnew->Branch("Btrk2PtErr",&Btrk2PtErrflatS,"Btrk2PtErr/F");

    background->SetBranchAddress("Bchi2cl",&Bchi2clB);
    backgroundnew->Branch("Bchi2cl",&Bchi2clflatB,"Bchi2cl/F");

    signal->SetBranchAddress("Bchi2cl",&Bchi2clS);
    signalnew->Branch("Bchi2cl",&Bchi2clflatS,"Bchi2cl/F");

    background->SetBranchAddress("Btrk1Eta",&Btrk1EtaB);
    backgroundnew->Branch("Btrk1Eta",&Btrk1EtaflatB,"Btrk1Eta/F");

    signal->SetBranchAddress("Btrk1Eta",&Btrk1EtaS);
    signalnew->Branch("Btrk1Eta",&Btrk1EtaflatS,"Btrk1Eta/F");

    background->SetBranchAddress("By",&ByB);
    backgroundnew->Branch("By",&ByflatB,"By/F");

    signal->SetBranchAddress("By",&ByS);
    signalnew->Branch("By",&ByflatS,"By/F");

    background->SetBranchAddress("Bmumumass",&BmumumassB);
    backgroundnew->Branch("Bmumumass",&BmumumassflatB,"Bmumumass/F");

    signal->SetBranchAddress("Bmumumass",&BmumumassS);
    signalnew->Branch("Bmumumass",&BmumumassflatS,"Bmumumass/F");

    background->SetBranchAddress("Bmu1eta",&Bmu1etaB);
    backgroundnew->Branch("Bmu1eta",&Bmu1etaflatB,"Bmu1eta/F");

    signal->SetBranchAddress("Bmu1eta",&Bmu1etaS);
    signalnew->Branch("Bmu1eta",&Bmu1etaflatS,"Bmu1eta/F");

    background->SetBranchAddress("Bmu1pt",&Bmu1ptB);
    backgroundnew->Branch("Bmu1pt",&Bmu1ptflatB,"Bmu1pt/F");

    signal->SetBranchAddress("Bmu1pt",&Bmu1ptS);
    signalnew->Branch("Bmu1pt",&Bmu1ptflatS,"Bmu1pt/F");

    background->SetBranchAddress("Bmu2eta",&Bmu2etaB);
    backgroundnew->Branch("Bmu2eta",&Bmu2etaflatB,"Bmu2eta/F");

    signal->SetBranchAddress("Bmu2eta",&Bmu2etaS);
    signalnew->Branch("Bmu2eta",&Bmu2etaflatS,"Bmu2eta/F");

    background->SetBranchAddress("Bmu2pt",&Bmu2ptB);
    backgroundnew->Branch("Bmu2pt",&Bmu2ptflatB,"Bmu2pt/F");

    signal->SetBranchAddress("Bmu2pt",&Bmu2ptS);
    signalnew->Branch("Bmu2pt",&Bmu2ptflatS,"Bmu2pt/F");

    background->SetBranchAddress("Bmu1dxyPV",&Bmu1dxyPVB);
    backgroundnew->Branch("Bmu1dxyPV",&Bmu1dxyPVflatB,"Bmu1dxyPV/F");

    signal->SetBranchAddress("Bmu1dxyPV",&Bmu1dxyPVS);
    signalnew->Branch("Bmu1dxyPV",&Bmu1dxyPVflatS,"Bmu1dxyPV/F");

    background->SetBranchAddress("Bmu2dxyPV",&Bmu2dxyPVB);
    backgroundnew->Branch("Bmu2dxyPV",&Bmu2dxyPVflatB,"Bmu2dxyPV/F");

    signal->SetBranchAddress("Bmu2dxyPV",&Bmu2dxyPVS);
    signalnew->Branch("Bmu2dxyPV",&Bmu2dxyPVflatS,"Bmu2dxyPV/F");

    background->SetBranchAddress("Bmu1dzPV",&Bmu1dzPVB);
    backgroundnew->Branch("Bmu1dzPV",&Bmu1dzPVflatB,"Bmu1dzPV/F");

    signal->SetBranchAddress("Bmu1dzPV",&Bmu1dzPVS);
    signalnew->Branch("Bmu1dzPV",&Bmu1dzPVflatS,"Bmu1dzPV/F");

    background->SetBranchAddress("Bmu2dzPV",&Bmu2dzPVB);
    backgroundnew->Branch("Bmu2dzPV",&Bmu2dzPVflatB,"Bmu2dzPV/F");

    signal->SetBranchAddress("Bmu2dzPV",&Bmu2dzPVS);
    signalnew->Branch("Bmu2dzPV",&Bmu2dzPVflatS,"Bmu2dzPV/F");

    background->SetBranchAddress("Btrk1PixelHit",&Btrk1PixelHitB);
    backgroundnew->Branch("Btrk1PixelHit",&Btrk1PixelHitflatB,"Btrk1PixelHit/F");

    signal->SetBranchAddress("Btrk1PixelHit",&Btrk1PixelHitS);
    signalnew->Branch("Btrk1PixelHit",&Btrk1PixelHitflatS,"Btrk1PixelHit/F");

    background->SetBranchAddress("Btrk2PixelHit",&Btrk2PixelHitB);
    backgroundnew->Branch("Btrk2PixelHit",&Btrk2PixelHitflatB,"Btrk2PixelHit/F");

    signal->SetBranchAddress("Btrk2PixelHit",&Btrk2PixelHitS);
    signalnew->Branch("Btrk2PixelHit",&Btrk2PixelHitflatS,"Btrk2PixelHit/F");

    background->SetBranchAddress("Btrk1StripHit",&Btrk1StripHitB);
    backgroundnew->Branch("Btrk1StripHit",&Btrk1StripHitflatB,"Btrk1StripHit/F");

    signal->SetBranchAddress("Btrk1StripHit",&Btrk1StripHitS);
    signalnew->Branch("Btrk1StripHit",&Btrk1StripHitflatS,"Btrk1StripHit/F");

    background->SetBranchAddress("Btrk2StripHit",&Btrk2StripHitB);
    backgroundnew->Branch("Btrk2StripHit",&Btrk2StripHitflatB,"Btrk2StripHit/F");

    signal->SetBranchAddress("Btrk2StripHit",&Btrk2StripHitS);
    signalnew->Branch("Btrk2StripHit",&Btrk2StripHitflatS,"Btrk2StripHit/F");

    background->SetBranchAddress("Btrk1Chi2ndf",&Btrk1Chi2ndfB);
    backgroundnew->Branch("Btrk1Chi2ndf",&Btrk1Chi2ndfflatB,"Btrk1Chi2ndf/F");

    signal->SetBranchAddress("Btrk1Chi2ndf",&Btrk1Chi2ndfS);
    signalnew->Branch("Btrk1Chi2ndf",&Btrk1Chi2ndfflatS,"Btrk1Chi2ndf/F");

    background->SetBranchAddress("Btrk2Chi2ndf",&Btrk2Chi2ndfB);
    backgroundnew->Branch("Btrk2Chi2ndf",&Btrk2Chi2ndfflatB,"Btrk2Chi2ndf/F");

    signal->SetBranchAddress("Btrk2Chi2ndf",&Btrk2Chi2ndfS);
    signalnew->Branch("Btrk2Chi2ndf",&Btrk2Chi2ndfflatS,"Btrk2Chi2ndf/F");

    background->SetBranchAddress("Btrk1nStripLayer",&Btrk1nStripLayerB);
    backgroundnew->Branch("Btrk1nStripLayer",&Btrk1nStripLayerflatB,"Btrk1nStripLayer/F");

    signal->SetBranchAddress("Btrk1nStripLayer",&Btrk1nStripLayerS);
    signalnew->Branch("Btrk1nStripLayer",&Btrk1nStripLayerflatS,"Btrk1nStripLayer/F");

    background->SetBranchAddress("Btrk2nStripLayer",&Btrk2nStripLayerB);
    backgroundnew->Branch("Btrk2nStripLayer",&Btrk2nStripLayerflatB,"Btrk2nStripLayer/F");

    signal->SetBranchAddress("Btrk2nStripLayer",&Btrk2nStripLayerS);
    signalnew->Branch("Btrk2nStripLayer",&Btrk2nStripLayerflatS,"Btrk2nStripLayer/F");

    background->SetBranchAddress("Btrk1nPixelLayer",&Btrk1nPixelLayerB);
    backgroundnew->Branch("Btrk1nPixelLayer",&Btrk1nPixelLayerflatB,"Btrk1nPixelLayer/F");

    signal->SetBranchAddress("Btrk1nPixelLayer",&Btrk1nPixelLayerS);
    signalnew->Branch("Btrk1nPixelLayer",&Btrk1nPixelLayerflatS,"Btrk1nPixelLayer/F");

    background->SetBranchAddress("Btrk2nPixelLayer",&Btrk2nPixelLayerB);
    backgroundnew->Branch("Btrk2nPixelLayer",&Btrk2nPixelLayerflatB,"Btrk2nPixelLayer/F");

    signal->SetBranchAddress("Btrk2nPixelLayer",&Btrk2nPixelLayerS);
    signalnew->Branch("Btrk2nPixelLayer",&Btrk2nPixelLayerflatS,"Btrk2nPixelLayer/F");

    background->SetBranchAddress("Btrk2Pt",&Btrk2PtB);
    backgroundnew->Branch("Btrk2Pt",&Btrk2PtflatB,"Btrk2Pt/F");

    signal->SetBranchAddress("Btrk2Pt",&Btrk2PtS);
    signalnew->Branch("Btrk2Pt",&Btrk2PtflatS,"Btrk2Pt/F");

    background->SetBranchAddress("Balpha",&BalphaB);
    backgroundnew->Branch("Balpha",&BalphaflatB,"Balpha/F");

    signal->SetBranchAddress("Balpha",&BalphaS);
    signalnew->Branch("Balpha",&BalphaflatS,"Balpha/F");

    background->SetBranchAddress("Btrk2Eta",&Btrk2EtaB);
    backgroundnew->Branch("Btrk2Eta",&Btrk2EtaflatB,"Btrk2Eta/F");

    signal->SetBranchAddress("Btrk2Eta",&Btrk2EtaS);
    signalnew->Branch("Btrk2Eta",&Btrk2EtaflatS,"Btrk2Eta/F");

    background->SetBranchAddress("Bmu1isTriggered",&Bmu1isTriggeredB);
    backgroundnew->Branch("Bmu1isTriggered",&Bmu1isTriggeredflatB,"Bmu1isTriggered/O");

    signal->SetBranchAddress("Bmu1isTriggered",&Bmu1isTriggeredS);
    signalnew->Branch("Bmu1isTriggered",&Bmu1isTriggeredflatS,"Bmu1isTriggered/O");

    background->SetBranchAddress("Bmu2isTriggered",&Bmu2isTriggeredB);
    backgroundnew->Branch("Bmu2isTriggered",&Bmu2isTriggeredflatB,"Bmu2isTriggered/O");

    signal->SetBranchAddress("Bmu2isTriggered",&Bmu2isTriggeredS);
    signalnew->Branch("Bmu2isTriggered",&Bmu2isTriggeredflatS,"Bmu2isTriggered/O");

    background->SetBranchAddress("Bmu1TMOneStationTight",&Bmu1TMOneStationTightB);
    backgroundnew->Branch("Bmu1TMOneStationTight",&Bmu1TMOneStationTightflatB,"Bmu1TMOneStationTight/O");

    signal->SetBranchAddress("Bmu1TMOneStationTight",&Bmu1TMOneStationTightS);
    signalnew->Branch("Bmu1TMOneStationTight",&Bmu1TMOneStationTightflatS,"Bmu1TMOneStationTight/O");

    background->SetBranchAddress("Bmu2TMOneStationTight",&Bmu2TMOneStationTightB);
    backgroundnew->Branch("Bmu2TMOneStationTight",&Bmu2TMOneStationTightflatB,"Bmu2TMOneStationTight/O");

    signal->SetBranchAddress("Bmu2TMOneStationTight",&Bmu2TMOneStationTightS);
    signalnew->Branch("Bmu2TMOneStationTight",&Bmu2TMOneStationTightflatS,"Bmu2TMOneStationTight/O");

    background->SetBranchAddress("Bmu1isTrackerMuon",&Bmu1isTrackerMuonB);
    backgroundnew->Branch("Bmu1isTrackerMuon",&Bmu1isTrackerMuonflatB,"Bmu1isTrackerMuon/O");

    signal->SetBranchAddress("Bmu1isTrackerMuon",&Bmu1isTrackerMuonS);
    signalnew->Branch("Bmu1isTrackerMuon",&Bmu1isTrackerMuonflatS,"Bmu1isTrackerMuon/O");

    background->SetBranchAddress("Bmu2isTrackerMuon",&Bmu2isTrackerMuonB);
    backgroundnew->Branch("Bmu2isTrackerMuon",&Bmu2isTrackerMuonflatB,"Bmu2isTrackerMuon/O");

    signal->SetBranchAddress("Bmu2isTrackerMuon",&Bmu2isTrackerMuonS);
    signalnew->Branch("Bmu2isTrackerMuon",&Bmu2isTrackerMuonflatS,"Bmu2isTrackerMuon/O");

    background->SetBranchAddress("Bmu1isGlobalMuon",&Bmu1isGlobalMuonB);
    backgroundnew->Branch("Bmu1isGlobalMuon",&Bmu1isGlobalMuonflatB,"Bmu1isGlobalMuon/O");

    signal->SetBranchAddress("Bmu1isGlobalMuon",&Bmu1isGlobalMuonS);
    signalnew->Branch("Bmu1isGlobalMuon",&Bmu1isGlobalMuonflatS,"Bmu1isGlobalMuon/O");

    background->SetBranchAddress("Bmu2isGlobalMuon",&Bmu2isGlobalMuonB);
    backgroundnew->Branch("Bmu2isGlobalMuon",&Bmu2isGlobalMuonflatB,"Bmu2isGlobalMuon/O");

    signal->SetBranchAddress("Bmu2isGlobalMuon",&Bmu2isGlobalMuonS);
    signalnew->Branch("Bmu2isGlobalMuon",&Bmu2isGlobalMuonflatS,"Bmu2isGlobalMuon/O");

    background->SetBranchAddress("Btrk1highPurity",&Btrk1highPurityB);
    backgroundnew->Branch("Btrk1highPurity",&Btrk1highPurityflatB,"Btrk1highPurity/O");

    signal->SetBranchAddress("Btrk1highPurity",&Btrk1highPurityS);
    signalnew->Branch("Btrk1highPurity",&Btrk1highPurityflatS,"Btrk1highPurity/O");

    background->SetBranchAddress("Btrk2highPurity",&Btrk2highPurityB);
    backgroundnew->Branch("Btrk2highPurity",&Btrk2highPurityflatB,"Btrk2highPurity/O");

    signal->SetBranchAddress("Btrk2highPurity",&Btrk2highPurityS);
    signalnew->Branch("Btrk2highPurity",&Btrk2highPurityflatS,"Btrk2highPurity/O");

    background->SetBranchAddress("Bmu1InPixelLayer",&Bmu1InPixelLayerB);
    backgroundnew->Branch("Bmu1InPixelLayer",&Bmu1InPixelLayerflatB,"Bmu1InPixelLayer/I");

    signal->SetBranchAddress("Bmu1InPixelLayer",&Bmu1InPixelLayerS);
    signalnew->Branch("Bmu1InPixelLayer",&Bmu1InPixelLayerflatS,"Bmu1InPixelLayer/I");

    background->SetBranchAddress("Bmu1InStripLayer",&Bmu1InStripLayerB);
    backgroundnew->Branch("Bmu1InStripLayer",&Bmu1InStripLayerflatB,"Bmu1InStripLayer/I");

    signal->SetBranchAddress("Bmu1InStripLayer",&Bmu1InStripLayerS);
    signalnew->Branch("Bmu1InStripLayer",&Bmu1InStripLayerflatS,"Bmu1InStripLayer/I");

    background->SetBranchAddress("Bmu2InPixelLayer",&Bmu2InPixelLayerB);
    backgroundnew->Branch("Bmu2InPixelLayer",&Bmu2InPixelLayerflatB,"Bmu2InPixelLayer/I");

    signal->SetBranchAddress("Bmu2InPixelLayer",&Bmu2InPixelLayerS);
    signalnew->Branch("Bmu2InPixelLayer",&Bmu2InPixelLayerflatS,"Bmu2InPixelLayer/I");

    background->SetBranchAddress("Bmu2InStripLayer",&Bmu2InStripLayerB);
    backgroundnew->Branch("Bmu2InStripLayer",&Bmu2InStripLayerflatB,"Bmu2InStripLayer/I");

    signal->SetBranchAddress("Bmu2InStripLayer",&Bmu2InStripLayerS);
    signalnew->Branch("Bmu2InStripLayer",&Bmu2InStripLayerflatS,"Bmu2InStripLayer/I");

    backgroundhl->SetBranchAddress("HLT_HIL1DoubleMu0_v1",&HLT_HIL1DoubleMu0_v1B);
    backgroundnew->Branch("HLT_HIL1DoubleMu0_v1",&HLT_HIL1DoubleMu0_v1B,"HLT_HIL1DoubleMu0_v1/I");

    signalhl->SetBranchAddress("HLT_HIL1DoubleMu0_v1",&HLT_HIL1DoubleMu0_v1S);
    signalnew->Branch("HLT_HIL1DoubleMu0_v1",&HLT_HIL1DoubleMu0_v1S,"HLT_HIL1DoubleMu0_v1/I");

    backgroundskim->SetBranchAddress("pPAprimaryVertexFilter",&pPAprimaryVertexFilterB);
    backgroundnew->Branch("pPAprimaryVertexFilter",&pPAprimaryVertexFilterB,"pPAprimaryVertexFilter/I");

    signalskim->SetBranchAddress("pPAprimaryVertexFilter",&pPAprimaryVertexFilterS);
    signalnew->Branch("pPAprimaryVertexFilter",&pPAprimaryVertexFilterS,"pPAprimaryVertexFilter/I");

    backgroundskim->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilterB);
    backgroundnew->Branch("pBeamScrapingFilter",&pBeamScrapingFilterB,"pBeamScrapingFilter/I");

    signalskim->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilterS);
    signalnew->Branch("pBeamScrapingFilter",&pBeamScrapingFilterS,"pBeamScrapingFilter/I");

    int nsig = signal->GetEntries();
    int nback = background->GetEntries();
    int ngen = gen->GetEntries();

    //int nsig =  100;
    //int nback = 100;
    //int ngen = 100;
    
    cout<<nsig<<endl;
    cout<<nback<<endl;
    cout<<ngen<<endl;
   for (int k=0;k<nsig;k++){

        signal->GetEntry(k);
        signalroot->GetEntry(k);
        signalhl->GetEntry(k);
        signalskim->GetEntry(k);
        signalhi->GetEntry(k);
        for (int s=0;s<BsizeS;s++){

            auto passTrackingBPS = [&] (Tracking cut) {
                return (Btrk1PtErrS[s]/Btrk1PtS[s] < ptErr[cut] &&
                Btrk1Chi2ndfS[s]/(Btrk1nStripLayerS[s]+Btrk1nPixelLayerS[s])
              < chi2Nlayer[cut]);
            };

	        auto passTrackingBsS = [&] (Tracking cut) {
                return (Btrk1PtErrS[s]/Btrk1PtS[s] < ptErr[cut] &&
                  Btrk1Chi2ndfS[s]/(Btrk1nStripLayerS[s]+Btrk1nPixelLayerS[s])
                  < chi2Nlayer[cut] &&
                  Btrk2PtErrS[s]/Btrk2PtS[s] < ptErr[cut] &&
                  Btrk2Chi2ndfS[s]/(Btrk2nStripLayerS[s]+Btrk2nPixelLayerS[s])
                  < chi2Nlayer[cut]);
             };
            
            if (meson_n == 0){
                if (passTrackingBPS(Tracking::tight)) {
                    trackS = Tracking::tight;
                } else if (passTrackingBPS(Tracking::standard)) {
                    trackS = Tracking::standard;
                } else {
                    trackS = Tracking::loose;
                }
	         }
	        else {
                if (passTrackingBsS(Tracking::tight)) {
                    trackS = Tracking::tight;
                } else if (passTrackingBsS(Tracking::standard)) {
                    trackS = Tracking::standard;
                } else {
                    trackS = Tracking::loose;
                }
	        }

            Btrk1Dz1flatS=Btrk1Dz1S[s];
            Btrk1DzError1flatS=Btrk1DzError1S[s];
            Btrk2Dz1flatS=Btrk2Dz1S[s];
            Btrk2DzError1flatS=Btrk2DzError1S[s];
            Btrk1Dxy1flatS=Btrk1Dxy1S[s];
            Btrk1DxyError1flatS=Btrk1DxyError1S[s];
            Btrk2Dxy1flatS=Btrk2Dxy1S[s];
            Btrk2DxyError1flatS=Btrk2DxyError1S[s];
            BtktkmassflatS=BtktkmassS[s];
            BsvpvDistanceflatS=BsvpvDistanceS[s];
            BsvpvDisErrflatS=BsvpvDisErrS[s];
            Bd0flatS=Bd0S[s];
            BmassflatS=BmassS[s];
            BptflatS=BptS[s];
            Bmu1isTriggeredflatS=Bmu1isTriggeredS[s];
            Bmu2isTriggeredflatS=Bmu2isTriggeredS[s];
            Btrk1PtflatS=Btrk1PtS[s];
            Btrk1PtErrflatS=Btrk1PtErrS[s];
            Btrk2PtErrflatS=Btrk2PtErrS[s];
            Bchi2clflatS=Bchi2clS[s];
            Btrk1EtaflatS=Btrk1EtaS[s];
            ByflatS=ByS[s];
            BmumumassflatS=BmumumassS[s];
            Bmu1etaflatS=Bmu1etaS[s];
            Bmu1ptflatS=Bmu1ptS[s];
            Bmu2etaflatS=Bmu2etaS[s];
            Bmu2ptflatS=Bmu2ptS[s];
            Bmu1TMOneStationTightflatS=Bmu1TMOneStationTightS[s];
            Bmu2TMOneStationTightflatS=Bmu2TMOneStationTightS[s];
            Bmu1InPixelLayerflatS=Bmu1InPixelLayerS[s];
            Bmu1InStripLayerflatS=Bmu1InStripLayerS[s];
            Bmu2InPixelLayerflatS=Bmu2InPixelLayerS[s];
            Bmu2InStripLayerflatS=Bmu2InStripLayerS[s];
            Bmu1dxyPVflatS=Bmu1dxyPVS[s];
            Bmu2dxyPVflatS=Bmu2dxyPVS[s];
            Bmu1dzPVflatS=Bmu1dzPVS[s];
            Bmu2dzPVflatS=Bmu2dzPVS[s];
            Bmu1isTrackerMuonflatS=Bmu1isTrackerMuonS[s];
            Bmu2isTrackerMuonflatS=Bmu2isTrackerMuonS[s];
            Bmu1isGlobalMuonflatS=Bmu1isGlobalMuonS[s];
            Bmu2isGlobalMuonflatS=Bmu2isGlobalMuonS[s];
            Btrk1highPurityflatS=Btrk1highPurityS[s];
            Btrk2highPurityflatS=Btrk2highPurityS[s];
            Btrk1PixelHitflatS=Btrk1PixelHitS[s];
            Btrk2PixelHitflatS=Btrk2PixelHitS[s];
            Btrk1StripHitflatS=Btrk1StripHitS[s];
            Btrk2StripHitflatS=Btrk2StripHitS[s];
            Btrk1Chi2ndfflatS=Btrk1Chi2ndfS[s];
            Btrk2Chi2ndfflatS=Btrk2Chi2ndfS[s];
            Btrk2PtflatS=Btrk2PtS[s];
            BalphaflatS=BalphaS[s];
            Btrk2EtaflatS=Btrk2EtaS[s];
            Btrk1nStripLayerflatS=Btrk1nStripLayerS[s];
            Btrk2nStripLayerflatS=Btrk2nStripLayerS[s];
            Btrk1nPixelLayerflatS=Btrk1nPixelLayerS[s];
            Btrk2nPixelLayerflatS=Btrk2nPixelLayerS[s];
            BgenflatS=BgenS[s];
            BgenptflatS=BgenptS[s];
            BgenyflatS=BgenyS[s];
            signalnew->Fill();
            
        }
        ;
    } 
    for (int k=0;k<ngen;k++){
        gen->GetEntry(k);
        signalhi->GetEntry(k);
        signal->GetEntry(k);
        gennew->Fill();
    }

   for (int k=0;k<nback;k++){
        background->GetEntry(k);
        backgroundroot->GetEntry(k);
        backgroundhl->GetEntry(k);
        backgroundskim->GetEntry(k);
        for (int s=0;s<BsizeB;s++){


            auto passTrackingBPB = [&] (Tracking cut) {
                return (Btrk1PtErrB[s]/Btrk1PtB[s] < ptErr[cut] &&
                Btrk1Chi2ndfB[s]/(Btrk1nStripLayerB[s]+Btrk1nPixelLayerB[s])
              < chi2Nlayer[cut]);
            };

	        auto passTrackingBsB = [&] (Tracking cut) {
                return (Btrk1PtErrB[s]/Btrk1PtB[s] < ptErr[cut] &&
                  Btrk1Chi2ndfB[s]/(Btrk1nStripLayerB[s]+Btrk1nPixelLayerB[s])
                  < chi2Nlayer[cut] &&
                  Btrk2PtErrB[s]/Btrk2PtB[s] < ptErr[cut] &&
                  Btrk2Chi2ndfB[s]/(Btrk2nStripLayerB[s]+Btrk2nPixelLayerB[s])
                  < chi2Nlayer[cut]);
             };
            
            if (meson_n == 0){
                if (passTrackingBPB(Tracking::tight)) {
                    trackB = Tracking::tight;
                } else if (passTrackingBPB(Tracking::standard)) {
                    trackB = Tracking::standard;
                } else {
                    trackB = Tracking::loose;
                }
	         }
	        else {
                if (passTrackingBsB(Tracking::tight)) {
                    trackB = Tracking::tight;
                } else if (passTrackingBsB(Tracking::standard)) {
                    trackB = Tracking::standard;
                } else {
                    trackB = Tracking::loose;
                }
	        }

            Btrk1Dz1flatB=Btrk1Dz1B[s];
            Btrk1DzError1flatB=Btrk1DzError1B[s];
            Btrk2Dz1flatB=Btrk2Dz1B[s];
            Btrk2DzError1flatB=Btrk2DzError1B[s];
            Btrk1Dxy1flatB=Btrk1Dxy1B[s];
            Btrk1DxyError1flatB=Btrk1DxyError1B[s];
            Btrk2Dxy1flatB=Btrk2Dxy1B[s];
            Btrk2DxyError1flatB=Btrk2DxyError1B[s];
            BtktkmassflatB=BtktkmassB[s];
            BsvpvDistanceflatB=BsvpvDistanceB[s];
            BsvpvDisErrflatB=BsvpvDisErrB[s];
            Bd0flatB=Bd0B[s];
            BmassflatB=BmassB[s];
            BptflatB=BptB[s];
            Bmu1isTriggeredflatB=Bmu1isTriggeredB[s];
            Bmu2isTriggeredflatB=Bmu2isTriggeredB[s];
            Btrk1PtflatB=Btrk1PtB[s];
            Btrk1PtErrflatB=Btrk1PtErrB[s];
            Btrk2PtErrflatB=Btrk2PtErrB[s];
            Bchi2clflatB=Bchi2clB[s];
            Btrk1EtaflatB=Btrk1EtaB[s];
            ByflatB=ByB[s];
            BmumumassflatB=BmumumassB[s];
            Bmu1etaflatB=Bmu1etaB[s];
            Bmu1ptflatB=Bmu1ptB[s];
            Bmu2etaflatB=Bmu2etaB[s];
            Bmu2ptflatB=Bmu2ptB[s];
            Bmu1TMOneStationTightflatB=Bmu1TMOneStationTightB[s];
            Bmu2TMOneStationTightflatB=Bmu2TMOneStationTightB[s];
            Bmu1InPixelLayerflatB=Bmu1InPixelLayerB[s];
            Bmu1InStripLayerflatB=Bmu1InStripLayerB[s];
            Bmu2InPixelLayerflatB=Bmu2InPixelLayerB[s];
            Bmu2InStripLayerflatB=Bmu2InStripLayerB[s];
            Bmu1dxyPVflatB=Bmu1dxyPVB[s];
            Bmu2dxyPVflatB=Bmu2dxyPVB[s];
            Bmu1dzPVflatB=Bmu1dzPVB[s];
            Bmu2dzPVflatB=Bmu2dzPVB[s];
            Bmu1isTrackerMuonflatB=Bmu1isTrackerMuonB[s];
            Bmu2isTrackerMuonflatB=Bmu2isTrackerMuonB[s];
            Bmu1isGlobalMuonflatB=Bmu1isGlobalMuonB[s];
            Bmu2isGlobalMuonflatB=Bmu2isGlobalMuonB[s];
            Btrk1highPurityflatB=Btrk1highPurityB[s];
            Btrk2highPurityflatB=Btrk2highPurityB[s];
            Btrk1PixelHitflatB=Btrk1PixelHitB[s];
            Btrk2PixelHitflatB=Btrk2PixelHitB[s];
            Btrk1StripHitflatB=Btrk1StripHitB[s];
            Btrk2StripHitflatB=Btrk2StripHitB[s];
            Btrk1Chi2ndfflatB=Btrk1Chi2ndfB[s];
            Btrk2Chi2ndfflatB=Btrk2Chi2ndfB[s];
            Btrk2PtflatB=Btrk2PtB[s];
            BalphaflatB=BalphaB[s];
            Btrk2EtaflatB=Btrk2EtaB[s];
            Btrk1nStripLayerflatB=Btrk1nStripLayerB[s];
            Btrk2nStripLayerflatB=Btrk2nStripLayerB[s];
            Btrk1nPixelLayerflatB=Btrk1nPixelLayerB[s];
            Btrk2nPixelLayerflatB=Btrk2nPixelLayerB[s];
            backgroundnew->Fill();
        }
        
    } 
    
    
    ROOT::RDataFrame dS(*signalnew);
    ROOT::RDataFrame dB(*backgroundnew);
    ROOT::RDataFrame dG(*gennew);
    cout<<"still fine"<<endl;
        
    auto cut = "";
    if (meson_n==0){
        cut = "(pPAprimaryVertexFilter == 1 && pBeamScrapingFilter == 1 && HLT_HIL1DoubleMu0_v1 == 1)  &&  (Bmu1isTriggered == 1 && Bmu2isTriggered == 1 ) && (Btrk1Pt > 0.5 && Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2.0 && Bpt > 2 && abs(Btrk1Eta-0.0) < 2.4  && (TMath::Abs(By)<2.4&&TMath::Abs(Bmumumass-3.096916)<0.15&&((abs(Bmu1eta)<1.2&&Bmu1pt>3.5)||(abs(Bmu1eta)>1.2&&abs(Bmu1eta)<2.1&&Bmu1pt>(5.47-1.89*abs(Bmu1eta)))||(abs(Bmu1eta)>2.1&&abs(Bmu1eta)<2.4&&Bmu1pt>1.5))&&((abs(Bmu2eta)<1.2&&Bmu2pt>3.5)||(abs(Bmu2eta)>1.2&&abs(Bmu2eta)<2.1&&Bmu2pt>(5.47-1.89*abs(Bmu2eta)))||(abs(Bmu2eta)>2.1&&abs(Bmu2eta)<2.4&&Bmu2pt>1.5))&&Bmu1TMOneStationTight&&Bmu2TMOneStationTight&&Bmu1InPixelLayer>0&&(Bmu1InPixelLayer+Bmu1InStripLayer)>5&&Bmu2InPixelLayer>0&&(Bmu2InPixelLayer+Bmu2InStripLayer)>5&&Bmu1dxyPV<0.3&&Bmu2dxyPV<0.3&&Bmu1dzPV<20&&Bmu2dzPV<20&&Bmu1isTrackerMuon&&Bmu2isTrackerMuon&&Bmu1isGlobalMuon&&Bmu2isGlobalMuon&&Btrk1highPurity&&abs(Btrk1Eta)<2.4&&Btrk1Pt>0.5)  && (Btrk1PixelHit + Btrk1StripHit > 10) &&  (Btrk1PtErr/Btrk1Pt < 0.1)&& Btrk1Chi2ndf/(Btrk1nStripLayer+Btrk1nPixelLayer) < 0.18   && (abs(PVz)<15))";
    }
    else {
        cut = "(abs(Btktkmass-1.019455)<0.015)&&(((((abs(Btktkmass-1.019455)<0.015)&& TMath::Abs(Bmumumass-3.096916)<0.15 && Bpt > 0 && Bpt < 5 && (abs(Btrk1Eta)<2.4 && abs(Btrk2Eta)<2.4 && Btrk1Pt>0.0 && Btrk2Pt>0.0) && Btrk1Pt > 0.5 && Btrk2Pt > 0.5  && Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2.0)  && ( (Bpt < 2 && Bpt > 0) || (Bpt < 3 && Bpt > 2) || (Bpt < 5 && Bpt > 3)  )))  ||  ( Bpt > 3 && ((pPAprimaryVertexFilter == 1 && pBeamScrapingFilter == 1 && HLT_HIL1DoubleMu0_v1 == 1 && (abs(PVz)<15))  &&  (Bmu1isTriggered == 1 && Bmu2isTriggered == 1 ) &&  (Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2.0)    && (TMath::Abs(By)<2.4&&TMath::Abs(Bmumumass-3.096916)<0.15&&((abs(Bmu1eta)<1.2&&Bmu1pt>3.5)||(abs(Bmu1eta)>1.2&&abs(Bmu1eta)<2.1&&Bmu1pt>(5.47-1.89*abs(Bmu1eta)))||(abs(Bmu1eta)>2.1&&abs(Bmu1eta)<2.4&&Bmu1pt>1.5))&&((abs(Bmu2eta)<1.2&&Bmu2pt>3.5)||(abs(Bmu2eta)>1.2&&abs(Bmu2eta)<2.1&&Bmu2pt>(5.47-1.89*abs(Bmu2eta)))||(abs(Bmu2eta)>2.1&&abs(Bmu2eta)<2.4&&Bmu2pt>1.5))&&Bmu1InPixelLayer>0&&(Bmu1InPixelLayer+Bmu1InStripLayer)>5&&Bmu2InPixelLayer>0&&(Bmu2InPixelLayer+Bmu2InStripLayer)>5&&Bmu1dxyPV<0.3&&Bmu2dxyPV<0.3&&Bmu1dzPV<20&&Bmu2dzPV<20&&Bmu1isTrackerMuon&&Bmu2isTrackerMuon&&Bmu1isGlobalMuon&&Bmu2isGlobalMuon)  && ( Btrk1Pt > 0.5 && Btrk2Pt > 0.5 && abs(Btrk1Eta-0.0) < 2.4 && abs(Btrk2Eta-0.0) < 2.4  && Btrk1highPurity  && Btrk2highPurity  && Btrk1PixelHit + Btrk1StripHit > 10  && Btrk2PixelHit + Btrk2StripHit > 10) &&  (Btrk1PtErr/Btrk1Pt < 0.1)  &&  (Btrk2PtErr/Btrk2Pt < 0.1)    && Btrk1Chi2ndf/(Btrk1nStripLayer+Btrk1nPixelLayer) < 0.18   && Btrk2Chi2ndf/(Btrk2nStripLayer+Btrk2nPixelLayer) < 0.18 )))";
    }
    cout<<"still fine"<<endl;
    auto dS_cut=dS.Filter(cut);
    cout<<"still fine"<<endl;
    auto dB_cut=dB.Filter(cut);
    cout<<"still fine"<<endl;
    if (meson_n==0){
        dS_cut.Snapshot("ntKp",outS);
        dG.Snapshot("ntGen",outG);
        dB_cut.Snapshot("ntKp",outB);
    }
    else {
        dS_cut.Snapshot("ntphi",outS);
        dG.Snapshot("ntGen",outG);
        dB_cut.Snapshot("ntphi",outB);
    }
    
    cout<<"still fine final"<<endl;
    cout<<*(dS_cut.Count())<<endl;
    cout<<*(dB_cut.Count())<<endl;
    outfS->Close();
    outfG->Close();
    outfB->Close();
}