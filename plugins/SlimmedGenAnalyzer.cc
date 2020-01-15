#include "SlimmedGenAnalyzer.h"


GenAnalyzer::GenAnalyzer(edm::ParameterSet& PSet, edm::ConsumesCollector&& CColl):
    // GenToken(CColl.consumes<GenEventInfoProduct>(PSet.getParameter<edm::InputTag>("genProduct"))),
    // LheToken(CColl.consumes<LHEEventProduct>(PSet.getParameter<edm::InputTag>("lheProduct"))),
    GenParticlesToken(CColl.consumes<std::vector<reco::GenParticle> >(PSet.getParameter<edm::InputTag>("genParticles"))),
    ParticleList(PSet.getParameter<std::vector<int> >("pdgId")),
    ParticleStatus(PSet.getParameter<std::vector<int> >("status"))
    // SampleDYJetsToLL(PSet.getParameter<std::vector<std::string> >("samplesDYJetsToLL")),
    // SampleZJetsToNuNu(PSet.getParameter<std::vector<std::string> >("samplesZJetsToNuNu")),
    // SampleWJetsToLNu(PSet.getParameter<std::vector<std::string> >("samplesWJetsToLNu")),
    // SampleDir(PSet.getParameter<std::string>("samplesDir")),
    // Sample(PSet.getParameter<std::string>("sample")),
    // EWKFileName(PSet.getParameter<std::string>("ewkFile")),
    // ApplyEWK(PSet.getParameter<bool>("applyEWK")),
    // ApplyTopPtReweigth(PSet.getParameter<bool>("applyTopPtReweigth")),
    // PythiaLOSample(PSet.getParameter<bool>("pythiaLOSample"))
{
    // for(unsigned int i = 0; i < SampleDYJetsToLL.size(); i++) {
    //     Files[SampleDYJetsToLL[i]] = new TFile((SampleDir+SampleDYJetsToLL[i]+".root").c_str(), "READ");
    //     hPartons[SampleDYJetsToLL[i]] = (TH1F*)Files[SampleDYJetsToLL[i]]->Get("counter/c_lhePartons");
    //     hBPartons[SampleDYJetsToLL[i]] = (TH1F*)Files[SampleDYJetsToLL[i]]->Get("counter/c_lheBPartons");
    //     hHT[SampleDYJetsToLL[i]] = (TH1F*)Files[SampleDYJetsToLL[i]]->Get("counter/c_lheHT");
    //     hPtV[SampleDYJetsToLL[i]] = (TH1F*)Files[SampleDYJetsToLL[i]]->Get("counter/c_lhePtZ");
    // }
    // for(unsigned int i = 0; i < SampleZJetsToNuNu.size(); i++) {
    //     Files[SampleZJetsToNuNu[i]] = new TFile((SampleDir+SampleZJetsToNuNu[i]+".root").c_str(), "READ");
    //     hPartons[SampleZJetsToNuNu[i]] = (TH1F*)Files[SampleZJetsToNuNu[i]]->Get("counter/c_lhePartons");
    //     hBPartons[SampleZJetsToNuNu[i]] = (TH1F*)Files[SampleZJetsToNuNu[i]]->Get("counter/c_lheBPartons");
    //     hHT[SampleZJetsToNuNu[i]] = (TH1F*)Files[SampleZJetsToNuNu[i]]->Get("counter/c_lheHT");
    //     hPtV[SampleZJetsToNuNu[i]] = (TH1F*)Files[SampleZJetsToNuNu[i]]->Get("counter/c_lhePtZ");
    // }
    // for(unsigned int i = 0; i < SampleWJetsToLNu.size(); i++) {
    //     Files[SampleWJetsToLNu[i]] = new TFile((SampleDir+SampleWJetsToLNu[i]+".root").c_str(), "READ");
    //     hPartons[SampleWJetsToLNu[i]] = (TH1F*)Files[SampleWJetsToLNu[i]]->Get("counter/c_lhePartons");
    //     hBPartons[SampleWJetsToLNu[i]] = (TH1F*)Files[SampleWJetsToLNu[i]]->Get("counter/c_lheBPartons");
    //     hHT[SampleWJetsToLNu[i]] = (TH1F*)Files[SampleWJetsToLNu[i]]->Get("counter/c_lheHT");
    //     hPtV[SampleWJetsToLNu[i]] = (TH1F*)Files[SampleWJetsToLNu[i]]->Get("counter/c_lhePtZ");
    // }
    
    // EWKFile = new TFile(EWKFileName.c_str(), "READ");
    
    // fZEWK = (TF1*)EWKFile->Get("z_ewkcorr/z_ewkcorr_func");
    // fWEWK = (TF1*)EWKFile->Get("w_ewkcorr/w_ewkcorr_func");
    
    /*
    Sample=sample;
    isDYFile=false;
    Sum=Num=NULL;
    // Read root file
    DYFile=new TFile("data/DYWeight.root", "READ");
    DYFile->cd();
    if(!DYFile->IsZombie()) {
      Num=(TH3F*)DYFile->Get(Sample.c_str());
      Sum=(TH3F*)DYFile->Get("Sum");
      if(Sum && Num && Sum->GetEntries()>0 && Num->GetEntries()>0) {
        npMax=Sum->GetXaxis()->GetBinCenter(Sum->GetNbinsX());
        ptMax=Sum->GetYaxis()->GetBinCenter(Sum->GetNbinsY());
        htMax=Sum->GetZaxis()->GetBinCenter(Sum->GetNbinsZ());
        isDYFile=true;
      }
      else std::cout << " - GenAnalyzer Warning: Drell-Yan initialization failed, check rootfile" << std::endl;
      
      
    }
    else std::cout << " - GenAnalyzer Warning: No Drell-Yan File" << std::endl;
    */
    
    // PU reweighting
//    LumiWeights=new edm::LumiReWeighting("data/MC_True.root", "data/Prod6.root", "S10", "pileup");
    
    std::cout << " --- GenAnalyzer initialization ---" << std::endl;
    // std::cout << "  sample            :\t" << Sample << std::endl;
    // if(ApplyEWK) std::cout << "  EWK file          :\t" << EWKFileName << std::endl;
    std::cout << std::endl;
}

GenAnalyzer::~GenAnalyzer() {
//    delete LumiWeights;
    // for(auto const &it : Files) it.second->Close();
    // EWKFile->Close();
}


std::vector<reco::GenParticle> GenAnalyzer::FillGenVectorByIdStatusAndMotherAndKin(const edm::Event& iEvent, int partid, int partstatus, int motherid, float pt, float eta) {

    std::vector<reco::GenParticle> Vect;

    // check if is real data
    // isRealData = iEvent.isRealData();
    // if(isRealData or PythiaLOSample) return Vect;
    // fill collection for this event 
    iEvent.getByToken(GenParticlesToken, GenCollection);
    // Loop on Gen Particles collection
    for(std::vector<reco::GenParticle>::const_iterator it = GenCollection->begin(); it != GenCollection->end(); ++it) {
      if(abs(it->pdgId()) == partid && (it->status()) == partstatus && fabs(it->mother()->pdgId()) == motherid && (it->pt())>pt && fabs(it->eta())<fabs(eta)) Vect.push_back(*it); // Fill vector
    }
    return Vect;
}


std::vector<reco::GenParticle> GenAnalyzer::FillGenVectorByIdAndStatusAndKin(const edm::Event& iEvent, int partid, int partstatus, float pt, float eta) {

    std::vector<reco::GenParticle> Vect;

    // check if is real data
    // isRealData = iEvent.isRealData();
    // if(isRealData or PythiaLOSample) return Vect;
    // fill collection for this event 
    iEvent.getByToken(GenParticlesToken, GenCollection);
    // Loop on Gen Particles collection
    for(std::vector<reco::GenParticle>::const_iterator it = GenCollection->begin(); it != GenCollection->end(); ++it) {
        if(abs(it->pdgId()) == partid && (it->status()) == partstatus && (it->pt())>pt && fabs(it->eta())<fabs(eta) ) Vect.push_back(*it); // Fill vector
    }
//    std::cout << "\n\n\n" << std::endl;
    return Vect;
}

std::vector<reco::GenParticle> GenAnalyzer::FillGenVectorByIdAndStatus(const edm::Event& iEvent, int partid, int partstatus) {

    std::vector<reco::GenParticle> Vect;

    // check if is real data
    // isRealData = iEvent.isRealData();
    // if(isRealData or PythiaLOSample) return Vect;
    // fill collection for this event 
    iEvent.getByToken(GenParticlesToken, GenCollection);
    // Loop on Gen Particles collection
    for(std::vector<reco::GenParticle>::const_iterator it = GenCollection->begin(); it != GenCollection->end(); ++it) {
      //for(unsigned int i = 0; i < ParticleList.size(); i++) {

      //for(unsigned int s = 0; s < ParticleStatus.size(); s++) {
            if(abs(it->pdgId()) == partid && (it->status()) == partstatus) Vect.push_back(*it); // Fill vector
	    //  }
	    //}
    }
//    std::cout << "\n\n\n" << std::endl;
    return Vect;
}

