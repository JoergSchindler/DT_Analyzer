#ifndef GENANALYZER_H
#define GENANALYZER_H

#include <iostream>
#include <cmath>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHECommonBlocks.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "TFile.h"
#include "TH3.h"
#include "TF1.h"
#include "TKey.h"

class GenAnalyzer {
    public:
        GenAnalyzer(edm::ParameterSet&, edm::ConsumesCollector&&);
        ~GenAnalyzer();
        virtual std::vector<reco::GenParticle> FillGenVectorByIdAndStatus(const edm::Event&, int, int);

        virtual std::vector<reco::GenParticle> FillGenVectorByIdStatusAndMotherAndKin(const edm::Event&, int, int, int, float, float);

        virtual std::vector<reco::GenParticle> FillGenVectorByIdAndStatusAndKin(const edm::Event&, int, int, float, float);

     private:

    edm::EDGetTokenT<std::vector<reco::GenParticle> > GenParticlesToken;
    std::vector<int> ParticleList;
    std::vector<int> ParticleStatus;
    edm::Handle<std::vector<reco::GenParticle> > GenCollection;
};

#endif