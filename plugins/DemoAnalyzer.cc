// -*- C++ -*-
//
// Package:    Demo/DemoAnalyzer
// Class:      DemoAnalyzer
// 
/**\class DemoAnalyzer DemoAnalyzer.cc Demo/DemoAnalyzer/plugins/DemoAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Joerg Christoph Schindler
//         Created:  Fri, 29 Nov 2019 12:53:21 GMT
//
//


// system include files
#include <memory>
#include <string>
#include <vector>
#include <tuple>
#include <fstream>
#include <TRandom3.h>
// user include files


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackCandidate/interface/TrackCandidate.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackReco/interface/TrackExtraBase.h"


#include "CommonTools/UtilAlgos/interface/TFileService.h"



#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"


#include "DataFormats/DetId/interface/DetId.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/DigiSimLinks/interface/DTDigiSimLink.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/StripDigiSimLink.h"
#include "DataFormats/MuonData/interface/MuonDigiCollection.h"

#include "TTree.h"


#include "SlimmedGenAnalyzer.h"
#include "Objects.h"
#include "ObjectsFormat.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

struct muonDtLayer{
  int id;
  int wheel;
  int sector;
  int station;
  int superlayer;
  int layer;
  int nSimHits;
  int nRecHits;


};

class DemoAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit DemoAnalyzer(const edm::ParameterSet&);
      ~DemoAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      // ----------member data ---------------------------
      edm::EDGetTokenT<std::vector<PSimHit>> MuonDTSimHitsToken_;
      edm::EDGetTokenT<MuonDigiCollection<DTLayerId,DTDigiSimLink>> simMuonDTDigisToken_;
      edm::EDGetTokenT<DTRecHitCollection> dtRechitInputToken_;
      
      edm::Handle<std::vector<PSimHit>> MuonDTSimHits;
      edm::Handle<MuonDigiCollection<DTLayerId,DTDigiSimLink>> MuonDTDigi;
      edm::Handle<DTRecHitCollection> dtRechits; 
      GenAnalyzer* theGenAnalyzer;
      edm::ParameterSet GenPSet;
   
      std::vector<GenPType> GenBquarks;
      std::vector<GenPType> GenLLPs;
      double MinGenBpt, MaxGenBeta;
      int nDtRecHits;
      int nDtSimHits;
      int nmatched;
      int nGenBquarks;
      long int nGenLL;
      int nGenLL_inside_DT;
      float gen_b_radius;
      int nDtDetLayer;
      
      int dtDetLayer[15000];
      
      float DtSimHitsPhi[15000];
      float DtSimHitsEta[15000];
      int DtSimHitsStation[15000];
      int DtSimHitsSector[15000];
      int DTSimHitsWheel[15000];
      int DTSimHitSuperLayer[15000];
      int DTSimHitLayer[15000];
      int DTSimHitWire[15000];
      int DTSimHitDetId[15000];
      int nDtSimHitsperLayer[15000];
      
      int DtRecHitsStation[15000];
      int DtRecHitsSector[15000];
      int DTRecHitsWheel[15000];
      int DTRecHitSuperLayer[15000];
      int DTRecHitLayer[15000];
      int DTRecHitWire[15000];
      int DTRecHitDetId[15000];
      int nDtRecHitsperLayer[15000];
      
      float DtSimHits_match_gParticle_minDeltaR[15000];
      float DtSimHits_match_gParticle_index[15000];
      float GenLL_r[2];
     
      TTree *displacedJetMuonTree;


      TH1F *NEvents;

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
DemoAnalyzer::DemoAnalyzer(const edm::ParameterSet& iConfig):
   MuonDTSimHitsToken_(consumes<std::vector<PSimHit>>(edm::InputTag("g4SimHits","MuonDTHits","SIM"))),
   dtRechitInputToken_(consumes<DTRecHitCollection>(edm::InputTag("dt1DRecHits"))),
   GenPSet(iConfig.getParameter<edm::ParameterSet>("genSet")),
   MinGenBpt(iConfig.getParameter<double>("minGenBpt")),
   MaxGenBeta(iConfig.getParameter<double>("maxGenBeta"))

{
   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::Service<TFileService> fs;

   displacedJetMuonTree = fs->make<TTree>("llp", "selected AOD information for llp analyses");
   NEvents = fs->make<TH1F>("NEvents",";;NEvents;",1,-0.5,0.5);
 
   displacedJetMuonTree->Branch("nmatched",&nmatched,"nmatched/I");  
   displacedJetMuonTree->Branch("nDtRecHits",&nDtRecHits,"nDtRecHits/I");
   displacedJetMuonTree->Branch("nDtSimHits",&nDtSimHits,"nDtSimHits/I");
   displacedJetMuonTree->Branch("nDtDetLayer",&nDtDetLayer,"nDtDetLayer/I");
   

   displacedJetMuonTree->Branch("DtSimHits_match_gParticle_minDeltaR",&DtSimHits_match_gParticle_minDeltaR,"DtSimHits_match_gParticle_minDeltaR[nDtSimHits]/I");
   displacedJetMuonTree->Branch("DtSimHits_match_gParticle_index",&DtSimHits_match_gParticle_index,"DtSimHits_match_gParticle_index[nDtSimHits]/F");
   displacedJetMuonTree->Branch("nGenLL_inside_DT",&nGenLL_inside_DT,"nGenLL_inside_DT/I");
   
   displacedJetMuonTree->Branch("DtRecHitsStation",&DtRecHitsStation,"DtRecHitsStation[nDtRecHits]/I");
   displacedJetMuonTree->Branch("DtRecHitsSector",&DtRecHitsSector,"DtRecHitsSector[nDtRecHits]/I");
   displacedJetMuonTree->Branch("DTRecHitsWheel",&DTRecHitsWheel,"DTRecHitsWheel[nDtRecHits]/I");
   displacedJetMuonTree->Branch("DTRecHitSuperLayer",&DTRecHitSuperLayer,"DTRecHitSuperLayer[nDtRecHits]/I");
   displacedJetMuonTree->Branch("DTRecHitLayer",&DTRecHitLayer,"DTRecHitLayer[nDtRecHits]/I");
   displacedJetMuonTree->Branch("DTRecHitWire",&DTRecHitWire,"DTRecHitWire[nDtRecHits]/I");
            
   displacedJetMuonTree->Branch("DtSimHitsStation",&DtSimHitsStation,"DtSimHitsStation[nDtSimHits]/I");
   displacedJetMuonTree->Branch("DtSimHitsSector",&DtSimHitsSector,"DtSimHitsSector[nDtSimHits]/I");
   displacedJetMuonTree->Branch("DTSimHitsWheel",&DTSimHitsWheel,"DTSimHitsWheel[nDtSimHits]/I");
   displacedJetMuonTree->Branch("DTSimHitSuperLayer",&DTSimHitSuperLayer,"DTSimHitSuperLayer[nDtSimHits]/I");
   displacedJetMuonTree->Branch("DTSimHitLayer",&DTSimHitLayer,"DTSimHitLayer[nDtSimHits]/I");
   displacedJetMuonTree->Branch("DTSimHitWire",&DTSimHitWire,"DTSimHitWire[nDtSimHits]/I");
   
   displacedJetMuonTree->Branch("nDtSimHitsperLayer",&nDtSimHitsperLayer,"nDtSimHitsperLayer[nDtDetLayer]/I");
   displacedJetMuonTree->Branch("nDtRecHitsperLayer",&nDtRecHitsperLayer,"nDtRecHitsperLayer[nDtDetLayer]/I");
   
   displacedJetMuonTree->Branch("GenLLPs", &GenLLPs);
   displacedJetMuonTree->Branch("GenBquarks", &GenBquarks);
   displacedJetMuonTree->Branch("nGenBquarks" , &nGenBquarks , "nGenBquarks/I");
   displacedJetMuonTree->Branch("nGenLL" , &nGenLL , "nGenLL/I");
   displacedJetMuonTree->Branch("gen_b_radius" , &gen_b_radius , "gen_b_radius/F");
   
   theGenAnalyzer         = new GenAnalyzer(GenPSet, consumesCollector());

}


DemoAnalyzer::~DemoAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DemoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   NEvents->Fill(0); 
   GenLLPs.clear();
   GenBquarks.clear();
   iEvent.getByToken(dtRechitInputToken_,dtRechits);
   iEvent.getByToken(MuonDTSimHitsToken_, MuonDTSimHits);
   edm::ESHandle<DTGeometry> dtG;
   iSetup.get<MuonGeometryRecord>().get(dtG);
   

  
   std::vector<reco::GenParticle> GenLongLivedVect = theGenAnalyzer->FillGenVectorByIdAndStatus(iEvent,9000006,22);
   // std::vector<reco::GenParticle> GenHiggsVect = theGenAnalyzer->FillGenVectorByIdAndStatus(iEvent,25,22);
   nGenLL = GenLongLivedVect.size();
   std::vector<reco::GenParticle> GenBquarksVect;
   if(nGenLL>0)
   {
   GenBquarksVect = theGenAnalyzer->FillGenVectorByIdStatusAndMotherAndKin(iEvent,5,23,9000006,float(MinGenBpt),float(MaxGenBeta));
   }
   else
   {
   GenBquarksVect = theGenAnalyzer->FillGenVectorByIdAndStatusAndKin(iEvent,5,23,float(MinGenBpt),float(MaxGenBeta));
   }
   nGenLL_inside_DT = 0;
   double radius_temp = -999.;
   for(unsigned int i =0; i<GenBquarksVect.size(); i++)
   {
       double radius = sqrt( pow(GenBquarksVect[i].vx(),2) + pow(GenBquarksVect[i].vy(),2) );
//        std::cout<< radius << std::endl;
       
       if (abs(GenBquarksVect[i].vz())<660 && radius < 738 && radius> 402)
       {
           if(radius!=radius_temp)
           {    
                radius_temp = radius;
                nGenLL_inside_DT++;
           }
       }
       
   }
   nGenBquarks = GenBquarksVect.size();

   for(unsigned int i = 0; i < GenLongLivedVect.size(); i++) GenLLPs.push_back( GenPType() );
   for(unsigned int i = 0; i < GenBquarksVect.size(); i++) GenBquarks.push_back( GenPType() );

   if(nGenBquarks>0) gen_b_radius = GenBquarksVect.at(0).mother()? sqrt(pow(GenBquarksVect.at(0).vx() - GenBquarksVect.at(0).mother()->vx(),2) + pow(GenBquarksVect.at(0).vy() - GenBquarksVect.at(0).mother()->vy(),2) + pow(GenBquarksVect.at(0).vz() - GenBquarksVect.at(0).mother()->vz(),2)) : -1.;

    
   
   std::vector<int> detLayers;
    //std::cout << "Number of dt rec hits: "<<dtRechits->size()<<std::endl;
   nDtRecHits=0;

   for (const DTRecHit1DPair dtRechit : *dtRechits) {
       
        // LocalError  segmentLocalDirectionError = iDT->localDirectionError();
        DetId geoid = dtRechit.geographicalId();
        DTChamberId dtdetid = DTChamberId(geoid);
        DTSuperLayerId slid = DTSuperLayerId(geoid);
        DTLayerId layerid = DTLayerId(geoid);
        DTWireId wireid = DTWireId(geoid);
        
        const DTChamber * dtchamber = dtG->chamber(dtdetid);
   
        if (dtchamber) {
            //DtRecHitsStation[nDtRecHits] = dtchamber->id().station();
            DTRecHitDetId[nDtRecHits]= dtdetid.rawId();
             if (!(std::find(detLayers.begin(), detLayers.end(), DTRecHitDetId[nDtRecHits]) != detLayers.end()))
             {

              detLayers.push_back(DTRecHitDetId[nDtRecHits]);

             }

            
            
            DtRecHitsStation[nDtRecHits] = dtdetid.station();
            DtRecHitsSector[nDtRecHits] = dtdetid.sector();
            DTRecHitsWheel[nDtRecHits] = dtdetid.wheel();
            DTRecHitSuperLayer[nDtRecHits] = slid.superLayer();
            DTRecHitLayer[nDtRecHits]= layerid.layer();
            DTRecHitWire[nDtRecHits] = wireid.wire();
            nDtRecHits++;
        }
  }
  nDtSimHits=0;

  //std::cout<<"number of DT simhits: "<<MuonDTSimHits->size()<<std::endl;
  for(size_t i=0; i<MuonDTSimHits->size();i++)
  {
	Local3DPoint DTSimHitLocalPosition = (*MuonDTSimHits)[i].localPosition();
	//DetId geoid = (*MuonDTSimHits)[i].geographicalId();
	DetId dtdetid = (DetId)(*MuonDTSimHits)[i].detUnitId();
       
        DTSuperLayerId slid = DTSuperLayerId(dtdetid);
        DTLayerId layerid = DTLayerId(dtdetid);
        DTWireId wireid = DTWireId(dtdetid);
        
	DTChamberId dtchamberdetid = DTChamberId(dtdetid);
	const DTChamber * dtchamber = dtG->chamber(dtdetid);
        
        
	if (dtchamber) {
	GlobalPoint globalPosition = dtchamber->toGlobal(DTSimHitLocalPosition);
        DTSimHitDetId[i]= dtchamberdetid.rawId();
             if (!(std::find(detLayers.begin(), detLayers.end(), DTRecHitDetId[nDtRecHits]) != detLayers.end()))
             {

              detLayers.push_back(DTRecHitDetId[nDtRecHits]);

             }
        
        DtSimHitsStation[i] = dtchamberdetid.station();
        DtSimHitsSector[i] = dtchamberdetid.sector();
        DTSimHitsWheel[i] = dtchamberdetid.wheel();
        DTSimHitSuperLayer[i] = slid.superLayer();
        DTSimHitLayer[i]= layerid.layer();
        DTSimHitWire[i] = wireid.wire();
	DtSimHitsPhi[i] = globalPosition.phi();
	DtSimHitsEta[i] = globalPosition.eta();
        nDtSimHits++;
	}
  }
  std::vector<muonDtLayer> dtLayers;
  nDtDetLayer = detLayers.size();
  
      for(unsigned int i = 0; i < GenLongLivedVect.size(); i++) ObjectsFormat::FillGenPType(GenLLPs[i], &GenLongLivedVect[i]);
      for(unsigned int i = 0; i < GenBquarksVect.size(); i++) ObjectsFormat::FillGenPType(GenBquarks[i], &GenBquarksVect[i]);

  for( int i = 0; i < nDtDetLayer; i++)
  {
    muonDtLayer tmpLayer;
    tmpLayer.nSimHits = 0;
    tmpLayer.nRecHits = 0;
    
    for (int j = 0; j < nDtRecHits; j++)
    {
      if (detLayers[i] == DTRecHitDetId[j])
      tmpLayer.nRecHits++;
    }
    for (int j = 0; j < nDtSimHits; j++)
    {
      if (detLayers[i] ==  DTSimHitDetId[j])
      {
        tmpLayer.nSimHits++;
      }

    }
    tmpLayer.id = detLayers[i];
    dtLayers.push_back(tmpLayer);
  }
  
  for(unsigned int i = 0; i < dtLayers.size();i++)
  {

    dtDetLayer[i] = dtLayers[i].id;
    nDtRecHitsperLayer[i]  = dtLayers[i].nRecHits;
    nDtSimHitsperLayer[i]  = dtLayers[i].nSimHits;
      
    }   
      
      
  displacedJetMuonTree->Fill();

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
DemoAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DemoAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DemoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DemoAnalyzer);
