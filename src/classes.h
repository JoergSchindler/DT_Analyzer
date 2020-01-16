#include <vector>
#include "Demo/DemoAnalyzer/plugins/RecoObjects.h"
#include "Demo/DemoAnalyzer/plugins/Objects.h"

namespace {
  struct dictionary {

    //Structures                                                                
    RecoJetType dummy0;
    CaloJetType dummy1;
    GenPType dummy2;
    JetType dummy3;
    PFCandidateType dummy4;
    DT4DSegmentType dummy5;
    CSCSegmentType dummy6;
    TrackType dummy7;
    muonDtChamber dummy8;
    muonDtLayer dummy9;

    //Vector of structures                                                      
    std::vector<RecoJetType> dummyVector0;
    std::vector<CaloJetType> dummyVector1;
    std::vector<GenPType> dummyVector2;
    std::vector<JetType> dummyVector3;
    std::vector<PFCandidateType> dummyVector4;
    std::vector<DT4DSegmentType> dummyVector5;
    std::vector<CSCSegmentType> dummyVector6;
    std::vector<TrackType> dummyVector7;
    std::vector<muonDtChamber> dummyVector8;
    std::vector<muonDtLayer> dummyVector9;
  };
}
