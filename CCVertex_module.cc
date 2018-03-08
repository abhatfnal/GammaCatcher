////////////////////////////////////////////////////////////////////////
// Class:       CCVertex
// Plugin Type: filter (art v2_09_06)
// File:        CCVertex_module.cc
//
// Generated at Mon Mar  5 14:14:25 2018 by David Caratelli using cetskelgen
// from cetlib version v3_01_03.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

#include "art/Framework/Services/Optional/TFileService.h"

#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardata/Utilities/AssociationUtil.h"

class CCVertex;


class CCVertex : public art::EDFilter {
public:
  explicit CCVertex(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CCVertex(CCVertex const &) = delete;
  CCVertex(CCVertex &&) = delete;
  CCVertex & operator = (CCVertex const &) = delete;
  CCVertex & operator = (CCVertex &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  // producers
  std::string fAssnProducer;

};


CCVertex::CCVertex(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  produces< std::vector< recob::Vertex > >();
  fAssnProducer = p.get<std::string>("AssnProducer");
}

bool CCVertex::filter(art::Event & e)
{

  // produce vertex
  std::unique_ptr< std::vector<recob::Vertex> > Vtx_v(new std::vector<recob::Vertex>);

  Double_t xyz[3] = {};

  art::Handle< art::Assns<recob::Vertex,recob::Track,void> > numuCCassn_h;
  e.getByLabel(fAssnProducer,numuCCassn_h);
  if (numuCCassn_h->size() != 1){
    std::cout << "Number of vertices != 1 -> ERROR ERROR ERROR" << std::endl;
    return false;
  }
  
  numuCCassn_h->at(0).first->XYZ(xyz);

  recob::Vertex vtx(xyz);
  Vtx_v->emplace_back(vtx);
  
  e.put(std::move(Vtx_v));

  return true;
}

void CCVertex::beginJob()
{
  // Implementation of optional member function here.
}

void CCVertex::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(CCVertex)
