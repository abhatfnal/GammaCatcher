////////////////////////////////////////////////////////////////////////
// Class:       SkimNeutrinos
// Plugin Type: analyzer (art v2_05_01)
// File:        SkimNeutrinos_module.cc
//
// Generated at Fri Feb  2 17:04:02 2018 by David Caratelli using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Services/Optional/TFileService.h"

#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardata/Utilities/AssociationUtil.h"

#include <TTree.h>

class SkimNeutrinos;


class SkimNeutrinos : public art::EDAnalyzer {
public:
  explicit SkimNeutrinos(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SkimNeutrinos(SkimNeutrinos const &) = delete;
  SkimNeutrinos(SkimNeutrinos &&) = delete;
  SkimNeutrinos & operator = (SkimNeutrinos const &) = delete;
  SkimNeutrinos & operator = (SkimNeutrinos &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  // producers
  std::string fAssnProducer;
  
  // TTree where to add event/vertex information
  TTree* _tree;
  int _evt;
  int _sub;
  int _run;
  float _xpos;
  float _ypos;
  float _zpos;

};


SkimNeutrinos::SkimNeutrinos(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  
{

  fAssnProducer    = p.get<std::string>("AssnProducer"   );

}

void SkimNeutrinos::analyze(art::Event const & e)
{

  Double_t xyz[3] = {};

  art::Handle< art::Assns<recob::Vertex,recob::Track,void> > numuCCassn_h;
  e.getByLabel(fAssnProducer,numuCCassn_h);
  if (numuCCassn_h->size() != 1){
    std::cout << "Number of vertices != 1 -> ERROR ERROR ERROR" << std::endl;
    return;
  }

  numuCCassn_h->at(0).first->XYZ(xyz);
  
  _run = e.run();
  _sub = e.subRun();
  _evt = e.event();

  _xpos = xyz[0];
  _ypos = xyz[1];
  _zpos = xyz[2];
  
  _tree->Fill();
}

void SkimNeutrinos::beginJob()
{

  art::ServiceHandle<art::TFileService> tfs;
  _tree = tfs->make<TTree>("_tree","Neutrino Info TTree");
  _tree->Branch("_evt",&_evt,"evt/I");
  _tree->Branch("_sub",&_sub,"sub/I");
  _tree->Branch("_run",&_run,"run/I");
  _tree->Branch("_xpos",&_xpos,"xpos/F");
  _tree->Branch("_ypos",&_ypos,"ypos/F");
  _tree->Branch("_zpos",&_zpos,"zpos/F");

}

void SkimNeutrinos::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(SkimNeutrinos)
