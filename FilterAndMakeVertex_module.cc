////////////////////////////////////////////////////////////////////////
// Class:       FilterAndMakeVertex
// Plugin Type: filter (art v2_05_01)
// File:        FilterAndMakeVertex_module.cc
//
// Generated at Sat Feb  3 12:29:04 2018 by David Caratelli using cetskelgen
// from cetlib version v1_21_00.
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

#include <TFile.h>
#include <TTree.h>

#include "lardataobj/RecoBase/Vertex.h"

class FilterAndMakeVertex;


class FilterAndMakeVertex : public art::EDFilter {
public:
  explicit FilterAndMakeVertex(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  FilterAndMakeVertex(FilterAndMakeVertex const &) = delete;
  FilterAndMakeVertex(FilterAndMakeVertex &&) = delete;
  FilterAndMakeVertex & operator = (FilterAndMakeVertex const &) = delete;
  FilterAndMakeVertex & operator = (FilterAndMakeVertex &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  TTree* _tree;

  std::string fTTreeName, fTFileName, fDirName;

  Int_t run, evt;
  Float_t x, y, z;

};


FilterAndMakeVertex::FilterAndMakeVertex(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{

  produces< std::vector< recob::Vertex > >();

  fTTreeName = p.get<std::string>("TTreeName");
  fTFileName = p.get<std::string>("TFileName");
  fDirName   = p.get<std::string>("DirName");

}

bool FilterAndMakeVertex::filter(art::Event & e)
{


  std::unique_ptr< std::vector<recob::Vertex> > Vertex_v(new std::vector<recob::Vertex>);

  auto erun = e.run();
  auto eevt = e.event();

  bool found = false;
  
  for (int n=0; n < _tree->GetEntries(); n++) {
    _tree->GetEntry(n);
    if ( ( ((unsigned int)run) == erun) && ( ((unsigned int)evt) == eevt) ) {
      found = true;
      std::cout << "FOUND EVENT! FILTERING" << std::endl;
      break;
    }
  }// for all TTree entries
  
  if (found) {
    Double_t xyz[3] = {};
    xyz[0] = (double)x;
    xyz[1] = (double)y;
    xyz[2] = (double)z;
    recob::Vertex vtx(xyz);
    Vertex_v->emplace_back(vtx);
  }

  e.put(std::move(Vertex_v));
  
  if (found)
    return true;
      
  return false;
}

void FilterAndMakeVertex::beginJob()
{

  // load ROOT file which contains events to filter
  TFile* f = new TFile(fTFileName.c_str());
  f->cd(fDirName.c_str());
  std::string treepath = fDirName+"/"+fTTreeName;
  _tree = (TTree*)f->Get(treepath.c_str());
  _tree->SetBranchAddress("run",&run);
  _tree->SetBranchAddress("evt",&evt);
  _tree->SetBranchAddress("x",&x);
  _tree->SetBranchAddress("y",&y);
  _tree->SetBranchAddress("z",&z);
  
}

void FilterAndMakeVertex::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(FilterAndMakeVertex)
