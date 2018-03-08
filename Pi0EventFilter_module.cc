////////////////////////////////////////////////////////////////////////
// Class:       Pi0EventFilter
// Plugin Type: filter (art v2_09_06)
// File:        Pi0EventFilter_module.cc
//
// Generated at Thu Mar  8 06:36:16 2018 by David Caratelli using cetskelgen
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

#include <TFile.h>
#include <TTree.h>


class Pi0EventFilter;


class Pi0EventFilter : public art::EDFilter {
public:
  explicit Pi0EventFilter(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  Pi0EventFilter(Pi0EventFilter const &) = delete;
  Pi0EventFilter(Pi0EventFilter &&) = delete;
  Pi0EventFilter & operator = (Pi0EventFilter const &) = delete;
  Pi0EventFilter & operator = (Pi0EventFilter &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  TTree* _tree;

  std::string fTTreeName, fTFileName;

  Int_t run, evt;


};


Pi0EventFilter::Pi0EventFilter(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  fTTreeName = p.get<std::string>("TTreeName");
  fTFileName = p.get<std::string>("TFileName");
}

bool Pi0EventFilter::filter(art::Event & e)
{



  for (int n=0; n < _tree->GetEntries(); n++) {
    _tree->GetEntry(n);
    if ( ( ((unsigned int)run) == e.run()) && ( ((unsigned int)evt) == e.event()) )
      return true;
  }// for all TTree entries

  return false;

}

void Pi0EventFilter::beginJob()
{
  // load ROOT file which contains events to filter
  TFile* f = new TFile(fTFileName.c_str());
  std::string treepath = fTTreeName;
  _tree = (TTree*)f->Get(treepath.c_str());
  _tree->SetBranchAddress("run",&run);
  _tree->SetBranchAddress("evt",&evt);
}

void Pi0EventFilter::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(Pi0EventFilter)
