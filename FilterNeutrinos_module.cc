////////////////////////////////////////////////////////////////////////
// Class:       FilterNeutrinos
// Plugin Type: filter (art v2_05_01)
// File:        FilterNeutrinos_module.cc
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

class FilterNeutrinos;


class FilterNeutrinos : public art::EDFilter {
public:
  explicit FilterNeutrinos(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  FilterNeutrinos(FilterNeutrinos const &) = delete;
  FilterNeutrinos(FilterNeutrinos &&) = delete;
  FilterNeutrinos & operator = (FilterNeutrinos const &) = delete;
  FilterNeutrinos & operator = (FilterNeutrinos &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  TTree* _tree;

  std::string fTTreeName, fTFileName, fDirName;

  Int_t _run, _evt;

};


FilterNeutrinos::FilterNeutrinos(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{

  fTTreeName = p.get<std::string>("TTreeName");
  fTFileName = p.get<std::string>("TFileName");
  fDirName   = p.get<std::string>("DirName");

}

bool FilterNeutrinos::filter(art::Event & e)
{

  auto run = e.run();
  auto evt = e.event();
  
  for (int n=0; n < _tree->GetEntries(); n++) {
    _tree->GetEntry(n);
    if ( ( ((unsigned int)_run) == run) && ( ((unsigned int)_evt) == evt) )
      return true;
  }// for all TTree entries

  return false;
}

void FilterNeutrinos::beginJob()
{

  // load ROOT file which contains events to filter
  TFile* f = new TFile(fTFileName.c_str());
  f->cd(fDirName.c_str());
  std::string treepath = fDirName+"/"+fTTreeName;
  _tree = (TTree*)f->Get(treepath.c_str());
  _tree->SetBranchAddress("_run",&_run);
  _tree->SetBranchAddress("_evt",&_evt);
  
}

void FilterNeutrinos::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(FilterNeutrinos)
