////////////////////////////////////////////////////////////////////////
// Class:       GammaCatcher
// Plugin Type: producer (art v2_05_01)
// File:        GammaCatcher_module.cc
//
// Generated at Sun Jan 28 22:25:13 2018 by David Caratelli using cetskelgen
// David Caratelli - davidc@fnal.gov
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Services
#include "art/Framework/Services/Optional/TFileService.h"

// Data Products
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RawData/RawDigit.h"

// C++
#include <memory>

// ROOT
#include <TTree.h>

class GammaCatcher;


class GammaCatcher : public art::EDProducer {
public:
  explicit GammaCatcher(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  GammaCatcher(GammaCatcher const &) = delete;
  GammaCatcher(GammaCatcher &&) = delete;
  GammaCatcher & operator = (GammaCatcher const &) = delete;
  GammaCatcher & operator = (GammaCatcher &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  /**
  Calculate baseline and RMS for a waveform's ADCs
  */
  std::pair<float,float> getBaselineRMS(const std::vector<short>& wf);

  // TTree where to store variables
  TTree* _chan_tree;
  int    _chan;
  float  _base;
  float  _rms;
  int    _run, _sub, _evt;

private:

  // Declare member data here.
  std::string fRawDigitProducer;

};


GammaCatcher::GammaCatcher(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{

  produces< std::vector< recob::Hit > >();
  
  // grab from fhicl file:
  fRawDigitProducer = p.get<std::string>("RawDigitProducer");

}

void GammaCatcher::produce(art::Event & e)
{

  _evt  = e.event();
  _sub  = e.subRun();
  _run  = e.run();

  // produce Hit objects
  std::unique_ptr< std::vector<recob::Hit> > Hit_v(new std::vector<recob::Hit>);

  // load RawDigits
  art::Handle<std::vector<raw::RawDigit> > rawdigit_h;
  e.getByLabel(fRawDigitProducer,rawdigit_h);

  // start looping through raw data
  for (auto const& rawdigit : *rawdigit_h) {

    _chan = rawdigit.Channel();

    // collection-plane only
    if (_chan < 4800) continue;
    
    auto channelvitals = getBaselineRMS(rawdigit.ADCs());
    _base = channelvitals.first;
    _rms  = channelvitals.second;
    
    _chan_tree->Fill();

    if (rawdigit.Channel() == 5000)
      std::cout << "DAVIDC REACHED CHANNEL 5000! DAVIDC" << std::endl;

  }// for all RawDigit objects
  
  e.put(std::move(Hit_v));

}

std::pair<float,float> GammaCatcher::getBaselineRMS(const std::vector<short>& wf) {

  // to do:
  // scan only portion of waveform
  // define acceptable noise level
  // if noise below -> return RMS and baseline
  // otherwise move to next segment in waveform.

  auto adc_v = wf;
  std::sort(adc_v.begin(), adc_v.end());

  // truncate top 10% of values [to remove real pulses]
  std::vector<short> truncated_adc(adc_v.begin(), adc_v.end() - adc_v.size()/10);

  float base  = 0.;
  float rms   = 0.;

  for (auto const& adc : truncated_adc)
    base += adc;
  base = base / truncated_adc.size();
  for (auto const& adc : truncated_adc)
    rms += (adc-base)*(adc-base);
  rms = sqrt( rms / ( truncated_adc.size() - 1) );

  return std::make_pair(base,rms);

}

void GammaCatcher::beginJob()
{

  // set TTree branches
  art::ServiceHandle<art::TFileService> tfs;
  _chan_tree = tfs->make<TTree>("_chan_tree","Channel Info TTree");
  _chan_tree->Branch("_chan",&_chan,"chan/I");
  _chan_tree->Branch("_base",&_base,"base/F");
  _chan_tree->Branch("_rms" ,&_rms ,"rms/F");
  _chan_tree->Branch("_evt" ,&_evt ,"evt/I");
  _chan_tree->Branch("_sub" ,&_sub ,"sub/I");
  _chan_tree->Branch("_run" ,&_run ,"run/I");
  

  
}

void GammaCatcher::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(GammaCatcher)
