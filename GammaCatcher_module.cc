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
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/GeometryCore.h"

// Services
#include "art/Framework/Services/Optional/TFileService.h"

// Data Products
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RawData/RawDigit.h"

// C++
#include <memory>

// ROOT
#include <TTree.h>

// Algorithms
#include "Algorithms/HitFinding.h"
#include "Algorithms/ProximityClusterer.h"

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

  // TTree where to store variables
  TTree* _chan_tree;
  int    _chan;
  float  _base;
  float  _rms;
  int    _run, _sub, _evt;

  // HitFinding class
  gammacatcher::HitFinding* _HitFinding;
  // Proximity clusterer class
  gammacatcher::ProximityClusterer* _ProximityClusterer;

private:

  // Declare member data here.
  std::string fRawDigitProducer;
  // NSigma for hit threshold on RMS noise
  double fNSigma;
  // minimum number of ticks above threshold to have a hit
  int fMinTickWidth;
  
  /**
     Make clusters
   */
  void MakeClusters(const std::unique_ptr< std::vector<recob::Hit> >& hits,
		    const std::vector<std::vector<unsigned int> >& cluster_idx_v,
		    const std::unique_ptr< std::vector<recob::Cluster> >& clusters);

};


GammaCatcher::GammaCatcher(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{

  produces< std::vector< recob::Hit > >();
  produces< std::vector< recob::Cluster > >();
  
  // grab from fhicl file:
  fRawDigitProducer = p.get<std::string>("RawDigitProducer");
  fNSigma           = p.get<double>     ("NSigma");
  fMinTickWidth     = p.get<int>        ("MinTickWidth");

}

void GammaCatcher::produce(art::Event & e)
{

  _evt  = e.event();
  _sub  = e.subRun();
  _run  = e.run();

  // geometry service
  art::ServiceHandle<geo::Geometry> geo;
  //geo::GeometryCore const* geom = lar::providerFrom<geo::Geometry>();

  // produce Hit objects
  std::unique_ptr< std::vector<recob::Hit> > Hit_v(new std::vector<recob::Hit>);
  // produce Cluster objects
  std::unique_ptr< std::vector<recob::Cluster> > Cluster_v(new std::vector<recob::Cluster>);


  // load RawDigits
  art::Handle<std::vector<raw::RawDigit> > rawdigit_h;
  e.getByLabel(fRawDigitProducer,rawdigit_h);

  // start looping through raw data
  for (auto const& rawdigit : *rawdigit_h) {

    _chan = rawdigit.Channel();

    // collection-plane only
    if (geo->View(_chan) != geo::kW) continue;
    
    auto channelvitals = _HitFinding->getBaselineRMS(rawdigit.ADCs());
    _base = channelvitals.first;
    _rms  = channelvitals.second;
    
    // get hits from channel
    auto hits = _HitFinding->getHits(rawdigit.ADCs(),_base,_rms);

    // loop through hits and create larsoft hits
    for (auto const& hit : hits) {
      recob::Hit arthit(_chan, (int)hit.tstart, (int)hit.tend, hit.time,
			0., _rms, hit.ampl, 0., hit.area, hit.area,
			0., 0., 0., 0., 0., 
			geo->View(_chan), geo->SignalType(_chan), geo->ChannelToWire(_chan)[0] );

      Hit_v->emplace_back(arthit);
    }// for all hits created

    _chan_tree->Fill();

  }// for all RawDigit objects

  // cluster index vectors will be stored here
  std::vector<std::vector<unsigned int> > cluster_v;
  // cluster hits together
  _ProximityClusterer->cluster(Hit_v,cluster_v);

  // go through indices and make clusters
  MakeClusters(Hit_v, cluster_v, Cluster_v);
	       

  std::cout << "DAVIDC created " << Hit_v->size() << " hits in event" << std::endl;
  std::cout << "DAVIDC created " << cluster_v.size() << " clusters in event" << std::endl;
  std::cout << "DAVIDC created " << Cluster_v->size() << " clusters in event" << std::endl;
  
  e.put(std::move(Hit_v));
  e.put(std::move(Cluster_v));

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
  
  _HitFinding = new gammacatcher::HitFinding();
  _HitFinding->setNSigma(fNSigma);
  _HitFinding->setMinTickWidth(fMinTickWidth);

  _ProximityClusterer = new gammacatcher::ProximityClusterer();
  _ProximityClusterer->initialize();
  
}

void GammaCatcher::endJob()
{
  // Implementation of optional member function here.
}


void GammaCatcher::MakeClusters(const std::unique_ptr< std::vector<recob::Hit> >& hits,
				const std::vector<std::vector<unsigned int> >& cluster_v,
				const std::unique_ptr< std::vector<recob::Cluster> >& clusters)
{

  clusters->clear();

  for (size_t n=0; n < cluster_v.size(); n++) {
    auto const& clus_idx_v = cluster_v.at(n);
    if (clus_idx_v.size() == 0) continue;
    float w_min = 9999;
    float t_min = 9999;
    float w_max = -1;
    float t_max = -1;
    float integral = 0;
    for (auto const& hit_idx : clus_idx_v) {
      auto const& hit = hits->at(hit_idx);
      float t = hit.PeakTime();
      float w = hit.WireID().Wire;
      if (t > t_max) t_max = t;
      if (t < t_min) t_min = t;
      if (w > w_max) w_max = w;
      if (w < w_min) w_min = w;
      integral += hit.Integral();
    }// for all hits
    recob::Cluster clus(w_min, 0., t_min, 0., 0., 0., 0., 
			w_max, 0., w_max, 0., 0., 0., 0., 
			integral, 0., integral, 0., 
			clus_idx_v.size(), 0., 0., n,
			hits->at(clus_idx_v[0]).View(),
			geo::PlaneID());

    clusters->emplace_back(clus);
    
  }// for all clusters

  return;
}

DEFINE_ART_MODULE(GammaCatcher)
