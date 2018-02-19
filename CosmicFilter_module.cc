////////////////////////////////////////////////////////////////////////
// Class:       CosmicFilter
// Plugin Type: producer (art v2_05_01)
// File:        CosmicFilter_module.cc
//
// Generated at Sun Feb 18 20:57:09 2018 by David Caratelli using cetskelgen
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

#include <memory>
#include <map>

#include "art/Framework/Services/Optional/TFileService.h"

#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/GeometryCore.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/Utilities/AssociationUtil.h"

class CosmicFilter;


class CosmicFilter : public art::EDProducer {
public:
  explicit CosmicFilter(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CosmicFilter(CosmicFilter const &) = delete;
  CosmicFilter(CosmicFilter &&) = delete;
  CosmicFilter & operator = (CosmicFilter const &) = delete;
  CosmicFilter & operator = (CosmicFilter &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;


private:

  // Declare member data here.

  /**
     Fill channel -> vector of SSNet hit indices map
   */
  void FillChannelMap(const art::ValidHandle<std::vector<::recob::Hit> > hit_h);

  // channel list for SSNet hits connecting channel number to list of hit indices
  std::map<unsigned int, std::vector<size_t> > _hitmap;

  Float_t _xpos, _ypos, _zpos; // xyz of vertex

  // producers
  std::string fVtxProducer, fTrkProducer, fHitProducer;
  // veto radius -> surrounds nu vtx and is used to identify crossing cosmic-rays
  double fVetoRadius;
  // square of radius for faster computation
  double fVetoRadiusSq;
  // minimum impact paramter for track to be considered cosmic
  double fIPmin;

  // vector of track-like hit indices
  std::vector<size_t> _trkhits;

  /**
     Return number of track - sphere intersection points and 
     minimum distance of track-points to sphere
   */
  std::pair<int,float> SphereIntersection(const recob::Track& trk);

  /**
     Square distance between point and reco'd vertex.
   */
  double SqDist(const TVector3& pt);

};


CosmicFilter::CosmicFilter(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  produces< std::vector< recob::Hit > >();

  fHitProducer = p.get<std::string>("HitProducer");
  fVtxProducer = p.get<std::string>("VtxProducer");
  fTrkProducer = p.get<std::string>("TrkProducer");
  fVetoRadius  = p.get<double>     ("VetoRadius" );
  fIPmin       = p.get<double>     ("IPmin"      );

  fVetoRadiusSq    = fVetoRadius * fVetoRadius;
}

void CosmicFilter::produce(art::Event & e)
{
  // Implementation of required member function here.

  // grab tracks
  auto const& trk_h = e.getValidHandle<std::vector<recob::Track>>(fTrkProducer);
  // grab ssnet hits
  auto const& hit_h = e.getValidHandle<std::vector<recob::Hit>>(fHitProducer);
  // grab vertex
  auto const& vtx_h = e.getValidHandle<std::vector<recob::Vertex>>(fVtxProducer);

  // grab hits associated to track
  art::FindManyP<recob::Hit> trk_hit_assn_v(trk_h, e, fTrkProducer);

  // produce hits
  std::unique_ptr< std::vector<recob::Hit> > Hit_v(new std::vector<recob::Hit>);

  // clear track-like hit vector
  _trkhits.clear();

  // BEGIN : LOAD VERTEX
  // check that only one vertex exists
  if (vtx_h->size() != 1) 
    std::cout << "\t\t DD \t\t ERROR" << std::endl;
  auto vtx = vtx_h->at(0);
  Double_t xyz[3] = {};
  vtx.XYZ(xyz);
  _xpos = xyz[0];
  _ypos = xyz[1];
  _zpos = xyz[2];
  // END : LOAD VERTEX

  // BEGIN : PERFORM HIT MATCHING
  // strategy:
  // fill map which links each channel with the vector of 
  // indices of the SSNet hits on that channel.
  FillChannelMap(hit_h);
  // END : PERFORM HIT MATCHING

  // BEGIN : IDENTIFY COSMIC TRACK HITS
  for (size_t t=0; t < trk_h->size(); t++) {

    auto const& trk = trk_h->at(t);
    auto trkdist = SphereIntersection(trk);
    
    // if no intersections -> skip
    if (trkdist.first == 0) continue;

    // if a single intersection -> check IP
    // if smaller then IP min -> neutrino track -> skip
    if ( (trkdist.first == 1) && (trkdist.second < fIPmin) ) continue;

    // in all other cases, track is cosmic-like
    // grab associated hits and compare to SSNet hits
    // if matched -> tag as one to be removed
    const std::vector<art::Ptr<recob::Hit> > hit_v = trk_hit_assn_v.at(t);
    for (size_t h=0; h < hit_v.size(); h++) {
      art::Ptr<recob::Hit> hit = hit_v.at(h);
      // if the hit channel is in the SSNet hit map:
      if (_hitmap.find( hit->Channel() ) != _hitmap.end() ){
	auto const& hitidx_v = _hitmap[ hit->Channel() ];
	for (auto const& idx : hitidx_v) {
	  // compare hit information
	  if ( hit_h->at(idx).PeakTime() == hit->PeakTime() )
	    _trkhits.push_back( idx ); // save idx of SSNet hit to be removed
	}// for all hit indices associated to this channel
      }// if the hit channel is in the SSNet hit map
    }// for all hits associated to track
  }// for all tracks
  // END : IDENTIFY COSMIC TRACK HITS

  // finally, save hits not identified as track-like
  for (size_t idx=0; idx < hit_h->size(); idx++) {
    // has this index been flagged?
    if (std::find(_trkhits.begin(),_trkhits.end(),idx) == _trkhits.end() )
      Hit_v->emplace_back(hit_h->at(idx));
  }// for all track hit indices

  e.put(std::move(Hit_v));

}

void CosmicFilter::beginJob()
{
  // Implementation of optional member function here.
}

void CosmicFilter::endJob()
{
  // Implementation of optional member function here.
}

void CosmicFilter::FillChannelMap(const art::ValidHandle<std::vector<::recob::Hit> > hit_h) {

  _hitmap.clear();

  for (size_t h=0; h < hit_h->size(); h++){

    unsigned int channel = hit_h->at(h).Channel();
    
    if (_hitmap.find(channel) == _hitmap.end() ){
      std::vector<size_t> chlist = {h};
      _hitmap[channel] = chlist;
    }// if entry did not exist
    else {
      _hitmap[channel].push_back( h );
    }// append to already existing list of hits
    
  }// for all SSNet hits

  return;
}

std::pair<int,float> CosmicFilter::SphereIntersection(const recob::Track& trk) {
  
  // loop through all points along the track
  // and calculate how many times the sphere radius is crossed
  // as well as the minimum distance to the vertex

  int ncross   = 0; // number of times the track crosses the bounding radius
  float dvtxSq = 1e6; // min distance of any point to the sphere [squared for computation]

  // are we inside or outside of the sphere? helps keep track
  // of whether we have stepped in/out
  bool insphere = false;
  // step sampled
  size_t nstep = 0;

  size_t ptn = 0; // loop trhough track points
  while (ptn < trk.NumberTrajectoryPoints()) {
    //auto const& vp = trk.NextValidPoint(validpoint);
    if (trk.HasValidPoint(ptn) == false) { ptn += 1; continue; }
    auto const& pt = trk.LocationAtPoint( ptn );
    ptn += 1;
    auto dSq = SqDist(pt);
    if (dSq < dvtxSq) { dvtxSq = dSq; }
    if (dSq < fVetoRadiusSq) {
      if (insphere == false) {
	insphere = true;
	if (nstep != 0) { ncross += 1; }
      }// if we just crossed in the sphere!
    }// if in the sphere
    else { // if we are outside of the sphere
      if (insphere == true) {
	ncross += 1;
	insphere = false;
      }// if we were inside the sphere and just crossed out
    }// if outside of the sphere
    nstep += 1;
  }// while looping through track points
  
  return std::make_pair(ncross,sqrt(dvtxSq));
}

// distance between point and vertex
double CosmicFilter::SqDist(const TVector3& pt) {

  return (_xpos - pt[0]) * (_xpos - pt[0]) + (_ypos - pt[1]) * (_ypos - pt[1]) + (_zpos - pt[2]) * (_zpos - pt[2]);
}

DEFINE_ART_MODULE(CosmicFilter)
