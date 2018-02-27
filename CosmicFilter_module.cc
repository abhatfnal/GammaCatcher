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

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Cluster.h"
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
  // map coonecting PFP associated track key to hit keys
  std::map<size_t, std::vector< art::Ptr<recob::Hit> > > _pfpmap;

  Float_t _xpos, _ypos, _zpos; // xyz of vertex

  // producers
  std::string fPFPProducer, fVtxProducer, fTrkProducer, fCluProducer, fHitProducer;
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
  fCluProducer = p.get<std::string>("CluProducer");
  fPFPProducer = p.get<std::string>("PFPProducer");
  fVetoRadius  = p.get<double>     ("VetoRadius" );
  fIPmin       = p.get<double>     ("IPmin"      );

  fVetoRadiusSq    = fVetoRadius * fVetoRadius;
}

void CosmicFilter::produce(art::Event & e)
{
  // Implementation of required member function here.

  // grab pfparticles
  auto const& pfp_h = e.getValidHandle<std::vector<recob::PFParticle>>(fPFPProducer);
  // grab tracks
  auto const& trk_h = e.getValidHandle<std::vector<recob::Track>>(fTrkProducer);
  // grab clusters
  auto const& clu_h = e.getValidHandle<std::vector<recob::Cluster>>(fCluProducer);
  // grab ssnet hits
  auto const& hit_h = e.getValidHandle<std::vector<recob::Hit>>(fHitProducer);
  // grab vertex
  auto const& vtx_h = e.getValidHandle<std::vector<recob::Vertex>>(fVtxProducer);

  // grab tracks and clusters associated to PFParticle
  art::FindManyP<recob::Track  > pfp_trk_assn_v(pfp_h, e, fPFPProducer);
  art::FindManyP<recob::Cluster> pfp_clu_assn_v(pfp_h, e, fPFPProducer);

  // grab hits associated to clusters
  art::FindManyP<recob::Hit> clu_hit_assn_v(clu_h, e, fCluProducer);

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
  
  // BEGIN : LOOP THROUGH ALL PFParticles
  for (size_t p=0; p < pfp_h->size(); p++) {

    auto const& pfp = pfp_h->at(p);

    // skip non-primary pfparticles
    if (pfp.IsPrimary() == false) continue;
    
    // is it muon-like? if not skip
    if (pfp.PdgCode() != 13) continue;
    
    // grab associated track ID
    const std::vector<art::Ptr<recob::Track> > pfp_trk_v = pfp_trk_assn_v.at(p);
    if (pfp_trk_v.size() != 1) 
      std::cout << "\t\t DD \t\t PFP associated to != 1 track" << std::endl;
    
    // if no tracks associated -> skip
    if (pfp_trk_v.size() == 0) continue;

    // get track key
    auto trkKey = pfp_trk_v.at(0).key();
    _pfpmap[trkKey] = std::vector<art::Ptr<recob::Hit>>{};

    // get daughters
    std::vector<size_t> daughters = pfp.Daughters();

    // find associated PFParticle daughters which are electron-like (delta-ray)
    for (size_t pp=0; pp < pfp_h->size(); pp++) {

      auto const& pfp2 = pfp_h->at(pp);

      // search through daughters
      for (auto const& daughterID : daughters) {
	if ( (pfp2.Self() == daughterID) && (pfp2.PdgCode() == 11) ) {

	  // grab associated clusters
	  const std::vector<art::Ptr<recob::Cluster> > pfp_clu_v = pfp_clu_assn_v.at(pp);
	  // for each cluster, find associated hits
	  for (size_t c=0; c < pfp_clu_v.size(); c++) {
	    // grab key and find hits
	    const std::vector<art::Ptr<recob::Hit> > clu_hit_v = clu_hit_assn_v.at( pfp_clu_v.at(c).key() );
	    // and add their indices to the pfp map
	    for (size_t h=0; h < clu_hit_v.size(); h++)
	      _pfpmap[trkKey].push_back( clu_hit_v.at(h) );
	  }// for all clusters associated to PFP
	}// if is daughter
      }// for all daughters
    }// second pfparticle loop
  }// for all PFParticles
  // END : PFPARTICLE MAP SCAN TO FIND DELTA-RAYS

  // BEGIN : IDENTIFY COSMIC TRACK HITS
  for (size_t t=0; t < trk_h->size(); t++) {

    auto const& trk = trk_h->at(t);
    auto trkdist = SphereIntersection(trk);
    
    std::cout << "Track has " << trkdist.first << " intersections w/ vertex ROI. IP min is : " << trkdist.second << std::endl;
    std::cout << "Track start [x,z] -> " << trk.Vertex().X() << ", " << trk.Vertex().Z() << std::endl;
    std::cout << "Track end   [x,z] -> " << trk.End().X()    << ", " << trk.End().Z()    << std::endl;

    // if no intersections -> check IP
    if ( (trkdist.first == 0) && (trkdist.second < fIPmin) ) continue;

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


    // for all deltay-rays associated to track, if they exist
    if ( _pfpmap.find(t) != _pfpmap.end() ) {
      auto const& pfp_hit_ptr_v = _pfpmap[t];
      for (size_t h=0; h < pfp_hit_ptr_v.size(); h++) {
	auto hitPtr = pfp_hit_ptr_v.at(h);// hit_h->at( pfp_hit_idx_v[h] );
	// if the hit channel is in the SSNet hit map:
	if (_hitmap.find( hitPtr->Channel() ) != _hitmap.end() ){
	  auto const& hitidx_v = _hitmap[ hitPtr->Channel() ];
	  for (auto const& idx : hitidx_v) {
	    // compare hit information
	    if ( hit_h->at(idx).PeakTime() == hitPtr->PeakTime() )
	      _trkhits.push_back( idx ); // save idx of SSNet hit to be removed
	  }// for all hit indices associated to this channel
	}// if the hit channel is in the SSNet hit map
      }// for all hits associated to track
    }// if delta-rays associated to track
    
    
  }// for all tracks
  // END : IDENTIFY COSMIC TRACK HITS
  
  // finally, save hits not identified as track-like
  for (size_t idx=0; idx < hit_h->size(); idx++) {
    // has this index been flagged?
    if (std::find(_trkhits.begin(),_trkhits.end(),idx) == _trkhits.end() )
      Hit_v->emplace_back(hit_h->at(idx));
  }// for all track hit indices

  std::cout << "input hits  : " << hit_h->size() << std::endl;
  std::cout << "output hits : " << Hit_v->size() << std::endl;
  
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
