////////////////////////////////////////////////////////////////////////
// Class:       NeutrinoGammas
// Plugin Type: producer (art v2_05_01)
// File:        NeutrinoGammas_module.cc
//
// Generated at Mon Feb  5 22:53:51 2018 by David Caratelli using cetskelgen
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

// Services
#include "art/Framework/Services/Optional/TFileService.h"

#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/GeometryCore.h"
// #include "larcorealg/GeoAlgo/GeoAlgo.h" CANNOT USE, CANNOT FIND IN DEPENDENCIES FOR MCC8...
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/Utilities/AssociationUtil.h"

// ROOT
#include <TFile.h>
#include <TTree.h>

class NeutrinoGammas;


class NeutrinoGammas : public art::EDProducer {
public:
  explicit NeutrinoGammas(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NeutrinoGammas(NeutrinoGammas const &) = delete;
  NeutrinoGammas(NeutrinoGammas &&) = delete;
  NeutrinoGammas & operator = (NeutrinoGammas const &) = delete;
  NeutrinoGammas & operator = (NeutrinoGammas &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  double _wire2cm, _time2cm;

  // buffer size
  double fBuff;
  // cluster size [cm2]
  double fClusAreaMax;
  // producers
  std::string fClusterProducer, fHitProducer, fVertexProducer, fAssnProducer, fTrackProducer;
  // debug mode
  bool fDebug;
  // veto radius -> surrounds nu vtx and is used to identify crossing cosmic-rays
  double fVetoRadius;
  // square of radius for faster computation
  double fVetoRadiusSq;
  // minimum impact paramter for track to be considered cosmic
  double fIPmin;

  // TTree variables
  TTree *_clus_tree;
  float _qbuffer; // charge in buffer region [ADC x TICK summed over all clusters]
  float _charge;  // cluster charge [ADC x TICK]
  float _vtxdist; // distance from reconstructed neutrino vertex, if being used

  // TTree for external data storing
  // event and vertex information
  TTree *_evt_tree;
  Int_t _run, _evt;
  Float_t _xpos, _ypos, _zpos; // xyz of vertex
  float _wirepos, _tickpos; // wire-tick coordinates of reco'd vertex on Y plane
  std::string fTTreeName, fTFileName, fDirName;

  /**
     Do cluster boxes overlap?
     Returns true of [wmin, tmin] -> [wmax, tmax]
     Of two clusters overlap
   */
  bool ClusterBoxOverlap(const recob::Cluster& c1,
			 const recob::Cluster& c2);

  /**
     Get cluster Center Of Mass
   */
  void ClusterCOM(const std::vector<const recob::Hit*>& hit_v,
		  double& COMw, double& COMt);

  /**
     Get cluster area in cm2
   */
  double ClusterArea(const recob::Cluster& c);

  /**
     2D distance between clusters:
     compare CM (w,t) of small cluster
     to all hits in larger cluster
   */
  double ClusterDistance(const double& CMw,
			 const double& CMt,
			 const std::vector<const recob::Hit*>& hit_v);

  /**
     Load vertex from TTree input information, if this is available
   */
  bool LoadVertex(const int& run, const int& event);

  /**
     Calculate number of tracks intersecting veto radius surrounding
     the neutrino interaction
   */
  int CosmicVeto(const  art::Handle<std::vector<recob::Track> >& tracks);

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


NeutrinoGammas::NeutrinoGammas(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{

  produces< std::vector< recob::Cluster > >();

  fBuff            = p.get<double>("Buffer"              );
  fClusAreaMax     = p.get<double>("ClusAreaMax"         );
  fClusterProducer = p.get<std::string>("ClusterProducer");
  fHitProducer     = p.get<std::string>("HitProducer"    );
  fTrackProducer   = p.get<std::string>("TrackProducer"  );
  fDebug           = p.get<bool>("Debug"                 );
  fTTreeName       = p.get<std::string>("TTreeName"      );
  fTFileName       = p.get<std::string>("TFileName"      );
  fDirName         = p.get<std::string>("DirName"        );
  fVetoRadius      = p.get<double     >("VetoRadius"     );
  fIPmin           = p.get<double     >("IPmin"          );
  fVetoRadiusSq    = fVetoRadius * fVetoRadius;

}

void NeutrinoGammas::produce(art::Event & e)
{

  // load clusters and hits already reconstructed by previous step
  art::Handle<std::vector<recob::Cluster> > cluster_h;
  e.getByLabel(fClusterProducer,cluster_h);
  art::FindMany<recob::Hit> clus_hit_assn_v(cluster_h, e, fHitProducer);

  // load tracks to do cosmic-removal
  art::Handle<std::vector<recob::Track> > track_h;
  e.getByLabel(fTrackProducer,track_h);

  // produce clusters
  std::unique_ptr< std::vector<recob::Cluster> > Cluster_v(new std::vector<recob::Cluster>);
    
  // if using vertex the vertex will be loaded by the below function:
  auto foundvtx = LoadVertex(e.run(),e.event());
  if (foundvtx == true) {

    if (fDebug) { std::cout << std::endl << "VERTEX @ " << _xpos << ", " << _ypos << ", " << _zpos << std::endl << std::endl; }

    // Identify how many tracks enter vertex buffer veto region
    auto const& ncosmics = CosmicVeto(track_h);
    
    std::cout << "Number of cosmics intersecting vertex veto : " << ncosmics << std::endl;
  }
  else
    std::cout << "No VERTEX FOUND!" << std::endl;

  e.put(std::move(Cluster_v));

}

void NeutrinoGammas::beginJob()
{

  // get detector specific properties
  auto const* geom = ::lar::providerFrom<geo::Geometry>();
  auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  _wire2cm = geom->WirePitch(0,1,0);
  _time2cm = detp->SamplingRate() / 1000.0 * detp->DriftVelocity( detp->Efield(), detp->Temperature() );

  art::ServiceHandle<art::TFileService> tfs;
  _clus_tree = tfs->make<TTree>("_clus_tree","Cluster Info TTree");
  _clus_tree->Branch("_charge" ,&_charge ,"charge/F" );
  _clus_tree->Branch("_qbuffer",&_qbuffer,"qbuffer/F");
  _clus_tree->Branch("_vtxdist",&_vtxdist,"vtxdist/F");

  // load TTree from disk where event/vertex info is stored
  TFile* f = new TFile(fTFileName.c_str());
  f->cd(fDirName.c_str());
  std::string treepath = fDirName+"/"+fTTreeName;
  _evt_tree = (TTree*)f->Get(treepath.c_str());
  _evt_tree->SetBranchAddress("_run",&_run);
  _evt_tree->SetBranchAddress("_evt",&_evt);
  _evt_tree->SetBranchAddress("_xpos",&_xpos);
  _evt_tree->SetBranchAddress("_ypos",&_ypos);
  _evt_tree->SetBranchAddress("_zpos",&_zpos);

}

void NeutrinoGammas::endJob()
{
  // Implementation of optional member function here.
}

bool NeutrinoGammas::ClusterBoxOverlap(const recob::Cluster& c1,
				  const recob::Cluster& c2)
{

  // assumption
  // end wire/tick always larger than start wire/tick

  // if not on the same plane -> no overlap!
  if (c1.Plane() != c2.Plane())
    return false;

  auto sw1 = c1.StartWire() * _wire2cm;
  auto ew1 = c1.EndWire() * _wire2cm;
  auto st1 = c1.StartTick() * _time2cm;
  auto et1 = c1.EndTick() * _time2cm;

  auto sw2 = c2.StartWire() * _wire2cm;
  auto ew2 = c2.EndWire() * _wire2cm;
  auto st2 = c2.StartTick() * _time2cm;
  auto et2 = c2.EndTick() * _time2cm;

  // if at least one point of one rectangle in the other
  // -> overlap
  if ( (sw1 > (sw2-fBuff)) && (sw1 < (ew2+fBuff)) && (st1 > (st2-fBuff)) && (st1 < (et2+fBuff)) )
    return true;
  if ( (sw1 > (sw2-fBuff)) && (sw1 < (ew2+fBuff)) && (et1 > (st2-fBuff)) && (et1 < (et2+fBuff)) )
    return true;
  if ( (ew1 > (sw2-fBuff)) && (ew1 < (ew2+fBuff)) && (st1 > (st2-fBuff)) && (st1 < (et2+fBuff)) )
    return true;
  if ( (ew1 > (sw2-fBuff)) && (ew1 < (ew2+fBuff)) && (et1 > (st2-fBuff)) && (et1 < (et2+fBuff)) )
    return true;

  if ( (sw2 > (sw1-fBuff)) && (sw2 < (ew1+fBuff)) && (st2 > (st1-fBuff)) && (st2 < (et1+fBuff)) )
    return true;
  if ( (sw2 > (sw1-fBuff)) && (sw2 < (ew1+fBuff)) && (et2 > (st1-fBuff)) && (et2 < (et1+fBuff)) )
    return true;
  if ( (ew2 > (sw1-fBuff)) && (ew2 < (ew1+fBuff)) && (st2 > (st1-fBuff)) && (st2 < (et1+fBuff)) )
    return true;
  if ( (ew2 > (sw1-fBuff)) && (ew2 < (ew1+fBuff)) && (et2 > (st1-fBuff)) && (et2 < (et1+fBuff)) )
    return true;

  // else -> no overlap
  return false;
}

void NeutrinoGammas::ClusterCOM(const std::vector<const recob::Hit*>& hit_v,
			   double& COMw, double& COMt) 
{

  COMw = 0.;
  COMt = 0.;
  double qtot = 0.;

  for (auto const& hit : hit_v){
    COMw += hit->WireID().Wire * _wire2cm * hit->Integral();
    COMt += hit->PeakTime()    * _time2cm * hit->Integral();
    qtot += hit->Integral();
  }
  
  COMw /= (hit_v.size() * qtot);
  COMt /= (hit_v.size() * qtot);

  return;
}

double NeutrinoGammas::ClusterArea(const recob::Cluster& c) 
{

  double dw = (c.EndWire()-c.StartWire());
  double dt = (c.EndTick()-c.StartTick());
  return fabs(dw * dt * _wire2cm * _time2cm);
}

double NeutrinoGammas::ClusterDistance(const double& CMw,
				  const double& CMt,
				  const std::vector<const recob::Hit*>& hit_v) 
{

  double dminSq = 99999.;
  
  for (auto const& hit : hit_v) {
    
    double dSq = ( pow(((hit->PeakTime()*_time2cm)-CMt),2) + 
		   pow(((hit->WireID().Wire*_wire2cm)-CMw),2) );
    
    if (dSq < dminSq) dminSq = dSq;
    
  }// for all hits

  return sqrt(dminSq);
}

bool NeutrinoGammas::LoadVertex(const int& run, const int& evt)
{

  bool found = false;
  
  for (int n=0; n < _evt_tree->GetEntries(); n++) {
    _evt_tree->GetEntry(n);
    if ( (_run == run) && (_evt == evt) ) {
      found = true;
      break; // exit loop with correct xyz pos information for this event
    }
  }// for all TTree entries
  // calculate wire-tick position of vertex
  // get 3D point to 2D:

  if (found == false) 
    return false;

  _wirepos = _zpos;
  _tickpos = _xpos;
  
  return true;
}

int NeutrinoGammas::CosmicVeto(const  art::Handle<std::vector<recob::Track> >& tracks) {

  int ntracks = 0;
  
  if (fDebug) 
    std::cout << "looping trhough " << tracks->size() << " cosmic tracks" << std::endl;

  // for each track in the event, determine if it intersects the veto region
  for (size_t t=0; t < tracks->size(); t++) {
    auto const& trk = tracks->at(t);
    auto SEdist = (trk.Vertex()-trk.End()).Mag(); // linear separation between vertex and end-point
    // measure if the track linear -> helps determine if really a cosmic muon or some nearby garbage
    if ( (trk.Length() > 50.) && (SEdist > 50) && (trk.Length() < 2* SEdist) ) {
      auto const& intersection = SphereIntersection(trk);
      if (fDebug) { std::cout << "track " << t << " has " << intersection.first 
			      << " intersections w/ Vertex Veto" << std::endl; }
      // if more then one intersection point -> cosmic track
      if (intersection.first > 1) {
	if (fDebug) { std::cout << "\t 2+ intersections -> cosmic!" << std::endl; }
	ntracks += 1;
      }
      else if ( (intersection.second > fIPmin) && (intersection.first == 1) ) {
	if (fDebug) { std::cout << "\t IP : " << intersection.second << " -> cosmic! " << std::endl; }
	ntracks += 1;
      }
    }// if track is straight and long
  }// for all tracks in the event
  
  return ntracks;
}

std::pair<int,float> NeutrinoGammas::SphereIntersection(const recob::Track& trk) {
  
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
double NeutrinoGammas::SqDist(const TVector3& pt) {

  return (_xpos - pt[0]) * (_xpos - pt[0]) + (_ypos - pt[1]) * (_ypos - pt[1]) + (_zpos - pt[2]) * (_zpos - pt[2]);
}

DEFINE_ART_MODULE(NeutrinoGammas)
