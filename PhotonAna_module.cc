////////////////////////////////////////////////////////////////////////
// Class:       PhotonAna
// Plugin Type: analyzer (art v2_05_01)
// File:        PhotonAna_module.cc
//
// Generated at Thu Feb  1 21:04:28 2018 by David Caratelli using cetskelgen
// contact : David Caratelli [davidc@fnal.gov]
// from cetlib version v1_21_00.
// The goal of this module is to take identified isolated clusters
// measure their spatial isolation and study their energy spectra
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

// Services
#include "art/Framework/Services/Optional/TFileService.h"

#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/GeometryCore.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/Utilities/AssociationUtil.h"

// ROOT
#include <TFile.h>
#include <TTree.h>

class PhotonAna;


class PhotonAna : public art::EDAnalyzer {
public:
  explicit PhotonAna(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PhotonAna(PhotonAna const &) = delete;
  PhotonAna(PhotonAna &&) = delete;
  PhotonAna & operator = (PhotonAna const &) = delete;
  PhotonAna & operator = (PhotonAna &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  double _wire2cm, _time2cm;

  // buffer size
  double fBuff;
  // cluster size [cm2]
  double fClusAreaMax;
  // producers
  std::string fClusterProducer, fHitProducer, fVertexProducer, fAssnProducer;
  // debug mode
  bool fDebug;
  // use vertex?
  bool fUseVertex;

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
  void LoadVertex(const int& run, const int& event);

};


PhotonAna::PhotonAna(fhicl::ParameterSet const & p)
  : EDAnalyzer(p)
{
  fBuff            = p.get<double>("Buffer"              );
  fClusAreaMax     = p.get<double>("ClusAreaMax"         );
  fClusterProducer = p.get<std::string>("ClusterProducer");
  fHitProducer     = p.get<std::string>("HitProducer"    );
  fDebug           = p.get<bool>("Debug"                 );
  fUseVertex       = p.get<bool>("UseVertex"             );
  fTTreeName       = p.get<std::string>("TTreeName"      );
  fTFileName       = p.get<std::string>("TFileName"      );
  fDirName         = p.get<std::string>("DirName"        );

}

void PhotonAna::analyze(art::Event const & e)
{

  art::Handle<std::vector<recob::Cluster> > cluster_h;
  e.getByLabel(fClusterProducer,cluster_h);

  art::FindMany<recob::Hit> clus_hit_assn_v(cluster_h, e, fHitProducer);

  // if using vertex the vertex will be loaded by the below function:
  LoadVertex(e.run(),e.event());

  // COM variables
  double COMw, COMt;

  // double loop through clusters
  for (size_t n1=0; n1 < cluster_h->size(); n1++) {

    auto const& c1 = cluster_h->at(n1);

    if (fDebug)
      std::cout << "Cluster  w/ bounds [ "
		<< c1.StartTick() << " , " << c1.StartWire()
		<< " ] -> [ "
		<< c1.EndTick() << " , " << c1.EndWire() << " ]"
		<< std::endl;
	

    _charge = c1.Integral();

    // if cluster is small enough
    if (ClusterArea(c1) > fClusAreaMax)
      continue;

    ClusterCOM(clus_hit_assn_v.at(n1),COMw,COMt);

    if (fDebug)
      std::cout << " w/ area " << ClusterArea(c1) 
		<< " & COM : [" << COMt << " , " << COMw << " ]"
		<< std::endl;

    // sum of charge in buffer region associated
    // to nearby clusters
    _qbuffer = 0.;

    int Noverlap = 0;
    int Nbuffer  = 0;
    
    for (size_t n2=0; n2 < cluster_h->size(); n2++) {

      if (n1 == n2) continue;

      auto const& c2 = cluster_h->at(n2);

      // do the clusters overlap? 
      if (ClusterBoxOverlap(c1,c2) == false)
	continue;

      Noverlap += 1;
      
      double d = ClusterDistance(COMw,COMt, clus_hit_assn_v.at(n2));
      
      if (d < fBuff){
	Nbuffer += 1;
	_qbuffer += c2.Integral();
      }

    }// 2nd cluster loop

    // vertex distance, if using vertex
    if (fUseVertex)
      _vtxdist = sqrt( (_wirepos - COMw) * (_wirepos - COMw) + (_tickpos - COMt) * (_tickpos - COMt) );
    
    if (fDebug)
      std::cout << "Noverlap : " << Noverlap 
		<< "\t Nbuffer : " << Nbuffer
		<< "\t QBuffer : " << _qbuffer
		<< std::endl << std::endl;

    _clus_tree->Fill();
      
  }// 1st cluster loop

}

void PhotonAna::beginJob()
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

void PhotonAna::endJob()
{
  // Implementation of optional member function here.
}

bool PhotonAna::ClusterBoxOverlap(const recob::Cluster& c1,
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

void PhotonAna::ClusterCOM(const std::vector<const recob::Hit*>& hit_v,
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
  
  COMw /= qtot;
  COMt /= qtot;

  return;
}

double PhotonAna::ClusterArea(const recob::Cluster& c) 
{

  double dw = (c.EndWire()-c.StartWire());
  double dt = (c.EndTick()-c.StartTick());
  return fabs(dw * dt * _wire2cm * _time2cm);
}

double PhotonAna::ClusterDistance(const double& CMw,
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

void PhotonAna::LoadVertex(const int& run, const int& evt)
{
  
  // if we are using a vertex, find the position in this event:
  if (fUseVertex) {
    
    for (int n=0; n < _evt_tree->GetEntries(); n++) {
      _evt_tree->GetEntry(n);
      if ( (_run == run) && (_evt == evt) )
	break; // exit loop with correct xyz pos information for this event
    }// for all TTree entries
    // calculate wire-tick position of vertex
    // get 3D point to 2D:
    _wirepos = _zpos;
    _tickpos = _xpos;
    /*
    auto const* geoM = ::lar::providerFrom<geo::Geometry>();
    Double_t origin[3] = {};
    geoM->PlaneOriginVtx(2,origin);
    Double_t xyz[3] = {_xpos,_ypos,_zpos};
    auto const* geoH = ::lar::providerFrom<util::GeometryUtilities>();
    auto vtx2D = _geoH->Get2DPointProjection(xyz,2);
    _wirepos = vtx2D.w / _wire2cm;
    _tickpos = vtx2D.t + 800 * _time2cm - origin[0]; 
    */
  }// if using a vertex

  return;
}

DEFINE_ART_MODULE(PhotonAna)
