////////////////////////////////////////////////////////////////////////
// Class:       ProximityClustering
// Plugin Type: producer (art v2_05_01)
// File:        ProximityClustering_module.cc
//
// Generated at Fri Feb 16 10:12:16 2018 by David Caratelli using cetskelgen
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

// Data Products
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/Utilities/AssociationUtil.h"

// this line good for "new" larsoft:
// #include "art/Persistency/Common/PtrMaker.h"
// for outdated versions (i.e MCC8) use this line:
#include "lardata/Utilities/PtrMaker.h"

// Algorithms
#include "Algorithms/ProximityClusterer.h"

// ROOT
#include <TTree.h>

class ProximityClustering;


class ProximityClustering : public art::EDProducer {
public:
  explicit ProximityClustering(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ProximityClustering(ProximityClustering const &) = delete;
  ProximityClustering(ProximityClustering &&) = delete;
  ProximityClustering & operator = (ProximityClustering const &) = delete;
  ProximityClustering & operator = (ProximityClustering &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  TTree* _clus_tree;
  float  _charge;
  int    _nhits;

  // clustering radius and cell size 
  std::string fHitProducer;
  // [radius is going to be the physically important quantity
  // cell-size is only for algorithm implementation]
  // cell-size should be ~x2 of the radius
  double fCellSize;
  double fClusterRadius;

  // Proximity clusterer class
  gammacatcher::ProximityClusterer* _ProximityClusterer;

  void MakeClusters(const art::ValidHandle<std::vector<recob::Hit> >& hits,
		    const std::vector<std::vector<unsigned int> >& cluster_idx_v,
		    std::unique_ptr< std::vector<recob::Cluster> >& clusters,
		    std::unique_ptr<art::Assns<recob::Cluster, recob::Hit> >& assns,
		    lar::PtrMaker<recob::Cluster> ClusPtrMaker);

};


ProximityClustering::ProximityClustering(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{

  produces< std::vector< recob::Cluster > >();
  produces< art::Assns<  recob::Cluster, recob::Hit> >();
  
  // grab from fhicl file:
  fHitProducer      = p.get<std::string>("HitProducer"  );
  fCellSize         = p.get<double>     ("CellSize"     );
  fClusterRadius    = p.get<double>     ("ClusterRadius");

}

void ProximityClustering::produce(art::Event & e)
{
  // Implementation of required member function here.

  // geometry service
  art::ServiceHandle<geo::Geometry> geo;

  // produce Cluster objects
  std::unique_ptr< std::vector<recob::Cluster> > Cluster_v(new std::vector<recob::Cluster>);
  // produce associations
  auto Cluster_Hit_Assn_v = std::make_unique< art::Assns<recob::Cluster, recob::Hit> >();

  // load hits
  // load RawDigits
  auto const& hit_h = e.getValidHandle<std::vector<recob::Hit> >(fHitProducer);

  // Art Pointer maker
  //lar::PtrMaker<recob::Hit>     makeHitPtr (e, *this);
  lar::PtrMaker<recob::Cluster> makeClusPtr(e, *this);

  // cluster index vectors will be stored here
  std::vector<std::vector<unsigned int> > cluster_v;
  // cluster hits together
  _ProximityClusterer->cluster(hit_h,cluster_v);

  // go through indices and make clusters
  MakeClusters(hit_h, cluster_v, Cluster_v, Cluster_Hit_Assn_v, makeClusPtr);

  e.put(std::move(Cluster_v));
  e.put(std::move(Cluster_Hit_Assn_v));

}

void ProximityClustering::beginJob()
{
  // Implementation of optional member function here.

  // set TTree branches
  art::ServiceHandle<art::TFileService> tfs;
  _clus_tree = tfs->make<TTree>("_clus_tree","Cluster Info TTree");
  _clus_tree->Branch("_charge",&_charge,"charge/F");
  _clus_tree->Branch("_nhits" ,&_nhits ,"nhits/I" );

  _ProximityClusterer = new gammacatcher::ProximityClusterer();
  _ProximityClusterer->initialize();
  _ProximityClusterer->setRadius(fClusterRadius);
  _ProximityClusterer->setCellSize(fCellSize);

  return;
}

void ProximityClustering::endJob()
{
  // Implementation of optional member function here.
  return;
}

void ProximityClustering::MakeClusters(const art::ValidHandle<std::vector<recob::Hit> >& hit_h,
				       const std::vector<std::vector<unsigned int> >& cluster_v,
				       std::unique_ptr< std::vector<recob::Cluster> >& clusters,
				       std::unique_ptr< art::Assns<recob::Cluster, recob::Hit> >& assns,
				       lar::PtrMaker<recob::Cluster> ClusPtrMaker)
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
      auto const& hit = hit_h->at(hit_idx);
      float t = hit.PeakTime();
      float w = hit.WireID().Wire;
      if (t > t_max) t_max = t;
      if (t < t_min) t_min = t;
      if (w > w_max) w_max = w;
      if (w < w_min) w_min = w;
      integral += hit.Integral();
    }// for all hits
    recob::Cluster clus(w_min, 0., t_min, 0., 0., 0., 0., 
			w_max, 0., t_max, 0., 0., 0., 0., 
			integral, 0., integral, 0., 
			clus_idx_v.size(), 0., 0., n,
			hit_h->at(clus_idx_v[0]).View(),
			geo::PlaneID(0,0,hit_h->at(clus_idx_v[0]).WireID().Plane));
    
    _charge = integral;
    _nhits  = clus_idx_v.size();
    _clus_tree->Fill();

    clusters->emplace_back(clus);

    // fill associations
    art::Ptr<recob::Cluster> const ClusPtr = ClusPtrMaker(clusters->size()-1);
    for (auto const& hit_idx : clus_idx_v) {
      //art::Ptr<recob::Hit> HitPtr(_hitptr.id(),key,_hitptr.productGetter());
      //art::Ptr<recob::Hit> const HitPtr = HitPtrMaker(hit_idx);
      const art::Ptr<recob::Hit> HitPtr(hit_h, hit_idx);
      assns->addSingle(ClusPtr,HitPtr);
    }// for all hits
    
  }// for all clusters

  return;
}

DEFINE_ART_MODULE(ProximityClustering)
