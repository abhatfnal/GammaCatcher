#ifndef GAMMACATCHER_PROXIMITYCLUSTERER_CXX
#define GAMMACATCHER_PROXIMITYCLUSTERER_CXX

#include "ProximityClusterer.h"

namespace gammacatcher {

  bool ProximityClusterer::initialize() {

    // get detector specific properties
    auto const* geom = ::lar::providerFrom<geo::Geometry>();
    auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    _wire2cm = geom->WirePitch(0,1,0);
    _time2cm = detp->SamplingRate() / 1000.0 * detp->DriftVelocity( detp->Efield(), detp->Temperature() );
    
    return true;
  }
  
  bool ProximityClusterer::cluster(const art::ValidHandle<std::vector<recob::Hit> >& hit_h,
				   std::vector<std::vector<unsigned int> >& _out_cluster_vector) {

    if (hit_h->size() == 0)
      return false;


    // a map to connect hit index wih a cluster index
    // each hit gets a cluster index
    // _clusterMap[hit_index] -> cluster_index
    std::map<size_t, size_t> _clusterMap;
    // a map to connect the cluster index with the vector of hit indices for that cluster
    // _clusters[index] -> vector of hit indices for that cluster
    std::map<size_t,std::vector<size_t> > _clusters;

    // keep track of largest cluster ID created
    size_t maxClusterID = 0;

    for (int pl=0; pl < 3; pl++){

      // hit map will only contain hits we want to use for clustering
      MakeHitMap(hit_h,pl);
      
      // iterator for hit cell map
      std::map<std::pair<int,int>, std::vector<size_t> >::iterator it;
      
      // loop through hits in each cell to find matches
      for (it = _hitMap.begin(); it != _hitMap.end(); it++){

	// pair = (i,j) indices of this cell in the _hitMap
	auto const& pair = it->first;
	
	// wire-space cell index
	// prepare a hit list of all neighboring cells
	// _________
	// |__|__|__|
	// |__|__|__|
	// |__|__|__|
	std::vector<size_t> cellhits = it->second;

	std::vector<size_t> neighborhits;
	getNeighboringHits(pair,neighborhits);

	for (size_t h1=0; h1 < cellhits.size(); h1++){

	  // has this hit been added to a cluster?
	  // if so not necessary to look at
	  auto const& hit1 = cellhits[h1];
	  // keep track if the hit will ever be matched to another
	  bool matched = false;
	  // if not find hits it should be clustered with and add it to the appropriate cluster
	  for (size_t h2=0; h2 < neighborhits.size(); h2++){
	    auto const& hit2 = neighborhits[h2];
	    if (hit1 == hit2) continue;
	    // are the hits compatible?
	    bool compat = HitsCompatible(hit_h->at(hit1),
					 hit_h->at(hit2));
	    // should the hits go in the same cluster?
	    if (compat){
	      matched = true;
	      // if both hits have already been assigned to a cluster then we can merge the cluster indices!
	      if ( (_clusterMap.find(hit1) != _clusterMap.end()) and
		   (_clusterMap.find(hit2) != _clusterMap.end()) ){
		// if in the same cluster -> do nothing
		// if they are in different clusters:
		if (_clusterMap[hit1] != _clusterMap[hit2]){
		  auto idx1 = _clusterMap[hit1];
		  auto idx2 = _clusterMap[hit2];
		  // hit indices for 1st cluster:
		  auto hits1 = _clusters[idx1];
		  auto hits2 = _clusters[idx2];
		  // append hits2 to hits1
		  for (auto h : hits2){
		    hits1.push_back(h);
		    // also change the index that the hit goes to (idx1 instead of idx2)
		    _clusterMap[h] = idx1;
		  }
		  _clusters[idx1] = hits1;
		  // erase cluster @ index2
		  _clusters.erase(idx2);
		}// if they are in different clusters
	      }
	      // if compatible and the 2nd hit has been added to a cluster
	      // add hit1 to the same cluster
	      else if ( (_clusterMap.find(hit2) != _clusterMap.end()) and
			(_clusterMap.find(hit1) == _clusterMap.end()) ){
		auto clusIdx = _clusterMap[hit2];
		_clusterMap[hit1] = clusIdx;
		_clusters[clusIdx].push_back(hit1);
	      }
	      // otherwise, add both to a new cluster
	      else if ( (_clusterMap.find(hit1) != _clusterMap.end()) and
			(_clusterMap.find(hit2) == _clusterMap.end()) ){
		auto clusIdx = _clusterMap[hit1];
		_clusterMap[hit2] = clusIdx;
		_clusters[clusIdx].push_back(hit2);
	      }
	      // if neither has a cluster yet
	      else{
		// create a new cluster for this match
		_clusterMap[hit1] = maxClusterID;
		_clusterMap[hit2] = maxClusterID;
		std::vector<size_t> cl = {hit1,hit2};
		_clusters[maxClusterID] = cl;
		maxClusterID += 1;
	      }
	    }// if the two hits are compatible
	  }// 2nd loop through hits in the cell
	  // has this hit been matched? if not we still need to add it as its own cluster
	  if (matched == false){
	    _clusterMap[hit1] = maxClusterID;
	    _clusters[maxClusterID] = {hit1};
	    maxClusterID += 1;
	  }
	}// 1st loop through hits in the cell
      }// loop through all cells

    }// loop through all planes

    // make a vector for the clusters
    for (auto it = _clusters.begin(); it != _clusters.end(); it++){
      auto indices = it->second;
      // if there are enough indices, make a cluster
      if (indices.size() >= 1){
	std::vector<unsigned int> clus;
	for (auto idx : indices)
	  clus.push_back(idx);
	_out_cluster_vector.push_back(clus);
      }// if there are 2 hits in cluster
    }
    
    return true;
  }

  // get all hits from neighboring cells
  void ProximityClusterer::getNeighboringHits(const std::pair<int,int>& pair, std::vector<size_t>& hitIndices){
   
    auto const& i       = pair.first;
    // time-space cell index
    auto const& j       = pair.second;

    // _________
    // |__|__|__|
    // |__|XX|__|
    // |__|__|__|
    if (_hitMap.find(std::make_pair(i,j)) != _hitMap.end()){
      for (auto &h : _hitMap[std::make_pair(i,j)])
	hitIndices.push_back(h);
    }

    // now look at neighboring cells, if they exist
    // _________
    // |__|__|__|
    // |XX|__|__|
    // |__|__|__|
    if (_hitMap.find(std::make_pair(i-1,j)) != _hitMap.end()){
      for (auto &h : _hitMap[std::make_pair(i-1,j)])
	hitIndices.push_back(h);
    }
    // _________
    // |__|__|__|
    // |__|__|__|
    // |__|XX|__|
    if (_hitMap.find(std::make_pair(i,j-1)) != _hitMap.end()){
      for (auto &h : _hitMap[std::make_pair(i,j-1)])
	hitIndices.push_back(h);
    }
    // _________
    // |__|__|__|
    // |__|__|__|
    // |XX|__|__|
    if ( _hitMap.find(std::make_pair(i-1,j-1)) != _hitMap.end() ){
      for (auto &h : _hitMap[std::make_pair(i-1,j-1)])
	hitIndices.push_back(h);
    }
    // _________
    // |__|XX|__|
    // |__|__|__|
    // |__|__|__|
    if ( _hitMap.find(std::make_pair(i,j+1)) != _hitMap.end() ){
      for (auto &h : _hitMap[std::make_pair(i,j+1)])
	hitIndices.push_back(h);
    }
    // _________
    // |__|__|__|
    // |__|__|XX|
    // |__|__|__|
    if ( _hitMap.find(std::make_pair(i+1,j)) != _hitMap.end() ){
      for (auto &h : _hitMap[std::make_pair(i+1,j)])
	hitIndices.push_back(h);
    }
    // _________
    // |__|__|XX|
    // |__|__|__|
    // |__|__|__|
    if ( _hitMap.find(std::make_pair(i+1,j+1)) != _hitMap.end() ){
      for (auto &h : _hitMap[std::make_pair(i+1,j+1)])
	hitIndices.push_back(h);
    }
    // _________
    // |XX|__|__|
    // |__|__|__|
    // |__|__|__|
    if ( _hitMap.find(std::make_pair(i-1,j+1)) != _hitMap.end() ){
      for (auto &h : _hitMap[std::make_pair(i-1,j+1)])
	hitIndices.push_back(h);
    }
    // _________
    // |__|__|__|
    // |__|__|__|
    // |__|__|XX|
    if ( _hitMap.find(std::make_pair(i+1,j-1)) != _hitMap.end() ){
      for (auto &h : _hitMap[std::make_pair(i+1,j-1)])
	hitIndices.push_back(h);
    }
  }

  // if two hits are further apart then the set distance -> not compatible
  bool ProximityClusterer::HitsCompatible(const recob::Hit& h1, const recob::Hit& h2){

    if (h1.WireID().Plane != h2.WireID().Plane)
      return false;

    double dt = ( h1.PeakTime() - h2.PeakTime() ) * _time2cm;
    //  if the hit time-ranges overlap, this distnce should be 0
    if (TimeOverlap(h1,h2,dt) == true)
      dt = 0;
    double dw = fabs(((double)h1.Channel()-(double)h2.Channel())*_wire2cm);
    if (dw >  0.3) dw -= 0.3;
    //if (dw < -0.3) dw  = 0.3;
    double d = dt*dt + dw*dw;

    if (d > (_radius*_radius))
      return false;

    return true;
  }

  bool ProximityClusterer::TimeOverlap(const recob::Hit& h1,
				       const recob::Hit& h2,
				       double& dmin) const
  {
    
    auto T1 = h1.PeakTime() * _time2cm; // time of first hit
    auto T2 = h2.PeakTime() * _time2cm;
    auto W1 = h1.RMS() * _time2cm;
    auto W2 = h2.RMS() * _time2cm;
    
    double d = dmin;
    
    if (T1 > T2) {
      
      if ( (T2+W2) > (T1-W1) ) return true;
      
      d = (T1-W1) - (T2+W2);
      if (d < dmin) dmin = d;
      
    }
    
    else {
      
      if ( (T1+W1) > (T2-W2) ) return true;
      
      d = (T2-W2) - (T1+W1);
      if (d < dmin) dmin = d;
      
    }

    return false;
  }
  
  void ProximityClusterer::MakeHitMap(const art::ValidHandle<std::vector<recob::Hit> >& hitlist, int plane){
    
    _hitMap.clear();
    // temporary pair
    std::pair<int,int> tmpPair;

    
    for (size_t h=0; h < hitlist->size(); h++){
      
      auto const& hit = hitlist->at(h);
      // skip if not of plane we want
      if (hit.View() != plane)
	continue;

      double t = hit.PeakTime()*_time2cm;
      double w = hit.WireID().Wire*_wire2cm;

      // if vertex is loaded, make sure within ROI
      if (_vertex) {
	// get 2D distance to vertex
	double d2D = ( (t-_vtx_t_cm[plane]) * (t-_vtx_t_cm[plane]) +
		       (w-_vtx_w_cm[plane]) * (w-_vtx_w_cm[plane]) );
	// ignore hit if out of ROI
	if (d2D > _ROISq) continue;
      }
	

      // map is (i,j) -> hit list
      // i : ith bin in wire of some width
      // j : jth bin in time of some width
      int i = int(w/_cellSize);
      int j = int(t/_cellSize);
      tmpPair = std::make_pair(i,j);
      // does this entry exist in the map?
      // if yes -> append to vector
      // if no create new vector and add to map
      if (_hitMap.find(tmpPair) == _hitMap.end()){
	std::vector<size_t> aaa = {h};
	_hitMap[tmpPair] = aaa;
      }
      else
	_hitMap[tmpPair].push_back(h);
    }// for all hits

    return;
  }

    bool ProximityClusterer::loadVertex(const art::ValidHandle<std::vector<::recob::Vertex> > vtx_h, const double& ROI) {

    _vtx_w_cm = {0., 0., 0.};
    _vtx_t_cm = {0., 0., 0.};

    _ROISq = ROI*ROI;

    // load required services to obtain offsets
    auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    auto const& geomH = ::util::GeometryUtilities::GetME();
    
    if (vtx_h->size() != 1) {
      _vertex = false;
      return false;
    }

    auto const& vtx = vtx_h->at(0);
    
    Double_t xyz[3] = {};
    vtx.XYZ(xyz);


    std::cout << "Vtx coordinates : [" << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << "]" << std::endl;

    for (size_t pl = 0; pl < 3; pl++) {

      auto const& pt = geomH->Get2DPointProjectionCM(xyz,pl);
      _vtx_w_cm[pl] = pt.w;
      _vtx_t_cm[pl] = pt.t + (detp->TriggerOffset() * _time2cm);

      std::cout << "trigger offset [cm] : " << (detp->TriggerOffset() * _time2cm) << std::endl;
      std::cout << "Vtx @ pl " << pl << " [" << _vtx_w_cm[pl] << ", " << _vtx_t_cm[pl] << " ]" << std::endl;

    }

    _vertex = true;

    return true;
  }


}
#endif
