/**
 * \file ProximityClusterer.h
 *
 * \ingroup Clusterer
 * 
 * \brief Class def header for a class ProximityClusterer
 *
 * @author david caratelli
 */

/** \addtogroup Clusterer

    @{*/

#ifndef GAMMACATCHER_PROXIMITYCLUSTERER_H
#define GAMMACATCHER_PROXIMITYCLUSTERER_H

#include <map>

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/GeometryCore.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

namespace gammacatcher {
  /**
     \class ProximityClusterer
     User custom analysis class made by david caratelli
   */
  class ProximityClusterer {
  
  public:

    /// Default constructor
    ProximityClusterer(){    
      _verbose     = false;
      _radius      = 2.0;
      _cellSize    = 2;
    }

    /// Default destructor
    virtual ~ProximityClusterer(){}

    /** IMPLEMENT in ProximityClusterer.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    bool initialize();

    /**
       Cluster function:
       @brief cluster hits based on proximity
       Input: pointer to event hit record.
       Output: vector of vector of hit indices which make up clusters
    */
    bool cluster(const art::ValidHandle<std::vector<recob::Hit> >& hit_h,
		 std::vector<std::vector<unsigned int> >& _out_cluster_vector);

    /**
       @brief Get total charge within a certain radius of a point
       Given the event hit record and the hits and COM for a specific cluster
       this function calculates the total charge depsoited within a 2D radius of the
       cluster on a specific plane, excluding any charge from hits already included in the cluster
       @input hit_h -> hit vector from the event record
       @input hit_v -> vector of hits associated to the cluster we are interested in
       @input plane  -> plane for which the COM info is given and for which to calculate the nearby hits
       @return total charge in the integration radius (in hit.Integral() units)
     */
    double nearbyCharge(const art::ValidHandle<std::vector<recob::Hit> >& hit_h,
			const std::vector<art::Ptr<recob::Hit> > hit_v,
			const int& plane, const double& radius);
			
    /**
       @brief Get nearest hit to a cluster, for hits in the plane not associated to the cluster itself.
       Given the event hit record and the hits and COM for a specific cluster
       this function calculates the distance to the nearest hit external to the cluster.
       @input hit_h -> hit vector from the event record
       @input hit_v -> vector of hits associated to the cluster we are interested in
       @input plane  -> plane for which the COM info is given and for which to calculate the nearby hits
       @return distance to nearest hit
     */
    double closestHit(const art::ValidHandle<std::vector<recob::Hit> >& hit_h,
		      const std::vector<art::Ptr<recob::Hit> > hit_v,
		      const int& plane);

			
    

    /// Set the size of each cell for hit-map
    void setCellSize(double d) { _cellSize = d; }
    /// Set the radius around which to search for hits
    /// if two hits are within this distance of each other
    /// then they go into the same cluster
    void setRadius(double d) { _radius = d; }
    /// Verbosity setter
    void setVerbose(bool on) { _verbose = on; }

    // vertex coordinates on each plane
    bool loadVertex(const art::ValidHandle<std::vector<::recob::Vertex> > vtx_h,
		    const double& ROI);
    
  protected:

    /// size of each cell [cm]
    double _cellSize;

    /// radius to count charge around [cm]
    double _radius;
    
    /// plane to select hits from
    int _plane;

    /// verbosity flag
    bool _verbose;

    /// conversion factors for hits
    double _wire2cm, _time2cm;

    /// Map making function
    void MakeHitMap(const art::ValidHandle<std::vector<recob::Hit> >& hit_h,
		    int plane);

    /// Functions to decide if two hits should belong to the same cluster or not
    bool HitsCompatible(const recob::Hit& h1, const recob::Hit& h2);

    /// Function to get neighboring hits with a variable cellspan (from self + neighoring cells)
    // if cellSpan = 0, get only from current cell
    // if == 1 -> get from 3x3 cell matrix, etc...
    void getNeighboringHits(const std::pair<int,int>& pair, const size_t& cellSpan,
			    std::vector<size_t>& hitIndices);

    /// Function to get neighboring hits (from self + neighoring cells)
    void getNeighboringHits(const std::pair<int,int>& pair, std::vector<size_t>& hitIndices);

    /// check if time overlaps
    bool TimeOverlap(const recob::Hit& h1, const recob::Hit& h2, double& dmin) const;
    
    /// map connecting coordinate index (i,j) to [h1,h2,h3] (hit index list)
    std::map<std::pair<int,int>, std::vector<size_t> > _hitMap;

    /// maximum i'th and j'th
    int _maxI;
    int _maxJ;

    // has the vertex been loaded?
    bool _vertex;
    // ROI squared distance max to vertex
    double _ROISq;

    /// vertex coordinates
    std::vector<double> _vtx_w_cm;
    std::vector<double> _vtx_t_cm;

  };
}
#endif

//**************************************************************************
// 
// For Analysis framework documentation, read Manual.pdf here:
//
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
//**************************************************************************

/** @} */ // end of doxygen group 
