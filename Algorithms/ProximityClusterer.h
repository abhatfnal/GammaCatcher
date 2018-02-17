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
