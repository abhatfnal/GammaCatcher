/**
 * \file HitFinding.h
 * \brief Class def header for a class HitFinding 
 * @author david caratelli [davidc@fnal.gov]
 */

#ifndef GAMMACATCHER_HITFINDING_H
#define GAMMACATCHER_HITFINDING_H

#include <iostream>
#include <vector>
#include <algorithm>

namespace gammacatcher{
  
  /**
     Struct RawHit which defines a simple hit
  */
  struct RawHit {
    int   chan;   // wire-channel
    float ampl;   // amplitude of pulse
    float area;   // integrated charge
    float time;   // time of pulse
    float tstart; // start time
    float tend;   // end time
  };// RawHit struct

  
  class HitFinding{
    
  public:
    
    /// Default constructor
    HitFinding(){}
    
    /// Default destructor 
    ~HitFinding(){}
    
    /**
       Calculate baseline and RMS for a waveform's ADCs
    */
    std::pair<float,float> getBaselineRMS(const std::vector<short>& wf);

    /**
       Identify hits on a channel given baseline and rms
     */
    std::vector<gammacatcher::RawHit> getHits(const std::vector<short> &wf,
					      const double &baseline,
					      const double &rms);
    

    /**
       set number of sigmas on RMS required to be above threshold
     */
    void setNSigma(const double& nsigma) { _nsigma = nsigma; }

    /**
       set minimum hit width in ticks
     */
    void setMinTickWidth(const int& nticks) { _mintickwidth = nticks; }

  private:
    
    double _nsigma; // threshold on sigma RMS for being in active hit region.
    int    _mintickwidth; // minimum tick number above threshold for a hit

  };
  
}// namespace

#endif
/** @} */ // end of doxygen group   
