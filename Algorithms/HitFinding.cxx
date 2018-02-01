#ifndef GAMMACATCHER_HITFINDING_CXX
#define GAMMACATCHER_HITFINDING_CXX

#include "HitFinding.h"

namespace gammacatcher {

  std::pair<float,float> HitFinding::getBaselineRMS(const std::vector<short>& wf) {
    
    // to do:
    // scan only portion of waveform
    // define acceptable noise level
    // if noise below -> return RMS and baseline
    // otherwise move to next segment in waveform.
    
    auto adc_v = wf;
    std::sort(adc_v.begin(), adc_v.end());
    
    // truncate top 10% of values [to remove real pulses]
    std::vector<short> truncated_adc(adc_v.begin(), adc_v.end() - adc_v.size()/10);
    
    float base  = 0.;
    float rms   = 0.;
    
    for (auto const& adc : truncated_adc)
      base += adc;
    base = base / truncated_adc.size();
    for (auto const& adc : truncated_adc)
      rms += (adc-base)*(adc-base);
    rms = sqrt( rms / ( truncated_adc.size() - 1) );
    
    return std::make_pair(base,rms);
    
  }

  std::vector<gammacatcher::RawHit> HitFinding::getHits(const std::vector<short> &wf,
							const double &baseline,
							const double &rms){

    // prepare a vector of hits to populate                        
    std::vector<gammacatcher::RawHit> hits;
    
    // keep track of:         
    double amp;   // amplitude of hit                                             
    double area;  // area of hit
    int    peak;  // peak-time of hit
    int    start; // start time of hit
    int    end;   // end time of hit
    // are we in an active hit?                                                                                                      
    bool active = false;

    // loop through waveform                                                                                                         
    for (size_t i=0; i < wf.size(); i++){

      double h = wf[i]-baseline;

      // is it above the cut we want to apply?                                                                                       
      if (h > rms*_nsigma){
        // if not active start the hit                                                                                               
	if (!active){
          amp = h;
          area = h;
          start = i;
          peak = i;
          active = true;
        }
        // if already in an active region keep adding to the hit                                                                     
        else{
          area += h;
          if (h > amp){
            amp = h;
            peak = i;
          }
        }
      }
      // else -> we are not in an active region                                                                                      
      else{
        // if we were in an active region we have reached the end of a hit                                                           
        if (active){
          end = i;
          // now make hit                                                                                                            
          // if passes minimum hit width condition                                                                                   
          if ( (end-start) >= _mintickwidth){
	    gammacatcher::RawHit hit;
	    hit.ampl   = amp;
	    hit.area   = area;
	    hit.time   = peak;
	    hit.tstart = start;
	    hit.tend   = end;
            hits.push_back(hit);
          }
          active = false;
          // clear hit attributes                                                                                                    
          peak = 0;
          area = 0;
          end = 0;
          start = 0;
          amp = 0;
        }
      }// if not in an active region                                                                                                 
    }// looping over waveform                                                                                                        

    return hits;
  }


}// namespace

#endif
