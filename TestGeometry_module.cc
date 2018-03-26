////////////////////////////////////////////////////////////////////////
// Class:       TestGeometry
// Plugin Type: analyzer (art v2_05_01)
// File:        TestGeometry_module.cc
//
// Generated at Mon Mar 26 10:47:14 2018 by David Caratelli using cetskelgen
// from cetlib version v1_21_00.
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

#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/GeometryCore.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"


class TestGeometry;


class TestGeometry : public art::EDAnalyzer {
public:
  explicit TestGeometry(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TestGeometry(TestGeometry const &) = delete;
  TestGeometry(TestGeometry &&) = delete;
  TestGeometry & operator = (TestGeometry const &) = delete;
  TestGeometry & operator = (TestGeometry &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

};


TestGeometry::TestGeometry(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{

  auto const* geom = ::lar::providerFrom<geo::Geometry>();
  auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();

  auto wire2cm = geom->WirePitch(0,1,0);
  auto time2cm = detp->SamplingRate() / 1000.0 * detp->DriftVelocity( detp->Efield(), detp->Temperature() );
  auto toffset = detp->TriggerOffset() * time2cm;

  std::cout << "wire2cm : " << wire2cm << std::endl;
  std::cout << "time2cm : " << time2cm << std::endl;
  std::cout << "trigger offset [cm] : " << toffset << std::endl;
  
  auto wirecm = geom->WireCoordinate(0.,500.,geo::PlaneID(0,0,2)) * wire2cm;
  std::cout << "Y,Z : [0, 500] cm -> " << " Y-plane wire "  << wirecm << std::endl;

  auto w = 1000. * wire2cm;
  auto t = (1000. - detp->TriggerOffset()) * time2cm;
  std::cout << "[w,t] = [1000,1000] goes to cm : [" << w << ", " << t << "] " << std::endl;



}

void TestGeometry::analyze(art::Event const & e)
{

  return;
}

void TestGeometry::beginJob()
{
  // Implementation of optional member function here.
}

void TestGeometry::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(TestGeometry)
