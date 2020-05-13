#include "constellationOpt.h"
#include "constellationLoader.h"
#include "toolsGeneral.h"
#include "toolsNavigation.h"
#include "toolsSimulation.h"
#include <iomanip>

using namespace tudat;
using namespace tudat::basic_mathematics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::unit_conversions;

constellationOptimiser::constellationOptimiser( std::vector< std::vector< double > > &bounds,
                                      const bool useParamaterX ) :
    problemBounds_( bounds ), useParamaterX_( useParamaterX )
{
    // Load Spice kernels.
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );
}


//! Descriptive name of the problem
std::string constellationOptimiser::get_name() const {
    return "Multi-constellation GDOP optimiser";
}

//! Get bounds
std::pair<std::vector<double>, std::vector<double> > constellationOptimiser::get_bounds() const {

    return { problemBounds_[0], problemBounds_[1] };
}

//! Implementation of the fitness function (return delta-v)
std::vector<double> constellationOptimiser::fitness( const std::vector<double> &xv ) const{

    std::vector<double> f;

    /// CONSTANTS
    // Set Environment
    // Set simulation start epoch.
    const double simulationStartEpoch = 0.0;

    // Set simulation end epoch.
    const double simulationEndEpoch = 1.0 * physical_constants::JULIAN_DAY;

    // Set numerical integration fixed step size.
    const double fixedStepSize = 30.0;


    /// VARIABLES
    // // Set PECMEO constellation
    // Set number of PECMEO satellites in constellation.
    const unsigned int numberOfSatellitesPEC = 9;

    // Set number of PECMEO planes in constellation.
    const unsigned int numberOfPlanesPEC = 3;

    //  Set orbital parameters of PECMEO constellation.
    const double semiMajorAxis = xv[0] + 6378.1e3;                                // [km]

    const double eccentricity = 0.0;                                                  // [-]
    const double inclination = convertDegreesToRadians( 0.0 );                        // [rad]
    const double argumentOfPeriapsis = 0.0;                                           // [rad]

    const double longitudeOfAscendingNodeSpacing = convertDegreesToRadians( xv[1] );    // [rad]
    const double argumentOfLattitudeShift1 = convertDegreesToRadians( xv[2] );
    const double argumentOfLattitudeShift2 = convertDegreesToRadians( xv[3] );
    const double argumentOfLattitudeShift3 = convertDegreesToRadians( xv[4] );
    std::vector< double > keplerElements = { semiMajorAxis, eccentricity, inclination, argumentOfPeriapsis, longitudeOfAscendingNodeSpacing, argumentOfLattitudeShift1, argumentOfLattitudeShift2, argumentOfLattitudeShift3 };

    // // Set GNSS constellations (Galileo, GLONASS, BeiDou, GPS)
    // Comment out if you don't want to use constellation.
    std::vector <constellationNames> constellationList;
    constellationList.push_back(Galileo);
//    constellationList.push_back(GLONASS);
//    constellationList.push_back(BeiDou);
//    constellationList.push_back(GPS);

    // // Set ISL true or false
    bool interSatelliteTracking = false; // Set intersatellite tracking on or off


    /// OBTAINABLES
    // Define and set satellite name list and initial conditions for GNSS constellations.
    std::vector< std::string > satelliteNames;
    unsigned int numberOfSatellitesTotal;
    unsigned int numberOfConstellations;
    std::string includedNames;

    std::vector< std::map< double, Eigen::VectorXd > > allSatellitesPropagationHistory;
    std::vector< std::map< double, Eigen::VectorXd > > allSatellitesPropagationHistoryKep;
    std::map< std::string, std::map< double, Eigen::Vector6d > > allSatellitesDOPHistory;

    std::cout << std::setprecision(20) << "These are the parameter: " << xv[0] << ", " << xv[1] << ", " << xv[2] << ", " << xv[3] << ", " << xv[4] << std::endl;

    /// OBTAIN RESULTS
    initiateAndObtainResults( simulationStartEpoch, simulationEndEpoch, fixedStepSize, numberOfSatellitesPEC, numberOfPlanesPEC, keplerElements, constellationList, interSatelliteTracking,
                              numberOfSatellitesTotal, satelliteNames, numberOfConstellations, includedNames, allSatellitesPropagationHistory, allSatellitesPropagationHistoryKep, allSatellitesDOPHistory );

    std::vector< double > results;
    results = obtainParamaterToBeOptimised( allSatellitesDOPHistory );
    double averageLinks = results.at( 0 );
    double averageGDOP = results.at( 1 );

    f.push_back( averageGDOP );

    if( useParamaterX_ )
    {
        f.push_back( xv[1] );
    }
    return f;
}
