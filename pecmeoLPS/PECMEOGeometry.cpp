/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <Tudat/SimulationSetup/tudatEstimationHeader.h>

#include "applicationOutput.h"
#include "constellationLoader.h"
#include "toolsGeneral.h"
#include "toolsNavigation.h"
#include "toolsSimulation.h"
#include <iomanip>

//! Execute simulation of Galileo constellation around the Earth.
int main( )
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             USING STATEMENTS             //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#define FOLDINGSTART {
    using namespace tudat;
    using namespace tudat::simulation_setup;
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::basic_mathematics;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::unit_conversions;
    using namespace tudat::input_output;
    using namespace tudat::observation_models;
    using namespace tudat::ephemerides;
#define FOLDINGEND }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE ENVIRONMENT AND VEHICLES             ///////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// CONSTANTS
    // Set Environment
    // Set simulation start epoch.
    const double simulationStartEpoch = 0.0;

    // Set simulation end epoch.
    const double simulationEndEpoch = 1.0 * physical_constants::JULIAN_DAY;

    // Set numerical integration fixed step size.
    const double fixedStepSize = 30.0;

    // Set number of PECMEO satellites in constellation.
    unsigned int numberOfSatellitesPEC = 9; // 9

    // Set number of PECMEO planes in constellation.
    const unsigned int numberOfPlanesPEC = 3; // 3




    std::stringstream characteristics;
    bool lighttimeCorrection = false;
    bool addNoise = true;
    double stddevCode = 1.0;
    double stddevPhase= 0.001;
    double codeWeight = pow( 1.0, -2.0 ); //sigma^-2 = 1e6 (sigma = stddev, sigma^2 = variance)
    double phaseWeight = pow( 0.001, -2.0 );
    if( addNoise ){ characteristics << "_Noise_" << stddevCode << "," << stddevPhase; }



//    bool writeCartesianResult = true;
//    bool writeKeplerianResult = true;





    // Load Spice kernels.
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );


    /// VARIABLES
    // // Set PECMEO constellation
    //  Set orbital parameters of PECMEO constellation.
    const double semiMajorAxis = 14000.0e3; //+ 6378.1e3;                                // [km]
    const double eccentricity = 0.0;                                                  // [-]
    const double inclination = convertDegreesToRadians( 0.0 );                        // [rad]
    const double argumentOfPeriapsis = 0.0;                                           // [rad]
    const double longitudeOfAscendingNodeSpacing = convertDegreesToRadians( 0.0 );    // [rad]

    const double argumentOfLattitudeShift1 = convertDegreesToRadians( 0.0 );
    const double argumentOfLattitudeShift2 = convertDegreesToRadians( 90.0 );
    const double argumentOfLattitudeShift3 = convertDegreesToRadians( 0.0 );
    std::vector< double > keplerElements = { semiMajorAxis, eccentricity, inclination, argumentOfPeriapsis, longitudeOfAscendingNodeSpacing, argumentOfLattitudeShift1, argumentOfLattitudeShift2, argumentOfLattitudeShift3 };

//    for( bool setGalileo : { false, true } )
//    {
////        std::cout << setGalileo << std::endl;

////        if( setGalileo )
////        {
////            std::cout << "panama" << std::endl;
////        }

//        for( bool setGLONASS : { false, true } )
//        {
//            for( bool setBeiDou : { false, true } )
//            {
//                for( bool setGPS : { false, true } )
//                {
//                    for( bool interSatelliteTracking : { false, true } )
//                    {
//                        std::cout << "New Run" << std::endl;


    // // Set GNSS constellations (Galileo, GLONASS, BeiDou, GPS)
    // Comment out if you don't want to use constellation.
    bool setGalileo = true;
    bool setGLONASS = true;
    bool setBeiDou = true;
    bool setGPS = true;
    bool interSatelliteTracking = true;

                        std::vector <constellationNames> constellationList;
                        std::stringstream constellations;
                        if( setGalileo ){ constellationList.push_back(Galileo); constellations << "Galileo"; }
                        if( setGLONASS ){ constellationList.push_back(GLONASS); constellations << "GLONASS"; }
                        if( setBeiDou ){ constellationList.push_back(BeiDou); constellations << "BeiDou"; }
                        if( setGPS ){ constellationList.push_back(GPS); constellations << "GPS"; }

                        // // Set ISL true or false
                        if( interSatelliteTracking ){ constellations << "ISL"; }

                        std::string outputSubFolder = "PECMEO_Geometry/" + constellations.str();


                        /// OBTAINABLES
                        // Define and set satellite name list and initial conditions for GNSS constellations.
                        std::vector< std::string > satelliteNames;
                        unsigned int numberOfSatellitesTotal;
                        unsigned int numberOfConstellations;
                        std::string includedNames;
                        std::vector< double > epochTimes;
                        simulation_setup::NamedBodyMap bodyMap;

                        // Create empty maps
                        std::vector< std::map< double, Eigen::VectorXd > > allSatellitesPropagationHistory;
                        std::vector< std::map< double, Eigen::VectorXd > > allSatellitesPropagationHistoryKep;
                        std::map< std::string, std::map< double, Eigen::VectorXd > > allUsersDOPResults;
                        std::map< std::string, std::map< double, Eigen::Vector4d > > estimationsPPSequencePerUser;
                        std::map< std::string, std::map< double, Eigen::Vector4d > > estimationsKinSmartPerUser;
                        std::map< std::string, std::map< double, Eigen::Vector4d > > estimationsPPSequencePerUserISL;
                        std::map< std::string, std::map< double, Eigen::Vector4d > > estimationsKinSmartPerUserISL;
                        std::map< std::string, std::map< double, Eigen::VectorXd > > allUsersPPIterations;
                        std::map< std::string, std::map< double, Eigen::VectorXd > > allUsersKinIterations;
                        std::map< int, Eigen::VectorXd > biasIterations;

                        /// PERFORM PROPAGATION AND OBSERVATION
                        std::cout.precision(16);
                        PodInputDataType observationsAndTimes = initiateAndObtainResults( simulationStartEpoch, simulationEndEpoch, fixedStepSize, numberOfSatellitesPEC, numberOfPlanesPEC,
                                                                                          keplerElements, constellationList, interSatelliteTracking, numberOfSatellitesTotal, satelliteNames,
                                                                                          numberOfConstellations, includedNames, epochTimes, bodyMap, allSatellitesPropagationHistory,
                                                                                          allSatellitesPropagationHistoryKep);
                    //    writeResults( numberOfSatellitesTotal, satelliteNames, outputSubFolder, allSatellitesPropagationHistory, allSatellitesPropagationHistoryKep, writeCartesianResult, writeKeplerianResult);

                        /// OBTAIN RESULTS ESTIMATION
                        // For only looking at one satellite
                    //    numberOfSatellitesPEC = 1;
                        // For set custom times for estimation
//                        unsigned long long numberOfEpochs = 10; // 151, 637, 842, 1109, 1578, 2067, 2555, 2576, 2881, 151, 462, 200, 262, 428, 428, 427, 20, 301
//                        double timeJump{0}; // 0, 5250, 19260, 25410, 34500, 49170, 63840, 76680, 77400
//                        epochTimes.resize(numberOfEpochs);
//                        std::iota(std::begin(epochTimes), std::end(epochTimes), 0); //0 is the starting number
//                        std::transform(epochTimes.begin(), epochTimes.end(), epochTimes.begin(), [fixedStepSize](double& c){return c*fixedStepSize;});
//                        std::transform(epochTimes.begin(), epochTimes.end(), epochTimes.begin(), [timeJump](double& c){return c+timeJump;});


                        // Set estimation settings
                        int maxIterationsPP = 1;
                        int maxIterationsKin = 0;
                        int maxIterationsPPISL = 0;
                        int maxIterationsKinISL = 0;
                    //    double correctionLimitPP = 0.001; // 1 random standard value
                    //    double correctionLimitKinematic = 1e-8;

                        std::vector< double > results;
                        retrieveOptimisationParameters( observationsAndTimes, epochTimes, bodyMap, satelliteNames, numberOfSatellitesPEC, allUsersDOPResults, results );
//                        writeEstimationResults( numberOfSatellitesPEC, numberOfSatellitesTotal, satelliteNames, includedNames, outputSubFolder, lighttimeCorrection, addNoise, maxIterationsPP, maxIterationsKin,
//                                                maxIterationsPPISL, maxIterationsKinISL, allUsersDOPResults, estimationsPPSequencePerUser, estimationsKinSmartPerUser, estimationsPPSequencePerUserISL,
//                                                estimationsKinSmartPerUserISL, allUsersPPIterations, allUsersKinIterations, biasIterations, true, false, false, false, false, false );

                        std::cout << "\nThis is a PECMEO Geometry run! With " << numberOfSatellitesTotal << " satellites in " << numberOfConstellations + 1 << " constellations!"<< std::endl;
//                    }
//                }
//            }
//        }
//    }

//    for( int SPP = 1; SPP < 5; SPP++ )
//    {
//        for( int Kin = 0; Kin < 2; Kin++ )
//        {
//            for( int SPPISL = 0; SPPISL < 2; SPPISL++ )
//            {
//                for( int KinISL = 0; KinISL < 2; KinISL++ )
//                {
//                    std::cout << "Deze combi: " << SPP << Kin << SPPISL << KinISL << std::endl;

//                    int maxIterationsPP = SPP;
//                    int maxIterationsKin = Kin;
//                    int maxIterationsPPISL = SPPISL;
//                    int maxIterationsKinISL = KinISL;

//                    performNavigationEstimation(observationsAndTimes, epochTimes, bodyMap, satelliteNames, interSatelliteTracking, numberOfSatellitesPEC, numberOfSatellitesTotal,
//                                                maxIterationsPP, maxIterationsKin, correctionLimitPP, correctionLimitKinematic, codeWeight, phaseWeight, lighttimeCorrection, addNoise, stddevCode, stddevPhase,
//                                                maxIterationsPPISL, maxIterationsKinISL, allUsersDOPResults, estimationsPPSequencePerUser, estimationsKinSmartPerUser,
//                                                estimationsPPSequencePerUserISL, estimationsKinSmartPerUserISL, allUsersPPIterations, allUsersKinIterations, biasIterations);

//                    writeEstimationResults( numberOfSatellitesPEC, numberOfSatellitesTotal, satelliteNames, includedNames, outputSubFolder, lighttimeCorrection, addNoise, maxIterationsPP, maxIterationsKin,
//                                            maxIterationsPPISL, maxIterationsKinISL, allUsersDOPResults, estimationsPPSequencePerUser, estimationsKinSmartPerUser, estimationsPPSequencePerUserISL,
//                                            estimationsKinSmartPerUserISL, allUsersPPIterations, allUsersKinIterations, biasIterations );
//                }
//            }
//        }
//    }

//    std::cout << "\nThis is a PECMEO Geometry run! With " << numberOfSatellitesTotal << " satellites in " << numberOfConstellations + 1 << " constellations!"<< std::endl;

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
