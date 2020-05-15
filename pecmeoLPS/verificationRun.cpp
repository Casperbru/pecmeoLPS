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
#include <memory>

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

    /// ENVIRONMENT/MODEL SETTINGS
/*
    // Set simulation start epoch.
    const double simulationStartEpoch = 0.0;

    // Set simulation end epoch.
    std::ostringstream timeSpan;
//    double days = 1;  const double simulationEndEpoch = days * physical_constants::JULIAN_DAY; timeSpan << "/" << days << "d/";
    double hours = 2; const double simulationEndEpoch = 60.0*60.0*hours; timeSpan << "/" << hours << "h/";
//    double minutes = 20; const double simulationEndEpoch = 60.0*minutes; timeSpan << "/" << minutes << "min/";

    // Set numerical integration fixed step size.
    const double fixedStepSize = 30.0;
    std::ostringstream timeStep; timeStep << "_" << fixedStepSize << "sec";

    // // PECMEO constellation settings
    #define FOLDINGSTART {
    // Set number of PECMEO satellites in constellation.
    unsigned int numberOfSatellitesPEC = 9; // 9

    // Set number of PECMEO planes in constellation.
    const unsigned int numberOfPlanesPEC = 3; // 3

    //  Set orbital parameters of PECMEO constellation.
    const double semiMajorAxis = 14000.0e3;                                           // [km]  + 6378.1e3
    const double eccentricity = 0.0;                                                  // [-]
    const double inclination = convertDegreesToRadians( 0.0 );                        // [rad]
    const double argumentOfPeriapsis = 0.0;                                           // [rad]
    const double longitudeOfAscendingNodeSpacing = convertDegreesToRadians( 0.0 );    // [rad]
    const double argumentOfLattitudeShift1 = convertDegreesToRadians( 0.0 );
//    const double argumentOfLattitudeShift2 = convertDegreesToRadians( 90.0 );
    const double argumentOfLattitudeShift2 = convertDegreesToRadians( 5.0 );
    const double argumentOfLattitudeShift3 = convertDegreesToRadians( 0.0 );
    std::vector< double > keplerElements = { semiMajorAxis, eccentricity, inclination, argumentOfPeriapsis, longitudeOfAscendingNodeSpacing, argumentOfLattitudeShift1, argumentOfLattitudeShift2, argumentOfLattitudeShift3 };
    #define FOLDINGEND }


    // // Set GNSS constellations (Galileo, GLONASS, BeiDou, GPS) and ISL
    bool setGalileo = true;
    bool setGLONASS = true;
    bool setBeiDou = true;
    bool setGPS = true;
    std::vector <constellationNames> constellationList;
    std::stringstream constellations; constellations << "/";
    if( setGalileo ){ constellationList.push_back(Galileo); constellations << "Galileo"; }
    if( setGLONASS ){ constellationList.push_back(GLONASS); constellations << "GLONASS"; }
    if( setBeiDou ){ constellationList.push_back(BeiDou); constellations << "BeiDou"; }
    if( setGPS ){ constellationList.push_back(GPS); constellations << "GPS"; }


    // // ISL Settings
    bool interSatelliteTracking = true;
    bool laserISL = false;
    bool KbandISL = false;
    bool clockPhase = false;
    std::stringstream ISLMap; ISLMap << "/";
    if( !interSatelliteTracking ){ ISLMap << "NoISL"; }
    double stddevCodeISL = 1.0e-0;  double codeWeightISL = pow( stddevCodeISL, -2.0 ); // sigma^-2 = weight (sigma = stddev, sigma^2 = variance)
    double stddevPhaseISL = 1.0e-3; double phaseWeightISL = pow( stddevPhaseISL, -2.0 );
//    codeWeightISL = 1; phaseWeightISL = 1;  // Setting for equal observation weights
    if( !laserISL && !KbandISL && interSatelliteTracking ) { ISLMap << "GNSSISL"; }
    if( laserISL ) {
        stddevCodeISL = 1.0e-7;  codeWeightISL = pow( stddevCodeISL, -2.0 );
        stddevPhaseISL = 1.0e-0; phaseWeightISL = 0;
        ISLMap << "LaserISL";
    }
    if( KbandISL ){
        stddevCodeISL = 5.0e-1;  codeWeightISL = 0;
        stddevPhaseISL = 3.0e-5; phaseWeightISL = pow( stddevPhaseISL, -2.0 );
        ISLMap << "KbandISL";
    }
    if( clockPhase ) { ISLMap << "_1C"; }


    // // Error Settings
    bool noiseError = false;
    std::stringstream errorMap; errorMap << "/error";
    double stddevCode  = 1.0e-0;  double codeWeight = pow( stddevCode, -2.0 ); // sigma^-2 = weight (sigma = stddev, sigma^2 = variance)
    double stddevPhase  = 1.0e-3; double phaseWeight = pow( stddevPhase, -2.0 );
//    codeWeight = 1; phaseWeight = 1;  // Setting for equal observation weights
    if( noiseError ){ errorMap << "_NE"; }
    bool ephemerisError = false;
    double arcTime = 1800;
    double beginErrorSMA = 0.1;
//    double beginErrorSMA = 2.0;
    if( ephemerisError ){ errorMap << "_EE"; }
    if( !noiseError && !ephemerisError ){ errorMap << "_none"; }


    // // Characteristics (LTC, room for other error sources)
    bool lighttimeCorrection = false;

    constellations << "/"; ISLMap << "/"; errorMap << "/"; //_eph01 _test
    std::string outputSubFolder = "PECMEO" + timeStep.str() + constellations.str() + timeSpan.str() + ISLMap.str() + errorMap.str();
    std::cout << "Output Subfolder: " << outputSubFolder << std::endl;


    // // Load Spice kernels.
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );
*/

    /// OBTAINABLES
/*
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


    /// PERFORM SIMULATION: PROPAGATION AND OBSERVATION
    std::cout.precision(16);
    PodInputDataType observationsAndTimes = initiateAndObtainResults( simulationStartEpoch, simulationEndEpoch, fixedStepSize, numberOfSatellitesPEC, numberOfPlanesPEC,
                                                                      keplerElements, constellationList, interSatelliteTracking, numberOfSatellitesTotal, satelliteNames,
                                                                      numberOfConstellations, includedNames, epochTimes, bodyMap, allSatellitesPropagationHistory,
                                                                      allSatellitesPropagationHistoryKep);
    bool writeCartesianResult = true;
    bool writeKeplerianResult = true;
    writeResults( numberOfSatellitesTotal, satelliteNames, outputSubFolder, allSatellitesPropagationHistory, allSatellitesPropagationHistoryKep, writeCartesianResult, writeKeplerianResult);
*/

    /// OBTAIN RESULTS ESTIMATION
/*
    // Convergence criteria
    double correctionLimitPP = 1e-5; // 1 random standard value
    double correctionLimitKinematic = 1e-8;

    // Loop to select number of max iteration per estimation scheme
    for( int SPP = 5; SPP < 6; SPP++ )
    {
        for( int Kin = 5; Kin < 6; Kin++ )
        {
            for( int SPPISL = 0; SPPISL < 1; SPPISL++ )
            {
                for( int KinISL = 10; KinISL < 11; KinISL++ )
                {
                    std::cout << "\nDeze combi: " << SPP << Kin << SPPISL << KinISL << std::endl;
//                    std::chrono::system_clock::time_point time1 = std::chrono::system_clock::now();
//                    getCurrentTime();

                    // Set max iterations
                    int maxIterationsPP = SPP;
                    int maxIterationsKin = Kin;
                    int maxIterationsPPISL = SPPISL;
                    int maxIterationsKinISL = KinISL;

                    // Create empty data maps
                    std::map< std::string, std::map< double, Eigen::VectorXd > > allUsersDOPResults;
                    std::map< std::string, std::map< double, Eigen::Vector4d > > estimationsPPSequencePerUser;
                    std::map< std::string, std::map< double, Eigen::Vector4d > > estimationsKinSmartPerUser;
                    std::map< std::string, std::map< double, Eigen::Vector4d > > estimationsPPSequencePerUserISL;
                    std::map< std::string, std::map< double, Eigen::Vector4d > > estimationsKinSmartPerUserISL;
                    std::map< std::string, std::map< double, Eigen::Vector4d > > estimationsFinal;
                    std::map< std::string, std::map< double, Eigen::VectorXd > > allUsersPPIterations;
                    std::map< std::string, std::map< double, Eigen::VectorXd > > allUsersKinIterations;
                    std::map< int, Eigen::VectorXd > biasIterations;
                    std::map< int, std::map< std::string, std::map< double, Eigen::Vector4d > > > estimationsKinPerUserPerIteration;

                    // Perform the navigation estimation algorithms
                    performNavigationEstimation(observationsAndTimes, epochTimes, bodyMap, satelliteNames, interSatelliteTracking, numberOfSatellitesPEC, numberOfSatellitesTotal,
                                                maxIterationsPP, maxIterationsKin, correctionLimitPP, correctionLimitKinematic, codeWeight, phaseWeight, codeWeightISL, phaseWeightISL, lighttimeCorrection, noiseError, clockPhase,
                                                ephemerisError, arcTime, beginErrorSMA, maxIterationsPPISL, maxIterationsKinISL, stddevCode, stddevPhase, stddevCodeISL, stddevPhaseISL, allUsersDOPResults,
                                                estimationsPPSequencePerUser, estimationsKinSmartPerUser, estimationsPPSequencePerUserISL, estimationsKinSmartPerUserISL, estimationsFinal, allUsersPPIterations, allUsersKinIterations, biasIterations, estimationsKinPerUserPerIteration);


//                    writeEstimationResults( numberOfSatellitesPEC, numberOfSatellitesTotal, satelliteNames, includedNames, outputSubFolder, lighttimeCorrection, noiseError, maxIterationsPP, maxIterationsKin,
//                                            maxIterationsPPISL, maxIterationsKinISL, allUsersDOPResults, estimationsPPSequencePerUser, estimationsKinSmartPerUser, estimationsPPSequencePerUserISL,
//                                            estimationsKinSmartPerUserISL, allUsersPPIterations, allUsersKinIterations, biasIterations, false, false, true, false, false, false );

//                    writeEstimationEndResults( numberOfSatellitesPEC, satelliteNames, outputSubFolder, maxIterationsPP, maxIterationsKin, maxIterationsPPISL, maxIterationsKinISL, estimationsKinPerUserPerIteration );
                    writeEstimationFinalResults( numberOfSatellitesPEC, satelliteNames, outputSubFolder, estimationsFinal );

                    allUsersDOPResults.clear();
                    estimationsPPSequencePerUser.clear();
                    estimationsKinSmartPerUser.clear();
                    estimationsPPSequencePerUserISL.clear();
                    estimationsKinSmartPerUserISL.clear();
                    allUsersPPIterations.clear();
                    allUsersKinIterations.clear();
                    biasIterations.clear();

                }
            }
        }
    }
*/

    /// GPS Book Example
    /// From: J. Spilker Jr,  P. Axelrad,  B. W. Parkinson,  P. Enge,  Global Positioning System:  Theory and Applications, Volume I, American Institute of Aeronautics and Astronautics, 1996. Page 414
#define FOLDINGSTART {
    std::cout.precision(16);
    Eigen::Vector4d trueState; trueState << 6378137.0, 0.0, 0.0, 85000.0;
    Eigen::Vector4d aPrioriState; aPrioriState << 6377000.0, 3000.0, 4000.0, 0.0;
    Eigen::Vector3d aPrioriPosition = aPrioriState.segment(0, 3);

    Eigen::Vector4d SV01; SV01 << 22808160.9, -12005866.6, -6609526.5, 0.0;
    Eigen::Vector4d SV02; SV02 << 21141179.5, -2355056.3, -15985716.1, 0.0;
    Eigen::Vector4d SV08; SV08 << 20438959.3, -4238967.1, 16502090.2, 0.0;
    Eigen::Vector4d SV14; SV14 << 18432296.2, -18613382.5, -4672400.8, 0.0;
    Eigen::Vector4d SV17; SV17 << 21772117.8, 13773269.7, 6656636.4, 0.0;
    Eigen::Vector4d SV23; SV23 << 15561523.9, 3469098.6, -21303596.2, 0.0;
    Eigen::Vector4d SV24; SV24 << 13773316.6, 15929331.4, -16266254.4, 0.0;

    Eigen::MatrixXd positionMatrixGNSSSatellites( 7, 3 );
    positionMatrixGNSSSatellites.block( 0, 0, 1, 3) << SV01.segment(0,3).transpose();
    positionMatrixGNSSSatellites.block( 1, 0, 1, 3) << SV02.segment(0,3).transpose();
    positionMatrixGNSSSatellites.block( 2, 0, 1, 3) << SV08.segment(0,3).transpose();
    positionMatrixGNSSSatellites.block( 3, 0, 1, 3) << SV14.segment(0,3).transpose();
    positionMatrixGNSSSatellites.block( 4, 0, 1, 3) << SV17.segment(0,3).transpose();
    positionMatrixGNSSSatellites.block( 5, 0, 1, 3) << SV23.segment(0,3).transpose();
    positionMatrixGNSSSatellites.block( 6, 0, 1, 3) << SV24.segment(0,3).transpose();
    std::cout << "Position matrix: \n" << positionMatrixGNSSSatellites << std::endl;


    Eigen::VectorXd measuredPseudorange(7);
    measuredPseudorange << 21480623.2, 21971919.2, 22175603.9, 22747561.5, 21787252.3, 23541613.4, 24022907.4;
    std::cout << "\nMeasured pseudorange: \n" << measuredPseudorange.transpose() << std::endl;


    Eigen::MatrixXd designMatrix = getDesignMatrix( aPrioriPosition, positionMatrixGNSSSatellites );
    std::cout << "\nDesign matrix: \n" << designMatrix << std::endl;


    Eigen::VectorXd computedPseudorange(7);

    Eigen::VectorXd truePseudorange(7);
    Eigen::VectorXd measuredPhaserange(7);
    Eigen::VectorXd aPrioriBias(7);
    for( int sat = 0; sat < 7; sat++ )
    {
        computedPseudorange( sat ) = (aPrioriPosition.transpose() - positionMatrixGNSSSatellites.block( sat, 0, 1, 3 )).norm();

        truePseudorange( sat ) = (trueState.segment(0,3).transpose() - positionMatrixGNSSSatellites.block( sat, 0, 1, 3 )).norm();
        measuredPhaserange( sat ) = truePseudorange( sat ) + 0.5/(sat+1) + (sat+1)*1000;
        aPrioriBias(sat) = measuredPseudorange(sat) - measuredPhaserange( sat );
    }
//    computedPseudorange(0) = 21399408.0;
    std::cout << "\nComputed pseudorange: \n" << computedPseudorange.transpose() << std::endl;


    Eigen::Vector6d DOPresults = getDOPValues( aPrioriPosition, positionMatrixGNSSSatellites );
//    std::cout << "\nDOP Results: \n" << DOPresults.transpose() << std::endl;

    Eigen::MatrixXd covarianceMatrix = getCovarianceMatrix( aPrioriPosition, positionMatrixGNSSSatellites );
    std::cout << "\nX DOP: " << sqrt( covarianceMatrix(0,0) ) << std::endl;
    std::cout << "Y DOP: " << sqrt( covarianceMatrix(1,1) ) << std::endl;
    std::cout << "Z DOP: " << sqrt( covarianceMatrix(2,2) ) << std::endl;
    std::cout << "Total DOP: " << DOPresults(1) << std::endl;

    std::cout << "\na priori state:\n " << aPrioriState.transpose() << std::endl;

    Eigen::VectorXd deltaPosition = ( designMatrix.transpose() * designMatrix ).inverse() * designMatrix.transpose() * ( measuredPseudorange - computedPseudorange );
//    Eigen::VectorXd deltaPosition = ( designMatrix.transpose() * designMatrix ).inverse() * designMatrix.transpose() * ( computedPseudorange - measuredPseudorange ); // Other way around - no difference
    std::cout << "\nState correction: \n" << deltaPosition << std::endl;

    Eigen::VectorXd correctedPosition = aPrioriState + deltaPosition;
//    Eigen::VectorXd correctedPosition = aPrioriState - deltaPosition; // Other way around - no difference
    std::cout << "\nCorrected position: \n" << correctedPosition << std::endl;

    // Next round
    std::cout << "\nRedo" << std::endl;
    designMatrix = getDesignMatrix( correctedPosition.segment(0,3), positionMatrixGNSSSatellites );
    std::cout << "\nDesign matrix: \n" << designMatrix << std::endl;
    for( int sat = 0; sat < 7; sat++ )
    {
        computedPseudorange( sat ) = (correctedPosition.segment(0,3).transpose() - positionMatrixGNSSSatellites.block( sat, 0, 1, 3 )).norm() + correctedPosition(3);
    }
    std::cout << "\nComputed pseudorange: \n" << computedPseudorange.transpose() << std::endl;

    deltaPosition = ( designMatrix.transpose() * designMatrix ).inverse() * designMatrix.transpose() * ( measuredPseudorange - computedPseudorange );
//    deltaPosition = ( designMatrix.transpose() * designMatrix ).inverse() * designMatrix.transpose() * ( computedPseudorange - measuredPseudorange ); // Other way around - no difference
    std::cout << "\nState correction: \n" << deltaPosition << std::endl;

    correctedPosition = correctedPosition + deltaPosition;
//    correctedPosition = correctedPosition - deltaPosition;  // Other way around - no difference
    std::cout << "\nCorrected position: \n" << correctedPosition << std::endl;

    std::cout << "\nActual error: \n" << correctedPosition - trueState << std::endl;


//    // Next round
//    std::cout << "\nRedo2" << std::endl;
//    designMatrix = getDesignMatrix( correctedPosition.segment(0,3), positionMatrixGNSSSatellites );
//    std::cout << "\nDesign matrix: \n" << designMatrix << std::endl;
//    for( int sat = 0; sat < 7; sat++ )
//    {
//        computedPseudorange( sat ) = (correctedPosition.segment(0,3).transpose() - positionMatrixGNSSSatellites.block( sat, 0, 1, 3 )).norm() + correctedPosition(3);
//    }
//    std::cout << "\nComputed pseudorange: \n" << computedPseudorange.transpose() << std::endl;

//    deltaPosition = ( designMatrix.transpose() * designMatrix ).inverse() * designMatrix.transpose() * ( measuredPseudorange - computedPseudorange );
//    std::cout << "\nState correction: \n" << deltaPosition << std::endl;

//    correctedPosition = correctedPosition + deltaPosition;
//    std::cout << "\nCorrected position: \n" << correctedPosition << std::endl;

//    std::cout << "\nActual error: \n" << correctedPosition - trueState << std::endl;
#define FOLDINGEND }

    std::cout << "\nThis is a PECMEO Verification run!"<< std::endl;

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
