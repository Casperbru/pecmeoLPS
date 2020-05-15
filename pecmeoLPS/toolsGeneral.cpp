#include "toolsGeneral.h"
#include <memory>
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include "constellationLoader.h"
#include "toolsSimulation.h"
#include "toolsNavigation.h"
#include "applicationOutput.h"
#include <random>

using namespace tudat;
using namespace tudat::input_output;


//! Function to initiate the propagation and simulation and obtain Cartesian and Keplerian
//! coordinates of all the satellites at the simulation times and the simulations self
PodInputDataType initiateAndObtainResults(
        const double simulationStartEpoch,
        const double simulationEndEpoch,
        const double fixedStepSize,
        unsigned int numberOfSatellitesPEC,
        unsigned int numberOfPlanesPEC,
        std::vector< double > keplerElements,
        std::vector <constellationNames> constellationList,
        bool interSatelliteTracking,
        unsigned int& numberOfSatellitesTotal,
        std::vector< std::string >& satelliteNames,
        unsigned int& numberOfConstellations,
        std::string& includedNames,
        std::vector< double >& epochTimes,
        simulation_setup::NamedBodyMap& bodyMap,
        std::vector< std::map< double, Eigen::VectorXd > >& allSatellitesPropagationHistory,
        std::vector< std::map< double, Eigen::VectorXd > >& allSatellitesPropagationHistoryKep)
{
    /// SET BODIES, VEHICLES AND BODYMAP AND CREATE SYSTEM INITIAL STATE VECTOR
    // Create bodymap
    double stepsInSimulation;
    getBodyMap( simulationStartEpoch, simulationEndEpoch, fixedStepSize, stepsInSimulation, bodyMap);

    // Set initial conditions of PECMEO constellation.
    Eigen::MatrixXd initialConditionsInKeplerianElementsPEC( 6, numberOfSatellitesPEC );
    setPECMEO( numberOfSatellitesPEC, numberOfPlanesPEC, keplerElements, satelliteNames, initialConditionsInKeplerianElementsPEC );

    // For only looking at one satellite
//    numberOfSatellitesPEC = 1;

    // Set initial conditions of GNSS constellation.
    Eigen::MatrixXd initialConditionsInKeplerianElementsTotal;
    setConstellation(numberOfSatellitesPEC, initialConditionsInKeplerianElementsPEC, constellationList,
                                     satelliteNames, initialConditionsInKeplerianElementsTotal, numberOfSatellitesTotal, numberOfConstellations, includedNames );

    // Set accelerations for each satellite and add to body map.
    finalizeBodyMap( numberOfSatellitesTotal, satelliteNames, initialConditionsInKeplerianElementsTotal, bodyMap );

    epochTimes = getEpochTimes( simulationStartEpoch, simulationEndEpoch, fixedStepSize );
    getStateHistory( numberOfSatellitesTotal, satelliteNames, bodyMap, epochTimes, allSatellitesPropagationHistory, allSatellitesPropagationHistoryKep );

    // OBSERVATION TOOLS
    PodInputDataType observationsAndTimes;
    observationsAndTimes = observationsAndTimesCalculator(numberOfSatellitesTotal, numberOfSatellitesPEC, satelliteNames, epochTimes, bodyMap, interSatelliteTracking);

    if( interSatelliteTracking )
    {
        std::cout << "Intersatellite tracking is turned on!" << std::endl;
        includedNames.append("_ISL");
    }

    return observationsAndTimes;
}


//! Perform navigation and estimation process for a multiple users
void performNavigationEstimation(
        PodInputDataType observationsAndTimes,
        std::vector< double > epochTimes,
        simulation_setup::NamedBodyMap bodyMap,
        std::vector< std::string > satelliteNames, bool interSatelliteTracking,
        unsigned int numberOfSatellitesPEC, unsigned int numberOfSatellitesTotal,
        int maxIterationsPP, int maxIterationsKin,
        double correctionLimitPP, double correctionLimitKinematic,
        double codeWeight, double phaseWeight,
        double codeWeightISL, double phaseWeightISL,
        bool lighttimeCorrection, bool addNoise, bool clockPhase,
        bool ephemerisError, double arcTime, double beginErrorSMA,
        int maxIterationsPPISL, int maxIterationsKinISL,
        double stddevCode, double stddevPhase,
        double stddevCodeISL, double stddevPhaseISL,
        std::map< std::string, std::map< double, Eigen::VectorXd > >& allUsersDOPResults,
        std::map< std::string, std::map< double, Eigen::Vector4d > >& estimationsPPSequencePerUser,
        std::map< std::string, std::map< double, Eigen::Vector4d > >& estimationsKinSmartPerUser,
        std::map< std::string, std::map< double, Eigen::Vector4d > >& estimationsPPSequencePerUserISL,
        std::map< std::string, std::map< double, Eigen::Vector4d > >& estimationsKinSmartPerUserISL,
        std::map< std::string, std::map< double, Eigen::Vector4d > >& estimationsFinal,
        std::map< std::string, std::map< double, Eigen::VectorXd > >& allUsersPPIterations,
        std::map< std::string, std::map< double, Eigen::VectorXd > >& allUsersKinIterations,
        std::map< int, Eigen::VectorXd >& biasIterations,
        std::map< int, std::map< std::string, std::map< double, Eigen::Vector4d > > >& estimationsKinPerUserPerIteration)
{
    // Get ephemeris and visible satellites
    std::map< std::string, std::shared_ptr< Ephemeris > > satelliteEphemeris = getEphemeridesOfSatellites( satelliteNames, bodyMap );
    std::vector< std::shared_ptr< Ephemeris > > satelliteEphemerisGNSS = getEphemeridesOfGNSS( satelliteNames, bodyMap, numberOfSatellitesPEC );

    /// GENERATE OBSERVATIONS
    // Create empty maps for the observations, visible satellites, availability matrix and travel times between the satellites and fill them
    std::map< std::string, Eigen::MatrixXi > availabilityMatrixAllUsers;
    std::map< std::string, std::map< double, std::vector< int > > > visibleSatelliteIDsPerEpochPerUser;
    std::map< std::string, std::map< double, std::map< int, double > > > correctedRangesPerUserPerEpoch;
    std::map< std::string, std::map< double, std::map< int, double > > > correctedPhasesPerUserPerEpoch;
    std::map< std::string, std::map< double, std::map< int, double > > > travelTimesPerUserPerEpoch;
    getMeasurementsAllUsers( numberOfSatellitesPEC, satelliteNames, observationsAndTimes, epochTimes, satelliteEphemeris, lighttimeCorrection, addNoise, stddevCode, stddevPhase, stddevCodeISL, stddevPhaseISL, availabilityMatrixAllUsers, visibleSatelliteIDsPerEpochPerUser, correctedRangesPerUserPerEpoch, correctedPhasesPerUserPerEpoch, travelTimesPerUserPerEpoch);
    std::map< std::string, std::map< double, std::map< int, double > > > rangeData = correctedRangesPerUserPerEpoch;
    std::map< std::string, std::map< double, std::map< int, double > > > phaseData = correctedPhasesPerUserPerEpoch;

    // Set user satellite position for all epoch
    std::map< std::string, std::map< double, Eigen::Vector3d > > userSatellitePositionPerUserPerEpoch;
    std::map< std::string, std::map< double, Eigen::Vector4d > > aPrioriEstimationsPerUserFirstEpoch;
    Eigen::Vector3d aPrioriDeviation; aPrioriDeviation << 1000000,-1000000,-1000000; // Initial guess of user deviates 1000 km from true position in every direction
    std::map< double, int > numberOfObservationsPerEpoch;
    for( unsigned long long i = 0; i < numberOfSatellitesPEC; i++ )
    {
        std::string userSatellite = satelliteNames.at(i);
        Eigen::MatrixXi availabilityMatrix = availabilityMatrixAllUsers[ userSatellite ];
//        writeMatrixToFile( availabilityMatrix, "availabilityMat_"+userSatellite+"_2.dat", 2, tudat_applications::getOutputPath( ) + "PECMEOAvailability"); // Save availability matrix option
        std::map< double, Eigen::Vector3d > userSatellitePositionPerEpoch;
        for( unsigned long long i = 0; i < epochTimes.size(); i++ )
        {
            double epochTime = epochTimes.at(i);
            Eigen::Vector3d userSatellitePosition = satelliteEphemeris[ userSatellite ]->getCartesianState( epochTime ).segment( 0, 3 );
            userSatellitePositionPerUserPerEpoch[ userSatellite ][ epochTime ] = userSatellitePosition;
            userSatellitePositionPerEpoch[ epochTime ] = userSatellitePosition;
            // For first epoch state, set deviation
            if(i==0){ aPrioriEstimationsPerUserFirstEpoch[ userSatellite ][ epochTime ] << userSatellitePosition + aPrioriDeviation, 0; }

            numberOfObservationsPerEpoch[ epochTime ] += availabilityMatrix.row(static_cast<int>(i)).sum(); // for remove loop
        }

        std::map< double, Eigen::MatrixXd > positionMatrixGNSSSatellitesAllEpochs = getPositionMatrixGNSSSatellites( availabilityMatrix, satelliteEphemerisGNSS, epochTimes );
        allUsersDOPResults[ userSatellite ] = getDOPValues( epochTimes, userSatellitePositionPerEpoch, positionMatrixGNSSSatellitesAllEpochs );
    }

    // Remove the epochs where less observations than unknowns are present and get new availability matrix
    // In the simulation of 4 GNSS constellations, this does not occur, but for less GNSS constellation it might
    unsigned int numberOfGNSS = numberOfSatellitesTotal - numberOfSatellitesPEC;
    unsigned int numberOfIndependentPECs = numberOfSatellitesPEC;
    std::map< double, std::vector< std::string > > usersPerEpoch; // For other method
    std::map< std::string, std::vector< double > > epochTimeCorrected;    
    std::map< double, double > enoughObservationsCounter;
    std::map< double, double > notEnoughObservationsCounter;
    std::vector<double> epochTimesAll;
    std::vector< int > rowsToRemoveAll;
    for( unsigned long long i = 0; i < numberOfSatellitesPEC; i++ )
    {
        std::string userSatellite = satelliteNames.at(i);
        std::vector< int > rowsToRemove;
        for( unsigned long long i = 0; i < epochTimes.size(); i++ )
        {
            double epochTime = epochTimes.at(i);
            if( interSatelliteTracking )
            {
                if( availabilityMatrixAllUsers[ userSatellite ].block(i,0,1,numberOfGNSS).sum() < 4 )
                {
                    notEnoughObservationsCounter[ i ] += 1;
                    allUsersDOPResults[ userSatellite ].erase( epochTime ); // This line can be removed for data processing purposes when, you want all results.
                }
                else
                {
                    enoughObservationsCounter[ epochTime ] += 1;
                }

                if( enoughObservationsCounter[ epochTime ] >= numberOfIndependentPECs ) // Number of PECMEOSatellites that are "independent" and have enough observations per epoch.
                {
                    epochTimesAll.push_back( epochTime ); // sort and remove
                }

                if( notEnoughObservationsCounter[ i ] > (numberOfSatellitesPEC - numberOfIndependentPECs) )
                {
                    rowsToRemoveAll.push_back( static_cast<int>(i) );
                }
            }
            else
            {
                if( availabilityMatrixAllUsers[ userSatellite ].block(i,0,1,numberOfGNSS).sum() < 4 )
                {
                    rowsToRemove.push_back( static_cast<int>(i) );
                    allUsersDOPResults[ userSatellite ].erase( epochTime ); // This line can be removed for data processing purposes when, you want all results.
                }
                else
                {
                    epochTimeCorrected[ userSatellite ].push_back( epochTime );
                }
            }
        }
        if( !interSatelliteTracking )
        {
            removeRows( availabilityMatrixAllUsers[ userSatellite ], rowsToRemove );
        }
    }
    if( interSatelliteTracking )
    {
        sort( rowsToRemoveAll.begin(), rowsToRemoveAll.end() );
        rowsToRemoveAll.erase( unique( rowsToRemoveAll.begin(), rowsToRemoveAll.end() ), rowsToRemoveAll.end() );

        sort( epochTimesAll.begin(), epochTimesAll.end() );
        epochTimesAll.erase( unique( epochTimesAll.begin(), epochTimesAll.end() ), epochTimesAll.end() );
        for( unsigned long long i = 0; i < numberOfSatellitesPEC; i++ )
        {
            std::string userSatellite = satelliteNames.at(i);
            removeRows( availabilityMatrixAllUsers[ userSatellite ], rowsToRemoveAll );
        }
    }


    // Create ephemeris error base matrix
    std::default_random_engine generator(123);
    std::uniform_real_distribution<double> distributionErrorSMA(-beginErrorSMA,beginErrorSMA);
    int numberOfArcs = epochTimes.back()/arcTime+1;
//    std::cout << "back: " << epochTimes.back() << " and no of arcs: " << numberOfArcs << " and no of gnss: " << numberOfGNSS << std::endl;
    Eigen::MatrixXd ephemerisErrorBase( numberOfArcs, numberOfGNSS );
    for( double arc = 0; arc < numberOfArcs; arc++ )
    {
        for( double prn = 0; prn < numberOfGNSS; prn++ )
        {
            ephemerisErrorBase( arc, prn ) = distributionErrorSMA(generator);
        }
    }
//    std::cout << ephemerisErrorBase << std::endl;



    /// PERFORM ESTIMATIONS
    std::map< std::string, std::map< double, Eigen::Vector4d > > estimationsPerUser;

    // // SPP
    std::cout << "\nSPP" << std::endl;
    for( unsigned long long i = 0; i < numberOfSatellitesPEC; i++ )
    {
        std::string userSatellite = satelliteNames.at(i);
        epochTimes.clear(); epochTimes = epochTimeCorrected[ userSatellite ];        
        std::map< std::string, Eigen::MatrixXi > availabilityMatrixOneUser;
        availabilityMatrixOneUser[ userSatellite ] = availabilityMatrixAllUsers[ userSatellite ];
        if( interSatelliteTracking ){ availabilityMatrixOneUser[ userSatellite ].conservativeResize( availabilityMatrixOneUser[ userSatellite ].rows(), availabilityMatrixOneUser[ userSatellite ].cols()-numberOfSatellitesPEC ); epochTimes.clear(); epochTimes = epochTimesAll; }

//        std::chrono::system_clock::time_point time1 = std::chrono::system_clock::now();
        std::cout << "\n" << userSatellite << " positions estimation" << std::endl;
        std::map< std::string, std::map< double, Eigen::Vector4d > > estimationsPPSequenceOneUser;
        estimationsPPSequenceOneUser = getPointPositionSequence( epochTimes, satelliteEphemerisGNSS, rangeData, availabilityMatrixOneUser, allUsersPPIterations, codeWeight, codeWeightISL,
                                                                 ephemerisError, arcTime, ephemerisErrorBase, aPrioriEstimationsPerUserFirstEpoch, maxIterationsPP, correctionLimitPP, lighttimeCorrection, travelTimesPerUserPerEpoch );
        estimationsPPSequencePerUser.insert(estimationsPPSequenceOneUser.begin(), estimationsPPSequenceOneUser.end());        
//        std::chrono::system_clock::time_point time2 = std::chrono::system_clock::now();
//        std::cout << userSatellite << " positions estimated with SPP with max " << maxIterationsPP << " iterations" << std::endl;
//        getCurrentTime( time1, time2 );
    }
    estimationsPerUser = estimationsPPSequencePerUser;


    // // SPP ISL
    if( interSatelliteTracking && maxIterationsPPISL > 0 ){
        std::cout << "\nSPPISL" << std::endl;
        epochTimes.clear(); epochTimes = epochTimesAll;

//        std::chrono::system_clock::time_point time1 = std::chrono::system_clock::now();
        estimationsPPSequencePerUserISL = getPointPositionSequence( epochTimes, satelliteEphemerisGNSS, rangeData, availabilityMatrixAllUsers, allUsersPPIterations, codeWeight, codeWeightISL,
                                                                 ephemerisError, arcTime, ephemerisErrorBase, estimationsPerUser, maxIterationsPPISL, correctionLimitPP, lighttimeCorrection, travelTimesPerUserPerEpoch ); //  estimationsKinSmartPerUser - estimationsPPSequencePerUser
        estimationsPerUser = estimationsPPSequencePerUserISL;
//        std::cout << "All users positions estimated with SPPISL with max " << maxIterationsPPISL << " iterations" << std::endl;
//        std::chrono::system_clock::time_point time2 = std::chrono::system_clock::now();
//        getCurrentTime( time1, time2 );
    }


    // // KOD
    if( maxIterationsKin > 0 ){
        std::cout << "\nKOD" << std::endl;
        for( unsigned long long i = 0; i < numberOfSatellitesPEC; i++ )
        {
            std::string userSatellite = satelliteNames.at(i);
            epochTimes.clear(); epochTimes = epochTimeCorrected[ userSatellite ];
            std::map< std::string, Eigen::MatrixXi > availabilityMatrixOneUser;
            availabilityMatrixOneUser[ userSatellite ] = availabilityMatrixAllUsers[ userSatellite ];
            if( interSatelliteTracking ){ availabilityMatrixOneUser[ userSatellite ].conservativeResize( availabilityMatrixOneUser[ userSatellite ].rows(), availabilityMatrixOneUser[ userSatellite ].cols()-numberOfSatellitesPEC ); epochTimes.clear(); epochTimes = epochTimesAll; }

//            std::chrono::system_clock::time_point time1 = std::chrono::system_clock::now();
            std::cout << "\n" << userSatellite << " positions estimation" << std::endl;
            std::map< std::string, std::map< double, Eigen::Vector4d > > estimationsKinSmartOneUser;
            estimationsKinSmartOneUser = getKinematicEstimation( epochTimes, satelliteEphemerisGNSS, rangeData, phaseData, availabilityMatrixOneUser, allUsersKinIterations, biasIterations, estimationsKinPerUserPerIteration,
                                                                 ephemerisError, arcTime, ephemerisErrorBase, estimationsPerUser, maxIterationsKin, correctionLimitKinematic, codeWeight, phaseWeight, codeWeightISL, phaseWeightISL,
                                                                 lighttimeCorrection, travelTimesPerUserPerEpoch, 2 );
            estimationsKinSmartPerUser.insert(estimationsKinSmartOneUser.begin(), estimationsKinSmartOneUser.end());

//            std::chrono::system_clock::time_point time2 = std::chrono::system_clock::now();
//            std::cout << userSatellite << " positions estimated with KOD with max " << maxIterationsKOD << " iterations" << std::endl;
//            getCurrentTime( time1, time2 );
        }
        estimationsPerUser = estimationsKinSmartPerUser;
    }


    // // KOD ISL
    estimationsKinPerUserPerIteration.clear();
    if( interSatelliteTracking && maxIterationsKinISL > 0 ){
        std::cout << "\nKODISL" << std::endl;
        epochTimes.clear(); epochTimes = epochTimesAll;

//        std::chrono::system_clock::time_point time1 = std::chrono::system_clock::now();
        estimationsKinSmartPerUserISL = getKinematicEstimation( epochTimes, satelliteEphemerisGNSS, rangeData, phaseData, availabilityMatrixAllUsers, allUsersKinIterations, biasIterations, estimationsKinPerUserPerIteration,
                                                             ephemerisError, arcTime, ephemerisErrorBase, estimationsPerUser, maxIterationsKinISL, correctionLimitKinematic, codeWeight, phaseWeight, codeWeightISL, phaseWeightISL, // - estimationsPPSequencePerUser
                                                             lighttimeCorrection, travelTimesPerUserPerEpoch, 2, clockPhase );
        estimationsPerUser = estimationsKinSmartPerUserISL;
//        std::cout << "All users positions estimated with KODISL with max " << maxIterationsKinISL << " iterations" << std::endl;
//        std::chrono::system_clock::time_point time2 = std::chrono::system_clock::now();
//        getCurrentTime( time1, time2 );
    }
     estimationsKinSmartPerUserISL = estimationsPerUser;
     estimationsFinal = estimationsPerUser;


    // Printing the results of the estimation and saving them in the created maps    
    for( unsigned i = 0; i < numberOfSatellitesPEC; i++ )
    {
        std::string userSatellite = satelliteNames.at(i);
        if( !interSatelliteTracking ){ epochTimes.clear(); epochTimes = epochTimeCorrected[ userSatellite ]; }
        for( unsigned long long i = 0; i < epochTimes.size(); i++ )
        {
            double epochTime = epochTimes.at(i);
            Eigen::Vector4d userPosition; userPosition << userSatellitePositionPerUserPerEpoch[ userSatellite ][ epochTime ], 0;

//            std::cout << "\n" << userSatellite << " with " << maxIterationsPP << maxIterationsKin << maxIterationsPPISL << maxIterationsKinISL << " at epoch " << epochTime << " has: " << std::endl;
//            std::cout << "Correct position : " << userPosition.transpose() << std::endl;

//            if(maxIterationsPP>0){std::cout << "PP estimated     : " << estimationsPPSequencePerUser[ userSatellite ][ epochTime ].transpose() << std::endl;}
//            if(maxIterationsKin>0){std::cout << "Kin estimated    : " << estimationsKinSmartPerUser[ userSatellite ][ epochTime ].transpose() << std::endl;}
//            if(maxIterationsPPISL>0){std::cout << "PP ISL estimated : " << estimationsPPSequencePerUserISL[ userSatellite ][ epochTime ].transpose() << std::endl;}
//            if(maxIterationsKinISL>0){std::cout << "Kin ISL estimated: " << estimationsKinSmartPerUserISL[ userSatellite ][ epochTime ].transpose() << std::endl;}

//            std::cout << "Errors in user x,y,z position:" << std::endl;
//            if(maxIterationsPP>0){    std::cout << "PP     : " << (estimationsPPSequencePerUser[ userSatellite ][ epochTime ].segment(0,3) - userPosition).transpose() << std::endl;}
//            if(maxIterationsKin>0){   std::cout << "Kin    : " << (estimationsKinSmartPerUser[ userSatellite ][ epochTime ].segment(0,3) - userPosition ).transpose() << std::endl;}
//            if(maxIterationsPPISL>0){ std::cout << "PP ISL : " << (estimationsPPSequencePerUserISL[ userSatellite ][ epochTime ].segment(0,3) - userPosition).transpose() << std::endl;}
//            if(maxIterationsKinISL>0){std::cout << "Kin ISL: " << (estimationsKinSmartPerUserISL[ userSatellite ][ epochTime ].segment(0,3) - userPosition ).transpose() << std::endl;}
        }
    }
}


//! Obtain the current time or measure time between two points in the code
void getCurrentTime(
        std::chrono::system_clock::time_point start,
        std::chrono::system_clock::time_point end)
{
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);
    std::time_t start_time = std::chrono::system_clock::to_time_t(start);

    if( start == end ){ std::cout << "current time: " << std::ctime(&end_time); }
    else
    {
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "started computation at " << std::ctime(&start_time)
                  << "finished computation at " << std::ctime(&end_time)
                  << "elapsed time: " << elapsed_seconds.count() << "s\n";
    }
}


//! Obtain the parameters to be optimised
//! To be used for optimisation problem
std::vector< double > obtainParamaterToBeOptimised(
        std::map< std::string, std::map< double, Eigen::VectorXd > > allUsersDOPResults)
{
    std::vector< double > epochVector, linkVector, GDOPVector;
    for( std::map< std::string, std::map< double, Eigen::VectorXd > >::iterator itrSatelliteMap = allUsersDOPResults.begin(); itrSatelliteMap != allUsersDOPResults.end(); itrSatelliteMap++ )
    {
        std::string currentSatellite = itrSatelliteMap->first;
        std::map< double, Eigen::VectorXd > epochMap = itrSatelliteMap->second;
        for( std::map< double, Eigen::VectorXd >::iterator itrEpochMap = epochMap.begin(); itrEpochMap != epochMap.end(); itrEpochMap++ )
        {
            double currentEpoch = itrEpochMap->first;
            Eigen::VectorXd resultsVector = itrEpochMap->second;
            epochVector.push_back( currentEpoch );
            linkVector.push_back( resultsVector( 0 ) );
            GDOPVector.push_back( resultsVector( 1 ) );
        }
    }
    double averageLink = std::accumulate( linkVector.begin(), linkVector.end(), 0.0)/linkVector.size();
    double averageGDOP = std::accumulate( GDOPVector.begin(), GDOPVector.end(), 0.0)/GDOPVector.size();
    std::vector< double > results;
    results = { averageLink, averageGDOP };
    return results;
}


//! Retrieve optimisation parameters for all users
//! Gives vector with average links, GDOP or PDOP over all users and epochs
//! Only used for basic computations and verification
void retrieveOptimisationParameters(
        PodInputDataType observationsAndTimes,
        std::vector< double > epochTimesTrue,
        simulation_setup::NamedBodyMap bodyMap,
        std::vector< std::string > satelliteNames,
        unsigned int numberOfSatellitesPEC,        
        std::map< std::string, std::map< double, Eigen::VectorXd > >& allUsersDOPResults,
        std::vector< double >& results)
{
    // Get ephemeris and visible satellites
    std::map< std::string, std::shared_ptr< Ephemeris > > satelliteEphemeris = getEphemeridesOfSatellites( satelliteNames, bodyMap );
    std::vector< std::shared_ptr< Ephemeris > > satelliteEphemerisGNSS = getEphemeridesOfGNSS( satelliteNames, bodyMap, numberOfSatellitesPEC );

    // Create empty maps for the observations, visible satellites, availability matrix and travel times between the satellites and fill them
    std::map< std::string, Eigen::MatrixXi > availabilityMatrixAllUsers;
    std::map< std::string, std::map< double, std::vector< int > > > visibleSatelliteIDsPerEpochPerUser;
    std::map< std::string, std::map< double, std::map< int, double > > > correctedRangesPerUserPerEpoch;
    std::map< std::string, std::map< double, std::map< int, double > > > correctedPhasesPerUserPerEpoch;
    std::map< std::string, std::map< double, std::map< int, double > > > travelTimesPerUserPerEpoch;

    unsigned long long numberOfEpochs = 0;
    std::vector< double > epochTimes;

    // Adding identifier to satellites for indexing
    std::vector< std::string > satelliteNamesGNSS = satelliteNames;
    std::rotate( satelliteNamesGNSS.begin(), satelliteNamesGNSS.begin() + numberOfSatellitesPEC, satelliteNamesGNSS.end() );
    unsigned long long numberOfSatellites = satelliteNamesGNSS.size();
    std::map< std::string, int > satelliteIDs;
    std::map< int, std::string > satelliteIDsNames;
    int satelliteID = 0;
    for( unsigned long long j = 0; j < numberOfSatellites; j++)
    {
        satelliteIDs[ satelliteNamesGNSS.at(j) ] = satelliteID;
        satelliteIDsNames[ satelliteID ] = satelliteNamesGNSS.at(j);
        satelliteID += 1;       
    }

    // Opening the simulated observations object and restructure the data
    for( PodInputDataType::iterator itrObservationsAndTimesMap = observationsAndTimes.begin(); itrObservationsAndTimesMap != observationsAndTimes.end(); itrObservationsAndTimesMap++ )
    {
        SingleObservablePodInputType singleObservablePodInputType = itrObservationsAndTimesMap->second;

        // Loop over different visible GNSS satellites
        for ( SingleObservablePodInputType::iterator itrSingleObservableMap = singleObservablePodInputType.begin( ) ; itrSingleObservableMap != singleObservablePodInputType.end( ) ;
              itrSingleObservableMap++ )
        {
            // Obtain name and ID of current satellites over all epochs where it is visible
            LinkEnds currentLinkEnds = itrSingleObservableMap->first;
            std::vector< std::string> satLinkNames;
            for ( LinkEnds::iterator itrLinkEnds = currentLinkEnds.begin(); itrLinkEnds != currentLinkEnds.end( ); itrLinkEnds++ )
            {
                satLinkNames.push_back(itrLinkEnds->second.first); // Gives body and point on body as link (for satellite just vehicle name) [->first is linkEndTyp nn]
            }
            int currentSatelliteID = satelliteIDs[ satLinkNames.at(0) ]; // GNSS satellite (.at(1) is user satellite)

            // Obtain list of epochs where a connection is made
            epochTimes = itrSingleObservableMap->second.second.first;
            numberOfEpochs = epochTimes.size();

            // Loop over the epochs and create list of visible satellites, range and phase data and travel times per epoch
            for ( unsigned int i = 0 ; i < numberOfEpochs ; i++ )
            {
                visibleSatelliteIDsPerEpochPerUser[ satLinkNames.at(1) ][ epochTimes[ i ] ].push_back( currentSatelliteID );
            }
        }
    }

    // Creating new epoch times list, which includes all of the epochs.
    epochTimes.clear();
    epochTimes = epochTimesTrue;
    numberOfEpochs = epochTimes.size();

    // Loop over different users
    for( unsigned long long h = 0; h < numberOfSatellitesPEC; h++ )
    {
        // Current user
        std::string userSatellite = satelliteNames.at( h );

        // Loop over all epochs
        std::map< double, Eigen::Vector3d > userSatellitePositionPerEpoch;
        for( unsigned int i = 0; i < numberOfEpochs; i++ )
        {
            double epochTime = epochTimes.at( i );
            std::vector< int > visibleSatelliteIDs = visibleSatelliteIDsPerEpochPerUser[ userSatellite ][ epochTime ];
            unsigned long long numberOfVisibleSatellites = visibleSatelliteIDs.size();
            Eigen::MatrixXd designMatrix( numberOfVisibleSatellites, 4 );
            Eigen::RowVector3d userSatellitePosition = satelliteEphemeris[ userSatellite ]->getCartesianState( epochTime ).segment( 0, 3 );

            for( unsigned long long j = 0; j < numberOfVisibleSatellites; j++ )
            {
                int currentSatelliteID = visibleSatelliteIDs.at(j);
                std::string currentSatelliteName = satelliteIDsNames[ currentSatelliteID ];
                Eigen::RowVector3d GNSSSatellitePosition = satelliteEphemeris[ currentSatelliteName ]->getCartesianState( epochTime ).segment( 0, 3 );
                Eigen::RowVector3d positionDifference = userSatellitePosition - GNSSSatellitePosition;
                designMatrix.block( j, 0, 1, 4 ) << positionDifference / positionDifference.norm(), 1;
            }

            Eigen::MatrixXd covarianceMatrix = ( designMatrix.transpose() * designMatrix ).inverse();
            double valueGDOP = sqrt( covarianceMatrix(0,0) + covarianceMatrix(1,1) + covarianceMatrix(2,2) + covarianceMatrix(3,3) );
            double valuePDOP = sqrt( covarianceMatrix(0,0) + covarianceMatrix(1,1) + covarianceMatrix(2,2) );
            double valueHDOP = sqrt( covarianceMatrix(0,0) + covarianceMatrix(1,1) );
            double valueVDOP = sqrt( covarianceMatrix(2,2) );
            double valueTDOP = sqrt( covarianceMatrix(3,3) );

            Eigen::VectorXd resultsVector(9);
            resultsVector << numberOfVisibleSatellites, valueGDOP, valuePDOP, valueHDOP, valueVDOP, valueTDOP, userSatellitePosition.transpose();
            allUsersDOPResults[ userSatellite ][ epochTime ] = resultsVector;
        }
    }

    std::vector< double > epochVector, linkVector, GDOPVector, PDOPVector;
    for( std::map< std::string, std::map< double, Eigen::VectorXd > >::iterator itrSatelliteMap = allUsersDOPResults.begin(); itrSatelliteMap != allUsersDOPResults.end(); itrSatelliteMap++ )
    {
        std::string currentSatellite = itrSatelliteMap->first;
        std::map< double, Eigen::VectorXd > epochMap = itrSatelliteMap->second;
        for( std::map< double, Eigen::VectorXd >::iterator itrEpochMap = epochMap.begin(); itrEpochMap != epochMap.end(); itrEpochMap++ )
        {
            double currentEpoch = itrEpochMap->first;
            Eigen::VectorXd resultsVector = itrEpochMap->second;
            epochVector.push_back( currentEpoch );
            linkVector.push_back( resultsVector( 0 ) );
            if( isnan( resultsVector( 1 ) ) ){ resultsVector( 1 ) = 99999999; }
            GDOPVector.push_back( resultsVector( 1 ) );
//            std::cout << "Op epoch " << currentEpoch << " deze gdop: " << resultsVector( 1 ) << " en zoveel links:" << resultsVector( 0 )  << std::endl;
            PDOPVector.push_back( resultsVector( 2 ) );
        }
    }
    double averageLink = std::accumulate( linkVector.begin(), linkVector.end(), 0.0)/linkVector.size();
    double averageGDOP = std::accumulate( GDOPVector.begin(), GDOPVector.end(), 0.0)/GDOPVector.size();
    double averagePDOP = std::accumulate( PDOPVector.begin(), PDOPVector.end(), 0.0)/PDOPVector.size();

    double maxGDOP = *std::max_element( GDOPVector.begin(), GDOPVector.end() );

//    auto mininumNumberOfLinks = std::min_element( linkVector.begin(), linkVector.end() );
//    std::cout << "\nThis is minimum number of links: " << *mininumNumberOfLinks << " at " << epochVector.at(std::distance(linkVector.begin(), mininumNumberOfLinks)) << " seconds" << std::endl;
//    auto mininumGDOP = std::min_element( GDOPVector.begin(), GDOPVector.end() );
//    std::cout << "\nThis is minimum number of links: " << *mininumGDOP << " at " << epochVector.at(std::distance(GDOPVector.begin(), mininumGDOP)) << " seconds" << std::endl;
//    std::cout << "Dit is GDOP op 210: " << GDOPVector.at(7) << std::endl;
//    std::cout << "Dit is average links: " << averageLink << std::endl;
    std::cout << "Dit is average GDOP: " << averageGDOP << std::endl;

    // Get ephemeris perturbed propagation history
    std::vector< std::map< double, Eigen::VectorXd > > allSatellitesPertubedHistory;
    std::vector< std::map< double, Eigen::VectorXd > > allSatellitesPertubedHistoryKep;
    allSatellitesPertubedHistory.resize( numberOfSatellites );
    allSatellitesPertubedHistoryKep.resize( numberOfSatellites );
    double arcTime = 1800;
    double beginErrorSMA = 0.1;

    // Create ephemeris error base matrix
    std::default_random_engine generator(123);
    std::uniform_real_distribution<double> distributionErrorSMA(-beginErrorSMA,beginErrorSMA);
    int numberOfArcs = epochTimes.back()/arcTime+1;
    Eigen::MatrixXd ephemerisErrorBase( numberOfArcs, numberOfSatellites );
    for( double arc = 0; arc < numberOfArcs; arc++ )
    {
        for( double prn = 0; prn < numberOfSatellites; prn++ )
        {
            ephemerisErrorBase( arc, prn ) = distributionErrorSMA(generator);
        }
    }

    // Loop over users
    for( unsigned long long sat = 0; sat < numberOfSatellites-9; sat++)
    {
        std::string currentSatelliteName = satelliteIDsNames[ sat ];
        auto currentSatelliteEphemeris = satelliteEphemeris[ currentSatelliteName ];
        for( unsigned int i = 0; i < numberOfEpochs; i++ )
        {
            double epochTime = epochTimes.at( i );
            int arcNumber = floor(epochTime/arcTime);
            double arcEpochTime = floor(epochTime/arcTime)*arcTime;
            double arcPropagationTime = epochTime - arcEpochTime;
            double maxErrorSMA = ephemerisErrorBase( arcNumber, sat );
            double errorSMA = maxErrorSMA;
            double earthGravitationalParameter = 3.986004418e14;

            Eigen::Vector6d arcBeginStateCart = satelliteEphemeris[ currentSatelliteName ]->getCartesianState( arcEpochTime );
            Eigen::Vector6d arcBeginStateKep = orbital_element_conversions::convertCartesianToKeplerianElements(arcBeginStateCart, earthGravitationalParameter);
            arcBeginStateKep(0) += errorSMA;
            Eigen::Vector6d arcPropagatedStateKep = orbital_element_conversions::propagateKeplerOrbit(arcBeginStateKep, arcPropagationTime, earthGravitationalParameter );
            Eigen::Vector6d arcPropagatedStateCart = orbital_element_conversions::convertKeplerianToCartesianElements(arcPropagatedStateKep, earthGravitationalParameter);

            allSatellitesPertubedHistory[ sat ][ epochTime ] = arcPropagatedStateCart;
            allSatellitesPertubedHistoryKep[ sat ][ epochTime ] = arcPropagatedStateKep;
        }

        // Set filename for output data.
        std::stringstream outputFilename;
        outputFilename << currentSatelliteName << ".dat";

        // Write propagation history to file.
        writeDataMapToTextFile( allSatellitesPertubedHistory.at( sat ),
                                outputFilename.str( ),
                                "C:/tudatBundle/tudatApplications/stepTwo/SimulationOutput/PECMEO_30sec/GalileoGLONASSBeiDouGPS/1d/propagationResultsEph",
                                "",
                                std::numeric_limits< double >::digits10,
                                std::numeric_limits< double >::digits10,
                                "," );

        // Set filename for output data.
        std::stringstream outputFilename2;
        outputFilename2 << currentSatelliteName << "Kep.dat";

        // Write propagation history to file.
        writeDataMapToTextFile( allSatellitesPertubedHistoryKep.at( sat ),
                                outputFilename2.str( ),
                                "C:/tudatBundle/tudatApplications/stepTwo/SimulationOutput/PECMEO_30sec/GalileoGLONASSBeiDouGPS/1d/propagationResultsEph",
                                "",
                                std::numeric_limits< double >::digits10,
                                std::numeric_limits< double >::digits10,
                                "," );
    }

    results = { averageLink, averageGDOP, averagePDOP, maxGDOP };
}


//! Function to write all the results into files
void writeResults(
        unsigned int numberOfSatellitesTotal,
        std::vector< std::string > satelliteNames,
        std::string outputSubFolder,
        std::vector< std::map< double, Eigen::VectorXd > > allSatellitesPropagationHistory,
        std::vector< std::map< double, Eigen::VectorXd > > allSatellitesPropagationHistoryKep,
        bool writeCartesianResult,
        bool writeKeplerianResult)
{
//    writeCartesianResult = false;
//    writeKeplerianResult = false;
    std::string propagationSubFolder = "/propagationResults";
    for ( unsigned int i = 0; i < numberOfSatellitesTotal; i++ )
    {        
        if( writeCartesianResult )
        {
            // Set filename for output data.
            std::stringstream outputFilename;
            outputFilename << satelliteNames.at( i ) << ".dat";

            // Write propagation history to file.
            writeDataMapToTextFile( allSatellitesPropagationHistory.at( i ),
                                    outputFilename.str( ),
                                    tudat_applications::getOutputPath( ) + outputSubFolder + propagationSubFolder,
                                    "",
                                    std::numeric_limits< double >::digits10,
                                    std::numeric_limits< double >::digits10,
                                    "," );
        }

        if( writeKeplerianResult )
        {
            // Set filename for output data. // FOR KEPLERIAN
            std::stringstream outputFilename;
            outputFilename << satelliteNames.at( i ) << "Kep.dat";

            // Write propagation history to file.
            writeDataMapToTextFile( allSatellitesPropagationHistoryKep.at( i ),
                                    outputFilename.str( ),
                                    tudat_applications::getOutputPath( ) + outputSubFolder + propagationSubFolder,
                                    "",
                                    std::numeric_limits< double >::digits10,
                                    std::numeric_limits< double >::digits10,
                                    "," );
        }
    }
}


//! Function to write all the estimation results into files
void writeEstimationResults(
        const unsigned int numberOfSatellitesPEC,
        unsigned int numberOfSatellitesTotal,
        std::vector< std::string > satelliteNames,
        std::string includedNames,
        std::string outputSubFolder,
        bool lighttimeCorrection, bool addNoise,
        int maxIterationsPP, int maxIterationsKin,
        int maxIterationsPPISL, int maxIterationsKinISL,
        std::map< std::string, std::map< double, Eigen::VectorXd > > allUsersDOPResults,
        std::map< std::string, std::map< double, Eigen::Vector4d > > estimationsPPSequencePerUser,
        std::map< std::string, std::map< double, Eigen::Vector4d > > estimationsKinSmartPerUser,
        std::map< std::string, std::map< double, Eigen::Vector4d > > estimationsPPSequencePerUserISL,
        std::map< std::string, std::map< double, Eigen::Vector4d > > estimationsKinSmartPerUserISL,
        std::map< std::string, std::map< double, Eigen::VectorXd > > allUsersPPIterations,
        std::map< std::string, std::map< double, Eigen::VectorXd > > allUsersKinIterations,
        std::map< int, Eigen::VectorXd > biasIterations,
        bool writeDOPResults,
        bool writePPResults,
        bool writeKinResults,
        bool writePPIterations,
        bool writeKinIterations,
        bool writeBiasIterations)
{
    std::stringstream characteristics;
    if(lighttimeCorrection){ characteristics << "LTC"; }
    if(addNoise){ characteristics << "Noise"; }
    std::string iterationMap = "/iterationMap";
    std::string estimationMap = "/estimationMap";
    std::string geometryMap = "/geometryMap";

//    writeDOPResults = false;
//    writeKinIterations = true;
//    writeBiasIterations = false;
//    writePPIterations = true;

    for ( unsigned int i = 0; i < numberOfSatellitesTotal; i++ )
    {
        if( writeDOPResults && i < numberOfSatellitesPEC )
        {
            // Set filename for output data.
            std::stringstream outputFilename;
            outputFilename << satelliteNames.at( i ) << ".dat";

            // Write propagation history to file.
            writeDataMapToTextFile( allUsersDOPResults[ satelliteNames.at( i ) ],
                                    outputFilename.str( ),
                                    tudat_applications::getOutputPath( ) + outputSubFolder + geometryMap,
                                    "",
                                    std::numeric_limits< double >::digits10,
                                    std::numeric_limits< double >::digits10,
                                    "," );
        }

        if( writePPResults && i < numberOfSatellitesPEC && maxIterationsPP > 0 +10 )
        {
            // Set filename for output data.
            std::stringstream outputFilename;
            outputFilename << satelliteNames.at( i ) << "_" << maxIterationsPP << ",0,0,0" << ".dat";

            // Write propagation history to file.
            writeDataMapToTextFile( estimationsPPSequencePerUser[ satelliteNames.at( i ) ],
                                    outputFilename.str( ),
                                    tudat_applications::getOutputPath( ) + outputSubFolder + estimationMap,
                                    "",
                                    std::numeric_limits< double >::digits10,
                                    std::numeric_limits< double >::digits10,
                                    "," );
        }

        if( writeKinResults && i < numberOfSatellitesPEC && maxIterationsKin > 0 +10 )
        {
            // Set filename for output data.
            std::stringstream outputFilename;
            outputFilename << satelliteNames.at( i ) << "_" << maxIterationsPP << "," << maxIterationsKin << ",0,0" << ".dat";

            // Write propagation history to file.
            writeDataMapToTextFile( estimationsKinSmartPerUser[ satelliteNames.at( i ) ],
                                    outputFilename.str( ),
                                    tudat_applications::getOutputPath( ) + outputSubFolder + estimationMap,
                                    "",
                                    std::numeric_limits< double >::digits10,
                                    std::numeric_limits< double >::digits10,
                                    "," );
        }

        if( writePPResults && i < numberOfSatellitesPEC && maxIterationsPPISL > 0 +10 )
        {
            // Set filename for output data.
            std::stringstream outputFilename;
            outputFilename << satelliteNames.at( i ) << "_" << maxIterationsPP << "," << maxIterationsKin << "," << maxIterationsPPISL << ",0" << ".dat";

            // Write propagation history to file.
            writeDataMapToTextFile( estimationsPPSequencePerUserISL[ satelliteNames.at( i ) ],
                                    outputFilename.str( ),
                                    tudat_applications::getOutputPath( ) + outputSubFolder + estimationMap,
                                    "",
                                    std::numeric_limits< double >::digits10,
                                    std::numeric_limits< double >::digits10,
                                    "," );
        }

        if( writeKinResults && i < numberOfSatellitesPEC )
        {
            // Set filename for output data.
            std::stringstream outputFilename;
            outputFilename << satelliteNames.at( i ) << "_" << maxIterationsPP << "," << maxIterationsKin << "," << maxIterationsPPISL << "," << maxIterationsKinISL << ".dat";

            // Write propagation history to file.
            writeDataMapToTextFile( estimationsKinSmartPerUserISL[ satelliteNames.at( i ) ],
                                    outputFilename.str( ),
                                    tudat_applications::getOutputPath( ) + outputSubFolder + estimationMap,
                                    "",
                                    std::numeric_limits< double >::digits10,
                                    std::numeric_limits< double >::digits10,
                                    "," );
        }

        if( writePPIterations && i < numberOfSatellitesPEC )
        {
            // Set filename for output data.
            std::stringstream outputFilename;
            outputFilename << satelliteNames.at( i ) << "PPIterations" << "_" << maxIterationsPP << "," << maxIterationsKin << "," << maxIterationsPPISL << "," << maxIterationsKinISL << ".dat";

            // Write propagation history to file.
            writeDataMapToTextFile( allUsersPPIterations[ satelliteNames.at( i ) ],
                                    outputFilename.str( ),
                                    tudat_applications::getOutputPath( ) + outputSubFolder + iterationMap,
                                    "",
                                    std::numeric_limits< double >::digits10,
                                    std::numeric_limits< double >::digits10,
                                    "," );
        }

        if( writeKinIterations && i < numberOfSatellitesPEC )
        {
            // Set filename for output data.
            std::stringstream outputFilename;
            outputFilename << satelliteNames.at( i ) << "KinIterations" << "_" << maxIterationsPP << "," << maxIterationsKin << "," << maxIterationsPPISL << "," << maxIterationsKinISL << ".dat";

            // Write propagation history to file.
            writeDataMapToTextFile( allUsersKinIterations[ satelliteNames.at( i ) ],
                                    outputFilename.str( ),
                                    tudat_applications::getOutputPath( ) + outputSubFolder + iterationMap,
                                    "",
                                    std::numeric_limits< double >::digits10,
                                    std::numeric_limits< double >::digits10,
                                    "," );
        }

        if( writeBiasIterations )
        {
            // Set filename for output data.
            std::stringstream outputFilename;
            outputFilename << "BiasIterations" << includedNames << characteristics.str() << ".dat";

            // Write propagation history to file.
            writeDataMapToTextFile( biasIterations,
                                    outputFilename.str( ),
                                    tudat_applications::getOutputPath( ) + outputSubFolder + iterationMap,
                                    "",
                                    std::numeric_limits< double >::digits10,
                                    std::numeric_limits< double >::digits10,
                                    "," );
        }
    }
}

//! Function to write all the final estimation results into files with iterations
void writeEstimationEndResults(
        const unsigned int numberOfSatellitesPEC,
        std::vector< std::string > satelliteNames,
        std::string outputSubFolder,
        int maxIterationsPP, int maxIterationsKin,
        int maxIterationsPPISL, int maxIterationsKinISL,
        std::map< int, std::map< std::string, std::map< double, Eigen::Vector4d > > > estimationsKinPerUserPerIteration)
{
    std::stringstream characteristics;
    std::string estimationMap = "/estimationMap";

    for( int iteration = 0; iteration < maxIterationsKinISL+1; iteration++ )
    {
        for( unsigned int i = 0; i < numberOfSatellitesPEC; i++ )
        {
            // Set filename for output data.
            std::stringstream outputFilename;
            outputFilename << satelliteNames.at( i ) << "_" << maxIterationsPP << "," << maxIterationsKin << "," << maxIterationsPPISL << "," << iteration << ".dat";

            if( estimationsKinPerUserPerIteration[ iteration ][ satelliteNames.at( i ) ].empty() )
            {
//                std::cout << "Iteration " << iteration << " is empty" << std::endl;
                continue;
            }

            // Write propagation history to file.
            writeDataMapToTextFile( estimationsKinPerUserPerIteration[ iteration ][ satelliteNames.at( i ) ],
                                    outputFilename.str( ),
                                    tudat_applications::getOutputPath( ) + outputSubFolder + estimationMap,
                                    "",
                                    std::numeric_limits< double >::digits10,
                                    std::numeric_limits< double >::digits10,
                                    "," );
        }
    }
}


//! Function to write all the final estimation results into files
void writeEstimationFinalResults(
        const unsigned int numberOfSatellitesPEC,
        std::vector< std::string > satelliteNames,
        std::string outputSubFolder,
        std::map< std::string, std::map< double, Eigen::Vector4d > > estimationsFinal)
{
    std::stringstream characteristics;
    std::string estimationMap = "/estimationMap";

    for( unsigned int i = 0; i < numberOfSatellitesPEC; i++ )
    {
        // Set filename for output data.
        std::stringstream outputFilename;
        outputFilename << satelliteNames.at( i ) << ".dat";

        // Write propagation history to file.
        writeDataMapToTextFile( estimationsFinal[ satelliteNames.at( i ) ],
                                outputFilename.str( ),
                                tudat_applications::getOutputPath( ) + outputSubFolder + estimationMap,
                                "",
                                std::numeric_limits< double >::digits10,
                                std::numeric_limits< double >::digits10,
                                "," );
    }
}


//! Function to remove rows from a matrix
void removeRows(Eigen::MatrixXi& matrix, std::vector<int> rowsToRemove)
{
    unsigned long long numberOfRowsToRemove = rowsToRemove.size();
    std::reverse( rowsToRemove.begin(), rowsToRemove.end() );
    for( unsigned long long i = 0; i < numberOfRowsToRemove; i++ )
    {
        int rowToRemove = rowsToRemove.at(i);
        long long numRows = matrix.rows()-1;
        long long numCols = matrix.cols();

        if( rowToRemove < numRows )
            matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.bottomRows(numRows-rowToRemove);

        matrix.conservativeResize(numRows,numCols);
    }
}
