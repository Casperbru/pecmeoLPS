#include "toolsOptimisation.h"
#include "constellationLoader.h"
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include "toolsSimulation.h"

using namespace tudat;


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
    finalizeBodyMap( numberOfSatellitesTotal, satelliteNames, bodyMap );

    // Set initial state in vector form
    Eigen::VectorXd systemInitialState;
    systemInitialState = getInitialStateVector( bodyMap, numberOfSatellitesTotal, initialConditionsInKeplerianElementsTotal );


    // CREATE ACCELERATIONS & PROPAGATION SETTINGS
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings;
    std::shared_ptr< IntegratorSettings< > > integratorSettings;
    getIntegratorPropagatorSettings( simulationStartEpoch, simulationEndEpoch, fixedStepSize, bodyMap, numberOfSatellitesTotal, satelliteNames, systemInitialState, propagatorSettings, integratorSettings );

    // PROPAGATE ORBIT
    std::vector< double > propagatedTimes;
    std::map< std::string, std::map< double, Eigen::Vector3d > > coordinateHistoryPerSat;
    propagatedTimes = propagationCalculator( numberOfSatellitesTotal, satelliteNames, bodyMap, integratorSettings, propagatorSettings,
                                             allSatellitesPropagationHistory, allSatellitesPropagationHistoryKep, coordinateHistoryPerSat);
    epochTimes = propagatedTimes;

    // OBSERVATION TOOLS
    PodInputDataType observationsAndTimes;
    observationsAndTimes = observationsAndTimesCalculator(numberOfSatellitesTotal, numberOfSatellitesPEC, satelliteNames, propagatedTimes, bodyMap, interSatelliteTracking);

    if( interSatelliteTracking )
    {
        std::cout << "Intersatellite tracking is turned on!" << std::endl;
        includedNames.append("_ISL");
    }

    return observationsAndTimes;
}


//! Retrieve optimisation parameters for all users
//! Gives vector with average links, GDOP or PDOP over all users and epochs
void retrieveOptimisationParameters(
        PodInputDataType observationsAndTimes,
        std::vector< double > epochTimesTrue,
        simulation_setup::NamedBodyMap bodyMap,
        std::vector< std::string > satelliteNames,
        unsigned int numberOfSatellitesPEC,
        bool lighttimeCorrection, bool addNoise,
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

//    std::cout << epochTimes.back() << std::endl;

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
//            std::cout << satLinkNames.at(0) << " has " << numberOfEpochs << " with " << satLinkNames.at(1) << std::endl;

            // Loop over the epochs and create list of visible satellites, range and phase data and travel times per epoch
            for ( unsigned int i = 0 ; i < numberOfEpochs ; i++ )
            {
                visibleSatelliteIDsPerEpochPerUser[ satLinkNames.at(1) ][ epochTimes[ i ] ].push_back( currentSatelliteID );
            }
        }
    }

    std::cout << "SWAT?" << std::endl;
    // Creating new epoch times list, which includes all of the epochs.
    epochTimes.clear();
    epochTimes = epochTimesTrue;
    numberOfEpochs = epochTimes.size();
    std::cout << numberOfEpochs << std::endl;

    // Loop over different users
    for( unsigned long long h = 0; h < numberOfSatellitesPEC; h++ )
    {
        // Current user
        std::string userSatellite = satelliteNames.at( h );
//        std::cout << "nr " << h << std::endl;
//        std::cout << userSatellite << std::endl;

        // Loop over all epochs
        std::map< double, Eigen::Vector3d > userSatellitePositionPerEpoch;
        for( unsigned int i = 0; i < numberOfEpochs; i++ )
        {
//            std::cout << "epoch " << epochTimes.at( i ) << std::endl;
            double epochTime = epochTimes.at( i );
            std::vector< int > visibleSatelliteIDs = visibleSatelliteIDsPerEpochPerUser[ userSatellite ][ epochTime ];
            unsigned long long numberOfVisibleSatellites = visibleSatelliteIDs.size();
            Eigen::MatrixXd designMatrix( numberOfVisibleSatellites, 4 );
            Eigen::RowVector3d userSatellitePosition = satelliteEphemeris[ userSatellite ]->getCartesianState( epochTime ).segment( 0, 3 );
//            std::cout << "nr of sat: " << numberOfVisibleSatellites << std::endl;

            for( unsigned long long j = 0; j < numberOfVisibleSatellites; j++ )
            {
//                std::cout << "sat: " << satelliteIDsNames[ visibleSatelliteIDs.at(j) ] << std::endl;
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
            resultsVector << numberOfVisibleSatellites, valueGDOP, valuePDOP, valueHDOP, valueVDOP, valueTDOP, userSatellitePosition;
            allUsersDOPResults[ userSatellite ][ epochTime ] = resultsVector;
        }
    }

//        std::cout << "MAK2" << std::endl;
////        availabilityMatrixAllUsers[ userSatellite ] = availabilityMatrix;
//        std::map< double, Eigen::MatrixXd > positionMatrixGNSSSatellitesAllEpochs = getPositionMatrixGNSSSatellites( availabilityMatrix, satelliteEphemerisGNSS, epochTimes );
//        std::cout << "MAK?" << std::endl;
//        allUsersDOPResults[ userSatellite ] = getDOPValues( epochTimes, userSatellitePositionPerEpoch, positionMatrixGNSSSatellitesAllEpochs );
//        std::cout << "MAK" << std::endl;
//    getMeasurementsAllUsers( numberOfSatellitesPEC, satelliteNames, observationsAndTimes, epochTimes, satelliteEphemeris, lighttimeCorrection, addNoise, availabilityMatrixAllUsers, visibleSatelliteIDsPerEpochPerUser, correctedRangesPerUserPerEpoch, correctedPhasesPerUserPerEpoch, travelTimesPerUserPerEpoch);
//    std::cout << "SWAY" << std::endl;

//    // Set user satellite position for all epoch
//    std::map< std::string, std::map< double, Eigen::Vector3d > > userSatellitePositionPerUserPerEpoch;
//    std::map< std::string, std::map< double, Eigen::Vector4d > > aPrioriEstimationsPerUserFirstEpoch;
//    Eigen::Vector3d aPrioriDeviation; aPrioriDeviation << 100,-100,-100; // Random value 1000000
//    std::map< double, int > numberOfObservationsPerEpoch;
//    for( unsigned long long i = 0; i < numberOfSatellitesPEC; i++ )
//    {
//        std::string userSatellite = satelliteNames.at(i);
//        Eigen::MatrixXi availabilityMatrix = availabilityMatrixAllUsers[ userSatellite ];
//        std::map< double, Eigen::Vector3d > userSatellitePositionPerEpoch;
//        for( unsigned long long i = 0; i < epochTimes.size(); i++ )
//        {
//            double epochTime = epochTimes.at(i);
//            Eigen::Vector3d userSatellitePosition = satelliteEphemeris[ userSatellite ]->getCartesianState( epochTime ).segment( 0, 3 );
//            userSatellitePositionPerUserPerEpoch[ userSatellite ][ epochTime ] = userSatellitePosition;
//            userSatellitePositionPerEpoch[ epochTime ] = userSatellitePosition;
//            if(i==0){ aPrioriEstimationsPerUserFirstEpoch[ userSatellite ][ epochTime ] << userSatellitePosition + aPrioriDeviation, 0; }

//            numberOfObservationsPerEpoch[ epochTime ] += availabilityMatrix.row(static_cast<int>(i)).sum(); // for remove loop
////            std::cout << "epoch " << epochTime << " has " << numberOfObservationsPerEpoch[ epochTime ] << std::endl;
//        }

//        std::map< double, Eigen::MatrixXd > positionMatrixGNSSSatellitesAllEpochs = getPositionMatrixGNSSSatellites( availabilityMatrix, satelliteEphemerisGNSS, epochTimes );
//        allUsersDOPResults[ userSatellite ] = getDOPValues( epochTimes, userSatellitePositionPerEpoch, positionMatrixGNSSSatellitesAllEpochs );
//    }

    std::cout << "SWAG" << std::endl;
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
            GDOPVector.push_back( resultsVector( 1 ) );
            PDOPVector.push_back( resultsVector( 2 ) );
        }
    }
    double averageLink = std::accumulate( linkVector.begin(), linkVector.end(), 0.0)/linkVector.size();
    double averageGDOP = std::accumulate( GDOPVector.begin(), GDOPVector.end(), 0.0)/GDOPVector.size();
    double averagePDOP = std::accumulate( PDOPVector.begin(), PDOPVector.end(), 0.0)/PDOPVector.size();

//    auto mininumNumberOfLinks = std::min_element( linkVector.begin(), linkVector.end() );
//    std::cout << "\nThis is minimum number of links: " << *mininumNumberOfLinks << " at " << epochVector.at(std::distance(linkVector.begin(), mininumNumberOfLinks)) << " seconds" << std::endl;
//    auto mininumGDOP = std::min_element( GDOPVector.begin(), GDOPVector.end() );
//    std::cout << "\nThis is minimum number of links: " << *mininumGDOP << " at " << epochVector.at(std::distance(GDOPVector.begin(), mininumGDOP)) << " seconds" << std::endl;
//    std::cout << "Dit is GDOP op 210: " << GDOPVector.at(7) << std::endl;
//    std::cout << "Dit is average links: " << averageLink << std::endl;
//    std::cout << "Dit is average GDOP: " << averageGDOP << std::endl;

    results = { averageLink, averageGDOP, averagePDOP };
}
