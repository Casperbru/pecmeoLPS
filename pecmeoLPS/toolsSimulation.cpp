#include "toolsSimulation.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Astrodynamics/ObservationModels/oneWayRangeObservationModel.h"
#include <random>
#include "Tudat/Astrodynamics/Ephemerides/keplerEphemeris.h"

using namespace tudat;
using namespace tudat::unit_conversions;
using namespace tudat::orbital_element_conversions;
using namespace tudat::ephemerides;
using namespace tudat::observation_models;


////////////////////////// SETTINGS //////////////////////////

//! Obtain the body map for the set simulation time
void getBodyMap(
        const double simulationStartEpoch,
        const double simulationEndEpoch,
        const double fixedStepSize,
        double& stepsInSimulation,
        simulation_setup::NamedBodyMap& bodyMap)
{
    // Set number of steps in simulation
    stepsInSimulation = ( ( simulationEndEpoch - simulationStartEpoch ) / fixedStepSize ) + 1;

    // Define environment settings
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( { "Earth" },
                                    simulationStartEpoch - 10.0 * fixedStepSize, simulationEndEpoch + 10.0 * fixedStepSize );
    bodySettings[ "Earth" ]->ephemerisSettings = std::make_shared< simulation_setup::ConstantEphemerisSettings >(
                Eigen::Vector6d::Zero( ), "SSB", "J2000" );
    bodySettings[ "Earth" ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    bodySettings[ "Earth" ]->atmosphereSettings = NULL;
    double earthRadius = 6378.1e3;
    bodySettings[ "Earth" ]->shapeModelSettings = std::make_shared< SphericalBodyShapeSettings >( earthRadius );

    // Create Earth object
    bodyMap = simulation_setup::createBodies( bodySettings );
}


//! Finalize body map by adding the satellites
void finalizeBodyMap(
        unsigned int numberOfSatellitesTotal,
        std::vector< std::string > satelliteNames,
        Eigen::MatrixXd initialConditionsInKeplerianElementsTotal,
        simulation_setup::NamedBodyMap& bodyMap)
{
    double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( ); // Kepler to cartesian that is going to be used to propagate

    std::string currentSatelliteName;
    Eigen::Vector6d currentInitialStateInKeplerianElements;
    for ( unsigned int i = 0; i < numberOfSatellitesTotal; i++ )
    {
        currentSatelliteName = satelliteNames.at( i );
        currentInitialStateInKeplerianElements = initialConditionsInKeplerianElementsTotal.block(0,i,6,1);

        bodyMap[ currentSatelliteName ] = std::make_shared< simulation_setup::Body >( );
        // Option for non Kepler propagation and ephemeris
//        bodyMap[ currentSatelliteName ]->setEphemeris( std::make_shared< TabulatedCartesianEphemeris< > >(
//                                                           std::shared_ptr< interpolators::OneDimensionalInterpolator
//                                                           < double, Eigen::Vector6d > >( ), "Earth", "J2000" ) );

        bodyMap[ currentSatelliteName ]->setEphemeris( std::make_shared< KeplerEphemeris >(
                                                           currentInitialStateInKeplerianElements,
                                                           0.0,
                                                           earthGravitationalParameter,
                                                           "Earth",
                                                           "J2000" ) );
    }

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );
}


////////////////////////// SIMULATION //////////////////////////

//! Get initial Cartesian state of all satellites in one vector
Eigen::VectorXd getInitialStateVector(
        simulation_setup::NamedBodyMap bodyMap,
        unsigned int numberOfSatellitesTotal,
        Eigen::MatrixXd initialConditionsInKeplerianElementsTotal)
{
    // Convert initial conditions to Cartesian elements of all constellations at once
    Eigen::MatrixXd initialConditions( 6, numberOfSatellitesTotal ); // For Cartesian initial condition set matrix (6 x # single satellites)
    double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( ); // Kepler to cartesian that is going to be used to propagate
    for ( unsigned int i = 0; i < numberOfSatellitesTotal; i++ )
    {
        Eigen::Vector6d initKepl = initialConditionsInKeplerianElementsTotal.col( i ).cast< double >();
        initialConditions.col( i ) = convertKeplerianToCartesianElements(
                    initKepl, static_cast< double >(earthGravitationalParameter) );
    }

    Eigen::VectorXd systemInitialState = Eigen::VectorXd( 6 * numberOfSatellitesTotal ); // For Cartesian initial condition set matrix (6 * # single satellites x 1)
    for ( unsigned int i = 0; i < numberOfSatellitesTotal; i++ )
    {
        systemInitialState.segment( i * 6, 6 ) = initialConditions.col( i );
    }
    return systemInitialState;
};


//! Obtain the ephemeris of all satellites
std::map< std::string, std::shared_ptr< Ephemeris > > getEphemeridesOfSatellites(
        std::vector< std::string > satelliteNames,
        simulation_setup::NamedBodyMap bodyMap)
{
    std::map< std::string, std::shared_ptr< Ephemeris > > satelliteEphemeris;
    for( unsigned int i = 0; i < static_cast<unsigned>( satelliteNames.size() ); i++ )
    {
        satelliteEphemeris[ satelliteNames.at(i) ] = bodyMap.at( satelliteNames.at(i) )->getEphemeris( );
    }
    return satelliteEphemeris;
};


//! Obtain ephemeris of GNSS satellites indexed by numbers
std::vector< std::shared_ptr< Ephemeris > > getEphemeridesOfGNSS(
        std::vector< std::string > satelliteNames,
        simulation_setup::NamedBodyMap bodyMap,
        unsigned int numberOfSatellitesPEC)
{
    std::vector< std::shared_ptr< Ephemeris > > satelliteEphemerisGNSS;
    auto numberOfSatellitesTotal = satelliteNames.size();
    for( auto i = numberOfSatellitesPEC; i < numberOfSatellitesTotal; i++ )
    {
        satelliteEphemerisGNSS.push_back( bodyMap.at( satelliteNames.at(i) )->getEphemeris( ) );
//        std::cout<< "Get ephemeris of ID " << i << " of " << satelliteNames.at(i)<<std::endl;
    }
    for( unsigned i = 0; i < numberOfSatellitesPEC; i++ )
    {
        satelliteEphemerisGNSS.push_back( bodyMap.at( satelliteNames.at(i) )->getEphemeris( ) );
//        std::cout<< "Get ephemeris of ID " << i << " of " << satelliteNames.at(i)<<std::endl;
    }
    return satelliteEphemerisGNSS;
}


//! Obtain the epochs over which the simulation is performed.
std::vector< double > getEpochTimes(
        const double simulationStartEpoch,
        const double simulationEndEpoch,
        const double fixedStepSize)
{
    double simulationPeriod = simulationEndEpoch-simulationStartEpoch;
    double numberOfEpochs = simulationPeriod/fixedStepSize;
    std::vector< double > epochTimes;
    for( int epoch = 0; epoch <= numberOfEpochs; epoch++ )
    {
        double currentEpochTime = epoch * fixedStepSize;
        epochTimes.push_back( currentEpochTime );
//        std::cout << currentEpochTime << std::endl;
    }

    return epochTimes;
}


//! Obtain the correct propagation history of the satellites.
void getStateHistory (
        unsigned int numberOfSatellitesTotal,
        std::vector< std::string > satelliteNames,
        simulation_setup::NamedBodyMap bodyMap,
        std::vector< double > epochTimes,
        std::vector< std::map< double, Eigen::VectorXd > >& allSatellitesPropagationHistory,
        std::vector< std::map< double, Eigen::VectorXd > >& allSatellitesPropagationHistoryKep)
{
    // Define Earth gravitational parameter for kepler writing!
    double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    allSatellitesPropagationHistory.resize( numberOfSatellitesTotal );
    allSatellitesPropagationHistoryKep.resize( numberOfSatellitesTotal );

    // Create map with all satellite states per epoch, both in Kepler and Cartesian elements
    for( unsigned int sat = 0; sat < numberOfSatellitesTotal; sat++ )
    {
        std::string currentSatellite = satelliteNames.at( sat );
        std::shared_ptr< Ephemeris > currentSatelliteEphemeris = bodyMap.at( satelliteNames.at(sat) )->getEphemeris( );

        for( unsigned int epoch = 0; epoch < epochTimes.size(); epoch++ )
        {
            double epochTime = epochTimes.at( epoch );
            Eigen::Vector6d currentSatelliteCartesian = currentSatelliteEphemeris->getCartesianState( epochTime );
            allSatellitesPropagationHistory[ sat ][ epochTime ] = currentSatelliteEphemeris->getCartesianState( epochTime );
            allSatellitesPropagationHistoryKep[ sat ][ epochTime ] = convertCartesianToKeplerianElements( currentSatelliteCartesian, earthGravitationalParameter );
        }
    }
}


//! Determine all the visible observations between the link ends over time
PodInputDataType observationsAndTimesCalculator(
        unsigned int numberOfSatellitesTotal,
        unsigned int numberOfSatellitesPEC,
        std::vector< std::string >  satelliteNames,
        std::vector< double > propagatedTimes,
        simulation_setup::NamedBodyMap bodyMap,
        bool interSatelliteTracking)
{
    unsigned int numberOfSatellitesGNSS = numberOfSatellitesTotal - numberOfSatellitesPEC;
    // // LINK ENDS
    // Create list of link ends in which station is receiver and in which station is transmitter (with spacecraft other link end).
    std::vector< LinkEnds > stationReceiverLinkEnds;
    std::vector< LinkEnds > stationTransmitterLinkEnds;
    unsigned int numberOfGNSSLinks = 0;
    unsigned int numberOfISLLinks = 0;
//    unsigned int numberOfLinks = numberOfGNSSLinks + numberOfISLLinks;

    for( unsigned int i = 0; i < numberOfSatellitesPEC; i++ )
    {
        for( unsigned int j = numberOfSatellitesPEC; j < numberOfSatellitesTotal; j++ )
        {
            LinkEnds linkEnds;
            linkEnds[ transmitter ] = std::make_pair( satelliteNames.at( j ), "" );
            linkEnds[ receiver ] = std::make_pair( satelliteNames.at( i ), "" );
            stationTransmitterLinkEnds.push_back( linkEnds );

            linkEnds.clear( );
            linkEnds[ receiver ] = std::make_pair( satelliteNames.at( j ), "" );
            linkEnds[ transmitter ] = std::make_pair( satelliteNames.at( i ), "" );
            stationReceiverLinkEnds.push_back( linkEnds );
            numberOfGNSSLinks += 1;
        }
    }

    // Add intersatellite tracking to the link ends
    if( interSatelliteTracking )
    {
        for( unsigned int i = 0; i < numberOfSatellitesPEC; i++ )
        {
            for( unsigned int j = 0; j < numberOfSatellitesPEC; j++ )
            {
                if( j != i )
                {
                    LinkEnds linkEnds;
                    linkEnds[ transmitter ] = std::make_pair( satelliteNames.at( j ), "" );
                    linkEnds[ receiver ] = std::make_pair( satelliteNames.at( i ), "" );
                    stationTransmitterLinkEnds.push_back( linkEnds );

                    linkEnds.clear( );
                    linkEnds[ receiver ] = std::make_pair( satelliteNames.at( j ), "" );
                    linkEnds[ transmitter ] = std::make_pair( satelliteNames.at( i ), "" );
                    stationReceiverLinkEnds.push_back( linkEnds );
                    numberOfISLLinks += 1;
                }
            }
        }
    }

    // Define (arbitrarily) link ends to be used for 1-way range, 1-way doppler and angular position observables
    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;
    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservableGNSS;
    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservableISL;
    for( unsigned int i = 0; i < numberOfSatellitesPEC; i++ )
    {
        for( unsigned int j = 0; j < numberOfSatellitesGNSS; j++ )
        {
            linkEndsPerObservable[ one_way_range ].push_back( stationTransmitterLinkEnds[ j + i*numberOfSatellitesGNSS] );
    //        linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ j + i*numberOfSatellitesGNSS] );
//            linkEndsPerObservable[ one_way_doppler ].push_back( stationTransmitterLinkEnds[ j + i*numberOfSatellitesGNSS] );

            linkEndsPerObservableGNSS[ one_way_range ].push_back( stationTransmitterLinkEnds[ j + i*numberOfSatellitesGNSS] );
//            linkEndsPerObservableGNSS[ one_way_doppler ].push_back( stationTransmitterLinkEnds[ j + i*numberOfSatellitesGNSS] );
        }
    }

    // Add intersatellite tracking to the link ends
    if( interSatelliteTracking )
    {
        for( unsigned int i = 0; i < numberOfISLLinks; i++)
        {
            linkEndsPerObservable[ one_way_range ].push_back( stationTransmitterLinkEnds[ i + numberOfGNSSLinks] );
    //        linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ j + i*(numberOfSatellitesPEC-1) + numberOfGNSSLinks] );

            linkEndsPerObservableISL[ one_way_range ].push_back( stationTransmitterLinkEnds[ i + numberOfGNSSLinks] );
        }
    }


    // // SETTINGS
    ObservationSettingsMap observationSettingsMap;
    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;

        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            // Define settings for observable, no light-time corrections, and biases for selected 1-way range links
            observationSettingsMap.insert( std::make_pair( currentLinkEndsList.at( i ),
                                                           std::make_shared< ObservationSettings >( currentObservable ) ) );
        }
    }


    // // SIMULATE OBSERVATIONS
    std::vector< double > baseTimeList;
    baseTimeList = propagatedTimes;

    // Create measurement simulation input
    std::map< ObservableType, std::map< LinkEnds, std::shared_ptr< ObservationSimulationTimeSettings< double > > > >
            measurementSimulationInput;
    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        // Define observable type and link ends
        ObservableType currentObservable = linkEndIterator->first;
        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;

        // Define observation times and reference link ends
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            measurementSimulationInput[ currentObservable ][ currentLinkEndsList.at( i ) ] =
                    std::make_shared< TabulatedObservationSimulationTimeSettings< double > >( receiver, baseTimeList );
        }
    }

    // Create observation viability settings and calculators
    std::vector< std::shared_ptr< ObservationViabilitySettings > > observationViabilitySettingsGNSS;
    for( unsigned int i = 0; i < numberOfSatellitesPEC; i++ )
    {
        observationViabilitySettingsGNSS.push_back( std::make_shared< ObservationViabilitySettings >(
                                                    body_occultation, std::make_pair( satelliteNames.at( i ), "" ), "Earth", TUDAT_NAN ) );
        observationViabilitySettingsGNSS.push_back( std::make_shared< ObservationViabilitySettings >(
                                                    body_avoidance_angle, std::make_pair( satelliteNames.at( i ), "" ), "Earth", convertDegreesToRadians( 90.0 ) ) ); // PECMEO can only look up
    }

    // Create special observation visibility settings for ISL
    std::vector< std::shared_ptr< ObservationViabilitySettings > > observationViabilitySettingsISL;
    for( unsigned int i = 0; i < numberOfSatellitesPEC; i++ )
    {
        observationViabilitySettingsISL.push_back( std::make_shared< ObservationViabilitySettings >(
                                                    body_occultation, std::make_pair( satelliteNames.at( i ), "" ), "Earth", TUDAT_NAN ) );
        observationViabilitySettingsGNSS.push_back( std::make_shared< ObservationViabilitySettings >(
                                                    body_avoidance_angle, std::make_pair( satelliteNames.at( i ), "" ), "Earth", convertDegreesToRadians( 0.0 ) ) ); // PECMEO ISL can look any direction
    }

    // Calculate visibilities
    PerObservableObservationViabilityCalculatorList viabilityCalculators;
    PerObservableObservationViabilityCalculatorList viabilityCalculatorsGNSS = createObservationViabilityCalculators(
                bodyMap, linkEndsPerObservableGNSS, observationViabilitySettingsGNSS );
    PerObservableObservationViabilityCalculatorList viabilityCalculatorsISL = createObservationViabilityCalculators(
                bodyMap, linkEndsPerObservableISL, observationViabilitySettingsISL );
    viabilityCalculators[ one_way_range ].insert(viabilityCalculatorsGNSS[ one_way_range ].begin(), viabilityCalculatorsGNSS[ one_way_range ].end());
    viabilityCalculators[ one_way_range ].insert(viabilityCalculatorsISL[ one_way_range ].begin(), viabilityCalculatorsISL[ one_way_range ].end());

    // Simulate observations
    std::map< ObservableType, std::shared_ptr< ObservationSimulatorBase< double, double > > > observationSimulators =
    createObservationSimulators( observationSettingsMap , bodyMap );

    PodInputDataType observationsAndTimes = simulateObservations< double, double >(
                measurementSimulationInput, observationSimulators, viabilityCalculators );

    return observationsAndTimes;
}


////////////////////////// PROCESSING OF SIMULATION RESULTS //////////////////////////

// OLD METHOD FOR MEASUREMENT SIMULATION
//! Obtain measurement data for visible GNSS satellites
//! for all epochs per user
std::map< int, std::map< double, Eigen::VectorXd > > getMeasurementData(
        std::vector< double > epochTimes,
        std::map< double, Eigen::Vector3d > userSatellitePositionPerEpoch,
        std::map< double, Eigen::MatrixXd > positionMatrixGNSSSatellitesPerEpoch)
{
    double frequency = 1575.42e6;
    double lightConstant = physical_constants::getSpeedOfLight< double >();
    double waveLength = lightConstant/frequency;

    std::map< int, std::map< double, Eigen::VectorXd > > measurementData;
    unsigned long long numberOfEpochs = epochTimes.size();
    for( unsigned long long epoch = 0; epoch < numberOfEpochs; epoch++)
    {
        double epochTime = epochTimes.at( epoch );
        Eigen::Vector3d userSatellitePosition = userSatellitePositionPerEpoch[ epochTime ];
        Eigen::MatrixXd positionMatrixGNSSSatellites = positionMatrixGNSSSatellitesPerEpoch[ epochTime ];

        long long numberOfNavigationSatellites = positionMatrixGNSSSatellites.rows();
        Eigen::VectorXd phaseData( numberOfNavigationSatellites );
        Eigen::VectorXd rangeData( numberOfNavigationSatellites );

        Eigen::RowVector3d gnssVector;
        Eigen::RowVector3d pecmeoVector = userSatellitePosition.segment( 0,3 );
        Eigen::RowVector3d positionDifference;

        for( int i = 0; i < numberOfNavigationSatellites; i++ )
        {
            gnssVector = positionMatrixGNSSSatellites.block( i, 0, 1, 3 );
            positionDifference = pecmeoVector - gnssVector;
            double distance = positionDifference.norm();
            double signalTravelTime = distance / lightConstant;
            double phaseRe0 = 0;
            double phaseTr0 = 1e7;
            double phaseRe = phaseRe0;
            double phaseTr = phaseTr0 - frequency * signalTravelTime;
            phaseData( i ) = (phaseRe - phaseTr)*waveLength;
            rangeData( i ) = distance;
//            if(epochTime == 30300.0){std::cout << "dist, sig, tr, data: " << distance << ", " << signalTravelTime << ", " << phaseTr0 - frequency * signalTravelTime << ", " << (phaseRe - phaseTr)*waveLength << std::endl;}
        }
        measurementData[ 1 ][ epochTime ] = rangeData;
        measurementData[ 2 ][ epochTime ] = phaseData;        
    }
    return measurementData;
}


//! Obtain all data for visible GNSS satellites for all epochs and one user
//! over all epochs per user
void getMeasurementsAllUsers(
        unsigned int numberOfSatellitesPEC,
        std::vector< std::string > satelliteNames,
        PodInputDataType observationsAndTimes, std::vector< double > epochTimesTrue,
        std::map< std::string, std::shared_ptr< Ephemeris > > satelliteEphemeris,
        bool lighttimeCorrection, bool addNoise,
        double stddevCode, double stddevPhase,
        double stddevCodeISL, double stddevPhaseISL,
        std::map< std::string, Eigen::MatrixXi >& availabilityMatrixAllUsers,
        std::map< std::string, std::map< double, std::vector< int > > >& visibleSatelliteIDsPerEpochPerUser,
        std::map< std::string, std::map< double, std::map< int, double > > >& correctedRangesPerUserPerEpoch,
        std::map< std::string, std::map< double, std::map< int, double > > >& correctedPhasesPerUserPerEpoch,
        std::map< std::string, std::map< double, std::map< int, double > > >& travelTimesPerUserPerEpoch)
{
    unsigned long long numberOfEpochs = 0;
    std::vector< double > epochTimes;
    double lightConstant = physical_constants::getSpeedOfLight< double >();

    // Adding identifier to satellites for indexing
    std::vector< std::string > satelliteNamesGNSS = satelliteNames;
    std::rotate( satelliteNamesGNSS.begin(), satelliteNamesGNSS.begin() + numberOfSatellitesPEC, satelliteNamesGNSS.end() );
    unsigned long long numberOfSatellites = satelliteNamesGNSS.size();
    std::map< std::string, int > satelliteIDs;
    int satelliteID = 0;
    for( unsigned long long j = 0; j < numberOfSatellites; j++)
    {
        satelliteIDs[ satelliteNamesGNSS.at(j) ] = satelliteID;
//        std::cout << "Sat name: " << satelliteNamesGNSS.at(j) << " with ID: " << satelliteID << std::endl;
        satelliteID += 1;
    }
    // Set index for ISL
    int firstPECID = numberOfSatellites - numberOfSatellitesPEC;

    // Opening the simulated observations object and restructure the data
    std::map< double, double > completeEpochList;
    std::map< std::string, int > satellitesUsedForNavigation;
    for ( PodInputDataType::iterator itrObservationsAndTimesMap = observationsAndTimes.begin( ) ; itrObservationsAndTimesMap != observationsAndTimes.end( ) ;
          itrObservationsAndTimesMap++ )
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

            ObservationVectorType observationVectorType = itrSingleObservableMap->second.first; // Gives observationVectorType which is a list of the measured value for al epochs
//            std::cout << "Current link is between: " << satLinkNames.at(0) << " & " << satLinkNames.at(1) << " satellites and they link from epoch " << epochTimes.front() << " to " << epochTimes.back() << ", " << numberOfEpochs << " times" << std::endl;

            // Loop over the epochs and create list of visible satellites, range and phase data and travel times per epoch
            for ( unsigned int i = 0 ; i < numberOfEpochs ; i++ )
            {
//                std::cout << "Current link is between: " << satLinkNames.at(0) << " & " << satLinkNames.at(1) << " satellites and they link at " << epochTimes[ i ] << std::endl;
                visibleSatelliteIDsPerEpochPerUser[ satLinkNames.at(1) ][ epochTimes[ i ] ].push_back( currentSatelliteID );
                correctedRangesPerUserPerEpoch[ satLinkNames.at(1) ][ epochTimes[ i ] ][ currentSatelliteID ] = observationVectorType(i);
                correctedPhasesPerUserPerEpoch[ satLinkNames.at(1) ][ epochTimes[ i ] ][ currentSatelliteID ] = observationVectorType(i);
                travelTimesPerUserPerEpoch[ satLinkNames.at(1) ][ epochTimes[ i ] ][ currentSatelliteID ] = observationVectorType(i)/lightConstant;
                completeEpochList[ epochTimes[ i ] ] = epochTimes[ i ];
                satellitesUsedForNavigation[ satLinkNames.at(0) ] = currentSatelliteID;
            }
        }
    }

    // Creating new epoch times list, which includes all of the epochs.
    epochTimes.clear();
    epochTimes = epochTimesTrue;
    numberOfEpochs = epochTimes.size();
    unsigned long long numberOfSatellitesGNSS = satellitesUsedForNavigation.size();

    // Loop over different users
    srand(456);
    std::default_random_engine generator(123);
    for( unsigned long long i = 0; i < numberOfSatellitesPEC; i++ )
    {
        // Current user
//        std::cout << i << std::endl;
        std::string userSatellite = satelliteNames.at( i );
//        std::cout << "\n" << i << ", " << userSatellite << std::endl;

        // Create availability matrix for current user
        Eigen::MatrixXi availabilityMatrix = Eigen::MatrixXi::Zero( numberOfEpochs, numberOfSatellitesGNSS );

        // Loop over all epochs to create availability matrix
        for( unsigned int i = 0; i < numberOfEpochs; i++ )
        {
            double epochTime = epochTimes.at( i );
            std::vector< int > visibleSatelliteIDs = visibleSatelliteIDsPerEpochPerUser[ userSatellite ][ epochTime ];
            unsigned long long numberOfVisibleSatellites = visibleSatelliteIDs.size();
            for( unsigned long long j = 0; j < numberOfVisibleSatellites; j++ )
            {
                int currentSatelliteID = visibleSatelliteIDs.at(j);
                availabilityMatrix( i, currentSatelliteID ) = 1;
            }
        }
        availabilityMatrixAllUsers[ userSatellite ] = availabilityMatrix;

        // Give errors and deviations
//        std::cout << "stddevCode: " << stddevCode << " & stddevPhase: " << stddevPhase << std::endl;
        std::normal_distribution<double> distributionCode(0.0, stddevCode);
        std::normal_distribution<double> distributionPhase(0.0, stddevPhase);
        std::normal_distribution<double> distributionCodeISL(0.0, stddevCodeISL);
        std::normal_distribution<double> distributionPhaseISL(0.0, stddevPhaseISL);

        std::map< int, int > currentCarrierAmbiguity;
//        double clockOffSet = (i+1); // Template for future option to simulate clock errors
        for( unsigned int i = 0; i < numberOfEpochs; i++ )
        {
            double epochTime = epochTimes.at( i );
            Eigen::Vector3d userSatellitePosition = satelliteEphemeris[ userSatellite ]->getCartesianState( epochTime ).segment( 0, 3 );
            std::vector< int > visibleSatelliteIDs = visibleSatelliteIDsPerEpochPerUser[ userSatellite ][ epochTime ];
            unsigned long long numberOfVisibleSatellites = visibleSatelliteIDs.size();
            for( unsigned long long j = 0; j < numberOfVisibleSatellites; j++ )
            {
                int currentSatelliteID = visibleSatelliteIDs.at(j);

                // Check if lighttime correction is used
                if(!lighttimeCorrection)
                {
                    std::string currentSatellite = satelliteNamesGNSS.at( static_cast<unsigned>(currentSatelliteID) );
//                    std::cout << "\n" << userSatellite << " and " << currentSatellite << " at " << epochTime << " and travel time: " << travelTimesPerUserPerEpoch[userSatellite][epochTime][currentSatelliteID] << std::endl;
                    Eigen::Vector3d currentSatellitePosition = satelliteEphemeris[ currentSatellite ]->getCartesianState( epochTime ).segment( 0, 3 );
                    correctedRangesPerUserPerEpoch[ userSatellite ][ epochTime ][ currentSatelliteID ] = (userSatellitePosition - currentSatellitePosition).norm();
                    correctedPhasesPerUserPerEpoch[ userSatellite ][ epochTime ][ currentSatelliteID ] = (userSatellitePosition - currentSatellitePosition).norm();
                }

                // Adjust code measurement
                double codeError = 0;
                if(addNoise){ codeError = distributionCode( generator ); }
                if( addNoise && currentSatelliteID >= firstPECID ){ codeError = distributionCodeISL( generator ); /*std::cout << codeError << ", ";*/ }
                correctedRangesPerUserPerEpoch[ userSatellite ][ epochTime ][ currentSatelliteID ] += codeError;
//                correctedRangesPerUserPerEpoch[ userSatellite ][ epochTime ][ currentSatelliteID ] += codeError + clockOffSet; // For future option to simulate clock errirs

                // Adjust phase measurement
                double phaseError = 0;
                if(addNoise){ phaseError = distributionPhase( generator ); }
                if( addNoise && currentSatelliteID >= firstPECID ){ phaseError = distributionPhaseISL( generator ); }
                double zeroPhaseDeviation = 0; // waveLength*(receiverZeroPhase - transmitterZeroPhase)
                if( i == 0 || availabilityMatrix( i-1, currentSatelliteID ) == 0 ) {
                    currentCarrierAmbiguity[ currentSatelliteID ] = rand() % 20000 + -10000;
//                    std::cout << currentCarrierAmbiguity[ currentSatelliteID ] << ", ";
                }
                double carrierAmbiguity = currentCarrierAmbiguity[ currentSatelliteID ]; //double carrierAmbiguity = 2e6;
                correctedPhasesPerUserPerEpoch[ userSatellite ][ epochTime ][ currentSatelliteID ] += ( phaseError + zeroPhaseDeviation - carrierAmbiguity );
            }
        }
    }
}
