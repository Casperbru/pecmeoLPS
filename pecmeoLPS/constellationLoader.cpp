#include "constellationLoader.h"
#include <memory>
#include <Eigen/Core>
#include <boost/make_shared.hpp>
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "toolsSimulation.h"

using namespace tudat;
using namespace tudat::unit_conversions;
using namespace tudat::orbital_element_conversions;

//! Obtain constellation characteristics (Galileo, GLONASS, BeiDou, GPS)
std::shared_ptr< constellationLoader > getConstellation( constellationNames constellationName )
{
    double nummer = 0.0;
    int numberOfSatellites;
    int numberOfPlanes;
    int numberOfSatellitesPerPlane;
    int sizeOfState;
    Eigen::Vector4i constellationConfig;
    Eigen::VectorXd orbitalElements(6);
    std::string constellationTitle;

    // Select body with prefined central gravity field.
    switch( constellationName )
    {
    case Galileo:
        // Set number of satellites in constellation.
        numberOfSatellites = 30;
        // Set number of planes in constellation.
        numberOfPlanes = 3;
        // Set number of satellites per plane in constellation.
        numberOfSatellitesPerPlane = numberOfSatellites / numberOfPlanes;
        // Set size of state
        sizeOfState = 6;
        // Vector with orbit configuration (number of satellites, number of planes, number of satellites per plane, size of state)
        constellationConfig << numberOfSatellites, numberOfPlanes, numberOfSatellitesPerPlane, sizeOfState;
        // Vector with orbital elements (SMA, Ecc, Inc, AoP, RAAN, TA)
        orbitalElements << 23222.0e3 + 6378.1e3, 0.0, convertDegreesToRadians( 56.0 ), 0.0, 2.0 * mathematical_constants::PI / numberOfPlanes, 2.0 * mathematical_constants::PI / numberOfSatellitesPerPlane;
        // String with the title of the constellation
        constellationTitle = "Galileo";

        break;

    case GLONASS:
        // Set number of satellites in constellation.
        numberOfSatellites = 24;
        // Set number of planes in constellation.
        numberOfPlanes = 3;
        // Set number of satellites per plane in constellation.
        numberOfSatellitesPerPlane = numberOfSatellites / numberOfPlanes;
        // Set size of state
        sizeOfState = 6;
        // Vector with orbit configuration (number of satellites, number of planes, number of satellites per plane, size of state)
        constellationConfig << numberOfSatellites, numberOfPlanes, numberOfSatellitesPerPlane, sizeOfState;
        // Vector with orbital elements (SMA, Ecc, Inc, AoP, RAAN, TA)
        orbitalElements << 19100.0e3 + 6378.1e3, 0.0, unit_conversions::convertDegreesToRadians( 64.8 ), 0.0, 2.0 * mathematical_constants::PI / numberOfPlanes, 2.0 * mathematical_constants::PI / numberOfSatellitesPerPlane;
        // String with the title of the constellation
        constellationTitle = "GLONASS";

        break;

    case BeiDou: // Note that for now only the MEO satellites are considered (there are also 3 IGSO and 5 GEO satellites)
        // Set number of satellites in constellation.
        numberOfSatellites = 27;
        // Set number of planes in constellation.
        numberOfPlanes = 3;
        // Set number of satellites per plane in constellation.
        numberOfSatellitesPerPlane = numberOfSatellites / numberOfPlanes;
        // Set size of state
        sizeOfState = 6;
        // Vector with orbit configuration (number of satellites, number of planes, number of satellites per plane, size of state)
        constellationConfig << numberOfSatellites, numberOfPlanes, numberOfSatellitesPerPlane, sizeOfState;
        // Vector with orbital elements (SMA, Ecc, Inc, AoP, RAAN, TA)
        orbitalElements << 21500.0e3 + 6378.1e3, 0.0, unit_conversions::convertDegreesToRadians( 55.8 ), 0.0, 2.0 * mathematical_constants::PI / numberOfPlanes, 2.0 * mathematical_constants::PI / numberOfSatellitesPerPlane;
        // String with the title of the constellation
        constellationTitle = "BeiDou";

        break;

    case GPS: // Note that the old constellation is used and not the current one with 30 satellites
        // Set number of satellites in constellation.
        numberOfSatellites = 24;
        // Set number of planes in constellation.
        numberOfPlanes = 6;
        // Set number of satellites per plane in constellation.
        numberOfSatellitesPerPlane = numberOfSatellites / numberOfPlanes;
        // Set size of state
        sizeOfState = 6;
        // Vector with orbit configuration (number of satellites, number of planes, number of satellites per plane, size of state)
        constellationConfig << numberOfSatellites, numberOfPlanes, numberOfSatellitesPerPlane, sizeOfState;
        // Vector with orbital elements (SMA, Ecc, Inc, AoP, RAAN, TA)
        orbitalElements << 20200.0e3 + 6378.1e3, 0.0, unit_conversions::convertDegreesToRadians( 55.0 ), 0.0, 2.0 * mathematical_constants::PI / numberOfPlanes, 2.0 * mathematical_constants::PI / numberOfSatellitesPerPlane;
        // String with the title of the constellation
        constellationTitle = "GPS";

        break;


    default:

        std::string errorMessage = "Desired predefined central gravity field " +
                std::to_string( constellationName ) + " does not exist";
        throw std::runtime_error( errorMessage );
    }
    return std::make_shared< constellationLoader >( nummer, constellationConfig, orbitalElements, constellationTitle );
}


//! Set constellations: Galileo, GLONASS, BeiDou, GPS
void setConstellation(
        unsigned int numberOfSatellitesPEC,
        Eigen::MatrixXd initialConditionsInKeplerianElementsPEC,
        std::vector <constellationNames> constellationList,
        std::vector< std::string >& satelliteNames,
        Eigen::MatrixXd& initialConditionsInKeplerianElementsTotal,
        unsigned int& numberOfSatellitesTotal,
        unsigned int& numberOfConstellations,
        std::string& includedNames)
{
    numberOfConstellations = static_cast< unsigned >( constellationList.size() );
    initialConditionsInKeplerianElementsTotal = initialConditionsInKeplerianElementsPEC;
    numberOfSatellitesTotal = numberOfSatellitesPEC;
    std::string currentName;

    // Load all the different constellations and set initial state of satellites
    for (unsigned int i = 0; i < numberOfConstellations; i++)
    {
        Eigen::VectorXd orbitalElements( 6 ) ;
        Eigen::Vector4i constellationConfig;
        Eigen::MatrixXd initialConditionsInKeplerianElements;
        Eigen::MatrixXd initialConditionsInKeplerianElementsTemp = initialConditionsInKeplerianElementsTotal;

        std::shared_ptr< constellationLoader > constellation = getConstellation( constellationList[ i ] );
        orbitalElements = constellation->getOrbitalElements();
        constellationConfig = constellation->getConstellationConfig();
        currentName = constellation->getConstellationTitle();
        includedNames.append( currentName );

        std::cout << currentName << " Constellation loaded!" <<std::endl;

        numberOfSatellitesTotal += static_cast< unsigned >(constellationConfig(0));

        initialConditionsInKeplerianElements = setInitialConds( orbitalElements, constellationConfig );

        if (currentName == "GPS"){ // Set different true anomaly spacing.

            Eigen::RowVectorXd trueAnomalySpacingIntegers( constellationConfig(2) );

            trueAnomalySpacingIntegers << 0.0, 30.0, 135.0, 255.0;
            trueAnomalySpacingIntegers = convertDegreesToRadians(trueAnomalySpacingIntegers);

            for ( unsigned int i = 0; i < static_cast< unsigned >(constellationConfig(1)); i++ )
            {
                initialConditionsInKeplerianElements.block( 5, i * static_cast< unsigned >(constellationConfig(2)), 1, constellationConfig(2) ) =
                        Eigen::MatrixXd::Constant( 1, constellationConfig(2), 1 ).array( ) * trueAnomalySpacingIntegers.array( );
            }
        }

        // Add initial conditions to complete matrix
        initialConditionsInKeplerianElementsTotal.resize(initialConditionsInKeplerianElementsTemp.rows(), initialConditionsInKeplerianElementsTemp.cols() + initialConditionsInKeplerianElements.cols() );
        initialConditionsInKeplerianElementsTotal << initialConditionsInKeplerianElementsTemp, initialConditionsInKeplerianElements;

        // Set satellite names
        for ( int i = 0; i < constellationConfig(0); i++)
        {
            std::string satDesignation;
            if( i < 9 )
            {satDesignation = "Sat0";}
            else
            {satDesignation = "Sat";}

            satelliteNames.push_back( currentName+satDesignation+std::to_string(i+1) );
        }
    }
}


//! Set PECMEO constellation initial conditions
void setPECMEO(
        const unsigned int numberOfSatellitesPEC,
        const unsigned int numberOfPlanesPEC,
        std::vector< double > keplerElements,
        std::vector< std::string >& satelliteNames,
        Eigen::MatrixXd& initialConditionsInKeplerianElementsPEC)
{
    // Set number of satellites per plane
    const unsigned int numberOfSatellitesPerPlanePEC = numberOfSatellitesPEC / numberOfPlanesPEC;

    // Split the Keplerian elements
    const double semiMajorAxisPEC = keplerElements.at(0);                    // [km]
    const double eccentricityPEC = keplerElements.at(1);                     // [-]
    const double inclinationPEC = keplerElements.at(2);                      // [rad]
    const double argumentOfPeriapsisPEC = keplerElements.at(3);              // [rad]
    const double longitudeOfAscendingNodeSpacingPEC = keplerElements.at(4);  // [rad]
//    const double trueAnomalySpacingPEC = keplerElements.at(5);               // [rad]
    std::vector< double > argumentOfLattitudeShift = { keplerElements.at(5), keplerElements.at(6), keplerElements.at(7) };
//    const double argumentOfLattitudeShift1 = keplerElements.at(5);               // [rad]
//    const double argumentOfLattitudeShift2 = keplerElements.at(6);               // [rad]
//    const double argumentOfLattitudeShift3 = keplerElements.at(7);               // [rad]

    const double trueAnomalySpacingPEC = (2.0 * mathematical_constants::PI / numberOfSatellitesPerPlanePEC);               // [rad]

    // Set orbital elements except true anomaly
    for (unsigned int i = 0; i < numberOfSatellitesPEC; i++)
    {
        double satNum = static_cast<double>(i+1);
        double satPartial = satNum/static_cast<double>(numberOfSatellitesPEC);

        initialConditionsInKeplerianElementsPEC(0, i) = semiMajorAxisPEC;
        initialConditionsInKeplerianElementsPEC(1, i) = eccentricityPEC;
        initialConditionsInKeplerianElementsPEC(2, i) = inclinationPEC;

        if (satPartial > 1.0/3.0)
        {
            initialConditionsInKeplerianElementsPEC(2, i) += convertDegreesToRadians( 90.0 );
        }
        initialConditionsInKeplerianElementsPEC(3, i) = argumentOfPeriapsisPEC;
        initialConditionsInKeplerianElementsPEC(4, i) = longitudeOfAscendingNodeSpacingPEC;

        if (satPartial > 2.0/3.0)
        {
            initialConditionsInKeplerianElementsPEC(4, i) += convertDegreesToRadians( 90.0 );
        }
    }

    // Set true anomaly
    Eigen::RowVectorXd trueAnomalySpacingIntegersPEC( numberOfSatellitesPerPlanePEC );

    for ( unsigned int i = 0; i < numberOfSatellitesPerPlanePEC; i++ )
    {
        trueAnomalySpacingIntegersPEC( i ) =  i * 1.0;
    }

    for ( unsigned int i = 0; i < numberOfPlanesPEC; i++ )
    {
        initialConditionsInKeplerianElementsPEC.block( 5, i * numberOfSatellitesPerPlanePEC, 1, numberOfSatellitesPerPlanePEC ) =
//                Eigen::MatrixXd::Constant( 1, numberOfSatellitesPerPlanePEC, trueAnomalySpacingPEC ).array( ) * trueAnomalySpacingIntegersPEC.array( );// + 0.0001*mathematical_constants::PI; // Manual offset so 1 and 4 wont "collide"
                Eigen::MatrixXd::Constant( 1, numberOfSatellitesPerPlanePEC, trueAnomalySpacingPEC ).array( ) * trueAnomalySpacingIntegersPEC.array( ) + argumentOfLattitudeShift.at(i); // Offset in argument of lattitude shift
    }

    // Set satellite names
    for (unsigned int i = 0; i <numberOfSatellitesPEC; i++)
    {
        satelliteNames.push_back( "PECMEOSat"+std::to_string(i+1) );
    }
}


//! Obtain the initial conditions matrix of costellation
Eigen::MatrixXd setInitialConds(
        Eigen::VectorXd orbitalElements,
        Eigen::Vector4i constellationConfig)
{
    // Split the constellation configuration in integers
    // With numberOfSatellites, numberOfPlanes, numberOfSatellitesPerPlane, sizeOfState
    int nos = constellationConfig(0);
    int nop = constellationConfig(1);
    int nospp = constellationConfig(2);
    int sos = constellationConfig(3);

    // Split the orbital elements vector in doubles
    double sma = orbitalElements(0);
    double ecc = orbitalElements(1);
    double inc = orbitalElements(2);
    double aop = orbitalElements(3);
    double raans = orbitalElements(4);
    double tas = orbitalElements(5);

    // Set Keplerian elements for Galileo satellites.
    Eigen::MatrixXd initialConditionsInKeplerianElements;
    initialConditionsInKeplerianElements.resize( sos, nos );

    // Set semiMajorAxis.
    initialConditionsInKeplerianElements.row( 0 ) = Eigen::MatrixXd::Constant( 1, nos, sma );

    // Set eccentricity.
    initialConditionsInKeplerianElements.row( 1 ) = Eigen::MatrixXd::Constant( 1, nos, ecc );

    // Set inclination.
    initialConditionsInKeplerianElements.row( 2 ) = Eigen::MatrixXd::Constant( 1, nos, inc );

    // Set argument of periapsis.
    initialConditionsInKeplerianElements.row( 3 ) = Eigen::MatrixXd::Constant( 1, nos, aop );

    // Set longitude of ascending node.
    for ( int i = 0; i < nop; i++ )
    {
        initialConditionsInKeplerianElements.block( 4, i * nospp, 1, nospp ) =
                Eigen::MatrixXd::Constant( 1, nospp, i * raans );
    }

    // Set true anomaly.
    Eigen::RowVectorXd trueAnomalySpacingIntegers( nospp );

    for ( int i = 0; i < nospp; i++ )
    {
        trueAnomalySpacingIntegers( i ) =  i * 1.0;
    }

    for ( int i = 0; i < nop; i++ )
    {
        initialConditionsInKeplerianElements.block( 5, i * nospp, 1, nospp ) =
                Eigen::MatrixXd::Constant( 1, nospp, tas ).array( ) * trueAnomalySpacingIntegers.array( );
    }

    return initialConditionsInKeplerianElements;
};


//! Get all the constellations (GNSS & PECMEO), initialize them in the body map and get initial state in vector form
void setAllConstellations(
        unsigned int numberOfSatellitesPEC,
        unsigned int numberOfPlanesPEC,
        std::vector< double > keplerElements,
        std::vector <constellationNames> constellationList,
        simulation_setup::NamedBodyMap& bodyMap,
        unsigned int& numberOfSatellitesTotal,
        std::vector< std::string >& satelliteNames,
        unsigned int& numberOfConstellations,
        std::string& includedNames,
        Eigen::VectorXd& systemInitialState)
{
    // Set initial conditions of PECMEO constellation.
    Eigen::MatrixXd initialConditionsInKeplerianElementsPEC( 6, numberOfSatellitesPEC );
    setPECMEO( numberOfSatellitesPEC, numberOfPlanesPEC, keplerElements, satelliteNames, initialConditionsInKeplerianElementsPEC );

    // Set initial conditions of GNSS constellation.
    Eigen::MatrixXd initialConditionsInKeplerianElementsTotal;
    setConstellation(numberOfSatellitesPEC, initialConditionsInKeplerianElementsPEC, constellationList,
                                     satelliteNames, initialConditionsInKeplerianElementsTotal, numberOfSatellitesTotal, numberOfConstellations, includedNames );

    // Set accelerations for each satellite and add to body map.
    finalizeBodyMap( numberOfSatellitesTotal, satelliteNames, initialConditionsInKeplerianElementsTotal, bodyMap );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    // Set initial state in vector form
    systemInitialState = getInitialStateVector( bodyMap, numberOfSatellitesTotal, initialConditionsInKeplerianElementsTotal );

}
