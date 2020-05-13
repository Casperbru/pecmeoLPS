#ifndef CONSTELLATIONLOADER_H
#define CONSTELLATIONLOADER_H
#include <memory>
#include <Eigen/Core>
#include <math.h>
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

using namespace tudat;
using namespace tudat::simulation_setup;

//! Load constellations: Galileo, GLONASS, BeiDou, GPS
class constellationLoader
{
public:
    constellationLoader( double nummer, Eigen::Vector4i constellationConfig, Eigen::VectorXd orbitalElements, std::string constellationTitle ) :
        nummer_( nummer ), constellationConfig_( constellationConfig), orbitalElements_( orbitalElements), constellationTitle_(constellationTitle)   {}

    double getNummer( )
    {
        return nummer_;
    }

    //! Returns constellation configuration vector (number of satellites, number of planes, number of satellites per plane, size of state
    Eigen::Vector4i getConstellationConfig( )
    {
        return constellationConfig_;
    }

    //! Returns orbital elements vector (SMA, Ecc, Inc, AoP, RAAN, TA)
    Eigen::VectorXd getOrbitalElements( )
    {
        return orbitalElements_;
    }

    //! Returns constellations title
    std::string getConstellationTitle( )
    {
        return constellationTitle_;
    }

protected:
    double nummer_;
    unsigned int numberOfSatellites_;
    unsigned int numberOfPlanes_;
    unsigned int numberOfSatellitesPerPlane_;
    unsigned int sizeOfState_;
    Eigen::Vector4i constellationConfig_;
    Eigen::VectorXd orbitalElements_;
    std::string constellationTitle_;

private:
};


enum constellationNames
{
    Galileo, GLONASS, BeiDou, GPS
};

typedef std::shared_ptr< constellationLoader > constellationLoaderPointer;

std::shared_ptr< constellationLoader > getConstellation( constellationNames constellationName );


//! Set constellations: Galileo, GLONASS, BeiDou, GPS
void setConstellation(
        unsigned int numberOfSatellitesPEC,
        Eigen::MatrixXd initialConditionsInKeplerianElementsPEC,
        std::vector <constellationNames> constellationList,
        std::vector< std::string >& satelliteNames,
        Eigen::MatrixXd& initialConditionsInKeplerianElementsTotal,
        unsigned int& numberOfSatellitesTotal,
        unsigned int& numberOfConstellations,
        std::string& includedNames);


//! Set PECMEO constellation initial conditions
void setPECMEO(const unsigned int numberOfSatellitesPEC,
        const unsigned int numberOfPlanesPEC,
        std::vector< double > keplerElements,
        std::vector< std::string >& satelliteNames,
        Eigen::MatrixXd& initialConditionsInKeplerianElementsPEC);


//! Obtain the initial conditions matrix of costellation
Eigen::MatrixXd setInitialConds(
        Eigen::VectorXd orbitalElements,
        Eigen::Vector4i constellationConfig);


//! Get all the constellations (GNSS & PECMEO) and get initialize them in the body map
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
        Eigen::VectorXd& systemInitialState);


#endif // CONSTELLATIONLOADER_H
