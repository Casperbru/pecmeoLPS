#ifndef TOOLSOPTIMISATION_H
#define TOOLSOPTIMISATION_H

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <Tudat/SimulationSetup/tudatEstimationHeader.h>
#include "constellationLoader.h"

using namespace tudat::observation_models;
using namespace tudat::ephemerides;


// Set typedefs for POD input (observation types, observation link ends, observation values, associated times with reference link ends.
typedef Eigen::Matrix< double, Eigen::Dynamic, 1 > ObservationVectorType;
typedef std::map< LinkEnds, std::pair< ObservationVectorType, std::pair< std::vector< double >, LinkEndType > > >
        SingleObservablePodInputType;
typedef std::map< ObservableType, SingleObservablePodInputType > PodInputDataType;


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
        std::vector< std::map< double, Eigen::VectorXd > >& allSatellitesPropagationHistoryKep);


//! Retrieve optimisation parameters for all users
//! Gives vector with average links, GDOP or PDOP over all users and epochs
void retrieveOptimisationParameters(
        PodInputDataType observationsAndTimes,
        std::vector< double > epochTimes,
        simulation_setup::NamedBodyMap bodyMap,
        std::vector< std::string > satelliteNames,
        unsigned int numberOfSatellitesPEC,
        bool lighttimeCorrection, bool addNoise,
        std::map< std::string, std::map< double, Eigen::VectorXd > >& allUsersDOPResults,
        std::vector< double >& results);

#endif // TOOLSOPTIMISATION_H
