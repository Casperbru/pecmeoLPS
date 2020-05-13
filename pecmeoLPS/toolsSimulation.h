#ifndef TOOLSSIMULATION_H
#define TOOLSSIMULATION_H
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <Tudat/SimulationSetup/tudatEstimationHeader.h>
#include <memory>
#include <Eigen/Core>

using namespace tudat;
using namespace tudat::simulation_setup;
using namespace tudat::observation_models;

using namespace tudat::ephemerides;


//! Obtain the body map for the set simulation time
void getBodyMap(
        const double simulationStartEpoch,
        const double simulationEndEpoch,
        const double fixedStepSize,
        double& stepsInSimulation,
        simulation_setup::NamedBodyMap& bodyMap);


//! Finalize body map by adding the satellites
void finalizeBodyMap(
        unsigned int numberOfSatellitesTotal,
        std::vector< std::string > satelliteNames,
        Eigen::MatrixXd initialConditionsInKeplerianElementsTotal,
        simulation_setup::NamedBodyMap& bodyMap);


//! Get initial Cartesian state of all satellites in one vector
Eigen::VectorXd getInitialStateVector(
        simulation_setup::NamedBodyMap bodyMap,
        unsigned int numberOfSatellitesTotal,
        Eigen::MatrixXd initialConditionsInKeplerianElementsTotal);


//! Obtain the ephemeris of the visible GNSS satellites
std::map< std::string, std::shared_ptr< Ephemeris > > getEphemeridesOfSatellites(
        std::vector< std::string > satelliteNames,
        simulation_setup::NamedBodyMap bodyMap);


//! Obtain ephemeris of GNSS satellites indexed by numbers
std::vector< std::shared_ptr< Ephemeris > > getEphemeridesOfGNSS(
        std::vector< std::string > satelliteNames,
        simulation_setup::NamedBodyMap bodyMap,
        unsigned int numberOfSatellitesPEC);


//! Obtain the epochs over which the simulation is performed.
std::vector< double > getEpochTimes(
        const double simulationStartEpoch,
        const double simulationEndEpoch,
        const double fixedStepSize);


//! Obtain the correct propagation history of the satellites.
void getStateHistory(
        unsigned int numberOfSatellitesTotal,
        std::vector< std::string > satelliteNames,
        simulation_setup::NamedBodyMap bodyMap,
        std::vector< double > epochTimes,
        std::vector< std::map< double, Eigen::VectorXd > >& allSatellitesPropagationHistory,
        std::vector< std::map< double, Eigen::VectorXd > >& allSatellitesPropagationHistoryKep);


// Set typedefs for POD input (observation types, observation link ends, observation values, associated times with reference link ends.
typedef Eigen::Matrix< double, Eigen::Dynamic, 1 > ObservationVectorType;
typedef std::map< LinkEnds, std::pair< ObservationVectorType, std::pair< std::vector< double >, LinkEndType > > >
        SingleObservablePodInputType;
typedef std::map< ObservableType, SingleObservablePodInputType > PodInputDataType;

//! Determine all the visible observations between the link ends over time
PodInputDataType observationsAndTimesCalculator(
        unsigned int numberOfSatellitesTotal,
        unsigned int numberOfSatellitesPEC,
        std::vector< std::string >  satelliteNames,
        std::vector<double> propagatedTimes,
        simulation_setup::NamedBodyMap bodyMap,
        bool interSatelliteTracking);


//! Obtain measurement data for visible GNSS satellites
//! for all epochs per user
std::map< int, std::map< double, Eigen::VectorXd > > getMeasurementData(
        std::vector< double > epochTimes,
        std::map< double, Eigen::Vector3d > userSatellitePositionPerEpoch,
        std::map< double, Eigen::MatrixXd > positionMatrixGNSSSatellitesPerEpoch);


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
        std::map< std::string, std::map< double, std::map< int, double > > >& travelTimesPerUserPerEpoch);


#endif // TOOLSSIMULATION_H
