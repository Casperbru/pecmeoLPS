#ifndef TOOLSGENERAL_H
#define TOOLSGENERAL_H
#include <memory>
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


//! Perform navigation and estimation process for all users and obtain results for DOP, PP, Kin and the iterations
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
        std::map< int, std::map< std::string, std::map< double, Eigen::Vector4d > > >& estimationsKinPerUserPerIteration);


//! Obtain the current time
void getCurrentTime(
        std::chrono::system_clock::time_point start = std::chrono::system_clock::now(),
        std::chrono::system_clock::time_point end = std::chrono::system_clock::now());


//! Obtain the parameters to be optimised
std::vector< double > obtainParamaterToBeOptimised(
        std::map< std::string, std::map< double, Eigen::VectorXd > > allSatellitesDOPHistory);


//! Retrieve optimisation parameters for all users
//! Gives vector with average links, GDOP or PDOP over all users and epochs
void retrieveOptimisationParameters(
        PodInputDataType observationsAndTimes,
        std::vector< double > epochTimes,
        simulation_setup::NamedBodyMap bodyMap,
        std::vector< std::string > satelliteNames,
        unsigned int numberOfSatellitesPEC,
        std::map< std::string, std::map< double, Eigen::VectorXd > >& allUsersDOPResults,
        std::vector< double >& results);


//! Function to write all the results into files
void writeResults(
        unsigned int numberOfSatellitesTotal,
        std::vector< std::string > satelliteNames,
        std::string outputSubFolder,
        std::vector< std::map< double, Eigen::VectorXd > > allSatellitesPropagationHistory,
        std::vector< std::map< double, Eigen::VectorXd > > allSatellitesPropagationHistoryKep,
        bool writeCartesianResult,
        bool writeKeplerianResult);


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
        bool writeDOPResults = true,
        bool writePPResults = true,
        bool writeKinResults = true,
        bool writePPIterations = true,
        bool writeKinIterations = true,
        bool writeBiasIterations = true);


//! Function to write all the final estimation results into files
void writeEstimationEndResults(
        const unsigned int numberOfSatellitesPEC,
        std::vector< std::string > satelliteNames,
        std::string outputSubFolder,
        int maxIterationsPP, int maxIterationsKin,
        int maxIterationsPPISL, int maxIterationsKinISL,
        std::map< int, std::map< std::string, std::map< double, Eigen::Vector4d > > > estimationsKinPerUserPerIteration);


//! Function to write all the final estimation results into files
void writeEstimationFinalResults(
        const unsigned int numberOfSatellitesPEC,
        std::vector< std::string > satelliteNames,
        std::string outputSubFolder,
        std::map< std::string, std::map< double, Eigen::Vector4d > > estimationsFinal);


//! Function to remove rows from a matrix
void removeRows(Eigen::MatrixXi& matrix, std::vector<int> rowsToRemove);


#endif // TOOLSGENERAL_H
