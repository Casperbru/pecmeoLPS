#ifndef TOOLSNAVIGATION_H
#define TOOLSNAVIGATION_H
#include <Tudat/SimulationSetup/tudatEstimationHeader.h>
#include <memory>
#include <Eigen/Core>

using namespace tudat::observation_models;
using namespace tudat::ephemerides;


// Set typedefs for POD input (observation types, observation link ends, observation values, associated times with reference link ends.
typedef Eigen::Matrix< double, Eigen::Dynamic, 1 > ObservationVectorType;
typedef std::map< LinkEnds, std::pair< ObservationVectorType, std::pair< std::vector< double >, LinkEndType > > >
        SingleObservablePodInputType;
typedef std::map< ObservableType, SingleObservablePodInputType > PodInputDataType;


// /// /// DOP CALCULATIONS
//! Find the point positioning design matrix for given point (user) and GNSS visible
Eigen::MatrixXd getDesignMatrix(
        Eigen::Vector3d userSatellitePosition,
        Eigen::MatrixXd positionMatrixGNSSSatellites);


//! Find the point positioning covariance matrix for given point (user) and GNSS visible
//! per epoch per user to be used for getDOPValues
Eigen::MatrixXd getCovarianceMatrix(
        Eigen::Vector3d userSatellitePosition,
        Eigen::MatrixXd positionMatrixGNSSSatellites);

//! Obtain the different DOP values and the number of links with navigation satellites
Eigen::Vector6d getDOPValues(
        Eigen::Vector3d userSatellitePosition,
        Eigen::MatrixXd positionMatrixGNSSSatellites);


//! Obtain the different DOP values and the number of links with navigation satellites
//! for all epochs per user
std::map< double, Eigen::VectorXd > getDOPValues(
        std::vector< double > epochTimes,
        std::map< double, Eigen::Vector3d > userSatellitePositionAllEpochs,
        std::map< double, Eigen::MatrixXd > positionMatrixGNSSSatellitesAllEpochs);


//! Obtain the state matrix of the visible GNSS satellites through ephemeris
//! for all epochs per user
std::map< double, Eigen::MatrixXd > getPositionMatrixGNSSSatellites(
        Eigen::MatrixXi availabilityMatrix,
        std::vector< std::shared_ptr< Ephemeris > > satelliteEphemerisGNSS,
        std::vector< double > epochTimes);



// ///  /// POSITION ESTIMATION
// ///  SINGLE POINT POSITIONING
//! Perform one single point positioning estimation iteration based on ephemeris of GNSS
//! a priori information is the input and the output are updated estimations using availability matrix
Eigen::VectorXd getPointPositionIteration(
        Eigen::VectorXd aPrioriEstimation,
        Eigen::VectorXd rangeData,
        std::map< std::string, Eigen::RowVectorXi > availabilityRowPerUser,
        std::vector< std::shared_ptr< Ephemeris > > satelliteEphemerisGNSS,
        double epochTime,
        double codeWeight, double codeWeightISL,
        bool ephemerisError, double arcTime, Eigen::MatrixXd ephemerisErrorBase,
        bool lighttimeCorrection = false,
        std::map< std::string, std::map< int, double > > travelTimesPerUser = {{ "0", {{ 0, 0.0 }} }} );


//! Iterate the single point positioning to obtain desired result controls iterations
//! and calls iterator using availability matrix
Eigen::VectorXd getPointPositionEstimation(
        Eigen::VectorXd rangeData,
        std::map< std::string, Eigen::RowVectorXi > availabilityRowPerUser,
        Eigen::MatrixXd& iterationCoordinateHistory,
        std::vector< std::shared_ptr< Ephemeris > > satelliteEphemerisGNSS,
        double epochTime,
        double codeWeight, double codeWeightISL,
        bool ephemerisError, double arcTime, Eigen::MatrixXd ephemerisErrorBase,
        int maxIterations = 20,
        double correctionLimit = 1e-4,
        Eigen::VectorXd aPrioriEstimation = { 0.0, 0.0, 0.0, 0.0 },
        bool lighttimeCorrection = false,
        std::map< std::string, std::map< int, double > > travelTimesPerUser = {{ "0", {{ 0, 0.0 }} }} );


//! Perform point positioning for one user over all epochs
std::map< std::string, std::map< double, Eigen::Vector4d > > getPointPositionSequence(
        std::vector< double > epochTimes,
        std::vector< std::shared_ptr< Ephemeris > > satelliteEphemerisGNSS,
        std::map< std::string, std::map< double, std::map< int, double > > > rangeData,
        std::map< std::string, Eigen::MatrixXi > availabilityMatrixPerUser,
        std::map< std::string, std::map< double, Eigen::VectorXd > >& allUsersPPIterations,
        double codeWeight, double codeWeightISL,
        bool ephemerisError, double arcTime, Eigen::MatrixXd ephemerisErrorBase,
        std::map< std::string, std::map< double, Eigen::Vector4d > > aPrioriEstimationsPerUser = {{ "0", {{ 0, { 0.0, 0.0, 0.0, 0.0 } }} }},
        int maxIterations = 20,
        double correctionLimit = 1e-4,
        bool lighttimeCorrection = false,
        std::map< std::string, std::map< double, std::map< int, double > > > travelTimesPerUserPerEpoch = {{ "0", {{ 0.0, {{ 0, 0.0 }} }} }} );


// /// KINEMATIC ORBIT DETERMINATION
//! Perform a single iteration of kinematic estimation of satellite position, clock errors and bias
//! Older version without block elimination and spars matrices, very slow
Eigen::VectorXd getKinematicIteration(
        std::vector< double > epochTimes,
        std::vector< std::shared_ptr< Ephemeris > > satelliteEphemerisGNSS,
        Eigen::VectorXd aPrioriEstimationsCurrent,
        std::map< std::string, std::map< double, std::map< int, double > > > rangeData,
        std::map< std::string, std::map< double, std::map< int, double > > > phaseData,
        std::map< std::string, Eigen::MatrixXi > availabilityMatrixPerUser,
        std::map< std::string, Eigen::MatrixXi > biasNumbersPerUser,
        bool ephemerisError, double arcTime, Eigen::MatrixXd ephemerisErrorBase,
        double codeWeight = 1.0,
        double phaseWeight = 0.01,
        double codeWeightISL = 1.0, double phaseWeightISL = 0.001,
        bool lighttimeCorrection = false,
        std::map< std::string, std::map< double, std::map< int, double > > > travelTimesPerUserPerEpoch = {{ "0", {{ 0.0, {{ 0, 0.0 }} }} }} );


//! Perform a single iteration of kinematic estimation of satellite position, clock errors and bias
//! Version that makes use of sparse matrices and block elimination, much quicker
Eigen::VectorXd getKinematicIterationSmart(
        std::vector< double > epochTimes,
        std::vector< std::shared_ptr< Ephemeris > > satelliteEphemerisGNSS,
        Eigen::VectorXd aPrioriEstimationsCurrent,
        std::map< std::string, std::map< double, std::map< int, double > > > rangeData,
        std::map< std::string, std::map< double, std::map< int, double > > > phaseData,
        std::map< std::string, Eigen::MatrixXi > availabilityMatrixPerUser,
        std::map< std::string, Eigen::MatrixXi > biasNumbersPerUser,
        bool ephemerisError, double arcTime, Eigen::MatrixXd ephemerisErrorBase,
        double codeWeight = 1.0,
        double phaseWeight = 0.001,
        double codeWeightISL = 1.0, double phaseWeightISL = 0.001,
        bool lighttimeCorrection = false,
        std::map< std::string, std::map< double, std::map< int, double > > > travelTimesPerUserPerEpoch = {{ "0", {{ 0.0, {{ 0, 0.0 }} }} }} );


//! Perform a single iteration of kinematic estimation of satellite position, clock errors and bias
//! Version that makes use of sparse matrices and block elimination, much quicker
//! Used when clocks of the navigation satellites are in phase and one clock is estimated
Eigen::VectorXd getKinematicIterationPhase(
        std::vector< double > epochTimes,
        std::vector< std::shared_ptr< Ephemeris > > satelliteEphemerisGNSS,
        Eigen::VectorXd aPrioriEstimationsCurrent,
        std::map< std::string, std::map< double, std::map< int, double > > > rangeData,
        std::map< std::string, std::map< double, std::map< int, double > > > phaseData,
        std::map< std::string, Eigen::MatrixXi > availabilityMatrixPerUser,
        bool ephemerisError, double arcTime, Eigen::MatrixXd ephemerisErrorBase,
        std::map< std::string, Eigen::MatrixXi > biasNumbersPerUser,
        double codeWeight = 1.0,
        double phaseWeight = 0.001,
        double codeWeightISL = 1.0, double phaseWeightISL = 0.001,
        bool lighttimeCorrection = false,
        std::map< std::string, std::map< double, std::map< int, double > > > travelTimesPerUserPerEpoch = {{ "0", {{ 0.0, {{ 0, 0.0 }} }} }} );


//! Perform entire kinematic estimation
//! Initiate the kinematic estimation
std::map< std::string, std::map< double, Eigen::Vector4d > > getKinematicEstimation(
        std::vector< double > epochTimes,
        std::vector< std::shared_ptr< Ephemeris > > satelliteEphemerisGNSS,
        std::map< std::string, std::map< double, std::map< int, double > > > rangeData,
        std::map< std::string, std::map< double, std::map< int, double > > > phaseData,
        std::map< std::string, Eigen::MatrixXi > availabilityMatrixPerUser,
        std::map< std::string, std::map< double, Eigen::VectorXd > >& allUsersKinIterations,
        std::map< int, Eigen::VectorXd >& biasIterations,
        std::map< int, std::map< std::string, std::map< double, Eigen::Vector4d > > >& estimationsKinPerUserPerIteration,
        bool ephemerisError, double arcTime, Eigen::MatrixXd ephemerisErrorBase,
        std::map< std::string, std::map< double, Eigen::Vector4d > > aPrioriPerUser = {{ "0", {{ 0, { 0.0, 0.0, 0.0, 0.0 } }} }},
        int maxIterationsKinematic = 5, //5
        double correctionLimitKinematic = 1e-8,
        double codeWeight = 1.0,
        double phaseWeight = 0.001,
        double codeWeightISL = 1.0, double phaseWeightISL = 0.001,
        bool lighttimeCorrection = false,
        std::map< std::string, std::map< double, std::map< int, double > > > travelTimesPerUserPerEpoch = {{ "0", {{ 0.0, {{ 0, 0.0 }} }} }},
        int solver = 2, bool clockPhase = false);


#endif // TOOLSNAVIGATION_H
