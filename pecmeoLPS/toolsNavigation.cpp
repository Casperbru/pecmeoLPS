#include "toolsNavigation.h"
#include "toolsGeneral.h"
#include <iostream>
#include "applicationOutput.h"

using namespace tudat::ephemerides;
using namespace tudat::input_output;

// ///  /// DOP CALCULATIONS
//! Find the point positioning design matrix for given point (user) and GNSS visible
//! per epoch per user
Eigen::MatrixXd getDesignMatrix(
        Eigen::Vector3d userSatellitePosition,
        Eigen::MatrixXd positionMatrixGNSSSatellites)
{
    auto numberOfNavigationSatellites = positionMatrixGNSSSatellites.rows();
    Eigen::MatrixXd designMatrix( numberOfNavigationSatellites, 4 );
    Eigen::RowVector3d gnssVector;
    Eigen::RowVector3d pecmeoVector = userSatellitePosition.segment( 0,3 );
    Eigen::RowVector3d positionDifference;

    for( int i = 0; i < numberOfNavigationSatellites; i++ )
    {
        gnssVector = positionMatrixGNSSSatellites.block( i, 0, 1, 3 );
        positionDifference = pecmeoVector - gnssVector;
        designMatrix.block( i, 0, 1, 4 ) << positionDifference / positionDifference.norm(), 1;
    }
//    std::cout << "Dit is de design mat: \n" << designMatrix << std::endl;
    return designMatrix;
};


//! Find the point positioning covariance matrix for given point (user) and GNSS visible
//! per epoch per user to be used for getDOPValues
Eigen::MatrixXd getCovarianceMatrix(
        Eigen::Vector3d userSatellitePosition,
        Eigen::MatrixXd positionMatrixGNSSSatellites)
{
    Eigen::MatrixXd designMatrix = getDesignMatrix( userSatellitePosition, positionMatrixGNSSSatellites );
    Eigen::MatrixXd covarianceMatrix = ( designMatrix.transpose() * designMatrix ).inverse();
//    std::cout << "Dit is de covariance mat: \n" << covarianceMatrix << std::endl;
    return covarianceMatrix;
};


//! Obtain the different DOP values and the number of links with navigation satellites ////////////////// kan weg
//! per epoch per user
Eigen::Vector6d getDOPValues(
        Eigen::Vector3d userSatellitePosition,
        Eigen::MatrixXd positionMatrixGNSSSatellites)
{
    double numberOfNavigationSatellites = positionMatrixGNSSSatellites.rows();
    Eigen::MatrixXd covarianceMatrix = getCovarianceMatrix( userSatellitePosition, positionMatrixGNSSSatellites );

    // Take the squares of the variances in the covariance matrix to obtain desired DOPs
    double valueGDOP = sqrt( covarianceMatrix(0,0) + covarianceMatrix(1,1) + covarianceMatrix(2,2) + covarianceMatrix(3,3) );
    double valuePDOP = sqrt( covarianceMatrix(0,0) + covarianceMatrix(1,1) + covarianceMatrix(2,2) );
    double valueHDOP = sqrt( covarianceMatrix(0,0) + covarianceMatrix(1,1) );
    double valueVDOP = sqrt( covarianceMatrix(2,2) );
    double valueTDOP = sqrt( covarianceMatrix(3,3) );

    Eigen::Vector6d resultsVector;
    resultsVector << numberOfNavigationSatellites, valueGDOP, valuePDOP, valueHDOP, valueVDOP, valueTDOP;

    return resultsVector;
};


//! Obtain the different DOP values and the number of links with navigation satellites
//! for all epochs per user used in toolsGeneral
std::map< double, Eigen::VectorXd > getDOPValues(
        std::vector< double > epochTimes,
        std::map< double, Eigen::Vector3d > userSatellitePositionAllEpochs,
        std::map< double, Eigen::MatrixXd > positionMatrixGNSSSatellitesAllEpochs)
{
    std::map< double, Eigen::VectorXd > resultsVectorMap;
    auto numberOfEpochs = epochTimes.size();
    for( unsigned long long i = 0; i < numberOfEpochs; i++ )
    {
        double epochTime = epochTimes.at(i);
        double numberOfNavigationSatellites = positionMatrixGNSSSatellitesAllEpochs[ epochTime ].rows();
        Eigen::MatrixXd covarianceMatrix = getCovarianceMatrix( userSatellitePositionAllEpochs[ epochTime ], positionMatrixGNSSSatellitesAllEpochs[ epochTime ] );

        // Take the squares of the variances in the covariance matrix to obtain desired DOPs
        double valueGDOP = sqrt( covarianceMatrix(0,0) + covarianceMatrix(1,1) + covarianceMatrix(2,2) + covarianceMatrix(3,3) );
        double valuePDOP = sqrt( covarianceMatrix(0,0) + covarianceMatrix(1,1) + covarianceMatrix(2,2) );
        double valueHDOP = sqrt( covarianceMatrix(0,0) + covarianceMatrix(1,1) );
        double valueVDOP = sqrt( covarianceMatrix(2,2) );
        double valueTDOP = sqrt( covarianceMatrix(3,3) );

        // Put all DOP results in vector
        Eigen::VectorXd resultsVector(9);
        resultsVector << numberOfNavigationSatellites, valueGDOP, valuePDOP, valueHDOP, valueVDOP, valueTDOP, userSatellitePositionAllEpochs[ epochTime ];
        resultsVectorMap[ epochTime ] = resultsVector;
    }
    return resultsVectorMap;
}


//! Obtain the state matrix of the visible GNSS satellites through ephemeris
//! for all epochs per user to be used with getDOPValues in toolsGeneral
std::map< double, Eigen::MatrixXd > getPositionMatrixGNSSSatellites(
        Eigen::MatrixXi availabilityMatrix,
        std::vector< std::shared_ptr< Ephemeris > > satelliteEphemerisGNSS,
        std::vector< double > epochTimes)
{
    std::map< double, Eigen::MatrixXd > positionMatrixGNSSSatellitesPerEpoch;
    auto numberOfEpochs = epochTimes.size();
    for( unsigned long long i = 0; i < numberOfEpochs; i++)
    {
        auto epochTime = epochTimes.at(i);
        auto numberOfNavigationSatellitesTotal = availabilityMatrix.cols();
        auto numberOfNavigationSatellitesAvailable =  availabilityMatrix.row(static_cast<int>(i)).sum();
        Eigen::MatrixXd positionMatrixGNSSSatellites( static_cast<int>(numberOfNavigationSatellitesAvailable), 3 );
        int rowIndex = 0;
        for( long long j = 0; j < numberOfNavigationSatellitesTotal; j++ )
        {
            if( availabilityMatrix( static_cast<unsigned>(i), j ) == 1 )
            {
                positionMatrixGNSSSatellites.block( rowIndex, 0, 1, 3) << ( satelliteEphemerisGNSS.at(static_cast<unsigned>(j))->getCartesianState( epochTime ).segment( 0, 3 ) ).transpose();
                rowIndex += 1;
            }
        }
        positionMatrixGNSSSatellitesPerEpoch[ epochTime ] = positionMatrixGNSSSatellites;
    }
    return positionMatrixGNSSSatellitesPerEpoch;
}



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
        bool lighttimeCorrection,
        std::map< std::string, std::map< int, double > > travelTimesPerUser)
{
    // Obtain constants
    int numberOfUsers = static_cast<int>( availabilityRowPerUser.size() );
    int numberOfSatellites = static_cast<int>( satelliteEphemerisGNSS.size() );
    int numberOfGNSS = numberOfSatellites - numberOfUsers;
    auto numberOfObservations = rangeData.size();
    double lightConstant = physical_constants::getSpeedOfLight< double >();

    // Create empty vectors and matrices
    std::vector< std::string > userNames;
    Eigen::MatrixXd largeDesignMatrix = Eigen::MatrixXd::Zero( numberOfObservations, numberOfUsers*4 );
    Eigen::MatrixXd weightMatrix = Eigen::MatrixXd::Identity( numberOfObservations, numberOfObservations );
    Eigen::VectorXd aPrioriRange= Eigen::VectorXd::Zero( numberOfObservations );
    int userIndex = 0;
    int observationIndex = 0;

    // Loop over users by looking at availilbity matrix map
    for( std::map< std::string, Eigen::RowVectorXi >::iterator itrUser = availabilityRowPerUser.begin(); itrUser != availabilityRowPerUser.end(); itrUser++ )
    {
        // Set user constants
        std::string userSatellite = itrUser->first;
        userNames.push_back( userSatellite );
        Eigen::RowVectorXi availabilityRow = itrUser->second;
        auto numberOfObservationsTotal = availabilityRowPerUser[ userSatellite ].size();
        std::map< int, double > travelTimes = travelTimesPerUser[ userSatellite ];

        // Find a priori position and clock offset
        Eigen::Vector3d aPrioriPosition = aPrioriEstimation.segment( userIndex*4, 3 );
        double clockError = aPrioriEstimation( userIndex*4+3 )/lightConstant;
//        double clockError = 0; // For no clock estimation

        // Loop over observations
        for(  int prn = 0; prn < numberOfObservationsTotal; prn++ )
        {
            if( availabilityRow( prn ) == 1 )
            {
                // Determine correct epoch time to receive correct state from ephemeris
                double correctedEpochTime = epochTime - clockError;
                if( lighttimeCorrection ){ correctedEpochTime = epochTime - clockError - travelTimes[ prn ]; }
                Eigen::Vector3d currentGNSSPosition = satelliteEphemerisGNSS.at(static_cast<unsigned>(prn))->getCartesianState( correctedEpochTime ).segment( 0, 3 );

                if( ephemerisError )
                {
                    // Determine ephemeris arc number
                    int arcNumber = floor(epochTime/arcTime);
                    double arcEpochTime = floor(epochTime/arcTime)*arcTime;
                    double arcPropagationTime = correctedEpochTime - arcEpochTime;
                    double maxErrorSMA = ephemerisErrorBase( arcNumber, prn );
                    double errorSMA = maxErrorSMA;
                    double earthGravitationalParameter = 3.986004418e14;

                    // Obtain state at start ephemeris error and propagate Kepler orbit to receive perturbed state at epoch
                    Eigen::Vector6d arcBeginStateCart = satelliteEphemerisGNSS.at(static_cast<unsigned>(prn))->getCartesianState( arcEpochTime );
                    Eigen::Vector6d arcBeginStateKep = orbital_element_conversions::convertCartesianToKeplerianElements(arcBeginStateCart, earthGravitationalParameter);
                    arcBeginStateKep(0) += errorSMA;
                    Eigen::Vector6d arcPropagatedStateKep = orbital_element_conversions::propagateKeplerOrbit(arcBeginStateKep, arcPropagationTime, earthGravitationalParameter );
                    Eigen::Vector6d arcPropagatedStateCart = orbital_element_conversions::convertKeplerianToCartesianElements(arcPropagatedStateKep, earthGravitationalParameter);

                    currentGNSSPosition = arcPropagatedStateCart.segment(0,3);
                }

                // Determine part of computed range
                Eigen::Vector3d positionDifference = aPrioriPosition - currentGNSSPosition;
                weightMatrix( observationIndex, observationIndex ) = codeWeight;
                double clockErrorLink = 0.0;
                if( prn >= numberOfGNSS )
                {
                    // If ISL, determine partial derivatives for other LPS satellite
                    int otherUserIndex = prn-numberOfGNSS;
                    positionDifference = aPrioriPosition - aPrioriEstimation.segment( otherUserIndex*4, 3 );                    
                    largeDesignMatrix.block( observationIndex, otherUserIndex*4, 1, 4 ) << -positionDifference.transpose() / positionDifference.norm(), -1;
                    clockErrorLink = -aPrioriEstimation( otherUserIndex*4+3 );
                    weightMatrix( observationIndex, observationIndex ) = codeWeightISL;
                }

                // Determine computed range
                aPrioriRange( observationIndex, 0 ) = positionDifference.norm() + aPrioriEstimation( userIndex*4+3 ) + clockErrorLink;
                // Store partial derivatives in large design matrix
                largeDesignMatrix.block( observationIndex, userIndex*4, 1, 4 ) << positionDifference.transpose() / positionDifference.norm(), 1;                
                observationIndex += 1;
            }
        }
        userIndex += 1;
    }

    // Create covariance and normal matrix
    Eigen::MatrixXd covarianceMatrix = ( largeDesignMatrix.transpose() * weightMatrix * largeDesignMatrix ).inverse();
    Eigen::MatrixXd normalMatrix = largeDesignMatrix.transpose() * weightMatrix * largeDesignMatrix;

    // Determine state correction
    Eigen::VectorXd deltaPosition = ( largeDesignMatrix.transpose() * weightMatrix * largeDesignMatrix ).inverse() * largeDesignMatrix.transpose() * weightMatrix * ( rangeData - aPrioriRange );

    // Some printing obtains for testing and verification
    if( numberOfUsers > 1 ) //&& (epochTime == 0.0||epochTime == 30.0 )
    {
//        Eigen::JacobiSVD<Eigen::MatrixXd> svd(normalMatrix);
//        double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
//        std::cout << "With code: " << pow( codeWeightISL, -0.5 ) << " and weight " << codeWeightISL << " op epoch: " << epochTime << " condition number: " << cond << std::endl;
//        if( cond > codeWeightISL/1.0e5 ){ std::cout << "Maybe" << std::endl; }

        //        std::cout << "\nAt epoch " << epochTime << ": "<< std::endl;
        //        std::cout << "\ndesign matrix: \n" << largeDesignMatrix << std::endl << std::endl;
        //        std::cout << "normal matrix: \n" << ( largeDesignMatrix.transpose() * largeDesignMatrix ) << std::endl << std::endl;
        //        std::cout << " rangeData   : " << rangeData.transpose() << "\n\n aPrioriRange: " << aPrioriRange.transpose() << "\n\n van elkaaraf: " << ( rangeData - aPrioriRange ).transpose() << std::endl;
        //        std::cout << "\n je delta  " << deltaPosition.size() << ": \n" << deltaPosition.transpose() << std::endl;
        //        std::cout << "\n Solution  " << (aPrioriEstimation + deltaPosition).size() << ": " << (aPrioriEstimation + deltaPosition).transpose() << std::endl;
    }

    return aPrioriEstimation + deltaPosition;
}


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
        int maxIterations,
        double correctionLimit,
        Eigen::VectorXd aPrioriEstimation,
        bool lighttimeCorrection,
        std::map< std::string, std::map< int, double > > travelTimesPerUser)
{
    // Set number of unknowns and current estimation
    long long numberOfUnknowns = aPrioriEstimation.size();
    Eigen::VectorXd currentEstimation = aPrioriEstimation;

    // Set first iteration in history vector
    iterationCoordinateHistory.resize( numberOfUnknowns, 1 );
    iterationCoordinateHistory.block( 0, 0, numberOfUnknowns, 1 ) = currentEstimation;

    // Set values to initiate while loop
    int numberOfIterations = 0;
    double correction = correctionLimit * 10;

    // Loop that corrects the estimation
    while( numberOfIterations < maxIterations && correction > correctionLimit )
    {
        numberOfIterations += 1;
        Eigen::VectorXd newEstimation = getPointPositionIteration( currentEstimation, rangeData, availabilityRowPerUser, satelliteEphemerisGNSS, epochTime, codeWeight, codeWeightISL, ephemerisError, arcTime, ephemerisErrorBase, lighttimeCorrection, travelTimesPerUser );
        Eigen::VectorXd estimationDifference = ( currentEstimation - newEstimation ).cwiseAbs();
        correction = estimationDifference.maxCoeff();
        currentEstimation = newEstimation;

        // Save in iteration history
        iterationCoordinateHistory.conservativeResize( iterationCoordinateHistory.rows(), iterationCoordinateHistory.cols() + 1 );
        iterationCoordinateHistory.block( 0, numberOfIterations, numberOfUnknowns, 1 ) = currentEstimation;
    }
//    std::cout << "At epoch " << epochTime << " SPP iteration " << numberOfIterations << " of 10" << std::endl;
    return currentEstimation;
}


//! Perform point positioning for one user over all epochs
std::map< std::string, std::map< double, Eigen::Vector4d > > getPointPositionSequence(
        std::vector< double > epochTimes,
        std::vector< std::shared_ptr< Ephemeris > > satelliteEphemerisGNSS,
        std::map< std::string, std::map< double, std::map< int, double > > > rangeData,
        std::map< std::string, Eigen::MatrixXi > availabilityMatrixPerUser,
        std::map< std::string, std::map< double, Eigen::VectorXd > >& allUsersPPIterations,
        double codeWeight, double codeWeightISL,
        bool ephemerisError, double arcTime, Eigen::MatrixXd ephemerisErrorBase,
        std::map< std::string, std::map< double, Eigen::Vector4d > > aPrioriEstimationsPerUser,
        int maxIterations,
        double correctionLimit,
        bool lighttimeCorrection,
        std::map< std::string, std::map< double, std::map< int, double > > > travelTimesPerUserPerEpoch)
{
    // Set constants
    unsigned long long numberOfEpochs = epochTimes.size();
    unsigned long long numberOfUsers = availabilityMatrixPerUser.size();

    // Start loop over all epochs and use getPointPositionEstimation
    std::map< std::string, std::map< double, Eigen::Vector4d > > estimationsPPPerUser;
    Eigen::VectorXd aPrioriEstimationAllUsers( numberOfUsers * 4 );
    for( unsigned long long epoch = 0; epoch < numberOfEpochs; epoch++ )
    {
        double epochTime = epochTimes.at( epoch );
        int userIndex = 0;
        Eigen::VectorXd rangeDataLocal(0);
        std::vector< std::string > userNames;
        std::map< std::string, std::map< int, double > > travelTimesPerUser;
        std::map< std::string, Eigen::RowVectorXi > availabilityRowPerUser;

        // Loop over users to obtain correct observation data
        for( std::map< std::string, Eigen::MatrixXi >::iterator itrUser = availabilityMatrixPerUser.begin(); itrUser != availabilityMatrixPerUser.end(); itrUser++ )
        {
            std::string userSatellite = itrUser->first;
            userNames.push_back( userSatellite );
            availabilityRowPerUser[ userSatellite ] = itrUser->second.row( static_cast<signed>(epoch) );
            travelTimesPerUser[ userSatellite ] = travelTimesPerUserPerEpoch[ userSatellite ][ epochTime ];
            aPrioriEstimationAllUsers.segment( userIndex*4, 4 ) = aPrioriEstimationsPerUser[ userSatellite ][ epochTime ];
            for ( std::map<int,double>::iterator itrRangeData = rangeData[ userSatellite ][ epochTime ].begin(); itrRangeData != rangeData[ userSatellite ][ epochTime ].end(); itrRangeData++ )
            {
                rangeDataLocal.conservativeResize( rangeDataLocal.size() + 1 );
                rangeDataLocal.tail(1)(0) = itrRangeData->second;
            }
            userIndex += 1;
        }

        Eigen::VectorXd aPrioriEstimation = aPrioriEstimationAllUsers;

        Eigen::MatrixXd iterationCoordinateHistory( aPrioriEstimation.size(), 1 );
        Eigen::VectorXd estimatedPositions = getPointPositionEstimation( rangeDataLocal, availabilityRowPerUser, iterationCoordinateHistory, satelliteEphemerisGNSS, epochTime, codeWeight, codeWeightISL, ephemerisError, arcTime, ephemerisErrorBase, maxIterations, correctionLimit, aPrioriEstimation, lighttimeCorrection, travelTimesPerUser );

        // Loop over users and save estimations and iteration history
        for( unsigned long long i = 0; i < numberOfUsers; i++ )
        {
            std::string userSatellite = userNames.at(i);
            estimationsPPPerUser[ userSatellite ][ epochTime ] = estimatedPositions.segment( static_cast<signed>(i)*4, 4 );
            unsigned long long numberOfaPrioriAvailable = aPrioriEstimationsPerUser[ userSatellite ].size();
            if( numberOfaPrioriAvailable < numberOfEpochs && epoch+1 < numberOfEpochs )
            {
                double nextEpochTime = epochTimes.at( epoch+1 );
                aPrioriEstimationsPerUser[ userSatellite ][ nextEpochTime ] = estimationsPPPerUser[ userSatellite ][ epochTime ];
            }

            Eigen::VectorXd iterationTemp( iterationCoordinateHistory.cols()*4 );
            for( auto it = 0; it < iterationCoordinateHistory.cols(); it++ )
            {
                iterationTemp.segment( it*4, 4 ) = iterationCoordinateHistory.block( i*4, it, 4, 1 );
            }
            allUsersPPIterations[ userSatellite ][ epochTime ] = iterationTemp;
        }
    }

//    std::cout << "\napriori: \n" << aPrioriEstimationsPerUser[ "PECMEOSat1" ][ 0.0 ].transpose() << ", " << aPrioriEstimationsPerUser[ "PECMEOSat1" ][ 30.0 ].transpose() << std::endl;
//    std::cout << "\nestimation: \n" << estimationsPPPerUser[ "PECMEOSat1" ][ 0.0 ].transpose() << ", " << estimationsPPPerUser[ "PECMEOSat1" ][ 30.0 ].transpose() << std::endl;
    return estimationsPPPerUser;
}


// /// KINEMATIC ORBIT DETERMINATION
//! Perform a single iteration of kinematic estimation of satellite position, clock errors and bias
//! Older version without block elimination and spars matrices, very slow, does not work for big
//! simulations, can not allocate the large matrices
Eigen::VectorXd getKinematicIteration(
        std::vector< double > epochTimes,
        std::vector< std::shared_ptr< Ephemeris > > satelliteEphemerisGNSS,
        Eigen::VectorXd aPrioriEstimationsCurrent,
        std::map< std::string, std::map< double, std::map< int, double > > > rangeData,
        std::map< std::string, std::map< double, std::map< int, double > > > phaseData,
        std::map< std::string, Eigen::MatrixXi > availabilityMatrixPerUser,
        std::map< std::string, Eigen::MatrixXi > biasNumbersPerUser,
        bool ephemerisError, double arcTime, Eigen::MatrixXd ephemerisErrorBase,
        double codeWeight, double phaseWeight,
        double codeWeightISL, double phaseWeightISL,
        bool lighttimeCorrection,
        std::map< std::string, std::map< double, std::map< int, double > > > travelTimesPerUserPerEpoch)
{
    // Determine constants
    double lightConstant = physical_constants::getSpeedOfLight< double >();
    long long numberOfEpochs = static_cast<signed>( epochTimes.size() );
    int numberOfUsers = static_cast<signed>( availabilityMatrixPerUser.size() );
    int numberOfBiases = biasNumbersPerUser.rbegin()->second.maxCoeff();
    int numberOfSatellites = static_cast<int>( satelliteEphemerisGNSS.size() );
    int numberOfGNSS = numberOfSatellites - numberOfUsers;
    int numberOfMeasurements = 2;

    // Dynamically grow apriori estimation vector
    Eigen::VectorXd aPrioriEstimations = aPrioriEstimationsCurrent;
    while( aPrioriEstimations.size() < numberOfEpochs*numberOfUsers*4 )
    {
        auto numberOfEstimationsCurrent = aPrioriEstimations.size();
        aPrioriEstimations.resize( numberOfEstimationsCurrent + 1 );
        aPrioriEstimations( numberOfEstimationsCurrent ) = 0.0;
    }

    // Dynamically grow apriori estimation vector by adding (estimated) bias
    while( aPrioriEstimations.size() < numberOfEpochs*numberOfUsers*4 + numberOfBiases )
    {        
        auto addedBiasNumber = aPrioriEstimations.size() - numberOfEpochs*numberOfUsers*4 + 1;
        int numberOfBiasNumbers = 0;
        double differenceRangePhaseSum = 0;
        numberOfMeasurements = 0;

        for( std::map< std::string, Eigen::MatrixXi >::iterator itrUser = biasNumbersPerUser.begin(); itrUser != biasNumbersPerUser.end(); itrUser++ )
        {
            std::string userSatellite = itrUser->first;
            Eigen::MatrixXi biasNumbers = itrUser->second;
            numberOfMeasurements += availabilityMatrixPerUser[ userSatellite ].sum()*2;
            for( unsigned i = 0; i < biasNumbers.rows(); i++)
            {
                double epochTime = epochTimes.at(i);
                for ( auto j = 0; j < biasNumbers.cols(); j++)
                {
                    if( biasNumbers( i, j ) == addedBiasNumber )
                    {
                        differenceRangePhaseSum += rangeData[ userSatellite ][ epochTime ][ j ] - phaseData[ userSatellite ][ epochTime ][ j ];
                        numberOfBiasNumbers += 1;
                    }
                }
            }
        }

        auto numberOfEstimationsCurrent = aPrioriEstimations.size();
        aPrioriEstimations.conservativeResize( numberOfEstimationsCurrent + 1 );
        double aPrioriBias = differenceRangePhaseSum / numberOfBiasNumbers;
        aPrioriEstimations( numberOfEstimationsCurrent ) = aPrioriBias;
    }

    // Pre-allocate vectors and matrices
    Eigen::MatrixXd designMatrix = Eigen::MatrixXd::Zero( numberOfMeasurements, numberOfEpochs*numberOfUsers*4 + numberOfBiases );
    Eigen::VectorXd measurements = Eigen::VectorXd::Zero( numberOfMeasurements );
    Eigen::VectorXd model = Eigen::VectorXd::Zero( numberOfMeasurements );
//    Eigen::MatrixXd weightMatrix = Eigen::MatrixXd::Identity( numberOfMeasurements, numberOfMeasurements );
//    Eigen::MatrixXd weightMatrix( numberOfMeasurements, numberOfMeasurements );
    Eigen::MatrixXd weightMatrix = Eigen::MatrixXd::Zero( numberOfMeasurements, numberOfMeasurements );
//    weightMatrix.setIdentity();
    Eigen::VectorXd correction = Eigen::VectorXd::Zero( numberOfEpochs*numberOfUsers*4 + numberOfBiases );

    // Loop over all the epochs
    int rowIterator = 0;
    for( int epoch = 0; epoch < numberOfEpochs; epoch++ )
    {
        double epochTime = epochTimes.at(static_cast<unsigned>(epoch));
        int userIndex = 0;
        for( std::map< std::string, Eigen::MatrixXi >::iterator itrUser = availabilityMatrixPerUser.begin(); itrUser != availabilityMatrixPerUser.end(); itrUser++ )
        {
            std::string userSatellite = itrUser->first;
            Eigen::RowVectorXi availabilityRow = itrUser->second.row( epoch );
            auto numberOfObservations = availabilityRow.size();

            Eigen::MatrixXi biasNumbers = biasNumbersPerUser[ userSatellite ];

            std::map< int, double > travelTimes = travelTimesPerUserPerEpoch[ userSatellite ][ epochTime ];

            std::cout << "userIndex: " << userIndex << std::endl;
            Eigen::VectorXd aPrioriEstimationsLocal = aPrioriEstimations.segment( epoch*numberOfUsers*4 + userIndex*4, 4 );
            Eigen::Vector3d aPrioriPosition = aPrioriEstimationsLocal.segment(0,3);
            double clockError = aPrioriEstimationsLocal(3)/lightConstant;

            for( int prn = 0; prn < numberOfObservations; prn++)
            {
                if( availabilityRow( prn ) == 1 )
                {
                    int biasNumber = biasNumbers( epoch, prn ) - 1;
                    // Determine correct epoch time to receive correct state from ephemeris
                    double correctedEpochTime = epochTime - clockError;
                    if( lighttimeCorrection ){ correctedEpochTime = epochTime - clockError - travelTimes[ prn ]; }
                    Eigen::Vector3d currentGNSSPosition = satelliteEphemerisGNSS.at( static_cast<unsigned>( prn ) )->getCartesianState( correctedEpochTime ).segment( 0, 3 );
                    if( ephemerisError )
                    {
                        // Determine ephemeris arc number
                        int arcNumber = floor(epochTime/arcTime);
                        double arcEpochTime = floor(epochTime/arcTime)*arcTime;
                        double arcPropagationTime = correctedEpochTime - arcEpochTime;
                        double maxErrorSMA = ephemerisErrorBase( arcNumber, prn );
                        double errorSMA = maxErrorSMA;
                        double earthGravitationalParameter = 3.986004418e14;

                        // Obtain state at start ephemeris error and propagate Kepler orbit to receive perturbed state at epoch
                        Eigen::Vector6d arcBeginStateCart = satelliteEphemerisGNSS.at(static_cast<unsigned>(prn))->getCartesianState( arcEpochTime );
                        Eigen::Vector6d arcBeginStateKep = orbital_element_conversions::convertCartesianToKeplerianElements(arcBeginStateCart, earthGravitationalParameter);
                        arcBeginStateKep(0) += errorSMA;
                        Eigen::Vector6d arcPropagatedStateKep = orbital_element_conversions::propagateKeplerOrbit(arcBeginStateKep, arcPropagationTime, earthGravitationalParameter );
                        Eigen::Vector6d arcPropagatedStateCart = orbital_element_conversions::convertKeplerianToCartesianElements(arcPropagatedStateKep, earthGravitationalParameter);

                        currentGNSSPosition = arcPropagatedStateCart.segment(0,3);
                    }

                    // Pre-allocate observation row of design matrix
                    Eigen::RowVectorXd observationRow = Eigen::RowVectorXd::Zero( numberOfUsers*4 );
                    Eigen::Vector3d positionDifference = aPrioriPosition - currentGNSSPosition;
                    double clockErrorLink = 0.0;
                    if( prn >= numberOfGNSS )
                    {
                        // If ISL, determine partial derivatives for other LPS satellite
                        codeWeight = codeWeightISL;
                        phaseWeight = phaseWeightISL;
                        int otherUserIndex = prn-numberOfGNSS;
                        positionDifference = aPrioriPosition - aPrioriEstimations.segment( epoch*numberOfUsers*4 + otherUserIndex*4, 3 );
                        observationRow.segment( otherUserIndex*4, 4 ) << -positionDifference.transpose() / positionDifference.norm(), -1;
                        clockErrorLink = aPrioriEstimations( epoch * numberOfUsers * 4 + otherUserIndex*4+3 );
                    }

                    // Determine observation row
                    observationRow.segment( userIndex*4, 4) << positionDifference.transpose()/positionDifference.norm(), 1;

                    // Add phase measurement to design matrix and weight matrix
                    designMatrix.block( rowIterator, epoch*numberOfUsers*4, 1, numberOfUsers*4 ) << observationRow;
                    designMatrix( rowIterator, numberOfEpochs*numberOfUsers*4 + biasNumber ) = -1.0; // Bias
                    weightMatrix( rowIterator, rowIterator ) = phaseWeight; // Weight for phase measurement = sigma^2

                    // Add code measurement to design matrix and weight matrix
                    designMatrix.block( rowIterator+1, epoch*numberOfUsers*4, 1, numberOfUsers*4 ) << observationRow;
                    weightMatrix( rowIterator+1, rowIterator+1 ) = codeWeight; // Weight for code measurement = sigma^2

                    // Add code and phase measurements to measurements vector
                    measurements( rowIterator ) = phaseData[ userSatellite ][ epochTime ][ prn ];
                    measurements( rowIterator+1 ) = rangeData[ userSatellite ][ epochTime ][ prn ];

                    // Find model code and phase data, and add to model vector
                    double model_phase = positionDifference.norm() + aPrioriEstimationsLocal(3) + clockErrorLink - aPrioriEstimations( numberOfEpochs*numberOfUsers*4 + biasNumber ); // model for phase
                    std::cout << "positionDifference.norm(): " << positionDifference.norm() << ", aPrioriEstimationsLocal(3): " << aPrioriEstimationsLocal(3) << ", clockErrorLink: " << clockErrorLink << ", aPrioriEstimations( numberOfEpochs*numberOfUsers*4 + biasNumber ): " << aPrioriEstimations( numberOfEpochs*numberOfUsers*4 + biasNumber ) << ", modelphase: " << model_phase << std::endl;
                    model( rowIterator) = model_phase;
                    double model_pseudorange = positionDifference.norm() + aPrioriEstimationsLocal(3) + clockErrorLink; // model for pseudo range
                    model( rowIterator+1 ) = model_pseudorange;

                    // Add 2 (1x phase + 1x code) to the row iteration number
                    rowIterator += 2;
                }
            }
            userIndex += 1;
        }
    }

    // Current error between measurements and model
    Eigen::VectorXd currentError = measurements - model;

    // Weighted least squares correction to estimation vector
    std::chrono::system_clock::time_point time1 = std::chrono::system_clock::now();
    correction = (designMatrix.transpose() * weightMatrix * designMatrix).ldlt().solve( designMatrix.transpose() * weightMatrix * currentError );
    std::chrono::system_clock::time_point time2 = std::chrono::system_clock::now();
    getCurrentTime( time1, time2 );

    // Determine normal matrix
    Eigen::MatrixXd normalMatrix = designMatrix.transpose() * weightMatrix * designMatrix;

    // Print options
    if( numberOfUsers > 1 ){
//        writeMatrixToFile( Eigen::MatrixXd(normalMatrix), "Ntot___.dat", 16, "C:/tudatBundle/tudatApplications/stepTwo/SimulationOutput/matrices");
    }

    // Determine new state
    Eigen::VectorXd estimation = aPrioriEstimations + correction;

    return estimation;
}


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
        double codeWeight, double phaseWeight,
        double codeWeightISL, double phaseWeightISL,
        bool lighttimeCorrection,
        std::map< std::string, std::map< double, std::map< int, double > > > travelTimesPerUserPerEpoch)
{
    // Determine constants
    double lightConstant = physical_constants::getSpeedOfLight< double >();
    long long numberOfEpochs = static_cast<signed>( epochTimes.size() );
    int numberOfUsers = static_cast<signed>( availabilityMatrixPerUser.size() );
    int numberOfBiases = biasNumbersPerUser.rbegin()->second.maxCoeff();
    int numberOfSatellites = static_cast<int>( satelliteEphemerisGNSS.size() );
    int numberOfGNSS = numberOfSatellites - numberOfUsers;
    int matrixSize = numberOfEpochs*numberOfUsers*4;
    int localMatrixSize = numberOfUsers*4;
    double codeWeightGNSS = codeWeight;
    double phaseWeightGNSS = phaseWeight;

    // Pre-allocate vectors and matrices
    Eigen::SparseMatrix<double> Nxx_inv( matrixSize, matrixSize );
    Eigen::SparseMatrix<double> Nxx_( matrixSize, matrixSize );
    Eigen::SparseMatrix<double> Nxb( matrixSize, numberOfBiases );
    Eigen::SparseMatrix<double> Nbb( numberOfBiases, numberOfBiases );
    Eigen::VectorXd nx = Eigen::VectorXd::Zero( matrixSize );
    Eigen::VectorXd nb = Eigen::VectorXd::Zero( numberOfBiases );    

    // Dynamically grow apriori estimation vector
    Eigen::VectorXd apriBias = Eigen::VectorXd::Zero( numberOfBiases );
    Eigen::VectorXd aPrioriEstimations = aPrioriEstimationsCurrent;
    while( aPrioriEstimations.size() < matrixSize )
    {
        long long numberOfEstimationsCurrent = aPrioriEstimations.size();
        aPrioriEstimations.resize( numberOfEstimationsCurrent + 1 );
        aPrioriEstimations( numberOfEstimationsCurrent ) = 0.0;
    }

    // Dynamically grow apriori estimation vector by adding (estimated) bias    
    while( aPrioriEstimations.size() < matrixSize + numberOfBiases )
    {
        auto addedBiasNumber = aPrioriEstimations.size() - matrixSize + 1;
        int numberOfBiasNumbers = 0;
        double differenceRangePhaseSum = 0;

        for( std::map< std::string, Eigen::MatrixXi >::iterator itrUser = biasNumbersPerUser.begin(); itrUser != biasNumbersPerUser.end(); itrUser++ )
        {
            std::string userSatellite = itrUser->first;
            Eigen::MatrixXi biasNumbers = itrUser->second;
            for( unsigned i = 0; i < biasNumbers.rows(); i++)
            {
                double epochTime = epochTimes.at(i);
                for ( auto j = 0; j < biasNumbers.cols(); j++)
                {
                    if( biasNumbers( i, j ) == addedBiasNumber )
                    {
                        differenceRangePhaseSum += rangeData[ userSatellite ][ epochTime ][ j ] - phaseData[ userSatellite ][ epochTime ][ j ];
                        numberOfBiasNumbers += 1;
                    }
                }
            }
        }

        auto numberOfEstimationsCurrent = aPrioriEstimations.size();
        aPrioriEstimations.conservativeResize( numberOfEstimationsCurrent + 1 );
        double aPrioriBias = differenceRangePhaseSum / numberOfBiasNumbers;
        aPrioriEstimations( numberOfEstimationsCurrent ) = aPrioriBias;
        apriBias( addedBiasNumber - 1 ) = aPrioriBias/numberOfBiasNumbers;
    }

    // Loop over all the epochs
    for( int epoch = 0; epoch < numberOfEpochs; epoch++ )
    {
        // Pre-allocate sub-matrix corresponding to current epoch.
        Eigen::MatrixXd NxxLocal = Eigen::MatrixXd::Zero( localMatrixSize, localMatrixSize );
        double epochTime = epochTimes.at(static_cast<unsigned>(epoch));
        int userIndex = 0;

        // Loop over different users
        for( std::map< std::string, Eigen::MatrixXi >::iterator itrUser = availabilityMatrixPerUser.begin(); itrUser != availabilityMatrixPerUser.end(); itrUser++ )
        {
            // Set current users, corresponing bias numbers and travel times
            std::string userSatellite = itrUser->first;
            Eigen::RowVectorXi availabilityRow = itrUser->second.row( epoch );
            auto numberOfObservations = availabilityRow.size();
//            if( numberOfUsers > 1 ){ numberOfObservations += -9; } // Option for no ISL, but simultaneous solving
            Eigen::MatrixXi biasNumbers = biasNumbersPerUser[ userSatellite ];
            std::map< int, double > travelTimes = travelTimesPerUserPerEpoch[ userSatellite ][ epochTime ];

            // Determine the a priori state for current user and epoch
            Eigen::VectorXd aPrioriEstimationsLocal = aPrioriEstimations.segment( epoch*localMatrixSize + userIndex*4, 4 );
            Eigen::Vector3d aPrioriPosition = aPrioriEstimationsLocal.segment(0,3);
            double clockError = aPrioriEstimationsLocal(3)/lightConstant;
//            double clockError = 0; // Option for estimating no clock
//            Eigen::MatrixXd matrixje = Eigen::MatrixXd::Zero( numberOfObservations, localMatrixSize ); // Some print verification function
//            Eigen::MatrixXd matrixje2 = Eigen::MatrixXd::Zero( numberOfObservations, numberOfBiases ); // Some print verification function

            // Loop over observations
            for( int prn = 0; prn < numberOfObservations; prn++)
            {
                if( availabilityRow( prn ) == 1 )
                {
                    // Determine current bias number
                    int biasNumberLocal = biasNumbers( epoch, prn );
                    // Determine correct epoch time to receive correct state from ephemeris
                    double correctedEpochTime = epochTime - clockError;
                    if( lighttimeCorrection ){ correctedEpochTime = epochTime - clockError - travelTimes[ prn ]; }
                    Eigen::Vector3d currentGNSSPosition = satelliteEphemerisGNSS.at( static_cast<unsigned>( prn ) )->getCartesianState( correctedEpochTime ).segment( 0, 3 );

                    if( ephemerisError )
                    {
                        // Determine ephemeris arc number
                        int arcNumber = floor(epochTime/arcTime);
                        double arcEpochTime = floor(epochTime/arcTime)*arcTime;
                        double arcPropagationTime = correctedEpochTime - arcEpochTime;
                        double maxErrorSMA = ephemerisErrorBase( arcNumber, prn );
                        double errorSMA = maxErrorSMA;
                        double earthGravitationalParameter = 3.986004418e14;

                        // Obtain state at start ephemeris error and propagate Kepler orbit to receive perturbed state at epoch
                        Eigen::Vector6d arcBeginStateCart = satelliteEphemerisGNSS.at(static_cast<unsigned>(prn))->getCartesianState( arcEpochTime );
                        Eigen::Vector6d arcBeginStateKep = orbital_element_conversions::convertCartesianToKeplerianElements(arcBeginStateCart, earthGravitationalParameter);
                        arcBeginStateKep(0) += errorSMA;
                        Eigen::Vector6d arcPropagatedStateKep = orbital_element_conversions::propagateKeplerOrbit(arcBeginStateKep, arcPropagationTime, earthGravitationalParameter );
                        Eigen::Vector6d arcPropagatedStateCart = orbital_element_conversions::convertKeplerianToCartesianElements(arcPropagatedStateKep, earthGravitationalParameter);

                        currentGNSSPosition = arcPropagatedStateCart.segment(0,3);
                    }

                    // Pre-allocate observation row for design matrix
                    Eigen::RowVectorXd observationRow = Eigen::RowVectorXd::Zero( localMatrixSize );
                    // Determine part of computed range
                    Eigen::Vector3d positionDifference = aPrioriPosition - currentGNSSPosition;
                    double clockErrorLink = 0.0;
                    codeWeight = codeWeightGNSS;
                    phaseWeight = phaseWeightGNSS;
                    if( prn >= numberOfGNSS )
                    {
                        // If ISL, determine partial derivatives for other LPS satellite and set correct weights
                        codeWeight = codeWeightISL;
                        phaseWeight = phaseWeightISL;
                        int otherUserIndex = prn-numberOfGNSS;
                        positionDifference = aPrioriPosition - aPrioriEstimations.segment( epoch*localMatrixSize + otherUserIndex*4, 3 );
                        // Set partial derivatives wrt other LPS sat in observation row
                        observationRow.segment( otherUserIndex*4, 4 ) << -positionDifference.transpose() / positionDifference.norm(), -1;
                        clockErrorLink = -aPrioriEstimations( epoch * numberOfUsers * 4 + otherUserIndex*4+3 );

                        // Set values in bias-state matrices
                        Nxb.coeffRef( epoch*localMatrixSize + otherUserIndex*4 + 0, biasNumberLocal - 1 ) += -observationRow( otherUserIndex*4 + 0) * phaseWeight;
                        Nxb.coeffRef( epoch*localMatrixSize + otherUserIndex*4 + 1, biasNumberLocal - 1 ) += -observationRow( otherUserIndex*4 + 1) * phaseWeight;
                        Nxb.coeffRef( epoch*localMatrixSize + otherUserIndex*4 + 2, biasNumberLocal - 1 ) += -observationRow( otherUserIndex*4 + 2) * phaseWeight;
                        Nxb.coeffRef( epoch*localMatrixSize + otherUserIndex*4 + 3, biasNumberLocal - 1 ) += -observationRow( otherUserIndex*4 + 3) * phaseWeight;
                    }

                    // Set partial derivatives wrt current LPS sat in observation row
                    observationRow.segment( userIndex*4, 4) << positionDifference.transpose()/positionDifference.norm(), 1;
//                    matrixje.block(prn, 0, 1, localMatrixSize) = observationRow; // Some print verification functions
//                    matrixje2(prn, biasNumberLocal-1) = -1; // Some print verification functions

                    // Find mode code and phase data, and add to model vector
                    double model_phase = positionDifference.norm() + aPrioriEstimationsLocal(3) + clockErrorLink - aPrioriEstimations( matrixSize + biasNumberLocal - 1 ); // model for phase
                    double model_pseudorange = positionDifference.norm() + aPrioriEstimationsLocal(3) + clockErrorLink; // model for pseudo range

                    // Build onto error matrices
                    nx.segment( epoch*localMatrixSize, localMatrixSize ) += observationRow * codeWeight * ( rangeData[ userSatellite ][ epochTime ][ prn ] - model_pseudorange );
                    nx.segment( epoch*localMatrixSize, localMatrixSize ) += observationRow * phaseWeight * ( phaseData[ userSatellite ][ epochTime ][ prn ] - model_phase );
                    nb( biasNumberLocal - 1 ) += -1 * phaseWeight * ( phaseData[ userSatellite ][ epochTime ][ prn ] - model_phase );

                    // Build onto normal equation matrices
                    NxxLocal += ( observationRow.transpose() * observationRow ) * ( phaseWeight + codeWeight );
                    Nbb.coeffRef( biasNumberLocal - 1, biasNumberLocal - 1 ) += phaseWeight;
                    Nxb.coeffRef( epoch*localMatrixSize + userIndex*4 + 0, biasNumberLocal - 1 ) += -observationRow( userIndex*4 + 0) * phaseWeight;
                    Nxb.coeffRef( epoch*localMatrixSize + userIndex*4 + 1, biasNumberLocal - 1 ) += -observationRow( userIndex*4 + 1) * phaseWeight;
                    Nxb.coeffRef( epoch*localMatrixSize + userIndex*4 + 2, biasNumberLocal - 1 ) += -observationRow( userIndex*4 + 2) * phaseWeight;
                    Nxb.coeffRef( epoch*localMatrixSize + userIndex*4 + 3, biasNumberLocal - 1 ) += -observationRow( userIndex*4 + 3) * phaseWeight;
                }
            }
//            writeMatrixToFile( matrixje, "Hxx_"+std::to_string(epoch)+"_"+std::to_string(userIndex)+".dat", 16, "C:/tudatBundle/tudatApplications/stepTwo/SimulationOutput/matrices"); // Some print verification functions
//            writeMatrixToFile( matrixje, "Hxx_"+std::to_string(userIndex)+"_"+std::to_string(epoch)+".dat", 16, "C:/tudatBundle/tudatApplications/pecmeoLPS/SimulationOutput/matrices"); // Some print verification functions
//            writeMatrixToFile( matrixje2, "Hbb_"+std::to_string(userIndex)+"_"+std::to_string(epoch)+".dat", 16, "C:/tudatBundle/tudatApplications/pecmeoLPS/SimulationOutput/matrices"); // Some print verification functions
            userIndex += 1;
        }

        // Inverse the local Nxx for faster inversion
        Eigen::MatrixXd NxxLocal_inv( localMatrixSize, localMatrixSize );
        NxxLocal_inv = NxxLocal.inverse(); // Inverse of local N_XX 4x4 or 36x36 matrix
        for( int i = 0; i < localMatrixSize; i++ )
        {
            for( int j = 0; j < localMatrixSize; j++ )
            {
                Nxx_inv.coeffRef( epoch*localMatrixSize + i, epoch*localMatrixSize + j ) = NxxLocal_inv( i, j );
                Nxx_.coeffRef( epoch*localMatrixSize + i, epoch*localMatrixSize + j ) = NxxLocal( i, j );
            }
        }
    }


    // Bias and Position updates
//    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver; // Other solver options
//    Eigen::SparseQR< Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver; // Other solver options
    Eigen::SparseQR< Eigen::SparseMatrix<double>, Eigen::AMDOrdering<int> > solver; // SparseQr
//    Eigen::SparseQR< Eigen::SparseMatrix<double>, Eigen::NaturalOrdering<int> > solver;/ // Other solver options
    solver.compute( Nbb - Nxb.transpose() * Nxx_inv * Nxb );
    Eigen::VectorXd deltaB = solver.solve( nb - Nxb.transpose() * Nxx_inv * nx );
    Eigen::VectorXd deltaX = Nxx_inv * ( nx - Nxb * deltaB );
    Eigen::VectorXd correction( deltaB.size() + deltaX.size() );
    correction << deltaX, deltaB;
    Eigen::VectorXd estimation = aPrioriEstimations + correction;

    #define FOLDINGSTART { // Create covariance matrices // This takes extremely long
//    Eigen::SparseMatrix<double> Qbb( numberOfBiases, numberOfBiases );
//    Eigen::MatrixXd Qbb( numberOfBiases, numberOfBiases );
//    Eigen::SparseMatrix<double> Ibb( numberOfBiases, numberOfBiases );
//    Ibb.setIdentity();
//    Qbb = solver.solve( Ibb );
//    Eigen::SparseMatrix<double> Qxx( matrixSize, matrixSize );
//    Eigen::MatrixXd Qxx( matrixSize, matrixSize );
//    Qxx = ( Nxx_inv * Nxb ) * Qbb * ( Nxx_inv * Nxb ).transpose();
//    Qxx += Nxx_inv;
//    Eigen::VectorXd QxxDiag = Qxx.diagonal();
//    std::cout << "covariance done" << std::endl; getCurrentTime();
    #define FOLDINGEND }

    if( numberOfUsers > 0 ){ // Options for printing and saving intermediate results
    #define FOLDINGSTART { // Saving results
//        writeMatrixToFile( aPrioriEstimations, "apri_"+setting+"_"+std::to_string(iteration)+".dat", 16, "C:/tudatBundle/tudatApplications/stepTwo/SimulationOutput/matrices");
//        writeMatrixToFile( Eigen::MatrixXd(Nxx_), "Nxx_20m.dat", 16, "C:/tudatBundle/tudatApplications/pecmeoLPS/SimulationOutput/matrices");
//        writeMatrixToFile( Eigen::MatrixXd(Nxb), "Nxb_10m.dat", 16, "C:/tudatBundle/tudatApplications/pecmeoLPS/SimulationOutput/matrices");
//        writeMatrixToFile( Eigen::MatrixXd(Nxb.transpose()), "NxbT_10m.dat", 16, "C:/tudatBundle/tudatApplications/pecmeoLPS/SimulationOutput/matrices");
//        writeMatrixToFile( Eigen::MatrixXd(Nbb), "Nbb_10m.dat", 16, "C:/tudatBundle/tudatApplications/pecmeoLPS/SimulationOutput/matrices");
//        writeMatrixToFile( nb, "nb_"+setting+"_"+std::to_string(iteration)+".dat", 16, "C:/tudatBundle/tudatApplications/stepTwo/SimulationOutput/matrices");
//        writeMatrixToFile( nx, "nx_"+setting+"_"+std::to_string(iteration)+".dat", 16, "C:/tudatBundle/tudatApplications/stepTwo/SimulationOutput/matrices");
//        writeMatrixToFile( deltaB, "deltaB_"+setting+"_"+std::to_string(iteration)+".dat", 16, "C:/tudatBundle/tudatApplications/stepTwo/SimulationOutput/matrices");
//        writeMatrixToFile( deltaX, "deltaX_"+setting+"_"+std::to_string(iteration)+".dat", 16, "C:/tudatBundle/tudatApplications/stepTwo/SimulationOutput/matrices");
//        writeMatrixToFile( estimation, "estimation_"+setting+"_"+std::to_string(iteration)+".dat", 16, "C:/tudatBundle/tudatApplications/stepTwo/SimulationOutput/matrices");
//        writeMatrixToFile( QxxDiag, "QxxDiag.dat", 16, "C:/tudatBundle/tudatApplications/stepTwo/SimulationOutput/PECMEO_30sec/GalileoGLONASSBeiDouGPS/1d/GNSSISL");
    #define FOLDINGEND }

    #define FOLDINGSTART { // Condition number
//        Eigen::JacobiSVD<Eigen::MatrixXd> svd(Nxx_);
//        double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
//        std::cout << "With code weight " << codeWeightISL << " condition number: " << cond << std::endl;
//        std::cout << "With code: " << pow( codeWeightISL, -0.5 ) << " and weight " << codeWeightISL << " op epoch: " << epochTime << " condition number: " << cond << std::endl;
//        if( cond > codeWeightISL/1.0e5 ){ std::cout << "Maybe" << std::endl; }
    #define FOLDINGEND }

    #define FOLDINGSTART { // Printing stuff
        std::cout << "\napriori: \n" << aPrioriEstimations.segment(0,36).transpose() << std::endl;
//        std::cout << "\nnb: \n" << nb.segment(0,16).transpose() << std::endl;
//        std::cout << "\nnx: \n" << nx.segment(0,16).transpose() << std::endl;
//        std::cout << "\nNbb: \n" << Nbb.block(0,0,16,16) << std::endl;
//        std::cout << "\nNxb: \n" << Nxb.block(0,0,16,16) << std::endl;
//        std::cout << "\nNxb Transpose: \n" << (Nxb.transpose()).block(0,0,16,16) << std::endl;
//        std::cout << "Nxx inv: \n" << Nxx_inv.block(0,0,16,16) << std::endl;
//        std::cout << "Nxx inv: \n" << Nxx_inv << std::endl;
//        std::cout << "Nxx inv: \n" << Eigen::MatrixXd(Nxx_inv) << std::endl;
//        std::cout << "\ndeltaB: \n" << deltaB.segment(0,16).transpose() << std::endl;
//        std::cout << "\ndeltaX: \n" << deltaX.segment(0,16).transpose() << std::endl;
        std::cout << "\nestimation: \n" << estimation.segment(0,36).transpose() << std::endl;
        std::cout << "\nbias: \n" << estimation.tail(numberOfBiases).segment(0,8).transpose() << std::endl;

//        std::cout << "Qxx: \n" << Qxx << std::endl;
//        std::cout << "Nxx_inv: \n" << Eigen::MatrixXd(Nxx_inv) << std::endl;
//        std::cout << "Qxx: \n" << Qxx << std::endl;
//        std::cout << "QxxDiag: \n" << QxxDiag.transpose() << std::endl;
    #define FOLDINGEND }
    }

    return estimation;
}


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
        std::map< std::string, Eigen::MatrixXi > biasNumbersPerUser,
        bool ephemerisError, double arcTime, Eigen::MatrixXd ephemerisErrorBase,
        double codeWeight, double phaseWeight,
        double codeWeightISL, double phaseWeightISL,
        bool lighttimeCorrection,
        std::map< std::string, std::map< double, std::map< int, double > > > travelTimesPerUserPerEpoch)
{
    // Determine constants
    double lightConstant = physical_constants::getSpeedOfLight< double >();
    long long numberOfEpochs = static_cast<signed>( epochTimes.size() );
    int numberOfUsers = static_cast<signed>( availabilityMatrixPerUser.size() );
    int numberOfBiases = biasNumbersPerUser.rbegin()->second.maxCoeff();
    int numberOfSatellites = static_cast<int>( satelliteEphemerisGNSS.size() );
    int numberOfGNSS = numberOfSatellites - numberOfUsers;
    int numberOfUnknowns = 3; //4;
    int localMatrixSize = numberOfUsers * numberOfUnknowns + 1; //numberOfUsers*4;
//    int localMatrixSize = numberOfUsers * numberOfUnknowns; // NO CLOCK OPTION
    int matrixSize = numberOfEpochs * localMatrixSize; //numberOfEpochs*numberOfUsers*4;
    double codeWeightGNSS = codeWeight;
    double phaseWeightGNSS = phaseWeight;

    // Pre-allocate vectors and matrices
    Eigen::SparseMatrix<double> Nxx_inv( matrixSize, matrixSize );
    Eigen::SparseMatrix<double> Nxb( matrixSize, numberOfBiases );
    Eigen::SparseMatrix<double> Nbb( numberOfBiases, numberOfBiases );
    Eigen::VectorXd nx = Eigen::VectorXd::Zero( matrixSize );
    Eigen::VectorXd nb = Eigen::VectorXd::Zero( numberOfBiases );

    // Dynamically grow apriori estimation vector
    Eigen::VectorXd aPrioriEstimations = aPrioriEstimationsCurrent;
    while( aPrioriEstimations.size() < matrixSize )
    {
        long long numberOfEstimationsCurrent = aPrioriEstimations.size();
        aPrioriEstimations.resize( numberOfEstimationsCurrent + 1 );
        aPrioriEstimations( numberOfEstimationsCurrent ) = 0.0;
    }

    // Dynamically grow apriori estimation vector by adding (estimated) bias
    while( aPrioriEstimations.size() < matrixSize + numberOfBiases )
    {
        auto addedBiasNumber = aPrioriEstimations.size() - matrixSize + 1;
        int numberOfBiasNumbers = 0;
        double differenceRangePhaseSum = 0;

        for( std::map< std::string, Eigen::MatrixXi >::iterator itrUser = biasNumbersPerUser.begin(); itrUser != biasNumbersPerUser.end(); itrUser++ )
        {
            std::string userSatellite = itrUser->first;
            Eigen::MatrixXi biasNumbers = itrUser->second;
            for( unsigned i = 0; i < biasNumbers.rows(); i++)
            {
                double epochTime = epochTimes.at(i);
                for ( auto j = 0; j < biasNumbers.cols(); j++)
                {
                    if( biasNumbers( i, j ) == addedBiasNumber )
                    {
                        differenceRangePhaseSum += rangeData[ userSatellite ][ epochTime ][ j ] - phaseData[ userSatellite ][ epochTime ][ j ];
                        numberOfBiasNumbers += 1;
                    }
                }
            }
        }

        auto numberOfEstimationsCurrent = aPrioriEstimations.size();
        aPrioriEstimations.conservativeResize( numberOfEstimationsCurrent + 1 );
        double aPrioriBias = differenceRangePhaseSum / numberOfBiasNumbers;
        aPrioriEstimations( numberOfEstimationsCurrent ) = aPrioriBias;
    }

    // Loop over all the epochs
    for( int epoch = 0; epoch < numberOfEpochs; epoch++ )
    {
        // Pre-allocate sub-matrix corresponding to current epoch.
        Eigen::MatrixXd NxxLocal = Eigen::MatrixXd::Zero( localMatrixSize, localMatrixSize );
        double epochTime = epochTimes.at(static_cast<unsigned>(epoch));
        int userIndex = 0;

        // Loop over different users
        for( std::map< std::string, Eigen::MatrixXi >::iterator itrUser = availabilityMatrixPerUser.begin(); itrUser != availabilityMatrixPerUser.end(); itrUser++ )
        {
            // Set current users, corresponing bias numbers and travel times
            std::string userSatellite = itrUser->first;
            Eigen::RowVectorXi availabilityRow = itrUser->second.row( epoch );
            auto numberOfObservations = availabilityRow.size();
            Eigen::MatrixXi biasNumbers = biasNumbersPerUser[ userSatellite ];
            std::map< int, double > travelTimes = travelTimesPerUserPerEpoch[ userSatellite ][ epochTime ];

            // Determine the a priori state for current user and epoch
            Eigen::VectorXd aPrioriEstimationsLocal = aPrioriEstimations.segment( epoch*localMatrixSize + userIndex*numberOfUnknowns, numberOfUnknowns );
            Eigen::Vector3d aPrioriPosition = aPrioriEstimationsLocal.segment(0,3);
            double clockError = aPrioriEstimations( epoch * localMatrixSize + localMatrixSize - 1 ) / lightConstant;
//            double clockError = 0; // NO CLOCK OPTION

            // Loop over observations
            for( int prn = 0; prn < numberOfObservations; prn++)
            {
                if( availabilityRow( prn ) == 1 )
                {
                    // Determine current bias number
                    int biasNumberLocal = biasNumbers( epoch, prn );
                    // Determine correct epoch time to receive correct state from ephemeris
                    double correctedEpochTime = epochTime - clockError;
                    if( lighttimeCorrection ){ correctedEpochTime = epochTime - clockError - travelTimes[ prn ]; }
                    Eigen::Vector3d currentGNSSPosition = satelliteEphemerisGNSS.at( static_cast<unsigned>( prn ) )->getCartesianState( correctedEpochTime ).segment( 0, 3 );
                    if( ephemerisError )
                    {
                        // Determine ephemeris arc number
                        int arcNumber = floor(epochTime/arcTime);
                        double arcEpochTime = floor(epochTime/arcTime)*arcTime;
                        double arcPropagationTime = correctedEpochTime - arcEpochTime;
                        double maxErrorSMA = ephemerisErrorBase( arcNumber, prn );
                        double errorSMA = maxErrorSMA;
                        double earthGravitationalParameter = 3.986004418e14;

                        // Obtain state at start ephemeris error and propagate Kepler orbit to receive perturbed state at epoch
                        Eigen::Vector6d arcBeginStateCart = satelliteEphemerisGNSS.at(static_cast<unsigned>(prn))->getCartesianState( arcEpochTime );
                        Eigen::Vector6d arcBeginStateKep = orbital_element_conversions::convertCartesianToKeplerianElements(arcBeginStateCart, earthGravitationalParameter);
                        arcBeginStateKep(0) += errorSMA;
                        Eigen::Vector6d arcPropagatedStateKep = orbital_element_conversions::propagateKeplerOrbit(arcBeginStateKep, arcPropagationTime, earthGravitationalParameter );
                        Eigen::Vector6d arcPropagatedStateCart = orbital_element_conversions::convertKeplerianToCartesianElements(arcPropagatedStateKep, earthGravitationalParameter);

                        currentGNSSPosition = arcPropagatedStateCart.segment(0,3);
                    }

                    // Pre-allocate observation row for design matrix
                    Eigen::RowVectorXd observationRow = Eigen::RowVectorXd::Zero( localMatrixSize );
                    observationRow( localMatrixSize - 1 ) = 1; // Comment out for NO CLOCK
                    // Determine part of computed range
                    Eigen::Vector3d positionDifference = aPrioriPosition - currentGNSSPosition;
                    double clockErrorLink = 0.0;
                    codeWeight = codeWeightGNSS;
                    phaseWeight = phaseWeightGNSS;
                    if( prn >= numberOfGNSS )
                    {
                        // If ISL, determine partial derivatives for other LPS satellite and set correct weights
                        codeWeight = codeWeightISL;
                        phaseWeight = phaseWeightISL;
                        int otherUserIndex = prn-numberOfGNSS;
                        positionDifference = aPrioriPosition - aPrioriEstimations.segment( epoch*localMatrixSize + otherUserIndex*numberOfUnknowns, 3 );
                        // Set partial derivatives wrt other LPS sat in observation row
                        observationRow.segment( otherUserIndex*numberOfUnknowns, numberOfUnknowns ) << -positionDifference.transpose() / positionDifference.norm();
                        observationRow( localMatrixSize - 1 ) = 0; // Comment out for NO CLOCK // testline
                        clockErrorLink = -aPrioriEstimations( epoch * localMatrixSize + localMatrixSize - 1 );

                        // Set values in bias-state matrices
                        Nxb.coeffRef( epoch*localMatrixSize + otherUserIndex*numberOfUnknowns + 0, biasNumberLocal - 1 ) += -observationRow( otherUserIndex*numberOfUnknowns + 0) * phaseWeight;
                        Nxb.coeffRef( epoch*localMatrixSize + otherUserIndex*numberOfUnknowns + 1, biasNumberLocal - 1 ) += -observationRow( otherUserIndex*numberOfUnknowns + 1) * phaseWeight;
                        Nxb.coeffRef( epoch*localMatrixSize + otherUserIndex*numberOfUnknowns + 2, biasNumberLocal - 1 ) += -observationRow( otherUserIndex*numberOfUnknowns + 2) * phaseWeight;
                    }

                    // Set partial derivatives wrt current LPS sat in observation row
                    observationRow.segment( userIndex*numberOfUnknowns, numberOfUnknowns ) << positionDifference.transpose()/positionDifference.norm();

                    // Find mode code and phase data, and add to model vector
                    double model_phase = positionDifference.norm() + aPrioriEstimations( epoch * localMatrixSize + localMatrixSize - 1 ) + clockErrorLink - aPrioriEstimations( matrixSize + biasNumberLocal - 1 ); // model for phase
                    double model_pseudorange = positionDifference.norm() + aPrioriEstimations( epoch * localMatrixSize + localMatrixSize - 1 ) + clockErrorLink; // model for pseudo range
//                    double model_phase = positionDifference.norm() - aPrioriEstimations( matrixSize + biasNumberLocal - 1 ); // model for phase NO CLOCK OPTION
//                    double model_pseudorange = positionDifference.norm(); // model for pseudo range NO CLOCK OPTION

                    // Build onto error matrices
                    nx.segment( epoch*localMatrixSize, localMatrixSize ) += observationRow * codeWeight * ( rangeData[ userSatellite ][ epochTime ][ prn ] - model_pseudorange );
                    nx.segment( epoch*localMatrixSize, localMatrixSize ) += observationRow * phaseWeight * ( phaseData[ userSatellite ][ epochTime ][ prn ] - model_phase );
                    nb( biasNumberLocal - 1 ) += -1 * phaseWeight * ( phaseData[ userSatellite ][ epochTime ][ prn ] - model_phase );

                    // Build onto normal equation matrices
                    NxxLocal += ( observationRow.transpose() * observationRow ) * ( phaseWeight + codeWeight );
                    Nbb.coeffRef( biasNumberLocal - 1, biasNumberLocal - 1 ) += phaseWeight;
                    Nxb.coeffRef( epoch*localMatrixSize + userIndex*numberOfUnknowns + 0, biasNumberLocal - 1 ) += -observationRow( userIndex*numberOfUnknowns + 0) * phaseWeight;
                    Nxb.coeffRef( epoch*localMatrixSize + userIndex*numberOfUnknowns + 1, biasNumberLocal - 1 ) += -observationRow( userIndex*numberOfUnknowns + 1) * phaseWeight;
                    Nxb.coeffRef( epoch*localMatrixSize + userIndex*numberOfUnknowns + 2, biasNumberLocal - 1 ) += -observationRow( userIndex*numberOfUnknowns + 2) * phaseWeight;
                    Nxb.coeffRef( epoch*localMatrixSize + localMatrixSize - 1, biasNumberLocal - 1 ) += -observationRow( localMatrixSize - 1 ) * phaseWeight; // Comment out for NO CLOCK
                }
            }
            userIndex += 1;
        }

        // Inverse the local Nxx for faster inversion
        Eigen::MatrixXd NxxLocal_inv( localMatrixSize, localMatrixSize );
        NxxLocal_inv = NxxLocal.inverse(); // Inverse of local N_XX 4x4 or 36x36 matrix
        for( int i = 0; i < localMatrixSize; i++ )
        {
            for( int j = 0; j < localMatrixSize; j++ )
            {
                Nxx_inv.coeffRef( epoch*localMatrixSize + i, epoch*localMatrixSize + j ) = NxxLocal_inv( i, j );
            }
        }
    }

    // Bias and Position updates
//    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver; // Other solver options
//    Eigen::SparseQR< Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver; // Other solver options
    Eigen::SparseQR< Eigen::SparseMatrix<double>, Eigen::AMDOrdering<int> > solver; // SparseQR
//    Eigen::SparseQR< Eigen::SparseMatrix<double>, Eigen::NaturalOrdering<int> > solver; // Other solver options
    solver.compute( Nbb - Nxb.transpose() * Nxx_inv * Nxb );
    Eigen::VectorXd deltaB = solver.solve( nb - Nxb.transpose() * Nxx_inv * nx );
    Eigen::VectorXd deltaX = Nxx_inv * ( nx - Nxb * deltaB );
    Eigen::VectorXd correction( deltaB.size() + deltaX.size() );
    correction << deltaX, deltaB;
    Eigen::VectorXd estimation = aPrioriEstimations + correction;

    #define FOLDINGSTART { // Create covariance matrices // This takes extremely long//    Eigen::SparseMatrix<double> Qbb( numberOfBiases, numberOfBiases );
//    Eigen::SparseMatrix<double> Ibb( numberOfBiases, numberOfBiases );
//    Ibb.setIdentity();
//    Qbb = solver.solve( Ibb );
//    Eigen::SparseMatrix<double> Qxx( matrixSize, matrixSize );
//    Qxx = ( Nxx_inv * Nxb ) * Qbb * ( Nxx_inv * Nxb ).transpose();
//    Qxx = Nxx_inv + Qxx;
#define FOLDINGEND }

    if( numberOfUsers > 0 ){ // Options for printing and saving intermediate results
        #define FOLDINGSTART { // Saving results
//        std::cout << "Nxx inv: \n" << Nxx_inv << std::endl;
//        std::cout << "Nxx inv: \n" << Eigen::MatrixXd(Nxx_inv) << std::endl;
//        writeMatrixToFile( aPrioriEstimations, "apri_"+setting+"_"+std::to_string(iteration)+".dat", 16, "C:/tudatBundle/tudatApplications/stepTwo/SimulationOutput/Nxx");
//        writeMatrixToFile( Eigen::MatrixXd(Nxx_inv), "Nxx_inv_"+setting+"_"+std::to_string(iteration)+".dat", 16, "C:/tudatBundle/tudatApplications/stepTwo/SimulationOutput/Nxx");
//        writeMatrixToFile( Eigen::MatrixXd(Nxb), "Nxb_"+setting+"_"+std::to_string(iteration)+".dat", 16, "C:/tudatBundle/tudatApplications/stepTwo/SimulationOutput/Nxx");
//        writeMatrixToFile( Eigen::MatrixXd(Nxb.transpose()), "NxbT_"+setting+"_"+std::to_string(iteration)+".dat", 16, "C:/tudatBundle/tudatApplications/stepTwo/SimulationOutput/Nxx");
//        writeMatrixToFile( nb, "nb_"+setting+"_"+std::to_string(iteration)+".dat", 16, "C:/tudatBundle/tudatApplications/stepTwo/SimulationOutput/Nxx");
//        writeMatrixToFile( nx, "nx_"+setting+"_"+std::to_string(iteration)+".dat", 16, "C:/tudatBundle/tudatApplications/stepTwo/SimulationOutput/Nxx");
//        writeMatrixToFile( deltaB, "deltaB_"+setting+"_"+std::to_string(iteration)+".dat", 16, "C:/tudatBundle/tudatApplications/stepTwo/SimulationOutput/Nxx");
//        writeMatrixToFile( deltaX, "deltaX_"+setting+"_"+std::to_string(iteration)+".dat", 16, "C:/tudatBundle/tudatApplications/stepTwo/SimulationOutput/Nxx");
//        writeMatrixToFile( estimation, "estimation_"+setting+"_"+std::to_string(iteration)+".dat", 16, "C:/tudatBundle/tudatApplications/stepTwo/SimulationOutput/Nxx");
    #define FOLDINGEND }

        #define FOLDINGSTART { // Condition number
//            Eigen::JacobiSVD<Eigen::MatrixXd> svd(normalMatrix);
//            double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
//            std::cout << "With code: " << pow( codeWeightISL, -0.5 ) << " and weight " << codeWeightISL << " op epoch: " << epochTime << " condition number: " << cond << std::endl;
//            if( cond > codeWeightISL/1.0e5 ){ std::cout << "Maybe" << std::endl; }
        #define FOLDINGEND }

        #define FOLDINGSTART { // Printing stuff
        std::cout << "\napriori: \n" << aPrioriEstimations.segment(0,28).transpose() << std::endl;
//        std::cout << "\nnb: \n" << nb.segment(0,16).transpose() << std::endl;
//        std::cout << "\nnx: \n" << nx.segment(0,16).transpose() << std::endl;
//        std::cout << "\nNbb: \n" << Nbb.block(0,0,16,16) << std::endl;
//        std::cout << "\nNxb: \n" << Nxb.block(0,0,16,16) << std::endl;
//        std::cout << "\nNxb Transpose: \n" << (Nxb.transpose()).block(0,0,16,16) << std::endl;
//        std::cout << "Nxx inv: \n" << Nxx_inv.block(0,0,16,16) << std::endl;
//        std::cout << "\ndeltaB: \n" << deltaB.segment(0,16).transpose() << std::endl;
//        std::cout << "\ndeltaX: \n" << deltaX.segment(0,16).transpose() << std::endl;
        std::cout << "\nestimation: \n" << estimation.segment(0,28).transpose() << std::endl;
        std::cout << "\nbias: \n" << estimation.tail(numberOfBiases).segment(0,12).transpose() << std::endl;
        #define FOLDINGEND }
    }

    return estimation;
}


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
        std::map< std::string, std::map< double, Eigen::Vector4d > > aPrioriPerUser,
        int maxIterationsKinematic,
        double correctionLimitKinematic,
        double codeWeight, double phaseWeight,
        double codeWeightISL, double phaseWeightISL,
        bool lighttimeCorrection,
        std::map< std::string, std::map< double, std::map< int, double > > > travelTimesPerUserPerEpoch,
        int solver, bool clockPhase)
{    
    // Determine constants
    long long numberOfEpochs = static_cast<signed>(epochTimes.size());
    long long numberOfUsers = static_cast<signed>(availabilityMatrixPerUser.size());
    std::vector< std::string > userNames;
    int numberOfUnknowns = 4;
    int localMatrixSize = numberOfUsers * numberOfUnknowns;

    if( numberOfUsers > 1 && clockPhase )
    {
        numberOfUnknowns = 3;
        localMatrixSize = numberOfUsers * numberOfUnknowns + 1;
//        localMatrixSize = numberOfUsers * numberOfUnknowns; // NO CLOCK ESTIMATION OPTION
    }
    int matrixSize = numberOfEpochs * localMatrixSize; //numberOfEpochs*numberOfUsers*numberOfUnknowns;

    // Create a vector for the a priori estimations and establish bias number matrices for all users
    Eigen::VectorXd aPrioriEstimations = Eigen::VectorXd::Zero( matrixSize );
    std::map< std::string, Eigen::MatrixXi > biasNumbersPerUser;
    int userIndex = 0;
    int previousBiasNumber = 0;
    for( std::map< std::string, Eigen::MatrixXi >::iterator itrUser = availabilityMatrixPerUser.begin(); itrUser != availabilityMatrixPerUser.end(); itrUser++ )
    {
        std::string userSatellite = itrUser->first;
        userNames.push_back( userSatellite );
        Eigen::MatrixXi availabilityMatrix = itrUser->second;
        int numberOfGNSS = static_cast<int>(availabilityMatrix.cols());
        Eigen::MatrixXi biasNumbers = Eigen::MatrixXi::Zero( numberOfEpochs, numberOfGNSS );
        for( long long epoch = 0; epoch < numberOfEpochs; epoch++ )
        {
            // Create apriori position and clock error estimation through point positioning
            double epochTime = epochTimes.at( static_cast<unsigned>(epoch) );
//            aPrioriEstimations.segment( epoch*numberOfUsers*4 + userIndex*4, 4 ) = aPrioriPerUser[ userSatellite ][ epochTime ];
//            aPrioriEstimations.segment( epoch*localMatrixSize + userIndex*numberOfUnknowns, 3 ) = aPrioriPerUser[ userSatellite ][ epochTime ].segment(0, 3);
            aPrioriEstimations.segment( epoch*localMatrixSize + userIndex*numberOfUnknowns, numberOfUnknowns ) = aPrioriPerUser[ userSatellite ][ epochTime ].segment(0, numberOfUnknowns);
            estimationsKinPerUserPerIteration[ 0 ][ userSatellite ][ epochTime ] = aPrioriPerUser[ userSatellite ][ epochTime ].segment(0, numberOfUnknowns);
            if( clockPhase ){
                aPrioriEstimations( epoch*localMatrixSize + localMatrixSize - 1 ) += aPrioriPerUser[ userSatellite ][ epochTime ]( numberOfUnknowns )/numberOfUsers; // Comment out for NO CLOCK
//                aPrioriEstimations( epoch*localMatrixSize + localMatrixSize - 1 ) += 0; // Uncomment for NO CLOCK
            }

            // Allocate bias numbers
            for( auto i = 0; i < numberOfGNSS; i++ )
            {
                if( availabilityMatrix( epoch, i ) == 0){continue;}
                else if( epoch > 0 && biasNumbers( epoch-1, i ) > 0)
                {
                    biasNumbers( epoch, i ) = biasNumbers( epoch-1, i );
                }
                else
                {
//                    biasNumbers( epoch, i ) = biasNumbers.maxCoeff() + 1;
                    biasNumbers( epoch, i ) = previousBiasNumber + 1;
                }
                previousBiasNumber = biasNumbers.maxCoeff();
            }
        }
        biasNumbersPerUser[ userSatellite ] = biasNumbers;
        userIndex += 1;
    }

    // Perform kinematic estimation iteration with while loop
    int iteration = 0;
    double correction = correctionLimitKinematic * 10;
    int numberOfBiases = previousBiasNumber;
    Eigen::VectorXd kinematicEstimation = aPrioriEstimations;
    Eigen::MatrixXd iterationHistory(matrixSize + numberOfBiases, 1);
    iterationHistory.block( 0, 0, matrixSize, 1 ) = kinematicEstimation;
    Eigen::VectorXd kinematicEstimationOld = iterationHistory;    
    while( iteration < maxIterationsKinematic && correction > correctionLimitKinematic )
    {
        // Select function for solver iteration
        if( solver == 1)
        {
            kinematicEstimation =      getKinematicIteration( epochTimes, satelliteEphemerisGNSS, kinematicEstimation, rangeData, phaseData, availabilityMatrixPerUser, biasNumbersPerUser, ephemerisError, arcTime, ephemerisErrorBase, codeWeight, phaseWeight, codeWeightISL, phaseWeightISL, lighttimeCorrection, travelTimesPerUserPerEpoch );
        }
        if( clockPhase )
        {
            if( numberOfUsers > 0 )
            {
                std::cout << "\nIteration " << iteration+1 << " of " << maxIterationsKinematic << " with clock in phase" << std::endl;
            }
            kinematicEstimation = getKinematicIterationPhase( epochTimes, satelliteEphemerisGNSS, kinematicEstimation, rangeData, phaseData, availabilityMatrixPerUser, biasNumbersPerUser, ephemerisError, arcTime, ephemerisErrorBase, codeWeight, phaseWeight, codeWeightISL, phaseWeightISL, lighttimeCorrection, travelTimesPerUserPerEpoch );
        }
        else if( solver == 2 )
        {
            if( numberOfUsers > 0 )
            {
                std::cout << "\nIteration " << iteration+1 << " of " << maxIterationsKinematic << std::endl;
            }
            kinematicEstimation = getKinematicIterationSmart( epochTimes, satelliteEphemerisGNSS, kinematicEstimation, rangeData, phaseData, availabilityMatrixPerUser, biasNumbersPerUser, ephemerisError, arcTime, ephemerisErrorBase, codeWeight, phaseWeight, codeWeightISL, phaseWeightISL, lighttimeCorrection, travelTimesPerUserPerEpoch );
        }

        Eigen::VectorXd estimationDifference = ( kinematicEstimationOld - kinematicEstimation ).cwiseAbs();
        correction = estimationDifference.maxCoeff();
        kinematicEstimationOld = kinematicEstimation;

        // Save the itereation history
        iterationHistory.conservativeResize(iterationHistory.rows(), iterationHistory.cols()+1);
        iterationHistory.col(iterationHistory.cols()-1) = kinematicEstimation;

        iteration += 1;
        if ( iteration == maxIterationsKinematic )
        {
            std::cout << "SKRT: " << correction << std::endl;
        }

        for( long long userIndex = 0; userIndex < numberOfUsers; userIndex++ )
        {
            std::string userSatellite = userNames.at( userIndex );
            for( long long epoch = 0; epoch < numberOfEpochs; epoch++ )
            {
                double epochTime = epochTimes.at( static_cast<unsigned>(epoch) );
                estimationsKinPerUserPerIteration[ iteration ][ userSatellite ][ epochTime ].segment(0, numberOfUnknowns) = kinematicEstimation.segment( epoch * localMatrixSize + userIndex * numberOfUnknowns, numberOfUnknowns );
                if( clockPhase ){ estimationsKinPerUserPerIteration[ iteration ][ userSatellite ][ epochTime ](3) = kinematicEstimation( epoch * localMatrixSize + localMatrixSize -1 ); }
            }
        }
    }

    // Process results to get estimation results in map and iteration history
    std::map< std::string, std::map< double, Eigen::Vector4d > > estimationsKinPerUser;
    for( long long userIndex = 0; userIndex < numberOfUsers; userIndex++ )
    {
        std::string userSatellite = userNames.at( userIndex );
        for( long long epoch = 0; epoch < numberOfEpochs; epoch++ )
        {
            double epochTime = epochTimes.at( static_cast<unsigned>(epoch) );
            estimationsKinPerUser[ userSatellite ][ epochTime ].segment(0, numberOfUnknowns) = kinematicEstimation.segment( epoch * localMatrixSize + userIndex * numberOfUnknowns, numberOfUnknowns );
            if( clockPhase ){ estimationsKinPerUser[ userSatellite ][ epochTime ](3) = kinematicEstimation( epoch * localMatrixSize + localMatrixSize - 1 ); }

            Eigen::VectorXd iterationsTemp = iterationHistory.block( epoch * localMatrixSize + userIndex * numberOfUnknowns, 0, numberOfUnknowns, 1 );
            biasIterations[ 0 ] = iterationHistory.col(0).tail(numberOfBiases);
            for( int it = 1; it < iteration+1; it++ )
            {
                iterationsTemp.conservativeResize( iterationsTemp.size() + numberOfUnknowns );
                iterationsTemp.segment( iterationsTemp.size()-numberOfUnknowns, numberOfUnknowns ) = iterationHistory.block( epoch*localMatrixSize + userIndex*numberOfUnknowns, it, numberOfUnknowns, 1 );

                biasIterations[ it ] = iterationHistory.col(it).tail(numberOfBiases);
            }
            allUsersKinIterations[ userSatellite ][ epochTime ] = iterationsTemp;
        }
    }

    return estimationsKinPerUser;
}
