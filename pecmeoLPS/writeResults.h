#ifndef WRITERESULTS_H
#define WRITERESULTS_H
#include <memory>
#include <Tudat/SimulationSetup/tudatEstimationHeader.h>


void writeResults(
        const unsigned int numberOfSatellitesPEC,
        unsigned int numberOfSatellitesTotal,
        std::vector< std::string > satelliteNames,
        std::string includedNames,
        std::string outputSubFolder,
        std::vector< std::map< double, Eigen::VectorXd > > allSatellitesPropagationHistory,
        std::vector< std::map< double, Eigen::VectorXd > > allSatellitesPropagationHistoryKep,
        std::map< std::string, std::map< double, Eigen::Vector6d > > allSatellitesDOPHistory,
        bool writeDOPResults,
        bool writeCartesianResult,
        bool writeKeplerianResult);


//class writeResults
//{
//public:
//    writeResults();
//};

#endif // WRITERESULTS_H
