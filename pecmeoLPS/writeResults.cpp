#include "writeResults.h"
#include <memory>
#include <Tudat/SimulationSetup/tudatEstimationHeader.h>
#include "applicationOutput.h"

using namespace tudat;
using namespace tudat::input_output;


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
        bool writeKeplerianResult)
{
    for ( unsigned int i = 0; i < numberOfSatellitesTotal; i++ )
    {
        if( writeDOPResults && i < numberOfSatellitesPEC )
        {
            // Set filename for output data.
            std::stringstream outputFilename;
            outputFilename << satelliteNames.at( i ) << "DOPResults" << includedNames << ".dat";

            // Write propagation history to file.
            writeDataMapToTextFile( allSatellitesDOPHistory[ satelliteNames.at( i ) ],
                                    outputFilename.str( ),
                                    tudat_applications::getOutputPath( ) + outputSubFolder,
                                    "",
                                    std::numeric_limits< double >::digits10,
                                    std::numeric_limits< double >::digits10,
                                    "," );
        }

        if( writeCartesianResult )
        {
            // Set filename for output data.
            std::stringstream outputFilename;
            outputFilename << satelliteNames.at( i ) << ".dat";

            // Write propagation history to file.
            writeDataMapToTextFile( allSatellitesPropagationHistory.at( i ),
                                    outputFilename.str( ),
                                    tudat_applications::getOutputPath( ) + outputSubFolder,
                                    "",
                                    std::numeric_limits< double >::digits10,
                                    std::numeric_limits< double >::digits10,
                                    "," );
        }

        if( writeKeplerianResult )
        {
            // Set filename for output data. // FOR KEPLERIAN
            std::stringstream outputFilename;
            outputFilename << satelliteNames.at( i ) << "Kep.dat";

            // Write propagation history to file.
            writeDataMapToTextFile( allSatellitesPropagationHistoryKep.at( i ),
                                    outputFilename.str( ),
                                    tudat_applications::getOutputPath( ) + outputSubFolder,
                                    "",
                                    std::numeric_limits< double >::digits10,
                                    std::numeric_limits< double >::digits10,
                                    "," );
        }
    }
}

//writeResults::writeResults()
//{

//}
