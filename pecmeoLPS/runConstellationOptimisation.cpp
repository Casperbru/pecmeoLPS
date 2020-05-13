/*    Copyright (c) 2010-2017, Delft University of Technology
*    All rigths reserved
*
*    This file is part of the Tudat. Redistribution and use in source and
*    binary forms, with or without modification, are permitted exclusively
*    under the terms of the Modified BSD license. You should have received
*    a copy of the license with this file. If not, please or visit:
*    http://tudat.tudelft.nl/LICENSE.
*/

#include <iostream>
#include <fstream>

#include <boost/filesystem.hpp>

#include "pagmo/algorithms/de1220.hpp"
#include "pagmo/algorithms/sade.hpp"
#include "pagmo/algorithms/cmaes.hpp"
#include "pagmo/algorithms/moead.hpp"
#include "pagmo/algorithms/nsga2.hpp"
#include "pagmo/algorithms/ihs.hpp"
#include "constellationOpt.h"
#include "C:\tudatBundle\tudatExampleApplications\libraryExamples\PaGMOEx\Problems\applicationOutput.h"
#include "C:\tudatBundle\tudatExampleApplications\libraryExamples\PaGMOEx\Problems\getAlgorithm.h"
#include "C:\tudatBundle\tudatExampleApplications\libraryExamples\PaGMOEx\Problems\saveOptimizationResults.h"
//#include "Problems/earthMarsTransfer.h"
//#include "Problems/applicationOutput.h"
//#include "Problems/getAlgorithm.h"
//#include "Problems/saveOptimizationResults.h"

using namespace tudat_pagmo_applications;

int main( )
{

    //Set seed for reproducible results
    pagmo::random_device::set_seed( 123 );

    // We have five decision variables each with a lower and upper
    // bound, create a vector of vectors that will contain these.
    unsigned long long numberOfParameters = 5;
    std::vector< std::vector< double > > bounds( 2, std::vector< double >( numberOfParameters, 0.0 ) );

    // Define bounds: Search between 2020 and 2025 for flight duration between 200 and 1000 days.
    bounds[ 0 ][ 0 ] = 10000.0e3;
    bounds[ 1 ][ 0 ] = 23000.0e3;
    bounds[ 0 ][ 1 ] = 0;
    bounds[ 1 ][ 1 ] = 359;
    bounds[ 0 ][ 2 ] = 0;
    bounds[ 1 ][ 2 ] = 120;
    bounds[ 0 ][ 3 ] = 0;
    bounds[ 1 ][ 3 ] = 120;
    bounds[ 0 ][ 4 ] = 0;
    bounds[ 1 ][ 4 ] = 120;

    // Create object to compute the problem fitness
    problem prob{constellationOptimiser ( bounds )};

    // Solve problem using x different optimizers
    for( int j = 0; j < 1; j++ )
    {
        // Retrieve MO algorithm
        int algorithmIndex = j;
        algorithm algo{getAlgorithm( algorithmIndex )};
//        algorithm algo{getMultiObjectiveAlgorithm( algorithmIndex )};

        // Create an island with 1024 individuals
        island isl{algo, prob, 1024};

        // Evolve for 25 generations
        for( int i = 0 ; i < 50; i++ )
        {
            isl.evolve( );
            while( isl.status( ) != pagmo::evolve_status::idle &&
                   isl.status( ) != pagmo::evolve_status::idle_error )
            {
                isl.wait( );
            }
            isl.wait_check( ); // Raises errors

            // Write current iteration results to file
            printPopulationToFile( isl.get_population( ).get_x( ), "PECMEOTest" + std::to_string( i ) + "_" + std::to_string( j ), false );
            printPopulationToFile( isl.get_population( ).get_f( ), "PECMEOTest" + std::to_string( i ) + "_" + std::to_string( j ), true );

            std::cout<<i<<" "<<j<<std::endl;


        }
    }

    return 0;

}
