#ifndef CONSTELLATIONOPT_H
#define CONSTELLATIONOPT_H

#include <vector>
#include <utility>
#include <limits>

#include <Eigen/Core>

#include <Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/convertMeanToEccentricAnomalies.h>
#include <Tudat/Astrodynamics/MissionSegments/multiRevolutionLambertTargeterIzzo.h>

#include "pagmo/island.hpp"
#include "pagmo/io.hpp"
#include "pagmo/serialization.hpp"
#include "pagmo/problem.hpp"

/*!
 *  The class defined in this file is to be used in a Pagmo optimization. It defines the objective function for an Earth-Mars
 *  two-burn impulsive transfer, using a Lambert targeter. The independent variables are:
 *
 *  1) Departure Julian day
 *  2) Time-of-flight of Earth-Mars transfer
 *
 *  The problem minimized the Delta V, and if requested, minimizes the time of flight
 */

using namespace pagmo;

//! Test function for a new interplanetary trajectory class in Tudat
struct constellationOptimiser
{

    typedef Eigen::Matrix< double, 6, 1 > StateType;

    //! Default constructor, required for Pagmo compatibility
    constellationOptimiser( ): useParamaterX_( false ){ }

    //! Constructor that sets boundaries of independent variables, and a boolean denoting whether the fitness is single-objective
    //! (GDOP), or dual objective (GDOP and ParameterX).
    constellationOptimiser( std::vector< std::vector< double > > &bounds, const bool useParamaterX = false );

    //! Calculate the fitness as a function of the parameter vector x
    std::vector< double > fitness( const std::vector< double > &x ) const;

    //! Retrieve the allowable limits of the parameter vector x: pair containing minima and maxima of parameter values
    std::pair< std::vector< double >, std::vector< double > > get_bounds() const;

    //! Retrieve the name of the problem
    std::string get_name( ) const;

    //! Serialization function for Pagmo compatibility
    template <typename Archive>
    void serialize(Archive &ar)
    {
        ar(problemBounds_);
    }

    //! Retrieve the number of objectives in problem, e.g. the size of the vector returned by the fitness function
    vector_double::size_type get_nobj() const
    {
        if(useParamaterX_ )
        {
            return 2u;
        }
        else
        {
            return 1u;
        }
    }

private:

    const std::vector< std::vector< double > > problemBounds_;

    bool useParamaterX_;
};

#endif // CONSTELLATIONOPT_H
