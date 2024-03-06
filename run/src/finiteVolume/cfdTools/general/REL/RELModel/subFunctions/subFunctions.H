/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::REL::subFuntions

Description
    

SourceFiles
    subFunctions.C

\*---------------------------------------------------------------------------*/

#ifndef subFunctions_H
#define subFunctions_H

#include "RELModel.H"

#include "autoPtr.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace REL
{

/*---------------------------------------------------------------------------*\
                            Class subFunctions Declaration
\*---------------------------------------------------------------------------*/

class subFunctions
:
    public RELModel
{
private:
    // Private Data

        //- Angular Velocity (rad/s)
        autoPtr<Function1<vector> > subAngularVelocity_;

        //- Translational Velocity (m/s)
        autoPtr<Function1<vector> > subTranslationalVelocity_;

public:

    //- Runtime type information
    TypeName("subFunctions");


    // Constructors

        //- Construct from components
        subFunctions(const volVectorField& U);

        //- Disallow default bitwise copy construction
        subFunctions(const subFunctions&) = delete;


    //- Destructor
    ~subFunctions();


    // Member Functions

        // I-O

            //- Read
            bool read();

            //- Source term component for momentum equation
            tmp<volVectorField::Internal> Su();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const subFunctions&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace REL
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
