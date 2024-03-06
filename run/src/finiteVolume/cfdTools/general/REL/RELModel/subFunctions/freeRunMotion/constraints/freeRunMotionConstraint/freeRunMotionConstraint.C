/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "freeRunMotionConstraint.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(freeRunMotionConstraint, 0);
    defineRunTimeSelectionTable(freeRunMotionConstraint, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::freeRunMotionConstraint::freeRunMotionConstraint
(
    const word& name,
    const dictionary& sDoFRBMCDict,
    const freeRunMotion& motion
)
:
    name_(name),
    sDoFRBMCCoeffs_(sDoFRBMCDict),
    motion_(motion)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::freeRunMotionConstraint::~freeRunMotionConstraint()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::freeRunMotionConstraint::read
(
    const dictionary& sDoFRBMCDict
)
{
    sDoFRBMCCoeffs_ = sDoFRBMCDict;

    return true;
}


void Foam::freeRunMotionConstraint::write(Ostream& os) const
{}


// ************************************************************************* //
