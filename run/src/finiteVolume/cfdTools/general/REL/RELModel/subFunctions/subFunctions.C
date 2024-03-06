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

\*---------------------------------------------------------------------------*/

#include "subFunctions.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

#include "fvCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace REL
    {
        defineTypeNameAndDebug(subFunctions, 0);

        addToRunTimeSelectionTable
        (
            RELModel,
            subFunctions,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::REL::subFunctions::subFunctions
(
    const volVectorField& U
)
:
    RELModel(typeName, U),
    subAngularVelocity_(Function1<vector>::New("angularVelocity", RELModelCoeffs_)),
    subTranslationalVelocity_(Function1<vector>::New("translationalVelocity", RELModelCoeffs_))
{
    // Initialise the angular velocity
    const vector subAngularVelocity0 = subAngularVelocity_->value(0);
    angularVelocity_.value() = subAngularVelocity0;

    const vector subTranslationalVelocity0 = subTranslationalVelocity_->value(0);
    translationalVelocity_.value() = subTranslationalVelocity0;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::REL::subFunctions::~subFunctions()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::REL::subFunctions::read()
{
    if (RELModel::read())
    {
        // Re-read angular velocity
        subAngularVelocity_.clear();
        subAngularVelocity_ = autoPtr<Function1<vector> >
        (
            Function1<vector>::New("angularVelocity", RELModelCoeffs_)
        );

        // Re-read translational velocity
        subTranslationalVelocity_.clear();
        subTranslationalVelocity_ = autoPtr<Function1<vector> >
        (
            Function1<vector>::New("translationalVelocity", RELModelCoeffs_)
        );
        return true;
    }
    else
    {
        return false;
    }
}

Foam::tmp<Foam::DimensionedField<Foam::vector, Foam::volMesh> >
Foam::REL::subFunctions::Su()
{
    scalar t  = db().time().value();
    scalar dt = db().time().deltaT().value();

    vector subAngularVelocity_t0 = angularVelocity_.value();
    vector subAngularVelocity_t1 = subAngularVelocity_->value(t);
    angularAcceleration_.value() = (subAngularVelocity_t1 - subAngularVelocity_t0) / dt;
    angularVelocity_.value() = subAngularVelocity_t1;
    Info<< "   Last angular velocity     value:          " << subAngularVelocity_t0 << endl;
    Info<< "Current angular velocity     value:          " << subAngularVelocity_t1 << endl;
    Info<< "Current angular Acceleration value:          " << angularAcceleration_.value() << endl;

    vector subTranslationalVelocity_t0 = translationalVelocity_.value();
    vector subTranslationalVelocity_t1 = subTranslationalVelocity_->value(t);
    translationalAcceleration_.value() = (subTranslationalVelocity_t1 - subTranslationalVelocity_t0) / dt;
    translationalVelocity_.value() = subTranslationalVelocity_t1;
    Info<< "   Last translational velocity     value:    " << subTranslationalVelocity_t0 << endl;
    Info<< "Current translational velocity     value:    " << subTranslationalVelocity_t1 << endl;
    Info<< "Current translational acceleration value:    " << translationalAcceleration_.value() << endl;

    if (mag(subAngularVelocity_t1) != 0)
    {
        axis_ = subAngularVelocity_t1;
        axis_ /= mag(axis_);
    }
    vector subAngle = angle().value();
    angle_.value() = subAngle + (subAngularVelocity_t0 + subAngularVelocity_t1) * 0.5 * dt;
    Info<< "Current axis value:     " << axis_ << endl;
    Info<< "Current angle value:    " << angle_.value() << endl;
    Info<< " "<< endl;

    return Fcoriolis() + Fcentrifugal() + Feular();
}

// ************************************************************************* //
