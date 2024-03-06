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

#include "RELModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
// Define type name and debugging options for the RELModel class
namespace Foam
{
    namespace REL
    {
        defineTypeNameAndDebug(RELModel, 0);
        defineRunTimeSelectionTable(RELModel, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
// Initialize the RELModel class with relevant parameters
Foam::REL::RELModel::RELModel
(
    const word& type,
    const volVectorField& Urel
)
:
    IOdictionary
    (
        IOobject
        (
            "RELProperties",                    // Object name
            Urel.time().constant(),             // Associated with time
            Urel.db(),                          // Database
            IOobject::MUST_READ_IF_MODIFIED,    // Must read if modified
            IOobject::NO_WRITE                  // Not writable
        )
    ),
    Urel_(Urel),                                                                                                        // Store the relative velocity field
    mesh_(Urel_.mesh()),                                                                                                // Store the mesh of the relative velocity field
    origin_("origin", dimLength, lookup("origin")),                                                                     // Store the rotation origin
    RELModelCoeffs_(optionalSubDict(type + "Coeffs")),                                                                  // Store the dictionary of model coefficients
    angle_(dimensionedVector("angle", dimless, Zero)),                                                                  // Store angle information
    angularVelocity_(dimensionedVector("angularVelocity", dimless/dimTime, Zero)),                                      // Store angular velocity
    angularAcceleration_(dimensionedVector("angularAcceleration", dimless / (dimTime * dimTime), Zero)),                // Store angular acceleration
    translationalVelocity_(dimensionedVector("translationalVelocity", dimLength / dimTime, Zero)),                      // Store translational velocity
    translationalAcceleration_(dimensionedVector("translationalAcceleration", dimLength / (dimTime * dimTime), Zero))   // Store translational acceleration
{
    // Normalise the axis
    axis_ = vector::zero;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::REL::RELModel::~RELModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::REL::RELModel::read()
{
    if (regIOobject::read())
    {
        // Re-read origin
        lookup("origin") >> origin_;

        // Re-read sub-model coeffs
        RELModelCoeffs_ = optionalSubDict(type() + "Coeffs");

        return true;
    }
    else
    {
        return false;
    }
}


const Foam::dimensionedVector& Foam::REL::RELModel::origin() const
{
    return origin_;
}


const Foam::vector& Foam::REL::RELModel::axis() const
{
    return axis_;
}

const Foam::dimensionedVector& Foam::REL::RELModel::angle() const
{
    return angle_;
}

const Foam::dimensionedVector& Foam::REL::RELModel::angularVelocity() const
{
    return angularVelocity_;
}

const Foam::dimensionedVector& Foam::REL::RELModel::angularAcceleration() const
{
    return angularAcceleration_;
}

const Foam::dimensionedVector& Foam::REL::RELModel::translationalVelocity() const
{
    return translationalVelocity_;
}

const Foam::dimensionedVector& Foam::REL::RELModel::translationalAcceleration() const
{
    return translationalAcceleration_;
}

Foam::tmp<Foam::DimensionedField<Foam::vector, Foam::volMesh>>
Foam::REL::RELModel::Fcoriolis() const
{
    return volVectorField::Internal::New
    (
        "Fcoriolis",
        2.0*angularVelocity_ ^ Urel_
    );
}


Foam::tmp<Foam::DimensionedField<Foam::vector, Foam::volMesh>>
Foam::REL::RELModel::Fcentrifugal() const
{
    return volVectorField::Internal::New
    (
        "Fcentrifugal",
        translationalAcceleration_ + (angularVelocity_ ^ (angularVelocity_ ^ (mesh_.C() - origin_)))
    );
}

Foam::tmp<Foam::DimensionedField<Foam::vector, Foam::volMesh>>
Foam::REL::RELModel::Feular() const
{
    return volVectorField::Internal::New
    (
        "Feular",
        angularAcceleration_ ^ (mesh_.C() - origin_)
    );
}


Foam::tmp<Foam::DimensionedField<Foam::vector, Foam::volMesh>>
Foam::REL::RELModel::Su()
{
    return Fcoriolis() + Fcentrifugal() + Feular();
}


Foam::vectorField Foam::REL::RELModel::velocity
(
    const vectorField& positions
) const
{
    tmp<vectorField> tfld =
        (angularVelocity_.value()
      ^ (
            (positions - origin_.value())
          - axis_*(axis_ & (positions - origin_.value()))
        ))
     + translationalVelocity_.value();

    return tfld();
}


Foam::tmp<Foam::volVectorField> Foam::REL::RELModel::U() const
{
    return volVectorField::New
    (
        "Urel",
        (angularVelocity_ ^ ((mesh_.C() - origin_) - axis_*(axis_ & (mesh_.C() - origin_)))) + translationalVelocity_
    );
}

// ************************************************************************* //
