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
#include "RELVelocityFvPatchVectorField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace REL
    {
        defineTypeNameAndDebug(RELModel, 0);
        defineRunTimeSelectionTable(RELModel, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

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
            "RELProperties",
            Urel.time().constant(),
            Urel.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    Urel_(Urel),
    mesh_(Urel_.mesh()),
    origin_("origin", dimLength, lookup("origin")),
    RELModelCoeffs_(optionalSubDict(type + "Coeffs")),
    angle_(dimensionedVector("angle", dimless, Zero)),
    angularVelocity_(dimensionedVector("angularVelocity", dimless/dimTime, Zero)),
    angularAcceleration_(dimensionedVector("angularAcceleration", dimless/(dimTime*dimTime), Zero)),
    translationalVelocity_(dimensionedVector("translationalVelocity", dimLength/dimTime, Zero)),
    translationalAcceleration_(dimensionedVector("translationalAcceleration", dimLength/(dimTime*dimTime), Zero))
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


Foam::tmp<Foam::volVectorField> Foam::REL::RELModel::Uabs() const
{
    tmp<volVectorField> Urel = U();

    tmp<volVectorField> tUabs
    (
        volVectorField::New("Uabs", Urel)
    );

    volVectorField& Uabs = tUabs.ref();

    // Add REL contribution to internal field
    Uabs.primitiveFieldRef() += Urel_.primitiveField();

    // Add Urel boundary contributions
    volVectorField::Boundary& Uabsbf = Uabs.boundaryFieldRef();
    const volVectorField::Boundary& bvf = Urel_.boundaryField();

    forAll(bvf, i)
    {
        if (isA<RELVelocityFvPatchVectorField>(bvf[i]))
        {
            // Only include relative contributions from
            // RELVelocityFvPatchVectorField's
            const RELVelocityFvPatchVectorField& UrelPatch =
                refCast<const RELVelocityFvPatchVectorField>(bvf[i]);
            if (UrelPatch.relative())
            {
                Uabsbf[i] += Urel_.boundaryField()[i];
            }
        }
        else
        {
            Uabsbf[i] += Urel_.boundaryField()[i];
        }
    }

    return tUabs;
}


// ************************************************************************* //
