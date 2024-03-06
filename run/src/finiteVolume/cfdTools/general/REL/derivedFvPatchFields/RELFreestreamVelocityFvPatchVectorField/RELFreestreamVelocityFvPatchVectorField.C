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

#include "RELFreestreamVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "RELModel.H"
#include "steadyStateDdtScheme.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RELFreestreamVelocityFvPatchVectorField::
RELFreestreamVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    inletOutletFvPatchVectorField(p, iF),
    relative_(false),
    UInf_(Zero)
{}


Foam::RELFreestreamVelocityFvPatchVectorField::
RELFreestreamVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchVectorField(p, iF),
    relative_(dict.lookupOrDefault("relative", false)),
    UInf_(dict.lookup("UInf"))
{
    this->phiName_ = dict.lookupOrDefault<word>("phi","phi");

    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}

Foam::RELFreestreamVelocityFvPatchVectorField::
RELFreestreamVelocityFvPatchVectorField
(
    const RELFreestreamVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchVectorField(ptf, p, iF, mapper),
    relative_(ptf.relative_),
    UInf_(ptf.UInf_)
{}


Foam::RELFreestreamVelocityFvPatchVectorField::
RELFreestreamVelocityFvPatchVectorField
(
    const RELFreestreamVelocityFvPatchVectorField& relvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    inletOutletFvPatchVectorField(relvpvf, iF),
    relative_(relvpvf.relative_),
    UInf_(relvpvf.UInf_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::RELFreestreamVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get reference to the REL model
    const REL::RELModel& REL =
        db().lookupObject<REL::RELModel>("RELProperties");

    word ddtScheme
    (
        this->internalField().mesh()
       .schemes().ddt(this->internalField().name())
    );

    if (ddtScheme == fv::steadyStateDdtScheme<scalar>::typeName)
    {
        // If not relative to the REL include the effect of the REL
        if (!relative_)
        {
            refValue() = UInf_ - REL.velocity(patch().Cf());
        }
        // If already relative to the REL simply supply the inlet value
        // as a fixed value
        else
        {
            refValue() = UInf_;
        }
    }
    else
    {
        scalar theta = mag(REL.angle().value());

        refValue() =
            cos(theta)*UInf_ + sin(theta)*(REL.axis() ^ UInf_)
            - REL.velocity(patch().Cf());
    }
  
    // Set the inlet-outlet choice based on the direction of the freestream
    valueFraction() = 1.0 - pos0(refValue() & patch().Sf());

    mixedFvPatchField<vector>::updateCoeffs();
}


void Foam::RELFreestreamVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "relative", relative_);
    writeEntry(os, "UInf", UInf_);
    writeEntry(os, "phi", this->phiName_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        RELFreestreamVelocityFvPatchVectorField
    );
}


// ************************************************************************* //
