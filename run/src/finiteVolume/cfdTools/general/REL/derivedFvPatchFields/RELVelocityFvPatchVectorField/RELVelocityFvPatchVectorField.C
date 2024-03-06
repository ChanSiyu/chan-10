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

#include "RELVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

#include "RELModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RELVelocityFvPatchVectorField::RELVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    relative_(0),
    inletValue_(p.size(), Zero)
{}


Foam::RELVelocityFvPatchVectorField::RELVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    relative_(dict.lookup("relative")),
    inletValue_("inletValue", dict, p.size())
{}


Foam::RELVelocityFvPatchVectorField::RELVelocityFvPatchVectorField
(
    const RELVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    relative_(ptf.relative_),
    inletValue_(mapper(ptf.inletValue_))
{}


Foam::RELVelocityFvPatchVectorField::RELVelocityFvPatchVectorField
(
    const RELVelocityFvPatchVectorField& relvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(relvpvf, iF),
    relative_(relvpvf.relative_),
    inletValue_(relvpvf.inletValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::RELVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    m(*this, *this);
    m(inletValue_, inletValue_);
}


void Foam::RELVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const RELVelocityFvPatchVectorField& tiptf =
        refCast<const RELVelocityFvPatchVectorField>(ptf);

    inletValue_.rmap(tiptf.inletValue_, addr);
}


void Foam::RELVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // If not relative to the REL include the effect of the REL
    if (!relative_)
    {
        // Get reference to the REL model
        const REL::RELModel& REL =
            db().lookupObject<REL::RELModel>("RELProperties");

        // Determine patch velocity due to REL
        const vectorField RELVelocity(REL.velocity(patch().Cf()));

        operator==(-RELVelocity + inletValue_);
    }
    // If already relative to the REL simply supply the inlet value as a fixed
    // value
    else
    {
        operator==(inletValue_);
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::RELVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "relative", relative_);
    writeEntry(os, "inletValue", inletValue_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        RELVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
