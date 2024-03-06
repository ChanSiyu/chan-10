/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
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
    Foam::RELWallVelocityFvPatchVectorField

Description
    Wall-velocity condition to be used in conjunction with the single rotating
    frame (REL) model (see: FOAM::RELModel)

    The condition applies the appropriate rotation transformation in time and
    space to determine the local REL velocity of the wall.

        \f[
            U_p = - U_{p,REL}
        \f]

    where
    \vartable
        U_p     = patch velocity [m/s]
        U_{p,REL} = REL velocity
    \endvartable

    The normal component of \f$ U_p \f$ is removed to ensure 0 wall-flux even
    if the wall patch faces are irregular.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            RELWallVelocity;
        value           uniform (0 0 0);    // Initial value
    }
    \endverbatim

See also
    Foam::RELModel
    Foam::RELVelocityFvPatchVectorField
    Foam::fixedValueFvPatchField

SourceFiles
    RELWallVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef RELWallVelocityFvPatchVectorField_H
#define RELWallVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class RELWallVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class RELWallVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

public:

    //- Runtime type information
    TypeName("RELWallVelocity");


    // Constructors

        //- Construct from patch and internal field
        RELWallVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        RELWallVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given RELWallVelocityFvPatchVectorField
        //  onto a new patch
        RELWallVelocityFvPatchVectorField
        (
            const RELWallVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        RELWallVelocityFvPatchVectorField
        (
            const RELWallVelocityFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        RELWallVelocityFvPatchVectorField
        (
            const RELWallVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new RELWallVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchVectorField&, const labelList&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
