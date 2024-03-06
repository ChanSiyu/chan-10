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
    Foam::RELFreestreamVelocityFvPatchVectorField

Description
    Freestream velocity condition to be used in conjunction with the single
    rotating frame (REL) model (see: RELModel class)

    Given the free stream velocity in the absolute frame, the condition
    applies the appropriate rotation transformation in time and space to
    determine the local velocity using:

        \f[
            U_p = cos(\theta)*U_{Inf} + sin(theta) (n^UInf) - U_{p,rel}
        \f]

    where
    \vartable
        U_p     = patch velocity [m/s]
        U_{Inf} = free stream velocity in the absolute frame [m/s]
        theta   = swept angle [rad]
        n       = axis direction of the REL
        U_{p,REL} = REL velocity of the patch
    \endvartable


Usage
    \table
        Property        | Description               | Required | Default value
        UInf            | freestream velocity       | yes      |
        relative        | UInf relative to the rel? | no       |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            RELFreestreamVelocity;
        UInf            uniform (0 0 0);
        relative        no;
        value           uniform (0 0 0);    // initial value
    }
    \endverbatim

See also
    Foam::freestreamFvPatchField
    Foam::RELVelocityFvPatchVectorField

SourceFiles
    RELFreestreamVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef RELFreestreamVelocityFvPatchVectorField_H
#define RELFreestreamVelocityFvPatchVectorField_H

#include "inletOutletFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class RELFreestreamVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class RELFreestreamVelocityFvPatchVectorField
:
    public inletOutletFvPatchVectorField
{
    // Private Data

        //- Is the supplied inlet value relative to the rel
        Switch relative_;

        //- Velocity of the free stream in the absolute frame [m/s]
        vector UInf_;


public:

    //- Runtime type information
    TypeName("RELFreestreamVelocity");


    // Constructors

        //- Construct from patch and internal field
        RELFreestreamVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        RELFreestreamVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given RELFreestreamVelocityFvPatchVectorField
        //  onto a new patch
        RELFreestreamVelocityFvPatchVectorField
        (
            const RELFreestreamVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        RELFreestreamVelocityFvPatchVectorField
        (
            const RELFreestreamVelocityFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        RELFreestreamVelocityFvPatchVectorField
        (
            const RELFreestreamVelocityFvPatchVectorField&,
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
                new RELFreestreamVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

        // Access

            //- Return the velocity at infinity
            const vector& UInf() const
            {
                return UInf_;
            }

            //- Return reference to the velocity at infinity to allow adjustment
            vector& UInf()
            {
                return UInf_;
            }


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
