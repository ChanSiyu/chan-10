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

Namespace
    Foam::REL

Description
    Namespace for single relative frame models

Class
    Foam::REL::RELModel

Description
    Top level model for single relative frame

SourceFiles
    RELModel.C

\*---------------------------------------------------------------------------*/

#ifndef RELModel_H
#define RELModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "fvMesh.H"
#include "volFields.H"
#include "vectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace REL
{

/*---------------------------------------------------------------------------*\
                           Class RELModel Declaration
\*---------------------------------------------------------------------------*/

class RELModel
:
    public IOdictionary
{

protected:

    // Protected data

        //- Reference to the relative velocity field
        const volVectorField& Urel_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Origin of the axis
        dimensionedVector origin_;

        //- Axis of rotation, a direction vector which passes through the origin
        vector axis_;

        //- REL model coefficients dictionary
        dictionary RELModelCoeffs_;

        //- Swept angle of the frame (rad)
        dimensionedVector angle_;

        //- Angular velocity of the frame (rad/s)
        dimensionedVector angularVelocity_;

        //- Angular acceleration of the frame (rad/s^2)
        dimensionedVector angularAcceleration_;

        //- Translational Velocity of the frame (m/s)
        dimensionedVector translationalVelocity_;

        //- Translational acceleration of the frame (rad/s^2)
        dimensionedVector translationalAcceleration_;

public:

    //- Runtime type information
    TypeName("RELModel");


    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             RELModel,
             dictionary,
             (
                 const volVectorField& Urel
             ),
             (Urel)
         );


    // Constructors

        //- Construct from components
        RELModel
        (
            const word& type,
            const volVectorField& Urel
        );

        //- Disallow default bitwise copy construction
        RELModel(const RELModel&) = delete;


    // Selectors

         //- Return a reference to the selected REL model
         static autoPtr<RELModel> New
         (
             const volVectorField& Urel
         );


    //- Destructor
    virtual ~RELModel();


    // Member Functions

        // Edit

            //- Read radiationProperties dictionary
            virtual bool read();


        // Access

            //- Return the origin of rotation
            const dimensionedVector& origin() const;

            //- Return the axis of rotation
            const vector& axis() const;

            //- Return the swept angle of rotation [rad]
            const dimensionedVector& angle() const;

            //- Return the angular velocity field [rad/s]
            const dimensionedVector& angularVelocity() const;

            //- Return the angular acceleration field [rad/s^2]
            const dimensionedVector& angularAcceleration() const;

            //- Return the translational velocity field [m/s]
            const dimensionedVector& translationalVelocity() const;

            //- Return the translational acceleration field [m/s^2]
            const dimensionedVector& translationalAcceleration() const;

            //- Return the acceleration field [m/s]
            const dimensionedVector& Uframe() const;

            //- Return the coriolis force
            tmp<volVectorField::Internal> Fcoriolis() const;

            //- Return the centrifugal force
            tmp<volVectorField::Internal> Fcentrifugal() const;

            //- Return the Eular force
            tmp<volVectorField::Internal> Feular() const;

            //- Source term component for momentum equation
            virtual tmp<volVectorField::Internal> Su();

            //- Return velocity vector from positions
            vectorField velocity(const vectorField& positions) const;

            //- Return velocity of REL for complete mesh
            tmp<volVectorField> U() const;

            //- Return absolute velocity for complete mesh
            tmp<volVectorField> Uabs() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const RELModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace REL
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
