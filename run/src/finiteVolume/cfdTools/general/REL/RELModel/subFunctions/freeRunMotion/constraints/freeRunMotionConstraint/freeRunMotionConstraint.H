/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::freeRunMotionConstraints

Description
    Namespace for six DoF motion constraints


Class
    Foam::freeRunMotionConstraint

Description
    Base class for defining constraints for sixDoF motions

SourceFiles
    freeRunMotionConstraint.C
    freeRunMotionConstraintNew.C

\*---------------------------------------------------------------------------*/

#ifndef freeRunMotionConstraint_H
#define freeRunMotionConstraint_H

#include "Time.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "point.H"
#include "pointConstraint.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace REL
{

// Forward declaration of classes
class freeRunMotion;


/*---------------------------------------------------------------------------*\
               Class freeRunMotionConstraint Declaration
\*---------------------------------------------------------------------------*/

class freeRunMotionConstraint
{

protected:

    // Protected data

        //- Name of the constraint
        word name_;

        //- Constraint model specific coefficient dictionary
        dictionary sDoFRBMCCoeffs_;

        //- Reference to the body motion
        const freeRunMotion& motion_;


public:

    //- Runtime type information
    TypeName("freeRunMotionConstraint");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            freeRunMotionConstraint,
            dictionary,
            (
                const word& name,
                const dictionary& sDoFRBMCDict,
                const freeRunMotion& motion
            ),
            (name, sDoFRBMCDict, motion)
        );


    // Constructors

        //- Construct from the sDoFRBMCDict dictionary and Time
        freeRunMotionConstraint
        (
            const word& name,
            const dictionary& sDoFRBMCDict,
            const freeRunMotion& motion
        );

        //- Construct and return a clone
        virtual autoPtr<freeRunMotionConstraint> clone() const = 0;


    // Selectors

        //- Select constructed from the sDoFRBMCDict dictionary and Time
        static autoPtr<freeRunMotionConstraint> New
        (
            const word& name,
            const dictionary& sDoFRBMCDict,
            const freeRunMotion& motion
        );


    //- Destructor
    virtual ~freeRunMotionConstraint();


    // Member Functions

        //- Return the name
        const word& name() const
        {
            return name_;
        }

        //- Set the centre of rotation if not the centre of mass
        virtual void setCentreOfRotation(point&) const
        {}

        //- Apply and accumulate translational constraints
        virtual void constrainTranslation(pointConstraint&) const = 0;

        //- Apply and accumulate rotational constraints
        virtual void constrainRotation(pointConstraint&) const = 0;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& sDoFRBMCDict);

        // Access

            // Return access to sDoFRBMCCoeffs
            inline const dictionary& coeffDict() const
            {
                return sDoFRBMCCoeffs_;
            }

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace REL
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
