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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::tensor Foam::REL::freeRunMotion::rotationTensorX
(
    scalar phi
) const
{
    return tensor
    (
        1, 0, 0,
        0, Foam::cos(phi), -Foam::sin(phi),
        0, Foam::sin(phi), Foam::cos(phi)
    );
}


inline Foam::tensor Foam::REL::freeRunMotion::rotationTensorY
(
    scalar phi
) const
{
    return tensor
    (
        Foam::cos(phi), 0, Foam::sin(phi),
        0, 1, 0,
        -Foam::sin(phi), 0, Foam::cos(phi)
    );
}


inline Foam::tensor Foam::REL::freeRunMotion::rotationTensorZ
(
    scalar phi
) const
{
    return tensor
    (
        Foam::cos(phi), -Foam::sin(phi), 0,
        Foam::sin(phi), Foam::cos(phi), 0,
        0, 0, 1
    );
}


inline Foam::Tuple2<Foam::tensor, Foam::vector>
Foam::REL::freeRunMotion::rotate
(
    const tensor& Q0,
    const vector& pi0,
    const scalar deltaT
) const
{
    Tuple2<tensor, vector> Qpi(Q0, pi0);
    tensor& Q = Qpi.first();
    vector& pi = Qpi.second();

    tensor R = rotationTensorX(0.5*deltaT*pi.x()/momentOfInertia_.xx());
    pi = pi & R;
    Q = Q & R;

    R = rotationTensorY(0.5*deltaT*pi.y()/momentOfInertia_.yy());
    pi = pi & R;
    Q = Q & R;

    R = rotationTensorZ(deltaT*pi.z()/momentOfInertia_.zz());
    pi = pi & R;
    Q = Q & R;

    R = rotationTensorY(0.5*deltaT*pi.y()/momentOfInertia_.yy());
    pi = pi & R;
    Q = Q & R;

    R = rotationTensorX(0.5*deltaT*pi.x()/momentOfInertia_.xx());
    pi = pi & R;
    Q = Q & R;

    return Qpi;
}


inline const Foam::PtrList<Foam::REL::freeRunMotionRestraint>&
Foam::REL::freeRunMotion::restraints() const
{
    return restraints_;
}


inline const Foam::PtrList<Foam::REL::freeRunMotionConstraint>&
Foam::REL::freeRunMotion::constraints() const
{
    return constraints_;
}


inline const Foam::point&
Foam::REL::freeRunMotion::initialCentreOfRotation() const
{
    return initialCentreOfRotation_;
}


inline const Foam::tensor&
Foam::REL::freeRunMotion::initialQ() const
{
    return initialQ_;
}


inline const Foam::tensor& Foam::REL::freeRunMotion::Q() const
{
    return motionState_.Q();
}


inline const Foam::vector& Foam::REL::freeRunMotion::v() const
{
    return motionState_.v();
}


inline const Foam::vector& Foam::REL::freeRunMotion::a() const
{
    return motionState_.a();
}


inline const Foam::vector& Foam::REL::freeRunMotion::pi() const
{
    return motionState_.pi();
}


inline const Foam::vector& Foam::REL::freeRunMotion::tau() const
{
    return motionState_.tau();
}


inline Foam::point& Foam::REL::freeRunMotion::initialCentreOfRotation()
{
    return initialCentreOfRotation_;
}


inline Foam::tensor& Foam::REL::freeRunMotion::initialQ()
{
    return initialQ_;
}


inline Foam::tensor& Foam::REL::freeRunMotion::Q()
{
    return motionState_.Q();
}


inline Foam::vector& Foam::REL::freeRunMotion::v()
{
    return motionState_.v();
}


inline Foam::vector& Foam::REL::freeRunMotion::a()
{
    return motionState_.a();
}


inline Foam::vector& Foam::REL::freeRunMotion::pi()
{
    return motionState_.pi();
}


inline Foam::vector& Foam::REL::freeRunMotion::tau()
{
    return motionState_.tau();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline Foam::scalar Foam::REL::freeRunMotion::mass() const
{
    return mass_;
}


inline const Foam::diagTensor&
Foam::REL::freeRunMotion::momentOfInertia() const
{
    return momentOfInertia_;
}


inline const Foam::REL::freeRunMotionState&
Foam::REL::freeRunMotion::state() const
{
    return motionState_;
}


inline const Foam::point& Foam::REL::freeRunMotion::centreOfRotation() const
{
    return motionState_.centreOfRotation();
}


inline const Foam::point&
Foam::REL::freeRunMotion::initialCentreOfMass() const
{
    return initialCentreOfMass_;
}


inline Foam::point Foam::REL::freeRunMotion::centreOfMass() const
{
    return transform(initialCentreOfMass_);
}


inline Foam::vector Foam::REL::freeRunMotion::momentArm() const
{
    return centreOfMass() - motionState_.centreOfRotation();
}


inline const Foam::tensor&
Foam::REL::freeRunMotion::orientation() const
{
    return Q();
}


inline Foam::vector Foam::REL::freeRunMotion::omega() const
{
    return  Q() & (inv(momentOfInertia_) & pi());
}


inline bool Foam::REL::freeRunMotion::report() const
{
    return report_;
}


inline void Foam::REL::freeRunMotion::newTime()
{
    motionState0_ = motionState_;
}


inline Foam::point& Foam::REL::freeRunMotion::centreOfRotation()
{
    return motionState_.centreOfRotation();
}


inline Foam::point Foam::REL::freeRunMotion::velocity
(
    const point& pt
) const
{
    return (omega() ^ (pt - centreOfRotation())) + v();
}


inline Foam::point Foam::REL::freeRunMotion::transform
(
    const point& initialPoint
) const
{
    return
    (
        centreOfRotation()
      + (Q() & initialQ_.T() & (initialPoint - initialCentreOfRotation_))
    );
}

// ************************************************************************* //
