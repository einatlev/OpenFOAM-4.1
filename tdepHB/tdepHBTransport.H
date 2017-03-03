/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::tdepHBTransport

Description
    Transport package using tdepHB's formula.

    Templated into a given thermodynamics package (needed for thermal
    conductivity).

    Dynamic viscosity [kg/m.s]
    \f[
        \mu = \pow{A_mu + \frac{B_mu}{T + C_mu },10}
    \f]

SourceFiles
    tdepHBTransportI.H
    tdepHBTransport.C

\*---------------------------------------------------------------------------*/

#ifndef tdepHBTransport_H
#define tdepHBTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class tdepHBTransport;

template<class Thermo>
inline tdepHBTransport<Thermo> operator+
(
    const tdepHBTransport<Thermo>&,
    const tdepHBTransport<Thermo>&
);

template<class Thermo>
inline tdepHBTransport<Thermo> operator-
(
    const tdepHBTransport<Thermo>&,
    const tdepHBTransport<Thermo>&
);

template<class Thermo>
inline tdepHBTransport<Thermo> operator*
(
    const scalar,
    const tdepHBTransport<Thermo>&
);

template<class Thermo>
inline tdepHBTransport<Thermo> operator==
(
    const tdepHBTransport<Thermo>&,
    const tdepHBTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const tdepHBTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                    Class tdepHBTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class tdepHBTransport
:
    public Thermo
{
    // Private data

        // tdepHB's coefficients
        scalar Amu_, Bmu_, Cmu_;


    // Private Member Functions

        //- Calculate the tdepHB coefficients
        //  given two viscosities and temperatures
        inline void calcCoeffs
        (
            const scalar mu1, const scalar T1,
            const scalar mu2, const scalar T2
        );


public:

    // Constructors

        //- Construct from components
        inline tdepHBTransport
        (
            const Thermo& t,
            const scalar Amu,
            const scalar Bmu,
            const scalar Cmu
        );

        //- Construct from two viscosities
        inline tdepHBTransport
        (
            const Thermo& t,
            const scalar mu1, const scalar T1,
            const scalar mu2, const scalar T2
        );

        //- Construct as named copy
        inline tdepHBTransport(const word&, const tdepHBTransport&);

        //- Construct from Istream
        tdepHBTransport(Istream&);

        //- Construct from dictionary
        tdepHBTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<tdepHBTransport> clone() const;

        // Selector from Istream
        inline static autoPtr<tdepHBTransport> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<tdepHBTransport> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "tdepHB<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;

        // Species diffusivity
        //inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline tdepHBTransport& operator=(const tdepHBTransport&);

        inline void operator+=(const tdepHBTransport&);

        inline void operator-=(const tdepHBTransport&);

        inline void operator*=(const scalar);


    // Friend operators

        friend tdepHBTransport operator+ <Thermo>
        (
            const tdepHBTransport&,
            const tdepHBTransport&
        );

        friend tdepHBTransport operator- <Thermo>
        (
            const tdepHBTransport&,
            const tdepHBTransport&
        );

        friend tdepHBTransport operator* <Thermo>
        (
            const scalar,
            const tdepHBTransport&
        );

        friend tdepHBTransport operator== <Thermo>
        (
            const tdepHBTransport&,
            const tdepHBTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const tdepHBTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tdepHBTransportI.H"

#ifdef NoRepository
#   include "tdepHBTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
