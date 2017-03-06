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
    Foam::powerLawTdepTransport

Description
    Transport package using powerLawTdep's formula.

    Templated into a given thermodynamics package (needed for thermal
    conductivity).

    Dynamic viscosity [kg/m.s]
    \f[
        \mu = \pow{A_mu + \frac{B_mu}{T + C_mu },10}
    \f]

SourceFiles
    powerLawTdepTransportI.H
    powerLawTdepTransport.C

\*---------------------------------------------------------------------------*/

#ifndef powerLawTdepTransport_H
#define powerLawTdepTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class powerLawTdepTransport;

template<class Thermo>
inline powerLawTdepTransport<Thermo> operator+
(
    const powerLawTdepTransport<Thermo>&,
    const powerLawTdepTransport<Thermo>&
);

template<class Thermo>
inline powerLawTdepTransport<Thermo> operator-
(
    const powerLawTdepTransport<Thermo>&,
    const powerLawTdepTransport<Thermo>&
);

template<class Thermo>
inline powerLawTdepTransport<Thermo> operator*
(
    const scalar,
    const powerLawTdepTransport<Thermo>&
);

template<class Thermo>
inline powerLawTdepTransport<Thermo> operator==
(
    const powerLawTdepTransport<Thermo>&,
    const powerLawTdepTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const powerLawTdepTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                    Class powerLawTdepTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class powerLawTdepTransport
:
    public Thermo
{
    // Private data

        // powerLawTdep's coefficients
        scalar K1_, K2_, n1_, n2_;


    // Private Member Functions

        //- Calculate the powerLawTdep coefficients
        //  given two viscosities and temperatures
        inline void calcCoeffs
        (
            const scalar mu1, const scalar T1,
            const scalar mu2, const scalar T2
        );


public:

    // Constructors

        //- Construct from components
        inline powerLawTdepTransport
        (
            const Thermo& t,
            const scalar K1,
            const scalar K2,
            const scalar n1,
            const scalar n2
        );

/*        //- Construct from two viscosities
        inline powerLawTdepTransport
        (
            const Thermo& t,
            const scalar mu1, const scalar T1,
            const scalar mu2, const scalar T2
        );
*/
        //- Construct as named copy
        inline powerLawTdepTransport(const word&, const powerLawTdepTransport&);

        //- Construct from Istream
        powerLawTdepTransport(Istream&);

        //- Construct from dictionary
        powerLawTdepTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<powerLawTdepTransport> clone() const;

        // Selector from Istream
        inline static autoPtr<powerLawTdepTransport> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<powerLawTdepTransport> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "powerLawTdep<" + Thermo::typeName() + '>';
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

        inline powerLawTdepTransport& operator=(const powerLawTdepTransport&);

        inline void operator+=(const powerLawTdepTransport&);

        inline void operator-=(const powerLawTdepTransport&);

        inline void operator*=(const scalar);


    // Friend operators

        friend powerLawTdepTransport operator+ <Thermo>
        (
            const powerLawTdepTransport&,
            const powerLawTdepTransport&
        );

        friend powerLawTdepTransport operator- <Thermo>
        (
            const powerLawTdepTransport&,
            const powerLawTdepTransport&
        );

        friend powerLawTdepTransport operator* <Thermo>
        (
            const scalar,
            const powerLawTdepTransport&
        );

        friend powerLawTdepTransport operator== <Thermo>
        (
            const powerLawTdepTransport&,
            const powerLawTdepTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const powerLawTdepTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "powerLawTdepTransportI.H"

#ifdef NoRepository
#   include "powerLawTdepTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
