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

\*---------------------------------------------------------------------------*/

#include "specie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo>
inline void Foam::powerLawTdepTransport<Thermo>::calcCoeffs
(
    const scalar mu1, const scalar T1,
    const scalar mu2, const scalar T2
)
{
    // These are ununsed right now
    scalar rootT1 = sqrt(T1);
    scalar mu1rootT2 = mu1*sqrt(T2);
    scalar mu2rootT1 = mu2*rootT1;

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::powerLawTdepTransport<Thermo>::powerLawTdepTransport
(
    const Thermo& t,
    const scalar K1,
    const scalar K2,
    const scalar n1,
    const scalar n2
)
:
    Thermo(t),
    K1_(K1),
    K2_(K2),
    n1_(n1),
    n2_(n2)
{}


template<class Thermo>
inline Foam::powerLawTdepTransport<Thermo>::powerLawTdepTransport
(
    const Thermo& t,
    const scalar mu1, const scalar T1,
    const scalar mu2, const scalar T2
)
:
    Thermo(t)
{
    calcCoeffs(mu1, T1, mu2, T2);
}


template<class Thermo>
inline Foam::powerLawTdepTransport<Thermo>::powerLawTdepTransport
(
    const word& name,
    const powerLawTdepTransport& st
)
:
    Thermo(name, st),
    K1_(st.K1_),
    K2_(st.K2_),
    n1_(st.n1_),
    n2_(st.n2_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::powerLawTdepTransport<Thermo> >
Foam::powerLawTdepTransport<Thermo>::clone() const
{
    return autoPtr<powerLawTdepTransport<Thermo> >
    (
        new powerLawTdepTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::powerLawTdepTransport<Thermo> >
Foam::powerLawTdepTransport<Thermo>::New
(
    Istream& is
)
{
    return autoPtr<powerLawTdepTransport<Thermo> >
    (
        new powerLawTdepTransport<Thermo>(is)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::powerLawTdepTransport<Thermo> >
Foam::powerLawTdepTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<powerLawTdepTransport<Thermo> >
    (
        new powerLawTdepTransport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::powerLawTdepTransport<Thermo>::mu
(
    const scalar p,
    const scalar T,
    const volVectorField& U //added U
) const
{
    // This the mu calculation
    tmp<volScalarField> sr(strainRate());

    scalar Kvalue_ = K1_ + K2_ * T; 
    scalar nvalue_ = n1_ + n2_ * T;
    return max
    (
        1.0,
	min
        (    
             10000000000,
             (Kvalue_*pow(sr(), (nValue_-1)))
	      /(max(sr(), dimensionedScalar ("VSMALL", dimless/dimTime, VSMALL)))
        )
    );


}


template<class Thermo>
inline Foam::scalar Foam::powerLawTdepTransport<Thermo>::kappa
(
    const scalar p, const scalar T
) const
{
    // This calculates kappa as a function of T
    //scalar Cv_ = this->Cv(p, T);
    //return mu(p, T)*Cv_*(1.32 + 1.77*this->R()/Cv_);
    return 0.2; // make kappa constant
}


template<class Thermo>
inline Foam::scalar Foam::powerLawTdepTransport<Thermo>::alphah
(
    const scalar p,
    const scalar T
) const
{

    return kappa(p, T)/this->Cpv(p, T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::powerLawTdepTransport<Thermo>&
Foam::powerLawTdepTransport<Thermo>::operator=
(
    const powerLawTdepTransport<Thermo>& st
)
{
    Thermo::operator=(st);

    K1_ = st.K1_;
    K2_ = st.K2_;
    n1_ = st.n1_;
    n2_ = st.n2_;

    return *this;
}


template<class Thermo>
inline void Foam::powerLawTdepTransport<Thermo>::operator+=
(
    const powerLawTdepTransport<Thermo>& st
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator+=(st);

    molr1 /= this->nMoles();
    scalar molr2 = st.nMoles()/this->nMoles();

    K1_ = molr1*K1_ + molr2*st.K1_;
    K2_ = molr1*K2_ + molr2*st.K2_;
    n1_ = molr1*n1_ + molr2*st.n1_;
    n2_ = molr1*n2_ + molr2*st.n2_;
}


template<class Thermo>
inline void Foam::powerLawTdepTransport<Thermo>::operator-=
(
    const powerLawTdepTransport<Thermo>& st
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator-=(st);

    molr1 /= this->nMoles();
    scalar molr2 = st.nMoles()/this->nMoles();

    K1_ = molr1*K1_ - molr2*st.K1_;
    K2_ = molr1*K2_ - molr2*st.K2_;
    n1_ = molr1*n1_ - molr2*st.n1_;
    n2_ = molr1*n2_ - molr2*st.n2_;
}


template<class Thermo>
inline void Foam::powerLawTdepTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::powerLawTdepTransport<Thermo> Foam::operator+
(
    const powerLawTdepTransport<Thermo>& st1,
    const powerLawTdepTransport<Thermo>& st2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(st1) + static_cast<const Thermo&>(st2)
    );

    scalar molr1 = st1.nMoles()/t.nMoles();
    scalar molr2 = st2.nMoles()/t.nMoles();

    return powerLawTdepTransport<Thermo>
    (
        t,
        molr1*st1.K1_ + molr2*st2.K1_,
        molr1*st1.K2_ + molr2*st2.K2_,
	molr1*st1.n1_ + molr2*st2.n1_,
        molr1*st1.n2_ + molr2*st2.n2_
    );
}


template<class Thermo>
inline Foam::powerLawTdepTransport<Thermo> Foam::operator-
(
    const powerLawTdepTransport<Thermo>& st1,
    const powerLawTdepTransport<Thermo>& st2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(st1) - static_cast<const Thermo&>(st2)
    );

    scalar molr1 = st1.nMoles()/t.nMoles();
    scalar molr2 = st2.nMoles()/t.nMoles();

    return powerLawTdepTransport<Thermo>
    (
        t,
        molr1*st1.K1_ - molr2*st2.K1_,
        molr1*st1.K2_ - molr2*st2.K2_,
        molr1*st1.n1_ - molr2*st2.n1_,
        molr1*st1.n2_ - molr2*st2.n2_
    );
}


template<class Thermo>
inline Foam::powerLawTdepTransport<Thermo> Foam::operator*
(
    const scalar s,
    const powerLawTdepTransport<Thermo>& st
)
{
    return powerLawTdepTransport<Thermo>
    (
        s*static_cast<const Thermo&>(st),
        st.K1_,
        st.K2_,
        st.n1_,
	st.n2_
    );
}


template<class Thermo>
inline Foam::powerLawTdepTransport<Thermo> Foam::operator==
(
    const powerLawTdepTransport<Thermo>& st1,
    const powerLawTdepTransport<Thermo>& st2
)
{
    return st2 - st1;
}


// ************************************************************************* //
