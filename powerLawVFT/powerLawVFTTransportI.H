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
inline void Foam::powerLawVFTTransport<Thermo>::calcCoeffs
(
    const scalar mu1, const scalar T1,
    const scalar mu2, const scalar T2
)
{
    scalar rootT1 = sqrt(T1);
    scalar mu1rootT2 = mu1*sqrt(T2);
    scalar mu2rootT1 = mu2*rootT1;

    Amu_ = -4.8;
    Bmu_ = 5997;
    Cmu_ = -603;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::powerLawVFTTransport<Thermo>::powerLawVFTTransport
(
    const Thermo& t,
    const scalar Amu,
    const scalar Bmu,
    const scalar Cmu
)
:
    Thermo(t),
    Amu_(Amu),
    Bmu_(Bmu),
    Cmu_(Cmu)
{}


template<class Thermo>
inline Foam::powerLawVFTTransport<Thermo>::powerLawVFTTransport
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
inline Foam::powerLawVFTTransport<Thermo>::powerLawVFTTransport
(
    const word& name,
    const powerLawVFTTransport& st
)
:
    Thermo(name, st),
    Amu_(st.Amu_),
    Bmu_(st.Bmu_),
    Cmu_(st.Cmu_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::powerLawVFTTransport<Thermo> >
Foam::powerLawVFTTransport<Thermo>::clone() const
{
    return autoPtr<powerLawVFTTransport<Thermo> >
    (
        new powerLawVFTTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::powerLawVFTTransport<Thermo> >
Foam::powerLawVFTTransport<Thermo>::New
(
    Istream& is
)
{
    return autoPtr<powerLawVFTTransport<Thermo> >
    (
        new powerLawVFTTransport<Thermo>(is)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::powerLawVFTTransport<Thermo> >
Foam::powerLawVFTTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<powerLawVFTTransport<Thermo> >
    (
        new powerLawVFTTransport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::powerLawVFTTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    scalar powerval_ = (Amu_ + Bmu_/max
    (
	T+Cmu_,
	620.0
    ));
    return max
    (
	1.0,
	min
	(
	     10000000000,
	     1.0*Foam::pow(10.0,powerval_)
	)
    );
}


template<class Thermo>
inline Foam::scalar Foam::powerLawVFTTransport<Thermo>::kappa
(
    const scalar p, const scalar T
) const
{
    //scalar Cv_ = this->Cv(p, T);
    //return mu(p, T)*Cv_*(1.32 + 1.77*this->R()/Cv_);
    return 0.2; // make kappa constant
}


template<class Thermo>
inline Foam::scalar Foam::powerLawVFTTransport<Thermo>::alphah
(
    const scalar p,
    const scalar T
) const
{

    return kappa(p, T)/this->Cpv(p, T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::powerLawVFTTransport<Thermo>&
Foam::powerLawVFTTransport<Thermo>::operator=
(
    const powerLawVFTTransport<Thermo>& st
)
{
    Thermo::operator=(st);

    Amu_ = st.Amu_;
    Bmu_ = st.Bmu_;
    Cmu_ = st.Cmu_;

    return *this;
}


template<class Thermo>
inline void Foam::powerLawVFTTransport<Thermo>::operator+=
(
    const powerLawVFTTransport<Thermo>& st
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator+=(st);

    molr1 /= this->nMoles();
    scalar molr2 = st.nMoles()/this->nMoles();

    Amu_ = molr1*Amu_ + molr2*st.Amu_;
    Bmu_ = molr1*Bmu_ + molr2*st.Bmu_;
    Cmu_ = molr1*Cmu_ + molr2*st.Cmu_;
}


template<class Thermo>
inline void Foam::powerLawVFTTransport<Thermo>::operator-=
(
    const powerLawVFTTransport<Thermo>& st
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator-=(st);

    molr1 /= this->nMoles();
    scalar molr2 = st.nMoles()/this->nMoles();

    Amu_ = molr1*Amu_ - molr2*st.Amu_;
    Bmu_ = molr1*Bmu_ - molr2*st.Bmu_;
    Cmu_ = molr1*Cmu_ - molr2*st.Cmu_;
}


template<class Thermo>
inline void Foam::powerLawVFTTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::powerLawVFTTransport<Thermo> Foam::operator+
(
    const powerLawVFTTransport<Thermo>& st1,
    const powerLawVFTTransport<Thermo>& st2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(st1) + static_cast<const Thermo&>(st2)
    );

    scalar molr1 = st1.nMoles()/t.nMoles();
    scalar molr2 = st2.nMoles()/t.nMoles();

    return powerLawVFTTransport<Thermo>
    (
        t,
        molr1*st1.Amu_ + molr2*st2.Amu_,
        molr1*st1.Bmu_ + molr2*st2.Bmu_,
        molr1*st1.Cmu_ + molr2*st2.Cmu_
    );
}


template<class Thermo>
inline Foam::powerLawVFTTransport<Thermo> Foam::operator-
(
    const powerLawVFTTransport<Thermo>& st1,
    const powerLawVFTTransport<Thermo>& st2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(st1) - static_cast<const Thermo&>(st2)
    );

    scalar molr1 = st1.nMoles()/t.nMoles();
    scalar molr2 = st2.nMoles()/t.nMoles();

    return powerLawVFTTransport<Thermo>
    (
        t,
        molr1*st1.Amu_ - molr2*st2.Amu_,
        molr1*st1.Bmu_ - molr2*st2.Bmu_,
        molr1*st1.Cmu_ - molr2*st2.Cmu_
    );
}


template<class Thermo>
inline Foam::powerLawVFTTransport<Thermo> Foam::operator*
(
    const scalar s,
    const powerLawVFTTransport<Thermo>& st
)
{
    return powerLawVFTTransport<Thermo>
    (
        s*static_cast<const Thermo&>(st),
        st.Amu_,
        st.Bmu_,
        st.Cmu_
    );
}


template<class Thermo>
inline Foam::powerLawVFTTransport<Thermo> Foam::operator==
(
    const powerLawVFTTransport<Thermo>& st1,
    const powerLawVFTTransport<Thermo>& st2
)
{
    return st2 - st1;
}


// ************************************************************************* //
