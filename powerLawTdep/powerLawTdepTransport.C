/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "powerLawTdepTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::powerLawTdepTransport<Thermo>::powerLawTdepTransport(Istream& is)
:
    Thermo(is),
    K1_(readScalar(is)),
    K2_(readScalar(is)),
    n1_(readScalar(is)),
    n2_(readScalar(is))
{
    is.check("powerLawTdepTransport<Thermo>::powerLawTdepTransport(Istream&)");
}


template<class Thermo>
Foam::powerLawTdepTransport<Thermo>::powerLawTdepTransport(const dictionary& dict)
:
    Thermo(dict),
    K1_(readScalar(dict.subDict("transport").lookup("K1"))),
    K2_(readScalar(dict.subDict("transport").lookup("K2"))),
    n1_(readScalar(dict.subDict("transport").lookup("n1"))),
    n2_(readScalar(dict.subDict("transport").lookup("n2")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::powerLawTdepTransport<Thermo>::write(Ostream& os) const
{
    os  << this->specie::name() << endl;
    os  << token::BEGIN_BLOCK  << incrIndent << nl;

    Thermo::write(os);

    dictionary dict("transport");
    dict.add("K1", K1_);
    dict.add("K2", K2_);
    dict.add("n1", n1_);
    dict.add("n2", n2_);
    os  << indent << dict.dictName() << dict;

    os  << decrIndent << token::END_BLOCK << nl;
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const powerLawTdepTransport<Thermo>& st
)
{
    os << static_cast<const Thermo&>(st) << tab << st.K1_ << st.K2_ << tab << st.n1_ << st.n2_;

    os.check
    (
        "Ostream& operator<<(Ostream&, const powerLawTdepTransport<Thermo>&)"
    );

    return os;
}


// ************************************************************************* //
