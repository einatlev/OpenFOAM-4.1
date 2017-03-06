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

#include "powerLawVFTTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::powerLawVFTTransport<Thermo>::powerLawVFTTransport(Istream& is)
:
    Thermo(is),
    Amu_(readScalar(is)),
    Bmu_(readScalar(is)),
    Cmu_(readScalar(is))
{
    is.check("powerLawVFTTransport<Thermo>::powerLawVFTTransport(Istream&)");
}


template<class Thermo>
Foam::powerLawVFTTransport<Thermo>::powerLawVFTTransport(const dictionary& dict)
:
    Thermo(dict),
    Amu_(readScalar(dict.subDict("transport").lookup("Amu"))),
    Bmu_(readScalar(dict.subDict("transport").lookup("Bmu"))),
    Cmu_(readScalar(dict.subDict("transport").lookup("Cmu")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::powerLawVFTTransport<Thermo>::write(Ostream& os) const
{
    os  << this->specie::name() << endl;
    os  << token::BEGIN_BLOCK  << incrIndent << nl;

    Thermo::write(os);

    dictionary dict("transport");
    dict.add("Amu", Amu_);
    dict.add("Bmu", Bmu_);
    dict.add("Cmu", Cmu_);
    os  << indent << dict.dictName() << dict;

    os  << decrIndent << token::END_BLOCK << nl;
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const powerLawVFTTransport<Thermo>& st
)
{
    os << static_cast<const Thermo&>(st) << tab << st.Amu_ << st.Bmu_ << tab << st.Cmu_;

    os.check
    (
        "Ostream& operator<<(Ostream&, const powerLawVFTTransport<Thermo>&)"
    );

    return os;
}


// ************************************************************************* //
