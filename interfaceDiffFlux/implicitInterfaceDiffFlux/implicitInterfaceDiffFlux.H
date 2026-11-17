/*---------------------------------------------------------------------------*\
            Copyright (c) 2017-2019, German Aerospace Center (DLR)
-------------------------------------------------------------------------------
License
    This file is part of the VoFLibrary source code library, which is an
	unofficial extension to OpenFOAM.
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
    Foam::implicitInterfaceDiffFlux

Description
    computes the diffusive flux at the interface and returns a implicit source
    term

    The interface is treated as a Dirichlet boundary (value boundaryValue)
    for a field phi with diffusivity gamma. For each marked (interface) cell,
    the normal gradient at the interface is estimated from the cells in the
    zoneDistribute stencil that lie on one side of the interface: those
    whose direction from the interface centre is within about 75 deg of the
    normal. Each of these neighbours gives (phi_nei - boundaryValue)/d, with
    d the normal distance to the interface, and they are averaged with weight
    cos^4 of the angle to the normal. otherSide = true flips the normal.

    In mulesGradFoam, phi is TL or TV, gamma is the phase conductivity and
    boundaryValue is Tsat.

SourceFiles
    implicitInterfaceDiffFlux.C

\*---------------------------------------------------------------------------*/

#ifndef implicitInterfaceDiffFlux_H
#define implicitInterfaceDiffFlux_H

#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"

#include "fvScalarMatrix.H"
#include "fvm.H"

// #include "reconstructionSchemes.H"
#include "zoneDistribute.H"
#include "tmp.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class implicitInterfaceDiffFlux Declaration
\*---------------------------------------------------------------------------*/

class implicitInterfaceDiffFlux
{

private:

    const fvMesh& mesh_;

    //- Map a local zoneDistribute index to a cell label. Indices >= nCells
    //  refer to boundary faces and are mapped to the cell next to that face.
    label getCellLabel
    (
        const label& localIndex
    );

    //- Collect the stencil cells of celli on the side faceNormal points to
    //  (cos(angle) > 0.25): their phi values, their distance vectors from
    //  the interface centre and their global indices
    void getNeigbours
    (
        const label celli,
        const vector faceCentre,
        const vector faceNormal,
        const volVectorField& centre,
        const volScalarField& phi,
        const Map<vector>& mapCC,
        const Map<scalar>& mapPhi,
        zoneDistribute& exchangeFields,
        DynamicField<scalar>& phiField,
        DynamicField<vector>& distField,
        DynamicField<label>& indicies
    );

    // Private Member Functions

        //- Disallow copy construct
        implicitInterfaceDiffFlux(const implicitInterfaceDiffFlux&);

        //- Disallow default bitwise assignment
        void operator=(const implicitInterfaceDiffFlux&);

public:

    // Constructors

        //- Construct from components
        implicitInterfaceDiffFlux
        (
            const fvMesh& mesh
        );

    //- Destructor
    virtual ~implicitInterfaceDiffFlux()
    {}

    // Member Functions

        // - return the diffusiveFlux
        //  Explicit flux per unit interface area, gamma*dphi/dn, in each
        //  marked cell (zero elsewhere). Units: [phi][gamma]/m
        tmp<volScalarField> diffusiveFlux
        (
            const boolList& markedCells,
            const volVectorField& normal,
            const volVectorField& centre,
            const volScalarField& phi,
            const volScalarField& gamma,
            const volScalarField& boundaryValue,
            const bool otherSide
        );

        // - return the diffusiveFlux
        //  The same flux as an implicit matrix with dimensions dim, to add to
        //  the phi equation. The flux is applied to the stencil neighbours
        //  (not to the marked cell): each neighbour gets c on the diagonal
        //  and c*boundaryValue in the source, with
        //  c = gamma*|interface area|*weight/(sum of weights*d).
        //  Contributions to cells on other processors are sent to them.
        tmp<fvScalarMatrix> diffusiveFlux
        (
            const boolList& markedCells,
            const volVectorField& normal,
            const volVectorField& centre,
            const volScalarField& phi,
            const volScalarField& gamma,
            const volScalarField& boundaryValue,
            const bool otherSide,
            const dimensionSet dim
        );

        // - return the diffusiveFlux
        //  Both of the above at once: diffFlux gets the explicit flux per
        //  unit volume (flux*interface area/cell volume), and the implicit
        //  coefficients are added to phiSource. Not used in mulesGradFoam.
        void diffusiveFlux
        (
            const boolList& markedCells,
            const volVectorField& normal,
            const volVectorField& centre,
            const volScalarField& phi,
            const volScalarField& gamma,
            const volScalarField& boundaryValue,
            const bool otherSide,
            volScalarField& diffFlux,
            fvScalarMatrix& phiSource
        );
	
        // tmp<volVectorField> intGrad
        // (
        //     const boolList& markedCells,
        //     const volVectorField& normal,
        //     const volVectorField& centre,
        //     const volScalarField& phi,
        //     const volScalarField& gamma,
        //     const volScalarField& boundaryValue,
        //     const bool otherSide
        // );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
