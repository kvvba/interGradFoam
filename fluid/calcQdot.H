// Interfacial heat flux, using the "implicit gradient" phase change model
// from the TwoPhaseFlow project (DLR-RY).
//
// The interface is reconstructed as the alpha1 = 0.5 iso-surface. Each cell
// is then classified as an interface, liquid or vapour cell. The interface
// is treated as a Dirichlet boundary at Tsat for both phase temperatures,
// and the conductive heat flux from each side is computed from the normal
// temperature gradient between the interface and nearby cells.

// Cell classification, filled in below
boolList interfaceCells(alpha1.mesh().nCells(), false);
labelList VapourCellsLabels;
labelList LiquidCellsLabels;
boolList vapourCells(alpha1.mesh().nCells(), false);
boolList liquidCells(alpha1.mesh().nCells(), false);

// alpha1 interpolated to the mesh points, needed by the iso-surface cutting
scalarField ap
(
 volPointInterpolation::New(mesh).interpolate(mixture.alpha1())
);

cutCellIso cutCell(mesh,ap);

forAll (alpha1, celli)
{
	// Cut the cell with the alpha1 = 0.5 iso-surface.
	// status: 0 = cut (interface cell), -1 = liquid, 1 = vapour
	label status = cutCell.calcSubCell(celli,0.5);
	if (status == 0) // cell is cut
		{
			// interfaceLabels.append(celli);
			interfaceCells[celli] = true;
			// Area vector and centre of the interface face in this cell
			normal[celli] = cutCell.faceArea();
			centre[celli] = cutCell.faceCentre();
			// Interface area per unit cell volume [1/m]
			SArea[celli] =
				mag(cutCell.faceArea())/mesh.V()[celli];
			// Volume-fraction weighted mixture properties in interface cells,
			// pure-phase values elsewhere (only used for output and
			// computeFluidHeatFluxes.H)
			kappa[celli] = k1.value() * alpha1[celli] + k2.value() * alpha2[celli];
			rhoCp[celli] = rho1.value() * cp1.value() * alpha1[celli] + rho2.value() * cp2.value() * alpha2[celli];
		}
	else
		{
			normal[celli] = Zero;
			centre[celli] = Zero;
			interfaceCells[celli] = false;
		}
	
	if (status == -1)							// liquid
		{
			liquidCells[celli] = true;
			LiquidCellsLabels.append(celli);
			kappa[celli] = k1.value();
			rhoCp[celli] = rho1.value() * cp1.value();
		}
	else if (status == 1)					// vapour
		{
			vapourCells[celli] = true;
			VapourCellsLabels.append(celli);
			kappa[celli] = k2.value();
			rhoCp[celli] = rho2.value() * cp2.value();
		}
}

implicitInterfaceDiffFlux impDiffFlux(mesh);

// Explicit conductive heat flux per unit interface area [W/m^2] in each
// interface cell: kappa*(T_nei - Tsat)/d, averaged over nearby cells on one
// side of the interface. The liquid side (TL) uses the cells the interface
// normal points towards; otherSide = true flips the normal for the vapour
// side (TV). A flux is positive when that phase is above Tsat, i.e. heat
// flows into the interface and evaporates liquid.
tmp<volScalarField> energyFlux1
(
 impDiffFlux.diffusiveFlux
 (
	interfaceCells,
	normal,
	centre,
	TL,
	kappaL,
	TsatField,
	false
 )
);

tmp<volScalarField> energyFlux2
(
 impDiffFlux.diffusiveFlux
 (
	interfaceCells,
	normal,
	centre,
	TV,
	kappaV,
	TsatField,
	true
 )
);

// The same interface heat flux, as implicit matrix contributions for the
// TL and TV equations (TEqn.H). For each stencil cell next to the interface,
// coefficient c is added to the diagonal and c*Tsat to the source. This
// gives the sink c*(T - Tsat), i.e. heat conducted from that cell to the
// interface held at Tsat.
tmp<fvScalarMatrix> TsourceL
(
 impDiffFlux.diffusiveFlux
 (
	interfaceCells,
	normal,
	centre,
	TL,
	kappaL,
	TsatField,
	false,
	dimPower
 )
);

tmp<fvScalarMatrix> TsourceV
(
 impDiffFlux.diffusiveFlux
 (
	interfaceCells,
	normal,
	centre,
	TV,
	kappaV,
	TsatField,
	true,
	dimPower
 )
);

// Convert to heat flux per unit volume [W/m^3]. The total qDot is the heat
// available for phase change, qDot > 0 means evaporation.
qL = energyFlux1 * SArea;
qV = energyFlux2 * SArea;

qDot = qL + qV;

Info << "max qDot = " << max(qDot).value() << nl;
