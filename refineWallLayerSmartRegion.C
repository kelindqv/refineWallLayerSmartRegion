/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Description
    Utility to refine cells next to patches.

    Takes a patchName and number of layers to refine. Works out cells within
    these layers and refines those in the wall-normal direction.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyTopoChange.H"
#include "polyTopoChanger.H"
#include "mapPolyMesh.H"
#include "polyMesh.H"
#include "cellCuts.H"
#include "cellSet.H"
#include "meshCutter.H"

#include "fvCFD.H"
#include "wallFvPatch.H"
#include "nearWallDist.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "A smart way to refine wall layer cells based on the y+/y* value.\n"
        "Please turn on the option \"basedOnYPlus\"to see its effect"
    );

    #include "addOverwriteOption.H"
    #include "addRegionOption.H"
    argList::noParallel();
    argList::validArgs.append("patchName");
    // argList::validArgs.append("edgeWeight");
    argList::validArgs.append("stretchingRatio");

    argList::addOption
    (
        "useSet",
        "name",
        "(obsoleted) restrict cells to refine based on specified cellSet name"
    );

    argList::addBoolOption
    (
        "basedOnYPlus",
        "Based on y+/y* field, run these utilities first. Target yPlus is unity"
    );

    argList::addOption
    (
        "targetYPlus",
        "value",
        "Manually set the target yPlus, e.g. 5 or 30"
    );

    argList::addOption
    (
        "refineLevel",
        "value",
        "Manually set the refineLevel, e.g. 3"
    );
    
    argList::addBoolOption
    (
        "adjacent",
        "Use stretchRatio as ratio between cell sizes"
    );
    
    argList::addBoolOption
    (
        "polyhedra",
        "Create polyhedra (hanging nodes) at patch edges, rather than cuts"
    );


#   include "setRootCase.H"
#   include "createTime.H"
    runTime.functionObjects().off();
// #   include "createPolyMesh.H"
#   include "createNamedMesh.H"
    const word oldInstance = mesh.pointsInstance();
    const bool basedOnYPlus = args.optionFound("basedOnYPlus");
    const bool adjacent = args.optionFound("adjacent");

    const word patchName = args[1];
    const scalar stretchRatio = args.argRead<scalar>(2);
    scalar weight = stretchRatio/(stretchRatio+1.);
    //scalar weight = 1./(args.argRead<scalar>(2)+1);
    if (stretchRatio < 1.)
    {
        FatalErrorIn(args.executable())
            << "stretchingRatio must >= 1.0 "
            << exit(FatalError);
    }
    
    const bool overwrite = args.optionFound("overwrite");

    label patchID = mesh.boundaryMesh().findPatchID(patchName);

    if (patchID == -1)
    {
        FatalErrorIn(args.executable())
            << "Cannot find patch " << patchName << endl
            << "Valid patches are " << mesh.boundaryMesh().names()
            << exit(FatalError);
    }
    const polyPatch& pp = mesh.boundaryMesh()[patchID];

    static label refineLevelValue;
    refineLevelValue = args.optionLookupOrDefault("refineLevel",1);

    if (basedOnYPlus)
    {
        volScalarField yPlus
        (
            IOobject
            (
                "yPlus",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("yPlus", dimless, 0.0)
        );

        if (!yPlus.headerOk())
        {
            volScalarField yPlus
            (
                IOobject
                (
                    "yStar",
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("yPlus", dimless, 0.0)
            );
        }
        else
        {
            Info<< "Note: this utility default is using yPlus, "
                << "if you have yStar file also in the folder, "
                << nl
                << "      please double check which one you want to use, "
                << "and delete the other. Thanks!\n"
                << endl;
        }

        if (!yPlus.headerOk())
        {
            FatalErrorIn
            (
             "smartRefineWallLayer"
            )
            << "No yStar and yPlus file to use!"
            << exit(FatalError);
        }

        forAll(yPlus.boundaryField()[patchID], faceI)
        {   
            label faceCellI = pp.faceCells()[faceI];
            yPlus[faceCellI]=yPlus.boundaryField()[patchID][faceI];
        }

        const scalarField& Yp = yPlus.boundaryField()[patchID];
        scalar YpMax = max(Yp);

        scalar targetYPlusValue;
        targetYPlusValue = args.optionLookupOrDefault("targetYPlus",1.0);
        if (args.optionFound("targetYPlus"))
        {
            if
            (
                targetYPlusValue >= YpMax
             && !args.optionFound("refineLevel")
            )
            {
                FatalErrorIn
                (
                 "smartRefineWallLayer"
                )
                << "targetYPlusValue is larger than max yPlus ("
                << YpMax
                << "), no need to refine the wall layer mesh."
                << nl
                << "Or you can use \"refineLevel\" option to force refine."
                << exit(FatalError);
            }
            refineLevelValue = ceil
            (
                std::log(YpMax/targetYPlusValue)/std::log(scalar(2.0))
            );
        }
        else
        {
            refineLevelValue = ceil
            (
                std::log(YpMax/targetYPlusValue)/std::log(scalar(2.0))
            );
            Info<< "Hiii ..." << endl;
            Info<< "refineLevelValue = " << refineLevelValue << endl;
        }

        if (args.optionFound("refineLevel"))
        {
            // word tmp = "refineLevel";
            // refineLevelValue = args.optionRead(tmp);
            // Info<< "Hi ..." << endl;
            refineLevelValue = args.optionLookupOrDefault("refineLevel",1);
        }


        Info<< "Using refineWallLayer utility in a smart way!"
            << nl << nl
            << "Max yPlus is "
            << YpMax
            << ", targetYPlusValue is "
            << targetYPlusValue
            << ", refineLevel is "
            << refineLevelValue
            << nl
            << endl;
    }
    else
    {
        Info<< "Using refineWallLayer utility in old way. Refine once!"
            << nl << endl;
    }

    volScalarField::GeometricBoundaryField dOrg = nearWallDist(mesh).y();
    Info<< "\tOriginally, patch " << patchID
        << " named " << pp.name() << nl
        << "\t d  : min: " << min(dOrg[patchID])
        << " max: " << max(dOrg[patchID])
        << " average: " << average(dOrg[patchID])
        << nl << endl;

    for (int i = 0; i < refineLevelValue; ++i)
    {
        //
        // Cells cut
        //
        labelHashSet cutCells(4*pp.size());
        
        if( adjacent )
        {
            scalar r(stretchRatio);
            int N(refineLevelValue-i+1);
            weight = (pow(r,N-1)-1.)/(pow(r,N)-1.);
        }

        const labelList& meshPoints = pp.meshPoints();
        
        if (args.optionFound("polyhedra"))
        {
            // Only select cells actually connected to patch
            forAll(pp, faceI)
            {   
                label faceCellI = pp.faceCells()[faceI];
                cutCells.insert(faceCellI);
            }
            
        }
        else
        {
            // Original behavior - cut everything nearby
            forAll(meshPoints, pointI)
            {
                label meshPointI = meshPoints[pointI];

                const labelList& pCells = mesh.pointCells()[meshPointI];

                forAll(pCells, pCellI)
                {
                    cutCells.insert(pCells[pCellI]);
                }
            }
        }
     


        Info<< "Selected " << cutCells.size()
            << " cells connected to patch " << pp.name() << endl << endl;

        //
        // List of cells to refine
        //

        word setName;
        if (args.optionReadIfPresent("useSet", setName))
        {
            Info<< "Subsetting cells to cut based on cellSet"
                << setName << nl << endl;

            cellSet cells(mesh, setName);

            Info<< "Read " << cells.size() << " cells from cellSet "
                << cells.instance()/cells.local()/cells.name()
                << nl << endl;

            forAllConstIter(cellSet, cells, iter)
            {
                cutCells.erase(iter.key());
            }
            Info<< "Removed from cells to cut all the ones not in set "
                << setName << nl << endl;
        }

        
        // Mark all meshpoints on patch

        boolList vertOnPatch(mesh.nPoints(), false);

        forAll(meshPoints, pointI)
        {
            const label meshPointI = meshPoints[pointI];

            vertOnPatch[meshPointI] = true;
        }
        



        // Mark cut edges.

        DynamicList<label> allCutEdges(pp.nEdges());

        DynamicList<scalar> allCutEdgeWeights(pp.nEdges());

        forAll(meshPoints, pointI)
        {
            label meshPointI = meshPoints[pointI];

            const labelList& pEdges = mesh.pointEdges()[meshPointI];

            forAll(pEdges, pEdgeI)
            {
                const label edgeI = pEdges[pEdgeI];
                const edge& e = mesh.edges()[edgeI];

                label otherPointI = e.otherVertex(meshPointI);

                if (!vertOnPatch[otherPointI])
                {
                    allCutEdges.append(edgeI);

                    if (e.start() == meshPointI)
                    {
                        allCutEdgeWeights.append(weight);
                    }
                    else
                    {
                        allCutEdgeWeights.append(1 - weight);
                    }
                }
            }
        }

        allCutEdges.shrink();
        allCutEdgeWeights.shrink();

        Info<< "Cutting:" << nl
            << "    cells:" << cutCells.size() << nl
            << "    edges:" << allCutEdges.size() << nl
            << endl;

        // Transfer DynamicLists to straight ones.
        scalarField cutEdgeWeights;
        cutEdgeWeights.transfer(allCutEdgeWeights);
        allCutEdgeWeights.clear();


        // Gets cuts across cells from cuts through edges.
        cellCuts cuts
        (
            mesh,
            cutCells.toc(),     // cells candidate for cutting
            labelList(0),       // cut vertices
            allCutEdges,        // cut edges
            cutEdgeWeights      // weight on cut edges
        );

        polyTopoChange meshMod(mesh);

        // Cutting engine
        meshCutter cutter(mesh);

        // Insert mesh refinement into polyTopoChange.
        cutter.setRefinement(cuts, meshMod);

        // Do all changes
        Info<< "Morphing ..." << endl;

        autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh, false);

        if (morphMap().hasMotionPoints())
        {
            mesh.movePoints(morphMap().preMotionPoints());
        }

        // Update stored labels on meshCutter.
        cutter.updateMesh(morphMap());

    }

    if (overwrite)
    {
        mesh.setInstance(oldInstance);
    }
    else
    {
        runTime++;
    }

    // Write resulting mesh
    Info<< "Writing refined morphMesh to time " << runTime.timeName() << endl;

    mesh.write();

    //
    // Output updated mesh information to the screen
    //
    Foam::fvMesh meshUpdated
    (
     Foam::IOobject
     (
         regionName,
         runTime.timeName(),
         runTime,
         Foam::IOobject::MUST_READ
     )
    );
    volScalarField::GeometricBoundaryField dNew = 
        nearWallDist(meshUpdated).y();
    Info<< "\n\tAfter refining, Patch " << patchID
        << " named " << pp.name() << nl
        << "\t d  : min: " << min(dNew[patchID])
        << " max: " << max(dNew[patchID])
        << " average: " << average(dNew[patchID])
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
