#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <SpatialDomains/MeshPartition.h>
#include <MultiRegions/DisContField1D.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr vSession;
    LibUtilities::CommSharedPtr vComm;
    string meshfile(argv[1]);
    string vCommModule("Serial");

    vSession = MemoryManager<LibUtilities::SessionReader>::AllocateSharedPtr(meshfile);

    if (vSession->DefinesSolverInfo("Communication"))
    {
        vCommModule = vSession->GetSolverInfo("Communication");
    }
    else if (LibUtilities::GetCommFactory().ModuleExists("ParallelMPI"))
    {
        vCommModule = "ParallelMPI";
    }

    vComm = LibUtilities::GetCommFactory().CreateInstance(vCommModule,argc,argv);

    if (vComm->GetSize() > 1)
    {
        if (vComm->GetRank() == 0)
        {
            LibUtilities::SessionReaderSharedPtr vSession = MemoryManager<LibUtilities::SessionReader>::AllocateSharedPtr(meshfile);
            SpatialDomains::MeshPartitionSharedPtr vPartitioner = MemoryManager<SpatialDomains::MeshPartition>::AllocateSharedPtr(vSession);
            vPartitioner->PartitionMesh(vComm->GetSize());
            vPartitioner->WritePartitions(vSession, meshfile);
        }

        vComm->Block();

        meshfile = meshfile + "." + boost::lexical_cast<std::string>(vComm->GetRank());
    }

    MultiRegions::DisContField1DSharedPtr Exp,Fce;
    int     i, nq,  coordim;
    Array<OneD,NekDouble>  fce;
    Array<OneD,NekDouble>  xc0,xc1,xc2;
    NekDouble  lambda;

    if(argc != 2)
    {
        fprintf(stderr,"Usage: Helmholtz1D  meshfile\n");
        exit(1);
    }

    //----------------------------------------------
    // Read in mesh from input file
    SpatialDomains::MeshGraph1D graph1D;
    graph1D.ReadGeometry(meshfile);
    graph1D.ReadExpansions(meshfile);
    //----------------------------------------------

    //----------------------------------------------
    // read the problem parameters from input file
    SpatialDomains::BoundaryConditions bcs(vSession, &graph1D);
    //----------------------------------------------

    //----------------------------------------------
    // Print summary of solution details
    lambda = vSession->GetParameter("Lambda");
    const SpatialDomains::CompositeMap domain = (graph1D.GetDomain());
    cout << "Solving 1D Helmholtz:"  << endl;
    cout << "         Lambda     : " << lambda << endl;
    //----------------------------------------------

    //----------------------------------------------
    // Define Expansion
    Exp = MemoryManager<MultiRegions::DisContField1D>::
        AllocateSharedPtr(vComm,graph1D,bcs);
    //----------------------------------------------

    //----------------------------------------------
    // Set up coordinates of mesh for Forcing function evaluation
    coordim = Exp->GetCoordim(0);
    nq      = Exp->GetTotPoints();

    xc0 = Array<OneD,NekDouble>(nq);
    xc1 = Array<OneD,NekDouble>(nq);
    xc2 = Array<OneD,NekDouble>(nq);

    switch(coordim)
    {
    case 1:
        Exp->GetCoords(xc0);
        Vmath::Zero(nq,&xc1[0],1);
        Vmath::Zero(nq,&xc2[0],1);
        break;
    case 2:
        Exp->GetCoords(xc0,xc1);
        Vmath::Zero(nq,&xc2[0],1);
        break;
    case 3:
        Exp->GetCoords(xc0,xc1,xc2);
        break;
    }
    //----------------------------------------------

    //----------------------------------------------
    // Define forcing function for first variable defined in file
    fce = Array<OneD,NekDouble>(nq);
    LibUtilities::EquationSharedPtr ffunc = vSession->GetFunction("Forcing", 0);
    for(i = 0; i < nq; ++i)
    {
        fce[i] = ffunc->Evaluate(xc0[i],xc1[i],xc2[i]);
    }
    //----------------------------------------------

    //----------------------------------------------
    // Setup expansion containing the  forcing function
    Fce = MemoryManager<MultiRegions::DisContField1D>::AllocateSharedPtr(*Exp);
    Fce->SetPhys(fce);
    //----------------------------------------------

    //----------------------------------------------
    // Helmholtz solution taking physical forcing
    Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateCoeffs(), lambda);
    //----------------------------------------------

    //----------------------------------------------
    // Backward Transform Solution to get solved values at
    Exp->BwdTrans(Exp->GetCoeffs(), Exp->UpdatePhys());
    //----------------------------------------------

    //----------------------------------------------
    // Write solution
    string   out(strtok(argv[1],"."));
    string   endfile(".fld");
    out += endfile;
    std::vector<SpatialDomains::FieldDefinitionsSharedPtr> FieldDef
        = Exp->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());
    for(i = 0; i < FieldDef.size(); ++i)
    {
        FieldDef[i]->m_fields.push_back("u");
        Exp->AppendFieldData(FieldDef[i], FieldData[i]);
    }
    graph1D.Write(out, FieldDef, FieldData);
    //----------------------------------------------

    //----------------------------------------------
    // See if there is an exact solution, if so
    // evaluate and plot errors
    LibUtilities::EquationSharedPtr ex_sol =
                                    vSession->GetFunction("ExactSolution", 0);

    if(ex_sol)
    {
        //----------------------------------------------
        // evaluate exact solution
        for(i = 0; i < nq; ++i)
        {
            fce[i] = ex_sol->Evaluate(xc0[i],xc1[i],xc2[i]);
        }
        //----------------------------------------------

        //--------------------------------------------
        // Calculate L_inf error
        Fce->SetPhys(fce);
        cout << "L infinity error: " << Exp->Linf(Fce->GetPhys()) << endl;
        cout << "L 2 error:        " << Exp->L2  (Fce->GetPhys()) << endl;
        //--------------------------------------------
    }
    //----------------------------------------------

    vComm->Finalise();

    return 0;
}

