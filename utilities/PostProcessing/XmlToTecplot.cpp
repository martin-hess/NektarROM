#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    int     i, j, nq,  coordim;
    Array<OneD,NekDouble>  fce; 
    Array<OneD,NekDouble>  xc0,xc1,xc2; 

    if(argc != 3)
    {
        fprintf(stderr,"Usage: FldToTecplot  meshfile fieldfile\n");
        exit(1);
    }

    //----------------------------------------------
    // Read in mesh from input file
    string meshfile(argv[argc-2]);
    SpatialDomains::MeshGraph graph; 
    SpatialDomains::MeshGraphSharedPtr graphShPt = graph.Read(meshfile);
    //----------------------------------------------

    //----------------------------------------------        
    // Define Expansion 
    graph.ReadExpansions(meshfile);
    int expdim   = graphShPt->GetMeshDimension();
    Array<OneD, MultiRegions::ExpListSharedPtr> Exp(1);

    switch(expdim)
    {
    case 1:
        {
            SpatialDomains::MeshGraph1DSharedPtr mesh;
            
            if(!(mesh = boost::dynamic_pointer_cast<SpatialDomains::MeshGraph1D>(graphShPt)))
            {
                ASSERTL0(false,"Dynamics cast failed");
            }
            
            MultiRegions::ExpList1DSharedPtr Exp1D;
            Exp1D = MemoryManager<MultiRegions::ExpList1D>::AllocateSharedPtr(*mesh);
            Exp[0] = Exp1D;
        }
        break;
    case 2:
        {
            SpatialDomains::MeshGraph2DSharedPtr mesh;
            
            if(!(mesh = boost::dynamic_pointer_cast<SpatialDomains::MeshGraph2D>(graphShPt)))
            {
                ASSERTL0(false,"Dynamics cast failed");
            }
            
            MultiRegions::ExpList2DSharedPtr Exp2D;
            Exp2D = MemoryManager<MultiRegions::ExpList2D>::AllocateSharedPtr(*mesh);
            Exp[0] =  Exp2D;

        }
        break;
    case 3:
        ASSERTL0(false,"3D not set up");
        break;
    default:
        ASSERTL0(false,"Expansion dimension not recognised");
        break;
    }

    //-----------------------------------------------
        
    //----------------------------------------------
    // Write solution  depending on #define
    string   outfile(strtok(argv[argc-1],"."));
    outfile +=  ".dat"; 
    ofstream outstrm(outfile.c_str());

    Exp[0]->WriteToFile(outstrm,eTecplot,"Dummy");
    outstrm.close();
    //----------------------------------------------
    return 0;
}
