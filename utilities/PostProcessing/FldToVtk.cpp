#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList3D.h>

using namespace Nektar;

#include <sys/stat.h>

int fexist( const char *filename ) {
  struct stat buffer ;
  if ( stat( filename, &buffer ) ) return 0 ;
  return 1 ;
}

int main(int argc, char *argv[])
{
    int i,j;

    if(argc < 3)
    {
        fprintf(stderr,"Usage: FldToVtk  meshfile  fieldfile(s)\n");
        exit(1);
    }

    //----------------------------------------------
    // Read in mesh from input file
    string meshfile(argv[1]);
    SpatialDomains::MeshGraph graph;
    SpatialDomains::MeshGraphSharedPtr graphShPt = graph.Read(meshfile);
    //----------------------------------------------

    for (int n = 2; n < argc; ++n)
    {
        string fname = std::string(argv[n]);
        fname = fname.substr(0,fname.find_last_of('.')) + ".vtu";
        if (argc > 3)
        {
            if (fexist(fname.c_str()))
            {
                cout << "Skipping converted file: " << argv[n] << endl;
                continue;
            }
            cout << "Processing " << argv[n] << endl;
        }

        //----------------------------------------------
        // Import field file.
        string fieldfile(argv[n]);
        vector<SpatialDomains::FieldDefinitionsSharedPtr> fielddef;
        vector<vector<NekDouble> > fielddata;
        graphShPt->Import(fieldfile,fielddef,fielddata);
        //----------------------------------------------

        //----------------------------------------------
        // Set up Expansion information
        vector< vector<LibUtilities::PointsType> > pointstype;
        for(i = 0; i < fielddef.size(); ++i)
        {
            vector<LibUtilities::PointsType> ptype;
            for(j = 0; j < 3; ++j)
            {
                ptype.push_back(LibUtilities::ePolyEvenlySpaced);
            }
            pointstype.push_back(ptype);
        }
        graphShPt->SetExpansions(fielddef,pointstype);
        //----------------------------------------------


        //----------------------------------------------
        // Define Expansion
        int expdim  = graphShPt->GetMeshDimension();
        int nfields = fielddef[0]->m_Fields.size();
        Array<OneD, MultiRegions::ExpListSharedPtr> Exp(nfields);

        switch(expdim)
        {
        case 1:
            {
                SpatialDomains::MeshGraph1DSharedPtr mesh;

                if(!(mesh = boost::dynamic_pointer_cast<
                                        SpatialDomains::MeshGraph1D>(graphShPt)))
                {
                    ASSERTL0(false,"Dynamic cast failed");
                }

                MultiRegions::ExpList1DSharedPtr Exp1D;
                Exp1D = MemoryManager<MultiRegions::ExpList1D>
                                                        ::AllocateSharedPtr(*mesh);
                Exp[0] = Exp1D;
                for(i = 1; i < nfields; ++i)
                {
                    Exp[i] = MemoryManager<MultiRegions::ExpList1D>
                                                        ::AllocateSharedPtr(*Exp1D);
                }
            }
            break;
        case 2:
            {
                SpatialDomains::MeshGraph2DSharedPtr mesh;

                if(!(mesh = boost::dynamic_pointer_cast<
                                        SpatialDomains::MeshGraph2D>(graphShPt)))
                {
                    ASSERTL0(false,"Dynamic cast failed");
                }

                MultiRegions::ExpList2DSharedPtr Exp2D;
                Exp2D = MemoryManager<MultiRegions::ExpList2D>
                                                        ::AllocateSharedPtr(*mesh);
                Exp[0] =  Exp2D;

                for(i = 1; i < nfields; ++i)
                {
                    Exp[i] = MemoryManager<MultiRegions::ExpList2D>
                                                        ::AllocateSharedPtr(*Exp2D);
                }
            }
            break;
        case 3:
            {
                SpatialDomains::MeshGraph3DSharedPtr mesh;

                if(!(mesh = boost::dynamic_pointer_cast<
                                        SpatialDomains::MeshGraph3D>(graphShPt)))
                {
                    ASSERTL0(false,"Dynamic cast failed");
                }

                MultiRegions::ExpList3DSharedPtr Exp3D;
                Exp3D = MemoryManager<MultiRegions::ExpList3D>
                                                        ::AllocateSharedPtr(*mesh);
                Exp[0] =  Exp3D;

                for(i = 1; i < nfields; ++i)
                {
                    Exp[i] = MemoryManager<MultiRegions::ExpList3D>
                                                        ::AllocateSharedPtr(*Exp3D);
                }
            }
            break;
        default:
            ASSERTL0(false,"Expansion dimension not recognised");
            break;
        }
        //----------------------------------------------

        //----------------------------------------------
        // Copy data from field file
        for(j = 0; j < nfields; ++j)
        {
            for(int i = 0; i < fielddata.size(); ++i)
            {
                Exp[j]->ExtractDataToCoeffs(fielddef [i],
                                            fielddata[i],
                                            fielddef [i]->m_Fields[j]);
            }
            Exp[j]->BwdTrans(Exp[j]->GetCoeffs(),Exp[j]->UpdatePhys());
        }
        //----------------------------------------------

        //----------------------------------------------
        // Write solution
        //string   outname(strtok(argv[n],"."));
        //outname += ".vtu";
        ofstream outfile(fname.c_str());

        Exp[0]->WriteVtkHeader(outfile);
        // For each field write out field data for each expansion.
        for(i = 0; i < Exp[0]->GetExpSize(); ++i)
        {
            Exp[0]->WriteVtkPieceHeader(outfile,i);
            // For this expansion, write out each field.
            for(j = 0; j < Exp.num_elements(); ++j)
            {
                Exp[j]->WriteVtkPieceData(outfile,i, fielddef[0]->m_Fields[j]);
            }
            Exp[0]->WriteVtkPieceFooter(outfile,i);
        }
        Exp[0]->WriteVtkFooter(outfile);
        //----------------------------------------------
    }
    return 0;
}

