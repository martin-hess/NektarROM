///////////////////////////////////////////////////////////////////////////////
//
// File: NodalDemo.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Demo for testing functionality of StdProject
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <boost/core/ignore_unused.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
//#include <LibUtilities/BasicUtils/Timer.h>

#include <StdRegions/StdPointExp.h>
#include <StdRegions/StdSegExp.h>
#include <StdRegions/StdTriExp.h>
#include <StdRegions/StdNodalTriExp.h>
#include <StdRegions/StdQuadExp.h>
#include <StdRegions/StdHexExp.h>
#include <StdRegions/StdPrismExp.h>
#include <StdRegions/StdNodalPrismExp.h>
#include <StdRegions/StdPyrExp.h>
#include <StdRegions/StdTetExp.h>
#include <StdRegions/StdNodalTetExp.h>

using namespace std;
using namespace Nektar;
using namespace Nektar::LibUtilities;
using namespace Nektar::StdRegions;

namespace po = boost::program_options;

template <class myType>
Array<OneD, NekDouble> commoncode(myType*E,
                                 Array<OneD, NekDouble> &physvals,
                             	  Array< OneD, Array<OneD, NekDouble> >evalPtsxy)
{
    int numevalvals = evalPtsxy[0].num_elements();
    int dim = evalPtsxy.num_elements();
    Array<OneD, NekDouble> ret(numevalvals);
    for(int i = 0; i < numevalvals; i++)
    {

        Array<OneD, NekDouble> coords(dim);
        for( int k = 0; k < dim; k++)
        {
            coords[k] = evalPtsxy[k][i];
        }

        NekDouble val1 = E->PhysEvaluate(coords,physvals);
        ret[i] = val1;

    }

    return ret;
}




// Evaluate polynomial for testing and save in ret (size same as pts[0])
// if tensorp = 0, we need tensorprod
// else just eval at pts
Array<OneD, NekDouble> EvalPoly(
                                Array< OneD, Array< OneD, NekDouble > >&pts
                                )
{
    Array<OneD, NekDouble> ret(pts[0].num_elements());
    
    // check if pts[0] and pts[1] have same size
    
    
    //polynomial = x^2 + y^2 - 3x - 4
    for(int i = 0; i < pts[0].num_elements(); i++)
    {
        ret[i] = pow(pts[0][i],2) + pow(pts[1][i],2) - 3*pts[0][i] - 4.0;
    }
    return ret;
    
}


int main(int argc, char *argv[])
{
    string shape, ntype;
    int baryinterp = 1;
    vector<string> basis(3, "NoBasisType"), pointstype(3, "NoPointsType");
    vector<int> order, points;
    po::options_description desc("Available options");
    desc.add_options()
        ("help,h",
         "Produce this help message and list basis and shape types.")
        ("nodal,n",
         po::value<string>(&ntype),
         "Optional nodal type, autofills shape and basis choices.")
        ("shape,s",
         po::value<string>(&shape),
         "Region shape to project function on.")
        ("basis,b",
         po::value<vector<string>>(&basis)->multitoken(),
         "Basis type, separate by spaces for higher dimensions.")
        ("order,o",
         po::value<vector<int>>(&order)->multitoken()->required(),
         "Order of basis sets, separate by spaces for higher dimensions.")
        ("points,p",
         po::value<vector<int>>(&points)->multitoken()->required(),
         "Number of quadrature points, separate by spaces for "
         "higher dimensions.")
        ("pointstype,P",
         po::value<vector<string>>(&pointstype)->multitoken(),
         "Optional points type, separate by spaces for higher dimensions.")
        ("diff,d",
         "Project derivative.");

    po::variables_map vm;
    try
    {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        if (vm.count("help"))
        {
            cout << desc;
            cout << endl << "All nodal types, -n [ --nodal ], are:" << endl;
            for (int i = 22; i < SIZE_PointsType; ++i)
            {
                cout << kPointsTypeStr[i] << endl;
            };
            cout << endl << "All shape types, -s [ --shape ], are:" << endl;
            for (int i = 1; i < SIZE_ShapeType; ++i)
            {
                cout << ShapeTypeMap[i] << endl;
            };
            cout << endl << "All basis types, -b [ --basis ], are:" << endl;
            for (int i = 1; i < SIZE_BasisType; ++i)
            {
                cout << BasisTypeMap[i] << endl;
            };
            cout << endl << "All points types, -P [ --pointstype ], are:"
                 << endl;
            for (int i = 1; i < SIZE_PointsType; ++i)
            {
                cout << kPointsTypeStr[i] << endl;
            };
            return 1;
        }
        po::notify(vm);
    }
    catch (const exception &e)
    {
        cerr << "Error: " << e.what() << endl << desc;
        return 0;
    }

    vector<PointsType> ptype;
    if (vm.count("pointstype"))
    {
        for (auto &p : pointstype)
        {
            PointsType tmp = eNoPointsType;
            for (int i = 1; i < SIZE_PointsType; ++i) // starts at nodal points
            {
                if (boost::iequals(kPointsTypeStr[i], p))
                {
                    tmp = static_cast<PointsType>(i);
                    break;
                }
                ASSERTL0(i != SIZE_PointsType - 1,
                         "The points type '" + p + "' does not exist");
            }

            ptype.push_back(tmp);
        }
    }

    //Convert string input argument to nodal type
    PointsType nodaltype = eNoPointsType;
    ShapeType stype = eNoShapeType;
    vector<BasisType> btype(3, eNoBasisType);
    if (vm.count("nodal"))
    {
        for (int i = 22; i < SIZE_PointsType; ++i) // starts at nodal points
        {
            if (boost::iequals(kPointsTypeStr[i], ntype))
            {
                nodaltype = static_cast<PointsType>(i);
                break;
            }
            ASSERTL0(i != SIZE_PointsType - 1, ("The nodal type '" + ntype +
                                                "' does not exist"))
                }
        switch (nodaltype)
        {
            case eNodalTriElec:
            case eNodalTriFekete:
            case eNodalTriSPI:
            case eNodalTriEvenlySpaced:
                btype[0] = eOrtho_A;
                btype[1] = eOrtho_B;
                stype = eTriangle;
                break;
            case eNodalQuadElec:
                btype[0] = eOrtho_A;
                btype[1] = eOrtho_B;
                stype = eQuadrilateral;
                break;
            case eNodalTetElec:
            case eNodalTetSPI:
            case eNodalTetEvenlySpaced:
                btype[0] = eOrtho_A;
                btype[1] = eOrtho_B;
                btype[2] = eOrtho_C;
                stype = eTetrahedron;
                break;
            case eNodalPrismElec:
            case eNodalPrismSPI:
            case eNodalPrismEvenlySpaced:
                btype[0] = eOrtho_A;
                btype[1] = eOrtho_A;
                btype[2] = eOrtho_B;
                stype = ePrism;
                break;
            case eNodalHexElec:
                btype[0] = eOrtho_A;
                btype[1] = eOrtho_A;
                btype[2] = eOrtho_A;
                stype = eHexahedron;
                break;
            default:
                ASSERTL0(!nodaltype, ("The nodal type '" + ntype +
                                      "' is invalid for StdProject."));
                break;
        }
    }

    //Convert string input argument to shape type
    if (stype == eNoShapeType)
    {
        for (int i = 1; i < SIZE_ShapeType; ++i)
        {
            if (boost::iequals(ShapeTypeMap[i], shape))
            {
                stype = static_cast<ShapeType>(i);
                break;
            }
            ASSERTL0(i != SIZE_ShapeType - 1, ("The shape type '" + shape +
                                               "' does not exist"))
                }
    }

    if (stype == ePoint && vm.count("diff"))
    {
        NEKERROR(ErrorUtil::efatal,
                 "It is not possible to run the diff version for shape: point")
            }

    //Check arguments supplied equals dimension
    const int dimension = (stype == ePoint) ? 1 : ShapeTypeDimMap[stype];
    ASSERTL0(order.size() == dimension,
             "Number of orders supplied should match shape dimension");
    ASSERTL0(points.size() == dimension,
             "Number of points supplied should match shape dimension");
    ASSERTL0(ptype.size() == dimension || ptype.size() == 0,
             "Number of points types should match shape dimension if "
             "supplied.");

    if (!vm.count("nodal"))
    {
        ASSERTL0(basis.size() == dimension,
                 "Number of bases supplied should match shape dimension");
        //Convert string input argument to basis types
        for (int i = 0; i < dimension; ++i)
        {
            for (int j = 1; j < SIZE_BasisType; ++j)
            {
                if (boost::iequals(BasisTypeMap[j], basis[i]))
                {
                    btype[i] = static_cast<BasisType>(j);
                    break;
                }
                ASSERTL0(j != SIZE_BasisType - 1, ("The basis type '" + basis[i]
                                                   + "' does not exist"))
                    }
        }
    }

    //check basis selection is permitted for chosen shape
    map<ShapeType, vector<vector<BasisType>>> allowableBasis;
    allowableBasis[ePoint] = {
        {eOrtho_A, eModified_A, eFourier, eGLL_Lagrange, eGauss_Lagrange,
         eLegendre, eChebyshev, eMonomial, eFourierSingleMode,
         eFourierHalfModeRe, eFourierHalfModeIm}
    };
    allowableBasis[eSegment] = { allowableBasis[ePoint][0] };
    allowableBasis[eTriangle] = {
        {eOrtho_A, eModified_A, eGLL_Lagrange, eGauss_Lagrange},
        {eOrtho_B, eModified_B, eGLL_Lagrange, eGauss_Lagrange}
    };
    allowableBasis[eQuadrilateral] = {
        allowableBasis[eSegment][0], allowableBasis[eSegment][0]
    };
    allowableBasis[eTetrahedron] = {
        {eOrtho_A, eModified_A, eGLL_Lagrange, eGauss_Lagrange},
        {eOrtho_B, eModified_B, eGLL_Lagrange, eGauss_Lagrange},
        {eOrtho_C, eModified_C, eGLL_Lagrange, eGauss_Lagrange}
    };
    allowableBasis[ePyramid] = {
        {eOrtho_A,    eModified_A,    eGLL_Lagrange, eGauss_Lagrange},
        {eOrtho_A,    eModified_A,    eGLL_Lagrange, eGauss_Lagrange},
        {eOrthoPyr_C, eModifiedPyr_C, eGLL_Lagrange, eGauss_Lagrange}
    };
    allowableBasis[ePrism] = {
        {eOrtho_A, eModified_A, eGLL_Lagrange, eGauss_Lagrange},
        {eOrtho_A, eModified_A, eGLL_Lagrange, eGauss_Lagrange},
        {eOrtho_B, eModified_B, eGLL_Lagrange, eGauss_Lagrange}
    };
    allowableBasis[eHexahedron] = {
        {eOrtho_A, eModified_A, eFourier, eGLL_Lagrange, eGauss_Lagrange,
         eLegendre, eChebyshev, eMonomial},
        {eOrtho_A, eModified_A, eFourier, eGLL_Lagrange, eGauss_Lagrange,
         eLegendre, eChebyshev, eMonomial},
        {eOrtho_A, eModified_A, eFourier, eGLL_Lagrange, eGauss_Lagrange,
         eLegendre, eChebyshev, eMonomial}
    };


    for (int i = 0; i < dimension; ++i)
    {
        const unsigned int basisListLength = allowableBasis[stype][i].size();
        for (int j = 0; j < basisListLength; ++j)
        {
            if (allowableBasis[stype][i][j] == btype[i])
            {
                break;
            }
            ASSERTL0(j != basisListLength - 1,
                     ("The basis type '" +
                      static_cast<string>(BasisTypeMap[btype[i]]) +
                      "' is invalid for basis argument " + to_string(i + 1) +
                      " for shape '" + ShapeTypeMap[stype] + "'."))
                }
    }

    //Declaration of other variables needed
    StdExpansion *E = nullptr;

    // Assign points type according to basis type selection, if not already
    // assigned.
    if (ptype.size() == 0)
    {
        ptype.resize(dimension);
        for (int i = 0; i < dimension; ++i)
        {
            if (btype[i] == eFourier)
            {
                ptype[i] = eFourierEvenlySpaced;
            }
            else if (btype[i] == eFourierSingleMode ||
                     btype[i] == eFourierHalfModeRe ||
                     btype[i] == eFourierHalfModeIm)
            {
                ptype[i] = eFourierSingleModeSpaced;
            }
            else
            {
                if (i == 1 && (stype == eTriangle || stype == eTetrahedron))
                {
                    ptype[i] = eGaussRadauMAlpha1Beta0;
                }
                else if (i == 2 && (stype == eTetrahedron || stype == ePyramid))
                {
                    ptype[i] = eGaussRadauMAlpha2Beta0;
                }
                else if (i == 2 && stype == ePrism)
                {
                    ptype[i] = eGaussRadauMAlpha1Beta0;
                }
                else
                {
                    ptype[i] = eGaussLobattoLegendre;
                }
            }
        }
    }

    vector<PointsKey> pkey;
    vector<BasisKey> bkey;
    for (int i = 0; i < dimension; ++i)
    {
        pkey.emplace_back(PointsKey(points[i], ptype[i]));
        bkey.emplace_back(BasisKey(btype[i], order[i], pkey[i]));
    }

    switch (stype)
    {
        cout<<"stype = "<<stype;
        case ePoint:
        {
            E = new StdPointExp(bkey[0]);
            break;
        }
        case eSegment:
        {
            E = new StdSegExp(bkey[0]);
            break;
        }
        case eTriangle:
        {
            E = nodaltype != eNoPointsType ? new StdNodalTriExp(bkey[0],
                                                                bkey[1],
                                                                nodaltype)
                : new StdTriExp(bkey[0], bkey[1]);
            break;
        }
        case eQuadrilateral:
        {
            E = new StdQuadExp(bkey[0], bkey[1]);
            break;
        }
        case eTetrahedron:
        {
            E = nodaltype != eNoPointsType ? new StdNodalTetExp(bkey[0],
                                                                bkey[1],
                                                                bkey[2],
                                                                nodaltype)
                : new StdTetExp(bkey[0], bkey[1],
                                bkey[2]);
            break;
        }
        case ePyramid:
        {
            E = new StdPyrExp(bkey[0], bkey[1], bkey[2]);
            break;
        }
        case ePrism:
        {
            E = nodaltype != eNoPointsType ? new StdNodalPrismExp(bkey[0],
                                                                  bkey[1],
                                                                  bkey[2],
                                                                  nodaltype)
                : new StdPrismExp(bkey[0], bkey[1],
                                  bkey[2]);
            break;
        }
        case eHexahedron:
        {
            E = new StdHexExp(bkey[0], bkey[1], bkey[2]);
            break;
        }
        default:
            break;
    }

    const auto totPoints = (unsigned) E->GetTotPoints();
    Array<OneD, NekDouble> x = Array<OneD, NekDouble>(totPoints);
    Array<OneD, NekDouble> y = Array<OneD, NekDouble>(totPoints);
    Array<OneD, NekDouble> z = Array<OneD, NekDouble>(totPoints);
    Array<OneD, NekDouble> dx = Array<OneD, NekDouble>(totPoints);
    Array<OneD, NekDouble> dy = Array<OneD, NekDouble>(totPoints);
    Array<OneD, NekDouble> dz = Array<OneD, NekDouble>(totPoints);

    switch (dimension)
    {
        case 1:
        {
            E->GetCoords(x);
            break;
        }

        case 2:
        {
            E->GetCoords(x, y);
            break;
        }

        case 3:
        {
            E->GetCoords(x, y, z);
            break;
        }
        default:
            break;
    }


    Array<OneD, PointsKey> pkey1(2);
    Array<OneD, PointsType> ptype1(2);
    vector<BasisKey> bkey1;
  
    int numpts = 3;
    
    ptype1[0] =  eFourierEvenlySpaced;;
    ptype1[1] =  eFourierEvenlySpaced;
    pkey1[0] = PointsKey(numpts, ptype[0]);
    pkey1[1] = PointsKey(numpts, ptype[1]);
    bkey1.emplace_back(BasisKey(btype[0], order[0], pkey1[0]));
    bkey1.emplace_back(BasisKey(btype[1], order[1], pkey1[1]));

    //evalPtsxy[0] = Array<OneD, NekDouble>(numpts);
    //evalPtsxy[1] = Array<OneD, NekDouble>(numpts);
   
    Array<OneD,NekDouble> sol(numpts);
    Array<OneD,NekDouble> phys(numpts);


    Array<OneD, NekDouble> ret;
    if(( strcmp(ShapeTypeMap[stype],"Triangle") == 0
	 || strcmp(ShapeTypeMap[stype],"Quadrilateral") == 0
	 || strcmp(ShapeTypeMap[stype],"Prism") == 0 
         || strcmp(ShapeTypeMap[stype],"Tetrahedron") == 0) &&( baryinterp == 1))
    {
        if( strcmp(ShapeTypeMap[stype], "Triangle") == 0 )
        {
            Array< OneD, Array<OneD, NekDouble> >evalPtsxy(2);
            Array<OneD, Array< OneD, NekDouble > > allQuadxy(2);
    

            StdTriExp exp1(bkey[0],bkey[1]);
            StdTriExp exp2(bkey1[0],bkey1[1]);

            allQuadxy[0] = x;
            allQuadxy[1] = y;
            evalPtsxy[0] = Array<OneD, NekDouble>(exp2.GetTotPoints());
            evalPtsxy[1] = Array<OneD, NekDouble>(exp2.GetTotPoints());
            exp2.GetCoords(evalPtsxy[0], evalPtsxy[1]);


            Array<OneD, NekDouble> temp = EvalPoly(allQuadxy);
            for(int ii = 0; ii<evalPtsxy[0].num_elements(); ii++)
            {
                Array<OneD, NekDouble> xy(2);
                xy[0] = evalPtsxy[0][ii];
                xy[1] = evalPtsxy[1][ii];
                Array<OneD, NekDouble> c(2);
                //exp2.LocCoordToLocCollapsed(xy,c);
                c = xy;
                evalPtsxy[0][ii] = c[0];
                evalPtsxy[1][ii] = c[1];

            }

            
            ret = commoncode(&exp1, temp, evalPtsxy);
            
            sol  = EvalPoly(evalPtsxy);

            phys = ret;
            
            cout << "\nL infinity error: \t" << exp2.Linf(ret, sol) << endl;
            cout << "L 2 error: \t \t" << exp2.L2(ret, sol) << endl;
    


        }
        else  if( strcmp(ShapeTypeMap[stype], "Quadrilateral") == 0 )
        {
            Array< OneD, Array<OneD, NekDouble> >evalPtsxy(2);
            Array<OneD, Array< OneD, NekDouble > > allQuadxy(2);
    

            StdQuadExp exp1(bkey[0],bkey[1]);
            StdQuadExp exp2(bkey1[0],bkey1[1]);

            allQuadxy[0] = x;
            allQuadxy[1] = y;
            evalPtsxy[0] = Array<OneD, NekDouble>(exp2.GetTotPoints());
            evalPtsxy[1] = Array<OneD, NekDouble>(exp2.GetTotPoints());
            exp2.GetCoords(evalPtsxy[0], evalPtsxy[1]);


            Array<OneD, NekDouble> temp = EvalPoly(allQuadxy);
            for(int ii = 0; ii<evalPtsxy[0].num_elements(); ii++)
            {
                Array<OneD, NekDouble> xy(2);
                xy[0] = evalPtsxy[0][ii];
                xy[1] = evalPtsxy[1][ii];
                Array<OneD, NekDouble> c(2);
                exp2.LocCoordToLocCollapsed(xy,c);
                evalPtsxy[0][ii] = c[0];
                evalPtsxy[1][ii] = c[1];
            }

            
            ret = commoncode(&exp1, temp, evalPtsxy);
            
            sol  = EvalPoly(evalPtsxy);

            phys = ret;
    
            cout << "\nL infinity error: \t" << exp2.Linf(ret, sol) << endl;
            cout << "L 2 error: \t \t" << exp2.L2(ret, sol) << endl;

        }
        else if( strcmp(ShapeTypeMap[stype], "Tetrahedron") == 0 )
        {

            Array<OneD, PointsKey> pkey11(3);
            Array<OneD, PointsType> ptype11(3);
            vector<BasisKey> bkey11;
                
            int numpts = 4;
            
            ptype11[0] =  eFourierEvenlySpaced;
            ptype11[1] =  eFourierEvenlySpaced;
            ptype11[2] =  eFourierEvenlySpaced;
            pkey11[0] = PointsKey(numpts, ptype11[0]);
            pkey11[1] = PointsKey(numpts, ptype11[1]);
            pkey11[2] = PointsKey(numpts, ptype11[2]);
            bkey11.emplace_back(BasisKey(btype[0], order[0], pkey11[0]));
            bkey11.emplace_back(BasisKey(btype[1], order[1], pkey11[1]));
            bkey11.emplace_back(BasisKey(btype[2], order[2], pkey11[2]));
            
            Array< OneD, Array<OneD, NekDouble> >evalPtsxy(3);
            Array<OneD, Array< OneD, NekDouble > > allQuadxy(3);
    
            StdTetExp exp1(bkey[0],bkey[1],bkey[2]);
            StdTetExp exp2(bkey11[0],bkey11[1],bkey11[2]);

            allQuadxy[0] = x;
            allQuadxy[1] = y;
            allQuadxy[2] = z;
            evalPtsxy[0] = Array<OneD, NekDouble>(exp2.GetTotPoints());
            evalPtsxy[1] = Array<OneD, NekDouble>(exp2.GetTotPoints());
            evalPtsxy[2] = Array<OneD, NekDouble>(exp2.GetTotPoints());
    
            exp2.GetCoords(evalPtsxy[0], evalPtsxy[1], evalPtsxy[2]);


            Array<OneD, NekDouble> temp = EvalPoly(allQuadxy);
            for(int ii = 0; ii<evalPtsxy[0].num_elements(); ii++)
            {
                Array<OneD, NekDouble> xyz(3);
                xyz[0] = evalPtsxy[0][ii];
                xyz[1] = evalPtsxy[1][ii];
                xyz[2] = evalPtsxy[2][ii];
                Array<OneD, NekDouble> c(3);
                //exp2.LocCoordToLocCollapsed(xy,c);
                c = xyz;
                evalPtsxy[0][ii] = c[0];
                evalPtsxy[1][ii] = c[1];
                evalPtsxy[2][ii] = c[2];

            }

            
            ret = commoncode(&exp1, temp, evalPtsxy);
            
            sol  = EvalPoly(evalPtsxy);

            phys = ret;
            
            cout << "\nL infinity error: \t" << exp2.Linf(ret, sol) << endl;
            cout << "L 2 error: \t \t" << exp2.L2(ret, sol) << endl;
    


        }
        else if( strcmp(ShapeTypeMap[stype], "Prism") == 0 )
        {
            Array<OneD, PointsKey> pkey11(3);
            Array<OneD, PointsType> ptype11(3);
            vector<BasisKey> bkey11;
                
            int numpts = 4;
            
            ptype11[0] =  eFourierEvenlySpaced;
            ptype11[1] =  eFourierEvenlySpaced;
            ptype11[2] =  eFourierEvenlySpaced;
            pkey11[0] = PointsKey(numpts, ptype11[0]);
            pkey11[1] = PointsKey(numpts, ptype11[1]);
            pkey11[2] = PointsKey(numpts, ptype11[2]);
            bkey11.emplace_back(BasisKey(btype[0], order[0], pkey11[0]));
            bkey11.emplace_back(BasisKey(btype[1], order[1], pkey11[1]));
            bkey11.emplace_back(BasisKey(btype[2], order[2], pkey11[2]));
            
            Array< OneD, Array<OneD, NekDouble> >evalPtsxy(3);
            Array<OneD, Array< OneD, NekDouble > > allQuadxy(3);
    
            StdPrismExp exp1(bkey[0],bkey[1],bkey[2]);
            StdPrismExp exp2(bkey11[0],bkey11[1],bkey11[2]);

            allQuadxy[0] = x;
            allQuadxy[1] = y;
            allQuadxy[2] = z;
            evalPtsxy[0] = Array<OneD, NekDouble>(exp2.GetTotPoints());
            evalPtsxy[1] = Array<OneD, NekDouble>(exp2.GetTotPoints());
            evalPtsxy[2] = Array<OneD, NekDouble>(exp2.GetTotPoints());
    
            exp2.GetCoords(evalPtsxy[0], evalPtsxy[1], evalPtsxy[2]);


            Array<OneD, NekDouble> temp = EvalPoly(allQuadxy);
            for(int ii = 0; ii<evalPtsxy[0].num_elements(); ii++)
            {
                Array<OneD, NekDouble> xyz(3);
                xyz[0] = evalPtsxy[0][ii];
                xyz[1] = evalPtsxy[1][ii];
                xyz[2] = evalPtsxy[2][ii];
                Array<OneD, NekDouble> c(3);
                //exp2.LocCoordToLocCollapsed(xy,c);
                c = xyz;
                evalPtsxy[0][ii] = c[0];
                evalPtsxy[1][ii] = c[1];
                evalPtsxy[2][ii] = c[2];

            }

            
            ret = commoncode(&exp1, temp, evalPtsxy);
            
            sol  = EvalPoly(evalPtsxy);

            phys = ret;
            
            cout << "\nL infinity error: \t" << exp2.Linf(ret, sol) << endl;
            cout << "L 2 error: \t \t" << exp2.L2(ret, sol) << endl;
    
            
        }

        else
        {
            cout<<"\n Error! Please enter Triangle or Quads only";
            exit(0);
        }

    }
    return 0;
}

