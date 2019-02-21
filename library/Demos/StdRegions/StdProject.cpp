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
// License for the specific language governing rights and limitations under
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
// Description: Demo for testing functionality of LocProject
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

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

NekDouble Shape_sol(NekDouble x, NekDouble y, NekDouble z, vector<int> order,
                    vector<BasisType> btype, ShapeType stype, bool diff);

//Modification to deal with exact solution for diff. Return 1 if integer < 0.
static double pow_loc(const double val, const int i)
{
    return (i < 0) ? 1.0 : pow(val, i);
}

int main(int argc, char *argv[])
{
    string shape, ntype;
    vector<string> basis(3, "NoBasisType");
    vector<int> order, points;
    po::options_description desc("Available options");
    desc.add_options()("help,h",
                       "Produce this help message and list "
                       "basis and shape types.")
            ("nodal,n", po::value<string>(&ntype),
             "Optional nodal type, autofills shape and basis choices.")
            ("shape,s", po::value<string>(&shape),
             "Region shape to project function on.")
            ("basis,b", po::value<vector<string>>(&basis)->multitoken(),
             "Basis type, separate by spaces for higher dimensions.")
            ("order,o", po::value<vector<int>>(&order)->multitoken()
                     ->required(),
             "Order of basis sets, separate by spaces for higher dimensions.")
            ("points,p", po::value<vector<int>>(&points)->multitoken()
                     ->required(),
             "Number of quadrature points, separate by spaces for "
             "higher dimensions.")
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
            return 1;
        }
        po::notify(vm);
    }
    catch (const exception &e)
    {
        cerr << "Error: " << e.what() << endl;
        cerr << desc;
        return 0;
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
                nodaltype = (PointsType) i;
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
                                      "' is invalid."));
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
                stype = (ShapeType) i;
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
    int dimension = ShapeTypeDimMap[stype];
    ASSERTL0(order.size() == dimension,
             "Number of orders supplied should match shape dimension");
    ASSERTL0(points.size() == dimension,
             "Number of points supplied should match shape dimension");

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
                    btype[i] = (BasisType) j;
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
    allowableBasis[eSegment] = {
            {eOrtho_A, eModified_A, eFourier, eGLL_Lagrange, eGauss_Lagrange,
                    eLegendre, eChebyshev, eMonomial, eFourierSingleMode,
                    eFourierHalfModeRe, eFourierHalfModeIm}
    };
    allowableBasis[eTriangle] = {
            {eOrtho_A, eModified_A, eGLL_Lagrange, eGauss_Lagrange},
            {eOrtho_B, eModified_B, eGLL_Lagrange, eGauss_Lagrange}
    };
    allowableBasis[eQuadrilateral] = {
            {eOrtho_A, eModified_A, eFourier, eGLL_Lagrange, eGauss_Lagrange,
                    eLegendre, eChebyshev, eMonomial, eFourierSingleMode,
                    eFourierHalfModeRe, eFourierHalfModeIm},
            {eOrtho_A, eModified_A, eFourier, eGLL_Lagrange, eGauss_Lagrange,
                    eLegendre, eChebyshev, eMonomial, eFourierSingleMode,
                    eFourierHalfModeRe, eFourierHalfModeIm}
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
        for (int j = 0; j < (allowableBasis[stype][i].size()); ++j)
        {
            if (allowableBasis[stype][i][j] == btype[i])
            {
                break;
            }
            ASSERTL0(j != allowableBasis[stype][i].size() - 1,
                     ("The basis type '" + (string) BasisTypeMap[btype[i]] +
                      "' is invalid for basis argument " + to_string(i + 1) +
                      " for shape '" + ShapeTypeMap[stype] + "'."))
        }
    }

    //Declaration of other variables needed
    vector<PointsType> pointstype(3, eNoPointsType);
    StdExpansion *E = nullptr;

    //Assign points type according to basis type selection
    for (int i = 0; i < dimension; ++i)
    {
        if (btype[i] == eFourier)
        {
            pointstype[i] = eFourierEvenlySpaced;
        }
        else if (btype[i] == eFourierSingleMode ||
                 btype[i] == eFourierHalfModeRe ||
                 btype[i] == eFourierHalfModeIm)
        {
            pointstype[i] = eFourierSingleModeSpaced;
        }
        else
        {
            if (i == 1 && (stype == eTriangle || stype == eTetrahedron))
            {
                pointstype[i] = eGaussRadauMAlpha1Beta0;
            }
            else if (i == 2 && (stype == eTetrahedron || stype == ePyramid))
            {
                pointstype[i] = eGaussRadauMAlpha2Beta0;
            }
            else if (i == 2 && stype == ePrism)
            {
                pointstype[i] = eGaussRadauMAlpha1Beta0;
            }
            else
            {
                pointstype[i] = eGaussLobattoLegendre;
            }
        }
    }

    vector<PointsKey> pkey;
    vector<BasisKey> bkey;
    for (int i = 0; i < dimension; ++i)
    {
        pkey.push_back(PointsKey(points[i], pointstype[i]));
        bkey.push_back(BasisKey(btype[i], order[i], pkey[i]));
    }

    switch (stype)
    {
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

    Array<OneD, NekDouble> x = Array<OneD, NekDouble>(
            (unsigned) E->GetTotPoints());
    Array<OneD, NekDouble> y = Array<OneD, NekDouble>(
            (unsigned) E->GetTotPoints());
    Array<OneD, NekDouble> z = Array<OneD, NekDouble>(
            (unsigned) E->GetTotPoints());
    Array<OneD, NekDouble> dx = Array<OneD, NekDouble>(
            (unsigned) E->GetTotPoints());
    Array<OneD, NekDouble> dy = Array<OneD, NekDouble>(
            (unsigned) E->GetTotPoints());
    Array<OneD, NekDouble> dz = Array<OneD, NekDouble>(
            (unsigned) E->GetTotPoints());
    Array<OneD, NekDouble> sol = Array<OneD, NekDouble>(
            (unsigned) E->GetTotPoints());

    switch (dimension)
    {
        case 0:
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

    //get solution array
    for (int i = 0; i < E->GetTotPoints(); ++i)
    {
        sol[i] = Shape_sol(x[i], y[i], z[i], order, btype, stype, 0);
    }

    Array<OneD, NekDouble> phys((unsigned) E->GetTotPoints());
    Array<OneD, NekDouble> coeffs((unsigned) E->GetNcoeffs());

    if (vm.count("diff"))
    {
        if (dimension == 1)
        {
            E->PhysDeriv(sol, sol);
        }
        else if (dimension == 2)
        {
            E->PhysDeriv(sol, dx, dy);
            Vmath::Vadd(E->GetTotPoints(), dx, 1, dy, 1, sol, 1);
        }
        else if (dimension == 3)
        {
            E->PhysDeriv(sol, dx, dy, dz);
            Vmath::Vadd(E->GetTotPoints(), dx, 1, dy, 1, sol, 1);
            Vmath::Vadd(E->GetTotPoints(), dz, 1, sol, 1, sol, 1);
        }
    }

    //Project onto expansion
    E->FwdTrans(sol, coeffs);

    //Backward transform solution to get projected values
    E->BwdTrans(coeffs, phys);

    if (vm.count("diff"))
    {
        for (int i = 0; i < E->GetTotPoints(); ++i)
        {
            sol[i] = Shape_sol(x[i], y[i], z[i], order, btype, stype, 1);
        }
    }

    //Calculate L_inf & L_2 error
    cout << "L infinity error: \t" << E->Linf(phys, sol) << endl;
    if (stype != ePoint)
    {
        cout << "L 2 error: \t \t \t" << E->L2(phys, sol) << endl;
    }

    if (!vm.count("diff"))
    {
        //Evaluate solution at x = y = 0 and print error
        Array<OneD, NekDouble> t = Array<OneD, NekDouble>(3);
        t[0] = -0.5;
        t[1] = -0.25;
        t[2] = -0.3;
        sol[0] = Shape_sol(t[0], t[1], t[2], order, btype, stype, 0);

        NekDouble nsol = E->PhysEvaluate(t, phys);

        cout << "Error at x = (";
        for (int i = 0; i < dimension; ++i)
        {
            cout << t[i] << ", ";
        }
        cout << "\b\b): " << nsol - sol[0] << endl;
    }

    return 0;
}

NekDouble Shape_sol(NekDouble x, NekDouble y, NekDouble z, vector<int> order,
                    vector<BasisType> btype, ShapeType stype, bool diff)
{
    map<ShapeType, function<int(int &, vector<int> &)>> shapeConstraint2;
    shapeConstraint2[eSegment] = [](int &k, vector<int> &order) { return 1; };
    shapeConstraint2[eTriangle] = [](int &k, vector<int> &order) {
        return order[1] - k; };
    shapeConstraint2[eQuadrilateral] = [](int &k, vector<int> &order) {
        return order[1]; };
    shapeConstraint2[eTetrahedron] = [](int &k, vector<int> &order) {
        return order[1] - k; };
    shapeConstraint2[ePyramid] = [](int &k, vector<int> &order) {
        return order[1] - k; };
    shapeConstraint2[ePrism] = [](int &k, vector<int> &order) {
        return order[1]; };
    shapeConstraint2[eHexahedron] = [](int &k, vector<int> &order) {
        return order[1]; };

    map<ShapeType, function<int(int &, int &,
                                vector<int> &order)>> shapeConstraint3;
    shapeConstraint3[eSegment] = [](int &k, int &l,
                                    vector<int> &order) { return 1; };
    shapeConstraint3[eTriangle] = [](int &k, int &l,
                                     vector<int> &order) { return 1; };
    shapeConstraint3[eQuadrilateral] = [](int &k, int &l,
                                          vector<int> &order) { return 1; };
    shapeConstraint3[eTetrahedron] = [](int &k, int &l, vector<int> &order) {
        return order[2] - k - l; };
    shapeConstraint3[ePyramid] = [](int &k, int &l, vector<int> &order) {
        return order[2] - k - l; };
    shapeConstraint3[ePrism] = [](int &k, int &l, vector<int> &order) {
        return order[2] - k; };
    shapeConstraint3[eHexahedron] = [](int &k, int &l, vector<int> &order) {
        return order[2]; };

    NekDouble sol = 0.0;
    if (!diff)
    {
        if (btype[0] == eFourier && stype == ePoint)
        {
            for (int k = 0; k < order[0] / 2 - 1; ++k)
            {
                sol += sin(k * M_PI * x) + cos(k * M_PI * x);
            }
        }
        else if (btype[0] == eFourierSingleMode && stype == ePoint)
        {
            sol += 0.45 * sin(M_PI * x) + 0.25 * cos(M_PI * x);
        }
        else if (btype[0] == eFourierHalfModeRe && stype == ePoint)
        {
            sol += 0.45 * cos(M_PI * x);
        }
        else if (btype[0] == eFourierHalfModeIm && stype == ePoint)
        {
            sol += 0.25 * sin(M_PI * x);
        }
        else if (btype[0] == eFourier && stype == eSegment)
        {
            for (int k = 0; k < order[0] / 2 - 1; ++k)
            {
                sol += sin(k * M_PI * x) + cos(k * M_PI * x);
            }
        }
        else if (btype[0] == eFourierSingleMode && stype == eSegment)
        {
            sol += 0.25 * sin(M_PI * x) + 0.25 * cos(M_PI * x);
        }
        else if (btype[0] == eFourier && stype == eQuadrilateral)
        {
            if (btype[1] == eFourier)
            {
                for (int k = 0; k < order[0] / 2; ++k)
                {
                    for (int l = 0; l < order[1] / 2; ++l)
                    {
                        sol += sin(k * M_PI * x) * sin(l * M_PI * y) +
                               sin(k * M_PI * x) * cos(l * M_PI * y) +
                               cos(k * M_PI * x) * sin(l * M_PI * y) +
                               cos(k * M_PI * x) * cos(l * M_PI * y);
                    }
                }
            }
            else if (btype[1] == eFourierSingleMode)
            {
                for (int k = 0; k < order[0] / 2; ++k)
                {
                    sol += sin(k * M_PI * x) * sin(M_PI * y) +
                           sin(k * M_PI * x) * cos(M_PI * y) +
                           cos(k * M_PI * x) * sin(M_PI * y) +
                           cos(k * M_PI * x) * cos(M_PI * y);
                }
            }
            else
            {
                for (int k = 0; k < order[0] / 2; ++k)
                {
                    for (int l = 0; l < order[1]; ++l)
                    {
                        sol += sin(k * M_PI * x) * pow_loc(y, l) +
                               cos(k * M_PI * x) * pow_loc(y, l);
                    }
                }
            }
        }
        else if (btype[0] == eFourierSingleMode && stype == eQuadrilateral)
        {
            if (btype[1] == eFourier)
            {
                for (int l = 0; l < order[1] / 2; ++l)
                {
                    sol += sin(M_PI * x) * sin(l * M_PI * y) +
                           sin(M_PI * x) * cos(l * M_PI * y) +
                           cos(M_PI * x) * sin(l * M_PI * y) +
                           cos(M_PI * x) * cos(l * M_PI * y);
                }

            }
            else if (btype[1] == eFourierSingleMode)
            {
                sol += sin(M_PI * x) * sin(M_PI * y) +
                       sin(M_PI * x) * cos(M_PI * y) +
                       cos(M_PI * x) * sin(M_PI * y) +
                       cos(M_PI * x) * cos(M_PI * y);
            }
            else
            {
                for (int l = 0; l < order[1]; ++l)
                {
                    sol += sin(M_PI * x) * pow_loc(y, l) +
                           cos(M_PI * x) * pow_loc(y, l);
                }
            }
        }
        else if (btype[1] == eFourier && stype == eQuadrilateral)
        {
            for (int k = 0; k < order[0]; ++k)
            {
                for (int l = 0; l < order[1] / 2; ++l)
                {
                    sol += sin(l * M_PI * y) * pow_loc(x, k) +
                           cos(l * M_PI * y) * pow_loc(x, k);
                }
            }
        }
        else if (btype[1] == eFourierSingleMode && stype == eQuadrilateral)
        {
            for (int k = 0; k < order[0]; ++k)
            {
                sol += sin(M_PI * y) * pow_loc(x, k) +
                       cos(M_PI * y) * pow_loc(x, k);
            }
        }
        else
        {
            for (int k = 0;
                 k < order[0]; ++k) //ShapeConstraint 1 is always < order1
            {
                for (int l = 0; l < shapeConstraint2[stype](k, order); ++l)
                {
                    for (int m = 0;
                         m < shapeConstraint3[stype](k, l, order); ++m)
                    {
                        sol += pow_loc(x, k) * pow_loc(y, l) * pow_loc(z, m);
                    }
                }
            }
        }
    }
    else if (diff)
    {
        if (btype[0] == eFourier && stype == eSegment)
        {
            for (int k = 0; k < order[0] / 2 - 1; ++k)
            {
                sol += k * M_PI * (cos(k * M_PI * z) - sin(k * M_PI * z));
            }
        }
        else if (btype[0] != eFourier && btype[1] == eFourier &&
                 stype == eQuadrilateral)
        {
            for (int k = 0; k < order[0]; ++k)
            {
                for (int l = 0; l < order[1] / 2; ++l)
                {
                    sol += k * pow_loc(x, k - 1) * sin(M_PI * l * y)
                           + M_PI * l * pow_loc(x, k) * cos(M_PI * l * y) +
                           +k * pow_loc(x, k - 1) * cos(M_PI * l * y)
                           - M_PI * l * pow_loc(x, k) * sin(M_PI * l * y);
                }
            }
        }
        else if (btype[0] == eFourier && btype[1] != eFourier &&
                 stype == eQuadrilateral)
        {
            for (int k = 0; k < order[0] / 2; ++k)
            {
                for (int l = 0; l < order[1]; ++l)
                {
                    sol += M_PI * k * cos(M_PI * k * x) * pow_loc(y, l)
                           + l * sin(M_PI * k * x) * pow_loc(y, l - 1) +
                           -M_PI * k * sin(M_PI * k * x) * pow_loc(y, l)
                           + l * sin(M_PI * k * x) * pow_loc(y, l - 1);
                }
            }
        }
        else if (btype[0] == eFourier && btype[1] == eFourier &&
                 stype == eQuadrilateral)
        {
            for (int k = 0; k < order[0] / 2; ++k)
            {
                for (int l = 0; l < order[1] / 2; ++l)
                {
                    sol += M_PI * k * cos(M_PI * k * x) * sin(M_PI * l * y)
                           + M_PI * l * sin(M_PI * k * x) * cos(M_PI * l * y)
                           + M_PI * k * cos(M_PI * k * x) * cos(M_PI * l * y)
                           - M_PI * l * sin(M_PI * k * x) * sin(M_PI * l * y)
                           - M_PI * k * sin(M_PI * k * x) * sin(M_PI * l * y)
                           + M_PI * l * cos(M_PI * k * x) * cos(M_PI * l * y)
                           - M_PI * k * sin(M_PI * k * x) * cos(M_PI * l * y)
                           - M_PI * l * cos(M_PI * k * x) * sin(M_PI * l * y);
                }
            }
        }
        else
        {
            NekDouble a;
            for (int k = 0;
                 k < order[0]; ++k) //ShapeConstraint 1 is always < order1
            {
                for (int l = 0; l < shapeConstraint2[stype](k, order); ++l)
                {
                    for (int m = 0;
                         m < shapeConstraint3[stype](k, l, order); ++m)
                    {
                        a = k * pow_loc(x, k - 1) * pow_loc(y, l) *
                            pow_loc(z, m);
                        sol += a;
                        a = l * pow_loc(x, k) * pow_loc(y, l - 1) *
                            pow_loc(z, m);
                        sol += a;
                        a = m * pow_loc(x, k) * pow_loc(y, l) *
                            pow_loc(z, m - 1);
                        sol += a;
                    }
                }
            }
        }
    }

    return sol;
}
