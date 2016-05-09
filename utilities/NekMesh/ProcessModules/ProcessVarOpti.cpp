////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessJac.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: Calculate Jacobians of elements.
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/ManagerAccess.h>

#include <NekMeshUtils/MeshElements/Element.h>
#include "ProcessVarOpti.h"

#include <boost/thread/mutex.hpp>

#include <StdRegions/StdTriExp.h>
#include <StdRegions/StdQuadExp.h>
#include <StdRegions/StdTetExp.h>
#include <StdRegions/StdPrismExp.h>

#include <LibUtilities/Foundations/NodalUtil.h>

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{

boost::mutex mtx;
NekMatrix<NekDouble> interp;
int ptsLow;
int ptsHigh;
int dim;
NekMatrix<NekDouble> VdmDx;
NekMatrix<NekDouble> VdmDy;
NekMatrix<NekDouble> VdmDz;
NekVector<NekDouble> quadW;
optimiser opti;

inline NekDouble GetElFunctional(ElDataSharedPtr d)
{
    vector<NodeSharedPtr> ns;
    d->el->GetCurvedNodes(ns);

    ASSERTL0(ptsHigh == d->maps.size(), "what");

    Array<OneD, Array<OneD, NekDouble> > jac(ptsHigh);

    if(dim == 2)
    {
        NekVector<NekDouble> X(ptsLow),Y(ptsLow),
                             x1(ptsLow),y1(ptsLow),
                             x2(ptsLow),y2(ptsLow);
        for(int i = 0; i < ptsLow; i++)
        {
            X(i) = ns[i]->m_x;
            Y(i) = ns[i]->m_y;
        }

        NekVector<NekDouble> x1i(ptsHigh),y1i(ptsHigh),
                             x2i(ptsHigh),y2i(ptsHigh);

        x1 = VdmDx*X;
        y1 = VdmDx*Y;
        x2 = VdmDy*X;
        y2 = VdmDy*Y;

        x1i = interp * x1;
        x2i = interp * x2;
        y1i = interp * y1;
        y2i = interp * y2;

        for(int i = 0; i < ptsHigh; i++)
        {
            Array<OneD, NekDouble> jaci(9,0.0);
            jaci[0] = x1i(i);
            jaci[1] = y1i(i);
            jaci[3] = x2i(i);
            jaci[4] = y2i(i);
            jac[i] = jaci;

        }

    }
    else
    {
        NekVector<NekDouble> X(ptsLow),Y(ptsLow),Z(ptsLow),
                             x1(ptsHigh),y1(ptsHigh),z1(ptsHigh),
                             x2(ptsHigh),y2(ptsHigh),z2(ptsHigh),
                             x3(ptsHigh),y3(ptsHigh),z3(ptsHigh);
        for(int i = 0; i < ptsLow; i++)
        {
            X(i) = ns[i]->m_x;
            Y(i) = ns[i]->m_y;
            Z(i) = ns[i]->m_z;
        }

        x1 = VdmDx*X;
        y1 = VdmDx*Y;
        z1 = VdmDx*Z;
        x2 = VdmDy*X;
        y2 = VdmDy*Y;
        z2 = VdmDy*Z;
        x3 = VdmDz*X;
        y3 = VdmDz*Y;
        z3 = VdmDz*Z;

        NekVector<NekDouble> x1i(ptsHigh),y1i(ptsHigh),z1i(ptsHigh),
                             x2i(ptsHigh),y2i(ptsHigh),z2i(ptsHigh),
                             x3i(ptsHigh),y3i(ptsHigh),z3i(ptsHigh);

        x1i = interp * x1;
        x2i = interp * x2;
        x3i = interp * x3;
        y1i = interp * y1;
        y2i = interp * y2;
        y3i = interp * y3;
        z1i = interp * z1;
        z2i = interp * z2;
        z3i = interp * z3;

        for(int i = 0; i < ptsHigh; i++)
        {
            Array<OneD, NekDouble> jaci(9,0.0);
            jaci[0] = x1i(i);
            jaci[1] = y1i(i);
            jaci[2] = z1i(i);
            jaci[3] = x2i(i);
            jaci[4] = y2i(i);
            jaci[5] = z2i(i);
            jaci[6] = x3i(i);
            jaci[7] = y3i(i);
            jaci[8] = z3i(i);
            jac[i] = jaci;

        }
    }

    Array<OneD, NekDouble> dW(ptsHigh);

    switch (opti)
    {
        case eLinEl:
        {
            for(int k = 0; k < ptsHigh; k++)
            {
                Array<OneD, NekDouble> jacIdeal(9,0.0);
                jacIdeal[0] = jac[k][0]*d->maps[k][0] + jac[k][3]*d->maps[k][1] + jac[k][6]*d->maps[k][2];
                jacIdeal[1] = jac[k][1]*d->maps[k][0] + jac[k][4]*d->maps[k][1] + jac[k][7]*d->maps[k][2];
                jacIdeal[2] = jac[k][2]*d->maps[k][0] + jac[k][5]*d->maps[k][1] + jac[k][8]*d->maps[k][2];

                jacIdeal[3] = jac[k][0]*d->maps[k][3] + jac[k][3]*d->maps[k][4] + jac[k][6]*d->maps[k][5];
                jacIdeal[4] = jac[k][1]*d->maps[k][3] + jac[k][4]*d->maps[k][4] + jac[k][7]*d->maps[k][5];
                jacIdeal[5] = jac[k][2]*d->maps[k][3] + jac[k][5]*d->maps[k][4] + jac[k][8]*d->maps[k][5];

                jacIdeal[6] = jac[k][0]*d->maps[k][6] + jac[k][3]*d->maps[k][7] + jac[k][6]*d->maps[k][8];
                jacIdeal[7] = jac[k][1]*d->maps[k][6] + jac[k][4]*d->maps[k][7] + jac[k][7]*d->maps[k][8];
                jacIdeal[8] = jac[k][2]*d->maps[k][6] + jac[k][5]*d->maps[k][7] + jac[k][8]*d->maps[k][8];

                NekDouble a = (jacIdeal[0]*jacIdeal[0]+jacIdeal[1]*jacIdeal[1]+jacIdeal[2]*jacIdeal[2]-1.0)*
                              (jacIdeal[0]*jacIdeal[0]+jacIdeal[1]*jacIdeal[1]+jacIdeal[2]*jacIdeal[2]-1.0);
                NekDouble b = (jacIdeal[3]*jacIdeal[3]+jacIdeal[4]*jacIdeal[4]+jacIdeal[5]*jacIdeal[5]-1.0)*
                              (jacIdeal[3]*jacIdeal[3]+jacIdeal[4]*jacIdeal[4]+jacIdeal[5]*jacIdeal[5]-1.0);
                NekDouble c = (jacIdeal[6]*jacIdeal[6]+jacIdeal[7]*jacIdeal[7]+jacIdeal[8]*jacIdeal[8]-1.0)*
                              (jacIdeal[6]*jacIdeal[6]+jacIdeal[7]*jacIdeal[7]+jacIdeal[8]*jacIdeal[8]-1.0);
                NekDouble D = (jacIdeal[0]*jacIdeal[6]+jacIdeal[1]*jacIdeal[7]+jacIdeal[2]*jacIdeal[8])*
                              (jacIdeal[0]*jacIdeal[6]+jacIdeal[1]*jacIdeal[7]+jacIdeal[2]*jacIdeal[8]);
                NekDouble e = (jacIdeal[3]*jacIdeal[6]+jacIdeal[4]*jacIdeal[7]+jacIdeal[5]*jacIdeal[8])*
                              (jacIdeal[3]*jacIdeal[6]+jacIdeal[4]*jacIdeal[7]+jacIdeal[5]*jacIdeal[8]);
                NekDouble f = (jacIdeal[0]*jacIdeal[3]+jacIdeal[1]*jacIdeal[4]+jacIdeal[3]*jacIdeal[5])*
                              (jacIdeal[0]*jacIdeal[3]+jacIdeal[1]*jacIdeal[4]+jacIdeal[3]*jacIdeal[5]);

                NekDouble trEtE = 0.25 * (a+b+c) + 0.5 * (D+e+f);

                NekDouble jacDet;
                if(dim == 2)
                {
                    jacDet = jacIdeal[0] * jacIdeal[4] - jacIdeal[3]*jacIdeal[1];
                }
                else
                {
                    jacDet = jacIdeal[0]*(jacIdeal[4]*jacIdeal[8]-jacIdeal[5]*jacIdeal[7])
                            -jacIdeal[3]*(jacIdeal[1]*jacIdeal[8]-jacIdeal[2]*jacIdeal[7])
                            +jacIdeal[6]*(jacIdeal[1]*jacIdeal[5]-jacIdeal[2]*jacIdeal[4]);
                }

                NekDouble ljacDet = log(jacDet);
                NekDouble nu = 0.45;
                NekDouble mu = 1.0/2.0/(1.0+nu);
                NekDouble K = 1.0 / 3.0 / (1.0 - 2.0 * nu);

                if(jacDet > 0)
                {
                    dW[k] = K *0.5 * ljacDet * ljacDet + mu * trEtE;
                }
                else
                {
                    NekDouble de = fabs(d->maps[k][9]) * sqrt(1E-3*1E-3 + 1E-3);
                    NekDouble sigma = 0.5*(jacDet + sqrt(jacDet*jacDet + 4.0*de*de));
                    NekDouble lsigma = log(sigma);
                    dW[k] = K *0.5 * lsigma * lsigma + mu * trEtE;
                }

            }
            break;
        }
        case eHypEl:
        {
            for(int k = 0; k < ptsHigh; k++)
            {
                Array<OneD, NekDouble> jacIdeal(9,0.0);
                jacIdeal[0] = jac[k][0]*d->maps[k][0] + jac[k][3]*d->maps[k][1] + jac[k][6]*d->maps[k][2];
                jacIdeal[1] = jac[k][1]*d->maps[k][0] + jac[k][4]*d->maps[k][1] + jac[k][7]*d->maps[k][2];
                jacIdeal[2] = jac[k][2]*d->maps[k][0] + jac[k][5]*d->maps[k][1] + jac[k][8]*d->maps[k][2];

                jacIdeal[3] = jac[k][0]*d->maps[k][3] + jac[k][3]*d->maps[k][4] + jac[k][6]*d->maps[k][5];
                jacIdeal[4] = jac[k][1]*d->maps[k][3] + jac[k][4]*d->maps[k][4] + jac[k][7]*d->maps[k][5];
                jacIdeal[5] = jac[k][2]*d->maps[k][3] + jac[k][5]*d->maps[k][4] + jac[k][8]*d->maps[k][5];

                jacIdeal[6] = jac[k][0]*d->maps[k][6] + jac[k][3]*d->maps[k][7] + jac[k][6]*d->maps[k][8];
                jacIdeal[7] = jac[k][1]*d->maps[k][6] + jac[k][4]*d->maps[k][7] + jac[k][7]*d->maps[k][8];
                jacIdeal[8] = jac[k][2]*d->maps[k][6] + jac[k][5]*d->maps[k][7] + jac[k][8]*d->maps[k][8];

                NekDouble I1 = jacIdeal[0]*jacIdeal[0] +
                               jacIdeal[1]*jacIdeal[1] +
                               jacIdeal[2]*jacIdeal[2] +
                               jacIdeal[3]*jacIdeal[3] +
                               jacIdeal[4]*jacIdeal[4] +
                               jacIdeal[5]*jacIdeal[5] +
                               jacIdeal[6]*jacIdeal[6] +
                               jacIdeal[7]*jacIdeal[7] +
                               jacIdeal[8]*jacIdeal[8];

                NekDouble jacDet;
                if(dim == 2)
                {
                    jacDet = jacIdeal[0] * jacIdeal[4] - jacIdeal[3]*jacIdeal[1];
                }
                else
                {
                    jacDet = jacIdeal[0]*(jacIdeal[4]*jacIdeal[8]-jacIdeal[5]*jacIdeal[7])
                            -jacIdeal[3]*(jacIdeal[1]*jacIdeal[8]-jacIdeal[2]*jacIdeal[7])
                            +jacIdeal[6]*(jacIdeal[1]*jacIdeal[5]-jacIdeal[2]*jacIdeal[4]);
                }

                NekDouble ljacDet = log(jacDet);
                NekDouble nu = 0.45;
                NekDouble mu = 1.0/2.0/(1.0+nu);
                NekDouble K = 1.0 / 3.0 / (1.0 - 2.0 * nu);

                if(jacDet > 0)
                {
                    dW[k] = 0.5 * mu * (I1 - 3.0 - 2.0*ljacDet) + 0.5 * K * ljacDet * ljacDet;
                }
                else
                {
                    NekDouble de = fabs(d->maps[k][9]) * sqrt(1E-3*1E-3 + 1E-3);
                    NekDouble sigma = 0.5*(jacDet + sqrt(jacDet*jacDet + 4.0*de*de));
                    NekDouble lsigma = log(sigma);
                    dW[k] = 0.5 * mu * (I1 - 3.0 - 2.0*lsigma) + 0.5 * K * lsigma * lsigma;
                }

            }
            break;
        }
        case eRoca:
        {
            for(int k = 0; k < ptsHigh; k++)
            {
                Array<OneD, NekDouble> jacIdeal(9,0.0);
                jacIdeal[0] = jac[k][0]*d->maps[k][0] + jac[k][3]*d->maps[k][1] + jac[k][6]*d->maps[k][2];
                jacIdeal[1] = jac[k][1]*d->maps[k][0] + jac[k][4]*d->maps[k][1] + jac[k][7]*d->maps[k][2];
                jacIdeal[2] = jac[k][2]*d->maps[k][0] + jac[k][5]*d->maps[k][1] + jac[k][8]*d->maps[k][2];

                jacIdeal[3] = jac[k][0]*d->maps[k][3] + jac[k][3]*d->maps[k][4] + jac[k][6]*d->maps[k][5];
                jacIdeal[4] = jac[k][1]*d->maps[k][3] + jac[k][4]*d->maps[k][4] + jac[k][7]*d->maps[k][5];
                jacIdeal[5] = jac[k][2]*d->maps[k][3] + jac[k][5]*d->maps[k][4] + jac[k][8]*d->maps[k][5];

                jacIdeal[6] = jac[k][0]*d->maps[k][6] + jac[k][3]*d->maps[k][7] + jac[k][6]*d->maps[k][8];
                jacIdeal[7] = jac[k][1]*d->maps[k][6] + jac[k][4]*d->maps[k][7] + jac[k][7]*d->maps[k][8];
                jacIdeal[8] = jac[k][2]*d->maps[k][6] + jac[k][5]*d->maps[k][7] + jac[k][8]*d->maps[k][8];

                NekDouble jacDet;
                if(dim == 2)
                {
                    jacDet = jacIdeal[0] * jacIdeal[4] - jacIdeal[3]*jacIdeal[1];
                }
                else
                {
                    jacDet = jacIdeal[0]*(jacIdeal[4]*jacIdeal[8]-jacIdeal[5]*jacIdeal[7])
                            -jacIdeal[3]*(jacIdeal[1]*jacIdeal[8]-jacIdeal[2]*jacIdeal[7])
                            +jacIdeal[6]*(jacIdeal[1]*jacIdeal[5]-jacIdeal[2]*jacIdeal[4]);
                }

                NekDouble frob = 0.0;

                frob += jacIdeal[0] * jacIdeal[0];
                frob += jacIdeal[1] * jacIdeal[1];
                frob += jacIdeal[2] * jacIdeal[2];
                frob += jacIdeal[3] * jacIdeal[3];
                frob += jacIdeal[4] * jacIdeal[4];
                frob += jacIdeal[5] * jacIdeal[5];
                frob += jacIdeal[6] * jacIdeal[6];
                frob += jacIdeal[7] * jacIdeal[7];
                frob += jacIdeal[8] * jacIdeal[8];

                if(jacDet > 0)
                {
                    dW[k] = frob / dim / pow(fabs(jacDet), 2.0/dim);
                }
                else
                {
                    NekDouble de = fabs(d->maps[k][9]) * sqrt(1E-3*1E-3 + 1E-3);
                    NekDouble sigma = 0.5*(jacDet + sqrt(jacDet*jacDet + 4.0*de*de));
                    dW[k] = frob / dim / pow(fabs(sigma), 2.0/dim) -1.0;
                }
            }
            break;
        }
        case eWins:
        {
            for(int k = 0; k < ptsHigh; k++)
            {
                Array<OneD, NekDouble> jacIdeal(9,0.0);
                jacIdeal[0] = jac[k][0]*d->maps[k][0] + jac[k][3]*d->maps[k][1] + jac[k][6]*d->maps[k][2];
                jacIdeal[1] = jac[k][1]*d->maps[k][0] + jac[k][4]*d->maps[k][1] + jac[k][7]*d->maps[k][2];
                jacIdeal[2] = jac[k][2]*d->maps[k][0] + jac[k][5]*d->maps[k][1] + jac[k][8]*d->maps[k][2];

                jacIdeal[3] = jac[k][0]*d->maps[k][3] + jac[k][3]*d->maps[k][4] + jac[k][6]*d->maps[k][5];
                jacIdeal[4] = jac[k][1]*d->maps[k][3] + jac[k][4]*d->maps[k][4] + jac[k][7]*d->maps[k][5];
                jacIdeal[5] = jac[k][2]*d->maps[k][3] + jac[k][5]*d->maps[k][4] + jac[k][8]*d->maps[k][5];

                jacIdeal[6] = jac[k][0]*d->maps[k][6] + jac[k][3]*d->maps[k][7] + jac[k][6]*d->maps[k][8];
                jacIdeal[7] = jac[k][1]*d->maps[k][6] + jac[k][4]*d->maps[k][7] + jac[k][7]*d->maps[k][8];
                jacIdeal[8] = jac[k][2]*d->maps[k][6] + jac[k][5]*d->maps[k][7] + jac[k][8]*d->maps[k][8];

                NekDouble jacDet;
                if(dim == 2)
                {
                    jacDet = jacIdeal[0] * jacIdeal[4] - jacIdeal[3]*jacIdeal[1];
                }
                else
                {
                    jacDet = jacIdeal[0]*(jacIdeal[4]*jacIdeal[8]-jacIdeal[5]*jacIdeal[7])
                            -jacIdeal[3]*(jacIdeal[1]*jacIdeal[8]-jacIdeal[2]*jacIdeal[7])
                            +jacIdeal[6]*(jacIdeal[1]*jacIdeal[5]-jacIdeal[2]*jacIdeal[4]);
                }

                NekDouble frob = 0.0;

                frob += jacIdeal[0] * jacIdeal[0];
                frob += jacIdeal[1] * jacIdeal[1];
                frob += jacIdeal[2] * jacIdeal[2];
                frob += jacIdeal[3] * jacIdeal[3];
                frob += jacIdeal[4] * jacIdeal[4];
                frob += jacIdeal[5] * jacIdeal[5];
                frob += jacIdeal[6] * jacIdeal[6];
                frob += jacIdeal[7] * jacIdeal[7];
                frob += jacIdeal[8] * jacIdeal[8];

                if(jacDet > 0)
                {
                    dW[k] = frob / jacDet;
                }
                else
                {
                    NekDouble de = fabs(d->maps[k][9]) * sqrt(1E-3*1E-3 + 1E-3);
                    NekDouble sigma = 0.5*(jacDet + sqrt(jacDet*jacDet + 4.0*de*de));
                    dW[k] = frob / sigma;
                }
            }
            break;
        }
    }

    NekDouble integral = 0.0;
    for(int i = 0; i < ptsHigh; i++)
    {
        integral += quadW[i] * dW[i];
    }
    return integral;
}

ModuleKey ProcessVarOpti::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "varopti"),
    ProcessVarOpti::create,
    "Optimise mesh locations.");

ProcessVarOpti::ProcessVarOpti(MeshSharedPtr m) : ProcessModule(m)
{
    m_config["linearelastic"] =
        ConfigOption(true, "", "Optimise for linear elasticity");
    m_config["winslow"] =
        ConfigOption(true, "", "Optimise for winslow");
    m_config["roca"] =
        ConfigOption(true, "", "Optimise for roca method");
    m_config["hyperelastic"] =
        ConfigOption(true, "", "Optimise for hyper elasticity");
    m_config["numthreads"] =
        ConfigOption(false, "1", "Number of threads");
    m_config["nq"] =
        ConfigOption(false, "0", "Number of quad points");
    m_config["stats"] =
        ConfigOption(false, "", "Write a file with list of scaled jacobians");
}

ProcessVarOpti::~ProcessVarOpti()
{
}

void ProcessVarOpti::Process()
{
    if (m_mesh->m_verbose)
    {
        cout << "ProcessVarOpti: Optimising... " << endl;
    }

    if(m_config["linearelastic"].beenSet)
    {
        opti = eLinEl;
    }
    else if(m_config["winslow"].beenSet)
    {
        opti = eWins;
    }
    else if(m_config["roca"].beenSet)
    {
        opti = eRoca;
    }
    else if(m_config["hyperelastic"].beenSet)
    {
        opti = eHypEl;
    }
    else
    {
        ASSERTL0(false,"not opti type set");
    }

    m_mesh->m_nummode = m_config["nq"].as<int>();

    ASSERTL0(m_mesh->m_nummode > 2,"not specified high-order");

    if(m_mesh->m_expDim == 2 && m_mesh->m_spaceDim == 3)
    {
        ASSERTL0(false,"cannot deal with manifolds");
    }

    res = boost::shared_ptr<Residual>(new Residual);
    res->val = 1.0;

    FillQuadPoints();

    //build Vandermonde information
    dim = m_mesh->m_spaceDim;
    switch (m_mesh->m_spaceDim)
    {
        case 2:
        {
            ptsLow  = m_mesh->m_nummode*(m_mesh->m_nummode+1)/2;

            LibUtilities::PointsKey pkey1(m_mesh->m_nummode,
                                          LibUtilities::eNodalTriElec);
            LibUtilities::PointsKey pkey2(m_mesh->m_nummode+2,
                                          LibUtilities::eNodalTriSPI);
            Array<OneD, NekDouble> u1, v1, u2, v2;

            LibUtilities::PointsManager()[pkey1]->GetPoints(u1, v1);
            LibUtilities::PointsManager()[pkey2]->GetPoints(u2, v2);
            NekVector<NekDouble> U1(u1), V1(v1);
            NekVector<NekDouble> U2(u2), V2(v2);
            ptsHigh = LibUtilities::PointsManager()[pkey2]->GetNumPointsAlt();

            interp = LibUtilities::GetInterpolationMatrix(U1, V1, U2, V2);

            NekMatrix<NekDouble> Vandermonde = LibUtilities::GetVandermonde(U1,V1);
            NekMatrix<NekDouble> VandermondeI = Vandermonde;
            VandermondeI.Invert();
            VdmDx = LibUtilities::GetVandermondeForXDerivative(U1,V1) *
                                                                VandermondeI;
            VdmDy = LibUtilities::GetVandermondeForYDerivative(U1,V1) *
                                                                VandermondeI;
            //quadW = LibUtilities::MakeQuadratureWeights(U2,V1);
            Array<OneD, NekDouble> qds = LibUtilities::PointsManager()[pkey2]->GetW();
            NekVector<NekDouble> quadWi(qds);
            quadW = quadWi;
        }
        break;
        case 3:
        {
            ptsLow  = m_mesh->m_nummode*(m_mesh->m_nummode+1)*(m_mesh->m_nummode+2)/6;
            LibUtilities::PointsKey pkey1(m_mesh->m_nummode,
                                          LibUtilities::eNodalTetElec);
            LibUtilities::PointsKey pkey2(m_mesh->m_nummode+2,
                                          LibUtilities::eNodalTetSPI);
            Array<OneD, NekDouble> u1, v1, u2, v2, w1, w2;
            LibUtilities::PointsManager()[pkey1]->GetPoints(u1, v1, w1);
            LibUtilities::PointsManager()[pkey2]->GetPoints(u2, v2, w2);
            NekVector<NekDouble> U1(u1), V1(v1), W1(w1);
            NekVector<NekDouble> U2(u2), V2(v2), W2(w2);
            ptsHigh = LibUtilities::PointsManager()[pkey2]->GetNumPointsAlt();

            interp = LibUtilities::GetTetInterpolationMatrix(U1, V1, W1,
                                                             U2, V2, W2);

            NekMatrix<NekDouble> Vandermonde =
                                LibUtilities::GetTetVandermonde(U1,V1,W1);
            NekMatrix<NekDouble> VandermondeI = Vandermonde;
            VandermondeI.Invert();
            VdmDx = LibUtilities::GetVandermondeForTetXDerivative(U1,V1,W1) *
                                                                VandermondeI;
            VdmDy = LibUtilities::GetVandermondeForTetYDerivative(U1,V1,W1) *
                                                                VandermondeI;
            VdmDz = LibUtilities::GetVandermondeForTetZDerivative(U1,V1,W1) *
                                                                VandermondeI;
            Array<OneD, NekDouble> qds = LibUtilities::PointsManager()[pkey2]->GetW();
            NekVector<NekDouble> quadWi(qds);
            quadW = quadWi;
        }
    }

    GetElementMap();

    vector<vector<NodeSharedPtr> > freenodes = GetColouredNodes();

    vector<vector<NodeOpti> > optiNodes;
    for(int i = 0; i < freenodes.size(); i++)
    {
        vector<NodeOpti> ns;
        for(int j = 0; j < freenodes[i].size(); j++)
        {
            NodeElMap::iterator it = nodeElMap.find(freenodes[i][j]->m_id);
            ASSERTL0(it != nodeElMap.end(),"could not find");
            ns.push_back(NodeOpti(freenodes[i][j],it->second,res));
        }
        optiNodes.push_back(ns);
    }

    res->startE = 0.0;
    for(int i = 0; i < m_mesh->m_element[2].size(); i++)
    {
        res->startE += GetElFunctional(dataSet[i]);
    }

    int nset = optiNodes.size();
    int p = 0;
    int mn = numeric_limits<int>::max();
    int mx = 0;
    for(int i = 0; i < nset; i++)
    {
        p += optiNodes[i].size();
        mn = min(mn, int(optiNodes[i].size()));
        mx = max(mx, int(optiNodes[i].size()));
    }

    cout << scientific << endl;
    cout << "N elements:\t\t" << m_mesh->m_element[m_mesh->m_expDim].size() << endl
         << "N elements invalid:\t" << res->startInv << endl
         << "N free nodes:\t\t" << res->n << endl
         << "N Dof:\t\t\t" << res->nDoF << endl
         << "N color sets:\t\t" << nset << endl
         << "Avg set colors:\t\t" << p/nset << endl
         << "Min set:\t\t" << mn << endl
         << "Max set:\t\t" << mx << endl;
    cout << "Starting energy:\t" << res->startE << endl;

    int nThreads = m_config["numthreads"].as<int>();

    int ctr = 0;
    Thread::ThreadMaster tms;
    tms.SetThreadingType("ThreadManagerBoost");
    Thread::ThreadManagerSharedPtr tm =
                tms.CreateInstance(Thread::ThreadMaster::SessionJob, nThreads);

    while (res->val > 1e-5)
    {
        ctr++;
        res->val = 0.0;
        for(int i = 0; i < optiNodes.size(); i++)
        {
            vector<Thread::ThreadJob*> jobs;
            for(int j = 0; j < optiNodes[i].size(); j++)
            {
                jobs.push_back(optiNodes[i][j].GetJob());
            }
            //cout << " -- inner loop " << i+1 << "/" << optiNodes.size()
            //     << " of size: " << jobs.size() << endl;

            tm->SetNumWorkers(0);
            tm->QueueJobs(jobs);
            tm->SetNumWorkers(nThreads);
            tm->Wait();

            /*for(int j = 0; j < optiNodes[i].size(); j++)
            {
                optiNodes[i][j].Run();
            }*/
        }

        res->val = sqrt(res->val / res->n);

        cout << ctr <<  "\tResidual: " << res->val << endl;
    }

    res->endE = 0.0;
    for(int i = 0; i < m_mesh->m_element[2].size(); i++)
    {
        res->endE += GetElFunctional(dataSet[i]);
    }
    cout << "end energy: " << res->endE << endl;

    if(m_config["stats"].beenSet)
    {
        string file = m_config["stats"].as<string>();
        cout << "writing stats to " << file.c_str() << endl;
        WriteStats(file);
    }
}

void ProcessVarOpti::NodeOpti::Optimise()
{
    //it doesnt matter at this point what the dim is so long that
    //in the 2d case z is left as zero

    Array<OneD, NekDouble> G = GetGrad();

    if(sqrt(G[0]*G[0] + G[1]*G[1] + G[2]*G[2]) > 1e-6)
    {
        //needs to optimise
        NekDouble currentW = GetFunctional();
        NekDouble xc = node->m_x;
        NekDouble yc = node->m_y;
        NekDouble zc = node->m_z;
        NekDouble alpha = 1.0;
        NekDouble delX=0.0;
        NekDouble delY=0.0;
        NekDouble delZ=0.0;
        if(dim == 2)
        {
             delX = 1.0/(G[3]*G[4]-G[6]*G[6])*(G[4]*G[0] - G[6]*G[1]);
             delY = 1.0/(G[3]*G[4]-G[6]*G[6])*(G[3]*G[1] - G[6]*G[0]);
        }
        else
        {
            DNekMat H(3,3,0.0);
            H(0,0) = G[3];
            H(1,1) = G[4];
            H(2,2) = G[5];
            H(0,1) = G[6];
            H(1,0) = G[6];
            H(0,2) = G[8];
            H(2,0) = G[8];
            H(2,1) = G[7];
            H(1,2) = G[7];
            H.Invert();
            NekVector<NekDouble> g(3);
            g[0] = G[0];
            g[1] = G[1];
            g[2] = G[2];
            NekVector<NekDouble> del = H * g;
            delX = del[0];
            delY = del[1];
            delZ = del[2];
        }

        bool found = false;
        while(alpha > 1e-6)
        {
            node->m_x = xc - alpha * delX;
            node->m_y = yc - alpha * delY;
            node->m_z = zc - alpha * delZ;
            if(GetFunctional() < currentW)
            {
                found = true;
                break;
            }

            alpha /= 2.0;
        }

        if(!found)
        {
            //reset the node
            node->m_x = xc;
            node->m_y = yc;
            node->m_z = zc;
        //    cout << "warning: had to reset node" << endl;
        //    cout << G[0] << " " << G[1] << " " << G[2] << " " << node->m_id << endl;
        }
        mtx.lock();
        res->val += (node->m_x-xc)*(node->m_x-xc)+(node->m_y-yc)*(node->m_y-yc)+
                         (node->m_z-zc)*(node->m_z-zc);
        mtx.unlock();
    }
}

NekDouble dir[19][3] = {{0.0,0.0,0.0},
                        {1.0,0.0,0.0},
                        {1.0,1.0,0.0},
                        {0.0,1.0,0.0},
                        {-1.0,1.0,0.0},
                        {-1.0,0.0,0.0},
                        {-1.0,-1.0,0.0},
                        {0.0,-1.0,0.0},
                        {1.0,-1.0,0.0},
                        {-1.0,0.0,-1.0},
                        {0.0,0.0,-1.0},
                        {1.0,0.0,-1.0},
                        {-1.0,0.0,1.0},
                        {0.0,0.0,1.0},
                        {1.0,0.0,1.0},
                        {0.0,1.0,-1.0},
                        {0.0,1.0,1.0},
                        {0.0,-1.0,-1.0},
                        {0.0,-1.0,1.0}};

Array<OneD, NekDouble> ProcessVarOpti::NodeOpti::GetGrad()
{
    NekDouble xc = node->m_x;
    NekDouble yc = node->m_y;
    NekDouble zc = node->m_z;
    NekDouble dx = 1e-3;

    vector<NekDouble> w;

    for(int i = 0; i < (dim == 2 ? 9 : 19); i++)
    {
        node->m_x = xc + dir[i][0] * dx;
        node->m_y = yc + dir[i][1] * dx;
        node->m_z = zc + dir[i][2] * dx;
        w.push_back(GetFunctional());
    }
    node->m_x = xc;
    node->m_y = yc;
    node->m_z = zc;

    Array<OneD, NekDouble> ret(9,0.0);

    //ret[0] d/dx
    //ret[1] d/dy
    //ret[2] d/dz

    //ret[3] d2/dx2
    //ret[4] d2/dy2
    //ret[5] d2/dz2
    //ret[6] d2/dxdy
    //ret[7] d2/dxdz
    //ret[8] d2/dydz


    ret[0] = (w[1] - w[5]) / 2.0 / dx;
    ret[1] = (w[3] - w[7]) / 2.0 / dx;
    ret[3] = (w[1] + w[5] - 2.0*w[0]) / dx / dx;
    ret[4] = (w[3] + w[7] - 2.0*w[0]) / dx / dx;
    ret[6] = (w[2] + w[6] - w[4] - w[8]) / 4.0 / dx /dx;
    if(dim == 3)
    {
        ret[2] = (w[13] - w[10]) / 2.0 / dx;
        ret[5] = (w[13] + w[10] - 2.0*w[0]) / dx / dx;
        ret[7] = (w[14] + w[9] - w[11] - w[12]) / 4.0 / dx /dx;
        ret[8] = (w[16] + w[17] - w[15] - w[18]) / 4.0 / dx /dx;
    }

    return ret;
}

NekDouble ProcessVarOpti::NodeOpti::GetFunctional()
{
    NekDouble r = 0.0;
    for(int i = 0; i < data.size(); i++)
    {
        r += GetElFunctional(data[i]);
    }
    return r;
}

vector<vector<NodeSharedPtr> > ProcessVarOpti::GetColouredNodes()
{
    //this figures out the dirclet nodes and colors the others into paralell sets
    NodeSet boundaryNodes;

    switch (m_mesh->m_spaceDim)
    {
        case 2:
        {
            EdgeSet::iterator it;
            for(it = m_mesh->m_edgeSet.begin(); it != m_mesh->m_edgeSet.end(); it++)
            {
                if((*it)->m_elLink.size() == 2)
                {
                    continue;
                }

                boundaryNodes.insert((*it)->m_n1);
                boundaryNodes.insert((*it)->m_n2);
                for(int i = 0; i < (*it)->m_edgeNodes.size(); i++)
                {
                    boundaryNodes.insert((*it)->m_edgeNodes[i]);
                }
            }
            break;
        }
        case 3:
        {
            FaceSet::iterator it;
            for(it = m_mesh->m_faceSet.begin(); it != m_mesh->m_faceSet.end(); it++)
            {
                if((*it)->m_elLink.size() == 2)
                {
                    continue;
                }

                vector<NodeSharedPtr> vs = (*it)->m_vertexList;
                for(int j = 0; j < vs.size(); j++)
                {
                    boundaryNodes.insert(vs[j]);
                }

                vector<EdgeSharedPtr> es = (*it)->m_edgeList;
                for(int j = 0; j < es.size(); j++)
                {
                    for(int k = 0; k < es[j]->m_edgeNodes.size(); k++)
                    {
                        boundaryNodes.insert(es[j]->m_edgeNodes[k]);
                    }
                }

                for(int i = 0; i < (*it)->m_faceNodes.size(); i++)
                {
                    boundaryNodes.insert((*it)->m_faceNodes[i]);
                }
            }
            break;
        }
        default:
            ASSERTL0(false,"space dim issue");
    }

    vector<NodeSharedPtr> remain;

    if(m_mesh->m_expDim == 2)
    {
        EdgeSet::iterator eit;
        for(eit = m_mesh->m_edgeSet.begin(); eit != m_mesh->m_edgeSet.end(); eit++)
        {
            vector<NodeSharedPtr> n = (*eit)->m_edgeNodes;
            n.push_back((*eit)->m_n1);
            n.push_back((*eit)->m_n2);
            for(int j = 0; j < n.size(); j++)
            {
                NodeSet::iterator nit = boundaryNodes.find(n[j]);
                if(nit == boundaryNodes.end())
                {
                    remain.push_back(n[j]);
                }
            }
        }

        for(int i = 0; i < m_mesh->m_element[2].size(); i++)
        {
            vector<NodeSharedPtr> ns = m_mesh->m_element[2][i]->GetVolumeNodes();
            for(int j = 0; j < ns.size(); j++)
            {
                NodeSet::iterator nit = boundaryNodes.find(ns[j]);
                if(nit == boundaryNodes.end())
                {
                    remain.push_back(ns[j]);
                }
            }
        }
    }
    else
    {
        FaceSet::iterator fit;
        for(fit = m_mesh->m_faceSet.begin(); fit != m_mesh->m_faceSet.end(); fit++)
        {
            vector<NodeSharedPtr> n;
            (*fit)->GetCurvedNodes(n);
            for(int j = 0; j < n.size(); j++)
            {
                NodeSet::iterator nit = boundaryNodes.find(n[j]);
                if(nit == boundaryNodes.end())
                {
                    remain.push_back(n[j]);
                }
            }
        }

        for(int i = 0; i < m_mesh->m_element[3].size(); i++)
        {
            vector<NodeSharedPtr> ns = m_mesh->m_element[3][i]->GetVolumeNodes();
            for(int j = 0; j < ns.size(); j++)
            {
                NodeSet::iterator nit = boundaryNodes.find(ns[j]);
                if(nit == boundaryNodes.end())
                {
                    remain.push_back(ns[j]);
                }
            }
        }
    }

    res->n = remain.size();
    res->nDoF = res->n * dim;

    vector<vector<NodeSharedPtr> > ret;

    while (remain.size() > 0)
    {
        vector<NodeSharedPtr> layer;
        set<int> locked;
        set<int> completed;
        for(int i = 0; i < remain.size(); i++)
        {
            NodeElMap::iterator it = nodeElMap.find(remain[i]->m_id);
            ASSERTL0(it != nodeElMap.end(),"could not find");
            bool islocked = false;
            for(int j = 0; j < it->second.size(); j++)
            {
                set<int>::iterator sit = locked.find(it->second[j]->el->GetId());
                if(sit != locked.end())
                {
                    islocked = true;
                    break;
                }
            }
            if(!islocked)
            {
                layer.push_back(remain[i]);
                completed.insert(remain[i]->m_id);
                for(int j = 0; j < it->second.size(); j++)
                {
                    locked.insert(it->second[j]->el->GetId());
                }
            }
        }

        vector<NodeSharedPtr> tmp = remain;
        remain.clear();
        for(int i = 0; i < tmp.size(); i++)
        {
            set<int>::iterator sit = completed.find(tmp[i]->m_id);
            if(sit == completed.end())
            {
                remain.push_back(tmp[i]);
            }
        }
        ret.push_back(layer);
    }
    return ret;
}

void ProcessVarOpti::GetElementMap()
{
    //build ideal maps and structs;
    for(int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); i++)
    {
        ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][i];
        ElDataSharedPtr d = boost::shared_ptr<ElData>(new ElData);
        d->el = el;

        d->maps = MappingIdealToRef(el);

        dataSet.push_back(d);
    }

    for(int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); i++)
    {
        ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][i];
        vector<NodeSharedPtr> ns;
        el->GetCurvedNodes(ns);

        for(int j = 0; j < ns.size(); j++)
        {
            nodeElMap[ns[j]->m_id].push_back(dataSet[i]);
        }
    }
}

vector<Array<OneD, NekDouble> > ProcessVarOpti::MappingIdealToRef(ElementSharedPtr el)
{
    //need to make ideal element out of old element
    ElmtConfig ec = el->GetConf();
    ec.m_order  = 1;
    ec.m_faceNodes = false;
    ec.m_volumeNodes = false;

    ElementSharedPtr E = GetElementFactory().CreateInstance(
                            ec.m_e, ec, el->GetVertexList(),
                            el->GetTagList());

    SpatialDomains::GeometrySharedPtr    geom = E->GetGeom(el->GetDim());
    geom->FillGeom();
    StdRegions::StdExpansionSharedPtr    chi  = geom->GetXmap();

    vector<Array<OneD, NekDouble> > ret;

    if(geom->GetShapeType() == LibUtilities::eQuadrilateral)
    {
        ASSERTL0(false,"Not coded");
        /*vector<Array<OneD, NekDouble> > xy;
        for(int i = 0; i < geom->GetNumVerts(); i++)
        {
            Array<OneD, NekDouble> loc(2);
            SpatialDomains::PointGeomSharedPtr p = geom->GetVertex(i);
            p->GetCoords(loc);
            xy.push_back(loc);
        }

        Array<OneD, const LibUtilities::BasisSharedPtr> b = chi->GetBase();
        Array<OneD, NekDouble> u = b[0]->GetZ();
        Array<OneD, NekDouble> v = b[1]->GetZ();

        for(int j = 0; j < b[1]->GetNumPoints(); j++)
        {
            for(int i = 0; i < b[0]->GetNumPoints(); i++)
            {
                NekDouble a1 = 0.5*(1.0-u[i]), a2 = 0.5*(1.0+u[i]);
                NekDouble b1 = 0.5*(1.0-v[j]), b2 = 0.5*(1.0+v[j]);
                DNekMat dxdz(2,2,1.0,eFULL);

                dxdz(0,0) = 0.5*(-b1*xy[0][0] + b1*xy[1][0] + b2*xy[2][0] - b2*xy[3][0]);
                dxdz(1,0) = 0.5*(-b1*xy[0][1] + b1*xy[1][1] + b2*xy[2][1] - b2*xy[3][1]);

                dxdz(0,1) = 0.5*(-a1*xy[0][0] - a2*xy[1][0] + a2*xy[2][0] + a1*xy[3][0]);
                dxdz(1,1) = 0.5*(-a1*xy[0][1] - a2*xy[1][1] + a2*xy[2][1] + a1*xy[3][1]);

                NekDouble det = 1.0/(dxdz(0,0)*dxdz(1,1) - dxdz(1,0)*dxdz(0,1));

                dxdz.Invert();
                Array<OneD, NekDouble> r(9,0.0);
                r[0] = dxdz(0,0);
                r[1] = dxdz(1,0);
                r[3] = dxdz(0,1);
                r[4] = dxdz(1,1);
                ret.push_back(r);
            }
        }*/
    }
    else if(geom->GetShapeType() == LibUtilities::eTriangle)
    {
        LibUtilities::PointsKey pkey(m_mesh->m_nummode,
                                     LibUtilities::eNodalTriElec);
        Array<OneD, NekDouble> u, v;
        LibUtilities::PointsManager()[pkey]->GetPoints(u, v);

        Array<OneD, NekDouble> xc(chi->GetTotPoints());
        Array<OneD, NekDouble> yc(chi->GetTotPoints());

        Array<OneD, NekDouble> coeffs0 = geom->GetCoeffs(0);
        Array<OneD, NekDouble> coeffs1 = geom->GetCoeffs(1);

        chi->BwdTrans(coeffs0,xc);
        chi->BwdTrans(coeffs1,yc);

        NekVector<NekDouble> X(ptsLow),Y(ptsLow),
                             x1(ptsLow),y1(ptsLow),
                             x2(ptsLow),y2(ptsLow);
        for(int j = 0; j < u.num_elements(); j++)
        {
            Array<OneD, NekDouble> xp(2);
            xp[0] = u[j];
            xp[1] = v[j];

            X(j) = chi->PhysEvaluate(xp, xc);
            Y(j) = chi->PhysEvaluate(xp, yc);
        }

        NekVector<NekDouble> x1i(ptsHigh),y1i(ptsHigh),
                             x2i(ptsHigh),y2i(ptsHigh);

        x1 = VdmDx*X;
        y1 = VdmDx*Y;
        x2 = VdmDy*X;
        y2 = VdmDy*Y;

        x1i = interp * x1;
        x2i = interp * x2;
        y1i = interp * y1;
        y2i = interp * y2;


        for(int i = 0 ; i < ptsHigh; i++)
        {
            DNekMat dxdz(2,2,1.0,eFULL);
            dxdz(0,0) = x1i(i);
            dxdz(0,1) = x2i(i);
            dxdz(1,0) = y1i(i);
            dxdz(1,1) = y2i(i);

            Array<OneD, NekDouble> r(10,0.0);
            r[9] = dxdz(0,0)*dxdz(1,1)-dxdz(1,0)*dxdz(0,1);

            dxdz.Invert();

            r[0] = dxdz(0,0);
            r[1] = dxdz(1,0);
            r[3] = dxdz(0,1);
            r[4] = dxdz(1,1);
            ret.push_back(r);
        }
    }
    else if(geom->GetShapeType() == LibUtilities::eTetrahedron)
    {
        LibUtilities::PointsKey pkey(m_mesh->m_nummode,
                                     LibUtilities::eNodalTetElec);
        Array<OneD, NekDouble> u, v, w;
        LibUtilities::PointsManager()[pkey]->GetPoints(u, v, w);

        Array<OneD, NekDouble> xc(chi->GetTotPoints());
        Array<OneD, NekDouble> yc(chi->GetTotPoints());
        Array<OneD, NekDouble> zc(chi->GetTotPoints());

        Array<OneD, NekDouble> coeffs0 = geom->GetCoeffs(0);
        Array<OneD, NekDouble> coeffs1 = geom->GetCoeffs(1);
        Array<OneD, NekDouble> coeffs2 = geom->GetCoeffs(2);

        chi->BwdTrans(coeffs0,xc);
        chi->BwdTrans(coeffs1,yc);
        chi->BwdTrans(coeffs2,zc);

        NekVector<NekDouble> X(ptsLow),Y(ptsLow),Z(ptsLow),
                             x1(ptsLow),y1(ptsLow),z1(ptsLow),
                             x2(ptsLow),y2(ptsLow),z2(ptsLow),
                             x3(ptsLow),y3(ptsLow),z3(ptsLow);
        for(int j = 0; j < u.num_elements(); j++)
        {
            Array<OneD, NekDouble> xp(3);
            xp[0] = u[j];
            xp[1] = v[j];
            xp[2] = w[j];

            X(j) = chi->PhysEvaluate(xp, xc);
            Y(j) = chi->PhysEvaluate(xp, yc);
            Z(j) = chi->PhysEvaluate(xp, zc);
        }

        NekVector<NekDouble> x1i(ptsHigh),y1i(ptsHigh),z1i(ptsHigh),
                             x2i(ptsHigh),y2i(ptsHigh),z2i(ptsHigh),
                             x3i(ptsHigh),y3i(ptsHigh),z3i(ptsHigh);

        x1 = VdmDx*X;
        y1 = VdmDx*Y;
        x2 = VdmDy*X;
        y2 = VdmDy*Y;
        z1 = VdmDx*Z;
        z2 = VdmDy*Z;
        x3 = VdmDz*X;
        y3 = VdmDz*Y;
        z3 = VdmDz*Z;

        x1i = interp * x1;
        x2i = interp * x2;
        x3i = interp * x3;
        y1i = interp * y1;
        y2i = interp * y2;
        y3i = interp * y3;
        z1i = interp * z1;
        z2i = interp * z2;
        z3i = interp * z3;

        for(int i = 0 ; i < ptsHigh; i++)
        {
            DNekMat dxdz(3,3,1.0,eFULL);
            dxdz(0,0) = x1i(i);
            dxdz(0,1) = x2i(i);
            dxdz(0,2) = x3i(i);
            dxdz(1,0) = y1i(i);
            dxdz(1,1) = y2i(i);
            dxdz(1,2) = y3i(i);
            dxdz(2,0) = z1i(i);
            dxdz(2,1) = z2i(i);
            dxdz(2,2) = z3i(i);

            Array<OneD, NekDouble> r(10,0.0); //store det in 10th entry

            r[9] = dxdz(0,0)*(dxdz(1,1)*dxdz(2,2)-dxdz(2,1)*dxdz(1,2))
                   -dxdz(0,1)*(dxdz(1,0)*dxdz(2,2)-dxdz(2,0)*dxdz(1,2))
                   +dxdz(0,2)*(dxdz(1,0)*dxdz(2,1)-dxdz(2,0)*dxdz(1,1));

            dxdz.Invert();

            r[0] = dxdz(0,0);
            r[1] = dxdz(1,0);
            r[2] = dxdz(2,0);
            r[3] = dxdz(0,1);
            r[4] = dxdz(1,1);
            r[5] = dxdz(2,1);
            r[6] = dxdz(0,2);
            r[7] = dxdz(1,2);
            r[8] = dxdz(2,2);
            ret.push_back(r);
        }
    }
    else if(geom->GetShapeType() == LibUtilities::ePrism)
    {
        ASSERTL0(false, "not coded");
        /*vector<Array<OneD, NekDouble> > xyz;
        for(int i = 0; i < geom->GetNumVerts(); i++)
        {
            Array<OneD, NekDouble> loc(3);
            SpatialDomains::PointGeomSharedPtr p = geom->GetVertex(i);
            p->GetCoords(loc);
            xyz.push_back(loc);
        }

        Array<OneD, const LibUtilities::BasisSharedPtr> b = chi->GetBase();
        Array<OneD, NekDouble> eta1 = b[0]->GetZ();
        Array<OneD, NekDouble> eta2 = b[1]->GetZ();
        Array<OneD, NekDouble> eta3 = b[2]->GetZ();

        for(int k = 0; k < b[2]->GetNumPoints(); k++)
        {

            for(int j = 0; j < b[1]->GetNumPoints(); j++)
            {
                for(int i = 0; i < b[0]->GetNumPoints(); i++)
                {
                    NekDouble xi1 = 0.5*(1+eta1[i])*(1-eta3[k])-1.0;
                    NekDouble a1 = 0.5*(1-xi1),     a2 = 0.5*(1+xi1);
                    NekDouble b1 = 0.5*(1-eta2[j]), b2 = 0.5*(1+eta2[j]);
                    NekDouble c1 = 0.5*(1-eta3[k]), c2 = 0.5*(1+eta3[k]);

                    DNekMat dxdz(3,3,1.0,eFULL);

                    dxdz(0,0) = 0.5*(-b1*xyz[0][0] + b1*xyz[1][0] + b2*xyz[2][0] - b2*xyz[3][0]);
                    dxdz(1,0) = 0.5*(-b1*xyz[0][1] + b1*xyz[1][1] + b2*xyz[2][1] - b2*xyz[3][1]);
                    dxdz(2,0) = 0.5*(-b1*xyz[0][2] + b1*xyz[1][2] + b2*xyz[2][2] - b2*xyz[3][2]);

                    dxdz(0,1) = 0.5*((a2-c1)*xyz[0][0] - a2*xyz[1][0] + a2*xyz[2][0] + (c1-a2)*xyz[3][0] - c2*xyz[4][0] + c2*xyz[5][0]);
                    dxdz(1,1) = 0.5*((a2-c1)*xyz[0][1] - a2*xyz[1][1] + a2*xyz[2][1] + (c1-a2)*xyz[3][1] - c2*xyz[4][1] + c2*xyz[5][1]);
                    dxdz(2,1) = 0.5*((a2-c1)*xyz[0][2] - a2*xyz[1][2] + a2*xyz[2][2] + (c1-a2)*xyz[3][2] - c2*xyz[4][2] + c2*xyz[5][2]);

                    dxdz(0,2) = 0.5*(-b1*xyz[0][0] - b2*xyz[3][0] + b1*xyz[4][0] + b2*xyz[5][0]);
                    dxdz(1,2) = 0.5*(-b1*xyz[0][1] - b2*xyz[3][1] + b1*xyz[4][1] + b2*xyz[5][1]);
                    dxdz(2,2) = 0.5*(-b1*xyz[0][2] - b2*xyz[3][2] + b1*xyz[4][2] + b2*xyz[5][2]);

                    dxdz.Invert();
                    Array<OneD, NekDouble> r(9,0.0);
                    r[0] = dxdz(0,0);
                    r[1] = dxdz(1,0);
                    r[3] = dxdz(0,1);
                    r[4] = dxdz(1,1);
                    r[2] = dxdz(2,0);
                    r[5] = dxdz(2,1);
                    r[6] = dxdz(0,2);
                    r[7] = dxdz(1,2);
                    r[8] = dxdz(2,2);
                    ret.push_back(r);
                }
            }
        }*/
    }
    else
    {
        ASSERTL0(false,"not coded");
    }

    return ret;
}

void ProcessVarOpti::FillQuadPoints()
{
    int nq = m_mesh->m_nummode;
    int id = m_mesh->m_vertexSet.size();

    LibUtilities::PointsKey ekey(m_mesh->m_nummode,
                                 LibUtilities::eGaussLobattoLegendre);
    Array<OneD, NekDouble> gll;
    LibUtilities::PointsManager()[ekey]->GetPoints(gll);

    EdgeSet::iterator eit;

    for(eit = m_mesh->m_edgeSet.begin(); eit != m_mesh->m_edgeSet.end(); eit++)
    {
        SpatialDomains::Geometry1DSharedPtr geom =
                                            (*eit)->GetGeom(m_mesh->m_spaceDim);
        geom->FillGeom();
        StdRegions::StdExpansionSharedPtr xmap = geom->GetXmap();

        vector<NodeSharedPtr> hons;

        switch (m_mesh->m_spaceDim)
        {
            case 2:
            {
                Array<OneD, NekDouble> coeffs0 = geom->GetCoeffs(0);
                Array<OneD, NekDouble> coeffs1 = geom->GetCoeffs(1);

                Array<OneD, NekDouble> xc(xmap->GetTotPoints());
                Array<OneD, NekDouble> yc(xmap->GetTotPoints());

                xmap->BwdTrans(coeffs0,xc);
                xmap->BwdTrans(coeffs1,yc);

                for(int j = 1; j < m_mesh->m_nummode - 1; j++)
                {
                    Array<OneD, NekDouble> xp(2);
                    xp[0] = gll[j];

                    hons.push_back(boost::shared_ptr<Node>(new Node(
                            id++,xmap->PhysEvaluate(xp,xc),
                                 xmap->PhysEvaluate(xp,yc),0.0)));
                }
                (*eit)->m_edgeNodes.clear();
                (*eit)->m_edgeNodes = hons;
                (*eit)->m_curveType = LibUtilities::eGaussLobattoLegendre;
            }
            break;
            case 3:
            {
                Array<OneD, NekDouble> coeffs0 = geom->GetCoeffs(0);
                Array<OneD, NekDouble> coeffs1 = geom->GetCoeffs(1);
                Array<OneD, NekDouble> coeffs2 = geom->GetCoeffs(2);

                Array<OneD, NekDouble> xc(xmap->GetTotPoints());
                Array<OneD, NekDouble> yc(xmap->GetTotPoints());
                Array<OneD, NekDouble> zc(xmap->GetTotPoints());

                xmap->BwdTrans(coeffs0,xc);
                xmap->BwdTrans(coeffs1,yc);
                xmap->BwdTrans(coeffs2,zc);

                for(int j = 1; j < m_mesh->m_nummode - 1; j++)
                {
                    Array<OneD, NekDouble> xp(2);
                    xp[0] = gll[j];

                    hons.push_back(boost::shared_ptr<Node>(new Node(
                            id++,xmap->PhysEvaluate(xp,xc),
                                 xmap->PhysEvaluate(xp,yc),
                                 xmap->PhysEvaluate(xp,zc))));
                }
                (*eit)->m_edgeNodes.clear();
                (*eit)->m_edgeNodes = hons;
                (*eit)->m_curveType = LibUtilities::eGaussLobattoLegendre;
            }
            break;
        }
    }

    if(m_mesh->m_expDim == 2)
    {
        //for faces need to do volume nodes
        for(int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); i++)
        {
            ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][i];

            SpatialDomains::GeometrySharedPtr geom =
                                            el->GetGeom(m_mesh->m_spaceDim);
            geom->FillGeom();
            StdRegions::StdExpansionSharedPtr xmap = geom->GetXmap();

            LibUtilities::PointsKey pkey(m_mesh->m_nummode,
                                         LibUtilities::eNodalTriElec);
            Array<OneD, NekDouble> u, v;
            LibUtilities::PointsManager()[pkey]->GetPoints(u, v);

            Array<OneD, NekDouble> coeffs0 = geom->GetCoeffs(0);
            Array<OneD, NekDouble> coeffs1 = geom->GetCoeffs(1);

            Array<OneD, NekDouble> xc(xmap->GetTotPoints());
            Array<OneD, NekDouble> yc(xmap->GetTotPoints());

            xmap->BwdTrans(coeffs0,xc);
            xmap->BwdTrans(coeffs1,yc);

            vector<NodeSharedPtr> hons;

            for(int j = nq * (nq + 1) / 2 - (nq-2)*(nq-3) / 2;
                                                    j < u.num_elements(); j++)
            {
                Array<OneD, NekDouble> xp(2);
                xp[0] = u[j];
                xp[1] = v[j];

                hons.push_back(boost::shared_ptr<Node>(new Node(
                        id++,xmap->PhysEvaluate(xp,xc),
                             xmap->PhysEvaluate(xp,yc),0.0)));
            }

            el->SetVolumeNodes(hons);
            el->SetCurveType(LibUtilities::eNodalTriElec);
        }
    }
    else
    {
        FaceSet::iterator it;
        for(it = m_mesh->m_faceSet.begin(); it != m_mesh->m_faceSet.end(); it++)
        {
            //this is a hack and needs to be fixed
            //it really should take the get geom of the whole element and
            //then pick the correct parts
            /////////
            (*it)->m_faceNodes.clear();
            ////////
            SpatialDomains::Geometry2DSharedPtr geom =
                                            (*it)->GetGeom(m_mesh->m_spaceDim);
            geom->FillGeom();
            StdRegions::StdExpansionSharedPtr xmap = geom->GetXmap();

            LibUtilities::PointsKey pkey(m_mesh->m_nummode,
                                         LibUtilities::eNodalTriElec);
            Array<OneD, NekDouble> u, v;
            LibUtilities::PointsManager()[pkey]->GetPoints(u, v);

            vector<NodeSharedPtr> hons;

            Array<OneD, NekDouble> coeffs0 = geom->GetCoeffs(0);
            Array<OneD, NekDouble> coeffs1 = geom->GetCoeffs(1);
            Array<OneD, NekDouble> coeffs2 = geom->GetCoeffs(2);

            Array<OneD, NekDouble> xc(xmap->GetTotPoints());
            Array<OneD, NekDouble> yc(xmap->GetTotPoints());
            Array<OneD, NekDouble> zc(xmap->GetTotPoints());

            xmap->BwdTrans(coeffs0,xc);
            xmap->BwdTrans(coeffs1,yc);
            xmap->BwdTrans(coeffs2,zc);

            for(int j = nq * (nq + 1) / 2 - (nq-2)*(nq-3) / 2;
                                                    j < u.num_elements(); j++)
            {
                Array<OneD, NekDouble> xp(2);
                xp[0] = u[j];
                xp[1] = v[j];

                hons.push_back(boost::shared_ptr<Node>(new Node(
                        id++,xmap->PhysEvaluate(xp,xc),
                             xmap->PhysEvaluate(xp,yc),
                             xmap->PhysEvaluate(xp,zc))));
            }
            (*it)->m_faceNodes.clear();
            (*it)->m_faceNodes = hons;
            (*it)->m_curveType = LibUtilities::eNodalTriElec;
        }
        for(int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); i++)
        {
            ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][i];

            SpatialDomains::GeometrySharedPtr geom =
                                            el->GetGeom(m_mesh->m_spaceDim);
            geom->FillGeom();
            StdRegions::StdExpansionSharedPtr xmap = geom->GetXmap();

            LibUtilities::PointsKey pkey(m_mesh->m_nummode,
                                         LibUtilities::eNodalTetElec);
            Array<OneD, NekDouble> u, v, w;
            LibUtilities::PointsManager()[pkey]->GetPoints(u, v, w);

            Array<OneD, NekDouble> coeffs0 = geom->GetCoeffs(0);
            Array<OneD, NekDouble> coeffs1 = geom->GetCoeffs(1);
            Array<OneD, NekDouble> coeffs2 = geom->GetCoeffs(2);

            Array<OneD, NekDouble> xc(xmap->GetTotPoints());
            Array<OneD, NekDouble> yc(xmap->GetTotPoints());
            Array<OneD, NekDouble> zc(xmap->GetTotPoints());

            xmap->BwdTrans(coeffs0,xc);
            xmap->BwdTrans(coeffs1,yc);
            xmap->BwdTrans(coeffs2,zc);

            vector<NodeSharedPtr> hons;

            ASSERTL0(4 + 6*(nq-2) + 4 * ((nq-2)*(nq-3) / 2) <= u.num_elements(),
                        "volume interior nodes in tet");

            /*//need to finish for tet
            for(int j = 4 + 6*(nq-2) + 4 * ((nq-2)*(nq-3) / 2);
                                                    j < u.num_elements(); j++)
            {
                Array<OneD, NekDouble> xp(3);
                xp[0] = u[j];
                xp[1] = v[j];
                xp[2] = w[j];

                hons.push_back(boost::shared_ptr<Node>(new Node(
                        id++,xmap->PhysEvaluate(xp,xc),
                             xmap->PhysEvaluate(xp,yc),
                             xmap->PhysEvaluate(xp,zc))));
            }

            el->SetVolumeNodes(hons);*/
            el->SetCurveType(LibUtilities::eNodalTetElec);
        }
    }

    res->startInv =0;
    res->worstJac = numeric_limits<double>::max();

    for(int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); i++)
    {
        ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][i];

        SpatialDomains::GeometrySharedPtr geom =
                                        el->GetGeom(m_mesh->m_spaceDim);
        SpatialDomains::GeomFactorsSharedPtr gfac = geom->GetGeomFactors();
        if(!gfac->IsValid())
        {
            res->startInv++;
        }
    }


    for(int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); i++)
    {
        ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][i];

        SpatialDomains::GeometrySharedPtr    geom = el->GetGeom(el->GetDim());
        geom->FillGeom();
        StdRegions::StdExpansionSharedPtr    chi  = geom->GetXmap();

        LibUtilities::PointsKey pkey(m_mesh->m_nummode,
                                     LibUtilities::eNodalTetElec);
        Array<OneD, NekDouble> u, v, w;
        LibUtilities::PointsManager()[pkey]->GetPoints(u, v, w);

        Array<OneD, NekDouble> xc(chi->GetTotPoints());
        Array<OneD, NekDouble> yc(chi->GetTotPoints());
        Array<OneD, NekDouble> zc(chi->GetTotPoints());

        Array<OneD, NekDouble> coeffs0 = geom->GetCoeffs(0);
        Array<OneD, NekDouble> coeffs1 = geom->GetCoeffs(1);
        Array<OneD, NekDouble> coeffs2 = geom->GetCoeffs(2);

        chi->BwdTrans(coeffs0,xc);
        chi->BwdTrans(coeffs1,yc);
        chi->BwdTrans(coeffs2,zc);

        for(int j = 0; j < u.num_elements(); j++)
        {
            Array<OneD, NekDouble> xp(3);
            xp[0] = u[j];
            xp[1] = v[j];
            xp[2] = w[j];
        }
    }
}

void ProcessVarOpti::WriteStats(string file)
{
    ASSERTL0(file != "", "no file name given");

    ofstream out;
    out.open(file.c_str());
    out << scientific;

    for(int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); i++)
    {
        vector<NodeSharedPtr> ns;
        m_mesh->m_element[m_mesh->m_expDim][i]->GetCurvedNodes(ns);
        int pts = ns.size();
        int dim = m_mesh->m_element[m_mesh->m_expDim][i]->GetDim();

        NekDouble mx = -1.0 * numeric_limits<double>::max();
        NekDouble mn = numeric_limits<double>::max();

        if(dim == 2)
        {
            NekVector<NekDouble> X(pts),Y(pts),Z(pts),
                                 x1(pts),y1(pts),
                                 x2(pts),y2(pts);
            for(int i = 0; i < pts; i++)
            {
                X(i) = ns[i]->m_x;
                Y(i) = ns[i]->m_y;
            }

            x1 = VdmDx*X;
            y1 = VdmDx*Y;
            x2 = VdmDy*X;
            y2 = VdmDy*Y;

            for(int i = 0; i < pts; i++)
            {
                mx = max(mx, x1(i)*y2(i)-y1(i)*x2(i));
                mn = min(mn, x1(i)*y2(i)-y1(i)*x2(i));
            }

        }
        else
        {
            ASSERTL0(false,"not coded");
        }
        out << mn / mx << endl;
    }
}

}
}
