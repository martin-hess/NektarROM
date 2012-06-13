///////////////////////////////////////////////////////////////////////////////
//
// File: UpwindSolver.cpp
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
// Description: Upwind Riemann solver.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/RiemannSolvers/UpwindSolver.h>

namespace Nektar
{
    namespace SolverUtils
    {
        std::string UpwindSolver::solverName = GetRiemannSolverFactory().
            RegisterCreatorFunction("Upwind", UpwindSolver::create, "Upwind solver");
        
        UpwindSolver::UpwindSolver() :
            RiemannSolver()
        {
        }
        
        void UpwindSolver::v_Solve(
            const Array<OneD, const Array<OneD, NekDouble> > &Fwd,
            const Array<OneD, const Array<OneD, NekDouble> > &Bwd,
                  Array<OneD,       Array<OneD, NekDouble> > &flux)
        {
            ASSERTL1(CheckScalars("Vn"), "Vn not defined.");
            
            Array<OneD, const NekDouble> &traceVel = m_scalars["Vn"]();
            
            for (int i = 0; i < Fwd.num_elements(); ++i)
            {
                for (int j = 0; j < traceVel.num_elements(); ++j)
                {
                    flux[i][j] = traceVel[j]*(traceVel[j] > 0.0 ? Fwd[i][j] : Bwd[i][j]);
                }
            }
        }
    }
}