////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessDeform.cpp
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
//  Description: Computes Q Criterion field.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
using namespace std;

#include "ProcessDeform.h"

#include <StdRegions/StdSegExp.h>
#include <StdRegions/StdQuadExp.h>
#include <StdRegions/StdTriExp.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <SolverUtils/Core/Deform.h>

namespace Nektar
{
    namespace Utilities
    {
        ModuleKey ProcessDeform::className =
        GetModuleFactory().RegisterCreatorFunction(
            ModuleKey(eProcessModule, "deform"), ProcessDeform::create,
            "Deform a mesh given an input field defining displacement");

        ProcessDeform::ProcessDeform(FieldSharedPtr f) :
            ProcessModule(f)
        {
        }

        ProcessDeform::~ProcessDeform()
        {
        }

        void ProcessDeform::Process(po::variables_map &vm)
        {
            if (m_f->m_verbose)
            {
                cout << "ProcessDeform: Deforming grid..." << endl;
            }

            SolverUtils::UpdateGeometry(m_f->m_graph, m_f->m_exp);
        }
    }
}
