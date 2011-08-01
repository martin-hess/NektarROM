///////////////////////////////////////////////////////////////////////////////
//
// File DriverArnoldi.hpp
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
// Description: Base Driver class for the stability solver
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_DRIVERARNOLDI_H
#define NEKTAR_SOLVERS_DRIVERARNOLDI_H

#include <Auxiliary/Driver.h>

namespace Nektar
{
    /// Base class for the development of solvers.
    class DriverArnoldi: public Driver
    {
    public:
        friend class MemoryManager<DriverArnoldi>;
		


    protected:
        int       m_kdim; /// Dimension of Krylov subspace
        int       m_nvec; /// Number of vectors to test 
        int       m_nits; /// Maxmum number of iterations
        NekDouble m_evtol;/// Tolerance of iteratiosn
        NekDouble m_period;/// Period of time stepping algorithm 
        bool      m_TimeSteppingAlgorithm; /// underlying operator is time stepping
        
        int m_nfields; 

        /// Constructor
        DriverArnoldi(LibUtilities::SessionReaderSharedPtr pSession)
        : Driver(pSession)
        {
        };
        
        //Destructor
        virtual ~DriverArnoldi()
        {
        };        

        /**
         * Copy Arnoldi array to field variables which depend from
         * either the m_fields or m_forces 
         */
        void v_CopyArnoldiArrayToField(Array<OneD, NekDouble> &array)
        {
            Array<OneD, MultiRegions::ExpListSharedPtr> fields;
            if(m_TimeSteppingAlgorithm)
            {
                fields = m_equ[0]->UpdateFields();
            }
            else
            {
                fields = m_equ[0]->UpdateForces();
            }            

            int nq = fields[0]->GetNpoints();
            
            for (int k = 0; k < m_nfields; ++k)
            {
                Vmath::Vcopy(nq, &array[k*nq], 1, &fields[k]->UpdatePhys()[0], 1);
                fields[k]->SetPhysState(true);
            }
        };


        /**
         * Copy field variables which depend from either the m_fields
         * or m_forces array the Arnoldi array
         */
        void v_CopyFieldToArnoldiArray(Array<OneD, NekDouble> &array)
        {
            Array<OneD, MultiRegions::ExpListSharedPtr> fields;
            fields = m_equ[0]->UpdateFields();
            int nq = fields[0]->GetNpoints();
            
            for (int k = 0; k < m_nfields; ++k)
            {
                if(!m_TimeSteppingAlgorithm)
                {
                    fields[k]->BwdTrans_IterPerExp(fields[k]->GetCoeffs(),fields[k]->UpdatePhys());
                }

                Vmath::Vcopy(nq,  &fields[k]->GetPhys()[0], 1, &array[k*nq], 1);
                fields[k]->SetPhysState(true);
            }
            
        };

    };
	
} //end of namespace

#endif //NEKTAR_SOLVERS_DRIVER_ARNOLDI_H

