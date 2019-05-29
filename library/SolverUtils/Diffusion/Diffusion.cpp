///////////////////////////////////////////////////////////////////////////////
//
// File: Diffusion.cpp
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
// Description: Abstract base class for diffusion.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Diffusion/Diffusion.h>

namespace Nektar
{
    namespace SolverUtils
    {
        DiffusionFactory& GetDiffusionFactory()
        {
            static DiffusionFactory instance;
            return instance;
        }

        void Diffusion::InitObject(
            const LibUtilities::SessionReaderSharedPtr        pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr>       pFields)
        {
            v_InitObject(pSession, pFields);

            // Div curl storage
            int nPts = pFields[0]->GetTotPoints();
            m_divVel        = Array<OneD, NekDouble>(nPts, 0.0);
            m_divVelSquare  = Array<OneD, NekDouble>(nPts, 0.0);
            m_curlVelSquare = Array<OneD, NekDouble>(nPts, 0.0);
        }

        void Diffusion::Diffuse(
            const int nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray,
            const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
        {
            v_Diffuse(nConvectiveFields, fields, inarray, outarray, pFwd, pBwd);
        }

        void Diffusion::Diffuse_coeff(
            const int nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray,
            const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
        {
            v_Diffuse_coeff(nConvectiveFields, fields, inarray, outarray, pFwd, pBwd);
        }

        void Diffusion::Diffuse(
            const int nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
            Array<OneD, Array<OneD, NekDouble> >              &outarray,
            NekDouble                                         time,
            const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
        {
            m_time  =    time;
            v_Diffuse(nConvectiveFields, fields, inarray, outarray, pFwd, pBwd);
        }

        void Diffusion::Diffuse_coeff(
            const int nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
            Array<OneD, Array<OneD, NekDouble> >              &outarray,
            NekDouble                                         time,
            const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
        {
            m_time  =    time;
            v_Diffuse_coeff(nConvectiveFields, fields, inarray, outarray, pFwd, pBwd);
        }
        void Diffusion::GetAVmu(
            const Array<OneD, MultiRegions::ExpListSharedPtr>           &fields,
            const Array<OneD, Array<OneD, NekDouble> >                  &inarray,
                  Array<OneD, NekDouble >                               &muvar,
                  Array<OneD, NekDouble >                               &MuVarTrace)
        {
            int nTracePts = fields[0]->GetTrace()->GetTotPoints();

            Array<OneD, NekDouble> Fwd(nTracePts,0.0);
            Array<OneD, NekDouble> Bwd(nTracePts,0.0);

            m_ArtificialDiffusionVector(inarray, muvar);

            // BwdMuvar is left to be 0.0 according to DiffusionLDG.cpp
            fields[0]->GetFwdBwdTracePhysNoBndFill(muvar,Fwd,Bwd);

            for(int k = 0; k < nTracePts; ++k)
            {
                MuVarTrace[k] = 0.5 * (Fwd[k] + Bwd[k]) ;
            }
        }

        void Diffusion::v_Diffuse_coeff(
            const int nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray,
            const Array<OneD, Array<OneD, NekDouble> > &pFwd,
            const Array<OneD, Array<OneD, NekDouble> > &pBwd)
        {
            ASSERTL0(false,"v_Diffuse_coeff not defined");
        }

        // Diffusion Calculate the physical derivatives
        void Diffusion::DiffuseCalculateDerivative(
                const int                                         nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble>>         &inarray,
                Array<OneD,Array<OneD, Array<OneD, NekDouble> > > &inarrayderivative,
                const Array<OneD, Array<OneD, NekDouble>>         &pFwd =NullNekDoubleArrayofArray,
                const Array<OneD, Array<OneD, NekDouble>>         &pBwd =NullNekDoubleArrayofArray)
            {
                v_DiffuseCalculateDerivative(nConvectiveFields, fields, inarray,inarrayderivative,pFwd, pBwd);

                GetDivCurl(fields, inarrayderivative);
            }

        const Array<OneD, const Array<OneD, NekDouble> > &Diffusion::v_GetTraceNormal()
        {
            ASSERTL0(false," not defined");
        }

        void Diffusion::v_ConsVarAveJump(
                const int                                           nConvectiveFields,
                const int                                           npnts,
                const Array<OneD, const Array<OneD, NekDouble> >    &vFwd,
                const Array<OneD, const Array<OneD, NekDouble> >    &vBwd,
                      Array<OneD,       Array<OneD, NekDouble> >    &aver,
                      Array<OneD,       Array<OneD, NekDouble> >    &jump)
        {
            ASSERTL0(false," not defined");
        }

        // No multiply(check if diffsionVolume difined)
        void Diffusion::v_DiffuseCalculateDerivative(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble>>         &inarray,
            Array<OneD,Array<OneD, Array<OneD, NekDouble> > > &inarrayderivative,
            const Array<OneD, Array<OneD, NekDouble>>         &pFwd,
            const Array<OneD, Array<OneD, NekDouble>>         &pBwd)
        {
            ASSERTL0(false, "Not defined for function DiffuseVolumeFLux.");
        }

        // No multiply(check if diffsionVolume difined)
        void Diffusion::v_DiffuseVolumeFlux(
            const int                                           nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
            const Array<OneD, Array<OneD, NekDouble>>           &inarray,
            Array<OneD,Array<OneD, Array<OneD, NekDouble> > >   &inarrayderivative,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &VolumeFlux,
            Array< OneD, int >                                  &nonZeroIndex)
        {
            ASSERTL0(false, "Not defined for function DiffuseVolumeFLux.");
        }

        // No multiply(check if diffsionTraceFlux difined)
        void Diffusion::v_DiffuseTraceFlux(
            const int                                           nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
            const Array<OneD, Array<OneD, NekDouble>>           &inarray,
            Array<OneD,Array<OneD, Array<OneD, NekDouble> > >   &inarrayderivative,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &VolumeFlux,
            Array<OneD, Array<OneD, NekDouble> >                &TraceFlux,
            const Array<OneD, Array<OneD, NekDouble>>           &pFwd,
            const Array<OneD, Array<OneD, NekDouble>>           &pBwd,
            Array< OneD, int >                                  &nonZeroIndex)
        {
            ASSERTL0(false, "Not defined function DiffuseTraceFLux.");
        }

        void Diffusion::v_DiffuseTraceFlux(
                const int                                                       nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>               &fields,
                const Array<OneD, Array<OneD, NekDouble>>                       &inarray,
                const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >  &qfield,
                Array<OneD, Array<OneD, NekDouble> >                            &TraceFlux,
                const Array<OneD, Array<OneD, NekDouble>>                       &pFwd,
                const Array<OneD, Array<OneD, NekDouble>>                       &pBwd,
                const Array<OneD, NekDouble>                                    &MuAVTrace,
                Array< OneD, int >                                              &nonZeroIndex  ,
                const Array<OneD, Array<OneD, NekDouble>>                       &Aver          ,
                const Array<OneD, Array<OneD, NekDouble>>                       &Jump          )
            {
                ASSERTL0(false, "Not defined function DiffuseTraceFLux.");
            }

#ifdef DEMO_IMPLICITSOLVER_JFNK_COEFF
        void Diffusion::v_MinusVolumDerivJacToMat(
            const int                                               nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr>       &pFields,
            const Array<OneD, const Array<OneD, DNekMatSharedPtr> > &ElmtJac,
            const int                                               nfluxDir,
            const int                                               nDervDir,
                  Array<OneD, Array<OneD, DNekBlkMatSharedPtr> >    &gmtxarray)
        {
            ASSERTL0(false," not defined");
        }
#endif

        /**
         * @brief Compute and store scalars needed for shock sensor
         *
         */
        void Diffusion::GetDivCurl(
            const Array<OneD, MultiRegions::ExpListSharedPtr>               &fields,
            const Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &derivativesO1)
        {
            int nDim      = fields[0]->GetCoordim(0);
            int nPts      = fields[0]->GetTotPoints();

            // Compute and store scalars needed for shock sensor
            Array<OneD, NekDouble>              tmp3(nPts, 0.0);
            Array<OneD, NekDouble>              tmp4(nPts, 0.0);

            // div vel
            for (int j = 0; j < nDim; ++j)
            {
                Vmath::Vadd(nPts, m_divVel, 1, derivativesO1[j][j], 1, m_divVel, 1);
            }
            // (div vel)**2
            Vmath::Vmul(nPts, m_divVel, 1, m_divVel, 1, m_divVelSquare, 1);

            // curl vel
            if ( nDim > 2 )
            {
                // curl[0] 3/2
                Vmath::Vadd(nPts, tmp4, 1, derivativesO1[2][1], 1, tmp4, 1);
                // curl[0]-2/3
                Vmath::Vcopy(nPts, derivativesO1[1][2], 1, tmp3, 1);
                Vmath::Neg(nPts, tmp3, 1);
                Vmath::Vadd(nPts, tmp4, 1, tmp3, 1, tmp4, 1);
                // square curl[0]
                Vmath::Vmul(nPts, tmp4, 1, tmp4, 1, tmp3, 1);
                Vmath::Vadd(nPts, m_curlVelSquare, 1, tmp3, 1, m_curlVelSquare, 1);

                tmp4 = Array<OneD, NekDouble>(nPts, 0.0);
                // curl[1] 3/1
                Vmath::Vadd(nPts, tmp4, 1, derivativesO1[2][0], 1, tmp4, 1);
                // curl[1]-1/3
                Vmath::Vcopy(nPts, derivativesO1[0][2], 1, tmp3, 1);
                Vmath::Neg(nPts, tmp3, 1);
                Vmath::Vadd(nPts, tmp4, 1, tmp3, 1, tmp4, 1);
                // square curl[1]
                Vmath::Vmul(nPts, tmp4, 1, tmp4, 1, tmp3, 1);
                Vmath::Vadd(nPts, m_curlVelSquare, 1, tmp3, 1, m_curlVelSquare, 1);

                tmp4 = Array<OneD, NekDouble>(nPts, 0.0);
                // curl[2] 1/2
                Vmath::Vadd(nPts, tmp4, 1, derivativesO1[0][1], 1, tmp4, 1);
                // curl[2]-2/1
                Vmath::Vcopy(nPts, derivativesO1[1][0], 1, tmp3, 1);
                Vmath::Neg(nPts, tmp3, 1);
                Vmath::Vadd(nPts, tmp4, 1, tmp3, 1, tmp4, 1);
                // square curl[2]
                Vmath::Vmul(nPts, tmp4, 1, tmp4, 1, tmp3, 1);
                Vmath::Vadd(nPts, m_curlVelSquare, 1, tmp3, 1, m_curlVelSquare, 1);

            }
            else if ( nDim > 1 )
            {
                // curl[2] 1/2
                Vmath::Vadd(nPts, tmp4, 1, derivativesO1[0][1], 1, tmp4, 1);
                // curl[2]-2/1
                Vmath::Vcopy(nPts, derivativesO1[1][0], 1, tmp3, 1);
                Vmath::Neg(nPts, tmp3, 1);
                Vmath::Vadd(nPts, tmp4, 1, tmp3, 1, tmp4, 1);
                // square curl[2]
                Vmath::Vmul(nPts, tmp4, 1, tmp4, 1, tmp3, 1);
                Vmath::Vadd(nPts, m_curlVelSquare, 1, tmp3, 1, m_curlVelSquare, 1);
            }
            else
            {
                Vmath::Fill(nPts, 0.0, m_curlVelSquare, 1);
            }
        }
    }
}
