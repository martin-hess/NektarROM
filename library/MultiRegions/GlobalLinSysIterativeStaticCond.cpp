///////////////////////////////////////////////////////////////////////////////
//
// File GlobalLinSysIterativeStaticCond.cpp
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
// Description: GlobalLinSysIterativeStaticCond definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/VDmathArray.hpp>
#include <MultiRegions/GlobalLinSysIterativeStaticCond.h>
#include <MultiRegions/LocalToGlobalC0ContMap.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class GlobalLinSysIterativeStaticCond
         *
         * Solves a linear system iteratively using single- or multi-level
         * static condensation.
         */

        /**
         * Registers the class with the Factory.
         */
        string GlobalLinSysIterativeStaticCond::className
                = GetGlobalLinSysFactory().RegisterCreatorFunction(
                    "IterativeStaticCond",
                    GlobalLinSysIterativeStaticCond::create,
                    "Iterative static condensation.");

        string GlobalLinSysIterativeStaticCond::className2
                = GetGlobalLinSysFactory().RegisterCreatorFunction(
                    "IterativeMultiLevelStaticCond",
                    GlobalLinSysIterativeStaticCond::create,
                    "Iterative multi-level static condensation.");

        /**
         * For a matrix system of the form @f[
         * \left[ \begin{array}{cc}
         * \boldsymbol{A} & \boldsymbol{B}\\
         * \boldsymbol{C} & \boldsymbol{D}
         * \end{array} \right]
         * \left[ \begin{array}{c} \boldsymbol{x_1}\\ \boldsymbol{x_2}
         * \end{array}\right]
         * = \left[ \begin{array}{c} \boldsymbol{y_1}\\ \boldsymbol{y_2}
         * \end{array}\right],
         * @f]
         * where @f$\boldsymbol{D}@f$ and
         * @f$(\boldsymbol{A-BD^{-1}C})@f$ are invertible, store and assemble
         * a static condensation system, according to a given local to global
         * mapping. #m_linSys is constructed by AssembleSchurComplement().
         * @param   mKey        Associated matrix key.
         * @param   pLocMatSys  LocalMatrixSystem
         * @param   locToGloMap Local to global mapping.
         */
        GlobalLinSysIterativeStaticCond::GlobalLinSysIterativeStaticCond(
                     const GlobalLinSysKey &pKey,
                     const boost::shared_ptr<ExpList> &pExpList,
                     const boost::shared_ptr<LocalToGlobalBaseMap>
                                                            &pLocToGloMap)
                : GlobalLinSysIterative(pKey, pExpList, pLocToGloMap)
        {
            ASSERTL1((pKey.GetGlobalSysSolnType()==eIterativeStaticCond)||
                     (pKey.GetGlobalSysSolnType()==eIterativeMultiLevelStaticCond),
                     "This constructor is only valid when using static "
                     "condensation");
            ASSERTL1(pKey.GetGlobalSysSolnType()
                        == pLocToGloMap->GetGlobalSysSolnType(),
                     "The local to global map is not set up for the requested "
                     "solution type");

            // Allocate memory for top-level structure
            SetupTopLevel(pLocToGloMap);

            // Construct this level
            Initialise(pLocToGloMap);
        }


        /**
         *
         */
        GlobalLinSysIterativeStaticCond::GlobalLinSysIterativeStaticCond(
                     const GlobalLinSysKey &pKey,
                     const boost::shared_ptr<ExpList> &pExpList,
                     const DNekScalBlkMatSharedPtr pSchurCompl,
                     const DNekScalBlkMatSharedPtr pBinvD,
                     const DNekScalBlkMatSharedPtr pC,
                     const DNekScalBlkMatSharedPtr pInvD,
                     const boost::shared_ptr<LocalToGlobalBaseMap>
                                                            &pLocToGloMap)
                : GlobalLinSysIterative(pKey, pExpList, pLocToGloMap),
                  m_schurCompl ( pSchurCompl ),
                  m_BinvD      ( pBinvD ),
                  m_C          ( pC ),
                  m_invD       ( pInvD )
        {
            // Construct this level
            Initialise(pLocToGloMap);
        }


        /**
         *
         */
        GlobalLinSysIterativeStaticCond::~GlobalLinSysIterativeStaticCond()
        {

        }


        /// Solve the linear system for given input and output vectors.
        void GlobalLinSysIterativeStaticCond::Solve( const Array<OneD,const NekDouble> &in,
                                                 Array<OneD, NekDouble> &out)
        {
            ASSERTL0(false, "NOT IMPLEMENTED");
        }


        /**
         * Solve a global linear system using the conjugate gradient method.
         * We solve only for the non-Dirichlet modes. The operator is evaluated
         * using the local-matrix representation. Distributed math routines are
         * used to support parallel execution of the solver.
         * @param       pInput      Input vector of non-Dirichlet DOFs.
         * @param       pOutput     Solution vector of non-Dirichlet DOFs.
         */
        void GlobalLinSysIterativeStaticCond::Solve(
                    const NekVector<NekDouble> &pInput,
                          NekVector<NekDouble> &pOutput)
        {
            ASSERTL1(pInput.GetRows()  == m_gmat->GetColumns(),
                     "Input rows do not match number of matrix columns.");
            ASSERTL1(pOutput.GetRows() == m_gmat->GetRows(),
                     "Output rows do not match number of matrix rows.");

            // Get the communicator for performing data exchanges
            LibUtilities::CommSharedPtr vComm = m_expList->GetComm();

            int nGlobHomBndDofs = pInput.GetRows();

            // Allocate array storage
            Array<OneD, NekDouble> d_A    (nGlobHomBndDofs, 0.0);
            Array<OneD, NekDouble> p_A    (nGlobHomBndDofs, 0.0);
            Array<OneD, NekDouble> z_A    (nGlobHomBndDofs, 0.0);
            Array<OneD, NekDouble> z_new_A(nGlobHomBndDofs, 0.0);
            Array<OneD, NekDouble> r_A    (nGlobHomBndDofs, 0.0);
            Array<OneD, NekDouble> r_new_A(nGlobHomBndDofs, 0.0);

            // Create NekVector wrappers for linear algebra operations
            NekVector<NekDouble> r(nGlobHomBndDofs,r_A,eWrapper);
            NekVector<NekDouble> r_new(nGlobHomBndDofs,r_new_A,eWrapper);
            NekVector<NekDouble> z(nGlobHomBndDofs,z_A,eWrapper);
            NekVector<NekDouble> z_new(nGlobHomBndDofs,z_new_A,eWrapper);
            NekVector<NekDouble> d(nGlobHomBndDofs,d_A, eWrapper);
            NekVector<NekDouble> p(nGlobHomBndDofs,p_A,eWrapper);

            int k;
            NekDouble alpha, beta, norm;

            // INVERSE of preconditioner matrix.
            DNekMat &M = m_preconditioner;

            // Initialise with zero as the initial guess.
            r = pInput;
            z = M * r;
            d = z;
            k = 0;

            //Array<OneD, int> map = m_locToGloMap->GetGlobalToUniversalMapUnique();
            Array<OneD, int> map(nGlobHomBndDofs, 1);

            // If input vector is zero, set zero output and skip solve.
            if (VDmath::Ddot2(vComm, nGlobHomBndDofs, r_A, 1, r_A, 1, map, 1)
                    < NekConstants::kNekZeroTol)
            {
                Vmath::Zero(nGlobHomBndDofs, pOutput.GetRawPtr(), 1);
                return;
            }

            // Continue until convergence
            while (true)
            {
                p = (*m_gmat)*d;

                // compute step length
                alpha = VDmath::Ddot2(vComm, nGlobHomBndDofs,
                                        d_A, 1,
                                        p_A, 1,
                                        map, 1);
                alpha = VDmath::Ddot2(vComm, nGlobHomBndDofs,
                                        z_A, 1,
                                        r_A, 1,
                                        map, 1) / alpha;

                // approximate solution
                pOutput   = pOutput + alpha*d;

                // compute residual
                r_new = r   - alpha*p;

                // Test if residual is small enough
                norm = VDmath::Ddot2(vComm, nGlobHomBndDofs,
                                        r_new_A, 1,
                                        r_new_A, 1,
                                        map, 1);

                if (sqrt(norm) < NekConstants::kNekIterativeTol)
                {
                    break;
                }

                // Apply preconditioner to new residual
                z_new = M * r_new;

                // Improvement achieved
                beta = VDmath::Ddot2(vComm, nGlobHomBndDofs,
                                        r_A, 1,
                                        z_A, 1,
                                        map, 1);
                beta = VDmath::Ddot2(vComm, nGlobHomBndDofs,
                                        r_new_A, 1,
                                        z_new_A, 1,
                                        map, 1) / beta;

                // Compute new search direction
                d = z_new + beta*d;

                // Next step
                r = r_new;
                z = z_new;
                k++;

                ASSERTL1(k < 20000,
                         "Exceeded maximum number of iterations (20000)");
            }
        }


        /**
         *
         */
        void GlobalLinSysIterativeStaticCond::Solve(
                    const Array<OneD, const NekDouble>  &in,
                          Array<OneD,       NekDouble>  &out,
                    const LocalToGlobalBaseMapSharedPtr &pLocToGloMap,
                    const Array<OneD, const NekDouble>  &dirForcing)
        {
            bool dirForcCalculated = (bool) dirForcing.num_elements();
            bool atLastLevel       = pLocToGloMap->AtLastLevel();

            int nGlobDofs          = pLocToGloMap->GetNumGlobalCoeffs();
            int nGlobBndDofs       = pLocToGloMap->GetNumGlobalBndCoeffs();
            int nDirBndDofs        = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            int nGlobHomBndDofs    = nGlobBndDofs - nDirBndDofs;
            int nLocBndDofs        = pLocToGloMap->GetNumLocalBndCoeffs();
            int nIntDofs           = pLocToGloMap->GetNumGlobalCoeffs()
                                                                - nGlobBndDofs;

            Array<OneD, NekDouble> F(nGlobDofs);
            if(nDirBndDofs && dirForcCalculated)
            {
                Vmath::Vsub(nGlobDofs,in.get(),1,dirForcing.get(),1,F.get(),1);
            }
            else
            {
                Vmath::Vcopy(nGlobDofs,in.get(),1,F.get(),1);
            }

            NekVector<NekDouble> F_HomBnd(nGlobHomBndDofs,F+nDirBndDofs,
                                          eWrapper);
            NekVector<NekDouble> F_Int(nIntDofs,F+nGlobBndDofs,eWrapper);

            NekVector<NekDouble> V_GlobBnd(nGlobBndDofs,out,eWrapper);
            NekVector<NekDouble> V_GlobHomBnd(nGlobHomBndDofs,out+nDirBndDofs,
                                              eWrapper);
            NekVector<NekDouble> V_Int(nIntDofs,out+nGlobBndDofs,eWrapper);
            NekVector<NekDouble> V_LocBnd(nLocBndDofs,0.0);

            NekVector<NekDouble> V_GlobHomBndTmp(nGlobHomBndDofs,0.0);

            if(nGlobHomBndDofs)
            {
                if(nIntDofs || ((nDirBndDofs) && (!dirForcCalculated)
                                              && (atLastLevel)) )
                {
                    // construct boundary forcing
                    if( nIntDofs  && ((nDirBndDofs) && (!dirForcCalculated)
                                                    && (atLastLevel)) )
                    {
                        //include dirichlet boundary forcing
                        DNekScalBlkMat &BinvD      = *m_BinvD;
                        DNekScalBlkMat &SchurCompl = *m_schurCompl;
                        pLocToGloMap->GlobalToLocalBnd(V_GlobBnd,V_LocBnd);
                        V_LocBnd = BinvD*F_Int + SchurCompl*V_LocBnd;
                    }
                    else if((nDirBndDofs) && (!dirForcCalculated)
                                          && (atLastLevel))
                    {
                        //include dirichlet boundary forcing
                        DNekScalBlkMat &SchurCompl = *m_schurCompl;
                        pLocToGloMap->GlobalToLocalBnd(V_GlobBnd,V_LocBnd);
                        V_LocBnd = SchurCompl*V_LocBnd;
                    }
                    else
                    {
                        DNekScalBlkMat &BinvD      = *m_BinvD;
                        V_LocBnd = BinvD*F_Int;
                    }
                    pLocToGloMap->AssembleBnd(V_LocBnd,V_GlobHomBndTmp,
                                             nDirBndDofs);
                    F_HomBnd = F_HomBnd - V_GlobHomBndTmp;
                }

                // solve boundary system
                if(atLastLevel)
                {
                    //Solve(F_HomBnd,V_GlobHomBnd);
                    Solve(F_HomBnd,V_GlobHomBnd);
                }
                else
                {
                    m_recursiveSchurCompl->Solve(F,
                                V_GlobBnd.GetPtr(),
                                pLocToGloMap->GetNextLevelLocalToGlobalMap());
                }
            }

            // solve interior system
            if(nIntDofs)
            {
                DNekScalBlkMat &invD  = *m_invD;

                if(nGlobHomBndDofs || nDirBndDofs)
                {
                    DNekScalBlkMat &C     = *m_C;

                    if(dirForcCalculated && nDirBndDofs)
                    {
                        pLocToGloMap->GlobalToLocalBnd(V_GlobHomBnd,V_LocBnd,
                                                      nDirBndDofs);
                    }
                    else
                    {
                        pLocToGloMap->GlobalToLocalBnd(V_GlobBnd,V_LocBnd);
                    }
                    F_Int = F_Int - C*V_LocBnd;
                }

                V_Int = invD*F_Int;
            }
        }



        /**
         * If at the last level of recursion (or the only level in the case of
         * single-level static condensation), assemble the Schur complement.
         * For other levels, in the case of multi-level static condensation,
         * the next level of the condensed system is computed.
         * @param   pLocToGloMap    Local to global mapping.
         */
        void GlobalLinSysIterativeStaticCond::Initialise(
                const boost::shared_ptr<LocalToGlobalBaseMap>& pLocToGloMap)
        {
            if(pLocToGloMap->AtLastLevel())
            {
                AssembleSchurComplement(pLocToGloMap);
            }
            else
            {
                ConstructNextLevelCondensedSystem(
                        pLocToGloMap->GetNextLevelLocalToGlobalMap());
            }
        }


        /**
         * For the first level in multi-level static condensation, or the only
         * level in the case of single-level static condensation, allocate the
         * condensed matrices and populate them with the local matrices
         * retrieved from the expansion list.
         * @param
         */
        void GlobalLinSysIterativeStaticCond::SetupTopLevel(
                const boost::shared_ptr<LocalToGlobalBaseMap>& pLocToGloMap)
        {
            int n;
            int n_exp = m_expList->GetNumElmts();

            const Array<OneD,const unsigned int>& nbdry_size
                    = pLocToGloMap->GetNumLocalBndCoeffsPerPatch();
            const Array<OneD,const unsigned int>& nint_size
                    = pLocToGloMap->GetNumLocalIntCoeffsPerPatch();

            // Setup Block Matrix systems
            MatrixStorage blkmatStorage = eDIAGONAL;
            m_schurCompl = MemoryManager<DNekScalBlkMat>
                    ::AllocateSharedPtr(nbdry_size, nbdry_size, blkmatStorage);
            m_BinvD      = MemoryManager<DNekScalBlkMat>
                    ::AllocateSharedPtr(nbdry_size, nint_size , blkmatStorage);
            m_C          = MemoryManager<DNekScalBlkMat>
                    ::AllocateSharedPtr(nint_size , nbdry_size, blkmatStorage);
            m_invD       = MemoryManager<DNekScalBlkMat>
                    ::AllocateSharedPtr(nint_size , nint_size , blkmatStorage);

            for(n = 0; n < n_exp; ++n)
            {
                if (m_linSysKey.GetMatrixType() == StdRegions::eHybridDGHelmBndLam)
                {
                    DNekScalMatSharedPtr loc_mat = GetBlock(n);
                    m_schurCompl->SetBlock(n,n,loc_mat);
                }
                else
                {
                    DNekScalBlkMatSharedPtr loc_mat = GetStaticCondBlock(n);
                    DNekScalMatSharedPtr tmp_mat;
                    m_schurCompl->SetBlock(n,n, tmp_mat = loc_mat->GetBlock(0,0));
                    m_BinvD     ->SetBlock(n,n, tmp_mat = loc_mat->GetBlock(0,1));
                    m_C         ->SetBlock(n,n, tmp_mat = loc_mat->GetBlock(1,0));
                    m_invD      ->SetBlock(n,n, tmp_mat = loc_mat->GetBlock(1,1));
                }
            }
        }


        /**
         * Assemble the schur complement matrix from the block matrices stored
         * in #m_blkMatrices and the given local to global mapping information.
         * @param   locToGloMap Local to global mapping information.
         */
        void GlobalLinSysIterativeStaticCond::AssembleSchurComplement(
                    const LocalToGlobalBaseMapSharedPtr &pLocToGloMap)
        {
            int i,j,n,cnt,gid1,gid2;
            NekDouble sign1,sign2,value;

            int nBndDofs  = pLocToGloMap->GetNumGlobalBndCoeffs();
            int NumDirBCs = pLocToGloMap->GetNumGlobalDirBndCoeffs();

            DNekScalBlkMatSharedPtr SchurCompl = m_schurCompl;
            DNekScalBlkMatSharedPtr BinvD      = m_BinvD;
            DNekScalBlkMatSharedPtr C          = m_C;
            DNekScalBlkMatSharedPtr invD       = m_invD;

            unsigned int rows = nBndDofs - NumDirBCs;
            unsigned int cols = nBndDofs - NumDirBCs;
            NekDouble zero = 0.0;


            int bwidth = pLocToGloMap->GetBndSystemBandWidth();
            MatrixStorage matStorage;

            switch(m_linSysKey.GetMatrixType())
            {
                // case for all symmetric matices
            case StdRegions::eMass:
            case StdRegions::eLaplacian:
            case StdRegions::eHelmholtz:
            case StdRegions::eHybridDGHelmBndLam:
                {
                    if( (2*(bwidth+1)) < rows)
                    {
                        try {
                            matStorage = ePOSITIVE_DEFINITE_SYMMETRIC_BANDED;
                            m_gmat = MemoryManager<DNekMat>
                                ::AllocateSharedPtr(rows, cols, zero,
                                                    matStorage,
                                                    bwidth, bwidth);
                        }
                        catch (...) {
                            NEKERROR(ErrorUtil::efatal,
                                     "Insufficient memory for GlobalLinSys.");
                        }
                    }
                    else
                    {
                        matStorage = ePOSITIVE_DEFINITE_SYMMETRIC;
                        m_gmat = MemoryManager<DNekMat>
                                        ::AllocateSharedPtr(rows, cols, zero,
                                                            matStorage);
                    }
                }
                break;
            case StdRegions::eLinearAdvectionReaction:
            case StdRegions::eLinearAdvectionDiffusionReaction:
                {
                    // Current inversion techniques do not seem to
                    // allow banded matrices to be used as a linear
                    // system
                    matStorage = eFULL;
                    m_gmat = MemoryManager<DNekMat>
                            ::AllocateSharedPtr(rows, cols, zero, matStorage);

                }
                break;
            default:
                {
                    NEKERROR(ErrorUtil::efatal, "Add MatrixType to switch "
                             "statement");
                }
            }

            // fill global matrix
            DNekScalMatSharedPtr loc_mat;
            int loc_lda;
            for(n = cnt = 0; n < SchurCompl->GetNumberOfBlockRows(); ++n)
            {
                loc_mat = SchurCompl->GetBlock(n,n);
                loc_lda = loc_mat->GetRows();

                // Set up  Matrix;
                for(i = 0; i < loc_lda; ++i)
                {
                    gid1  = pLocToGloMap->GetLocalToGlobalBndMap (cnt + i)
                                                                    - NumDirBCs;
                    sign1 = pLocToGloMap->GetLocalToGlobalBndSign(cnt + i);

                    if(gid1 >= 0)
                    {
                        for(j = 0; j < loc_lda; ++j)
                        {
                            gid2  = pLocToGloMap->GetLocalToGlobalBndMap(cnt+j)
                                                                 - NumDirBCs;
                            sign2 = pLocToGloMap->GetLocalToGlobalBndSign(cnt+j);

                            if(gid2 >= 0)
                            {
                                // As the global matrix should be
                                // symmetric, only add the value for
                                // the upper triangular part in order
                                // to avoid entries to be entered
                                // twice
                                if((matStorage == eFULL)||(gid2 >= gid1))
                                {
                                    value = m_gmat->GetValue(gid1,gid2)
                                                + sign1*sign2*(*loc_mat)(i,j);
                                    m_gmat->SetValue(gid1,gid2,value);
                                }
                            }
                        }
                    }
                }
                cnt += loc_lda;
            }
            ComputeDiagonalPreconditioner(pLocToGloMap);
        }


        /**
         *
         */
        void GlobalLinSysIterativeStaticCond::ConstructNextLevelCondensedSystem(
                        const LocalToGlobalBaseMapSharedPtr& pLocToGloMap)
        {
            int i,j,n,cnt;
            NekDouble one  = 1.0;
            NekDouble zero = 0.0;
            DNekScalBlkMatSharedPtr blkMatrices[4];

            // Create temporary matrices within an inner-local scope to ensure
            // any references to the intermediate storage is lost before
            // the recursive step, rather than at the end of the routine.
            // This allows the schur complement matrix from this level to be
            // disposed of in the next level after use without being retained
            // due to lingering shared pointers.
            {

                const Array<OneD,const unsigned int>& nBndDofsPerPatch
                                = pLocToGloMap->GetNumLocalBndCoeffsPerPatch();
                const Array<OneD,const unsigned int>& nIntDofsPerPatch
                                = pLocToGloMap->GetNumLocalIntCoeffsPerPatch();

                // STEP 1:
                // Based upon the schur complement of the the current level we
                // will substructure this matrix in the form
                //      --     --
                //      | A   B |
                //      | C   D |
                //      --     --
                // All matrices A,B,C and D are (diagonal) blockmatrices.
                // However, as a start we will use an array of DNekMatrices as
                // it is too hard to change the individual entries of a
                // DNekScalBlkMatSharedPtr.

                // In addition, we will also try to ensure that the memory of
                // the blockmatrices will be contiguous. This will probably
                // enhance the efficiency
                // - Calculate the total number of entries in the blockmatrices
                int nPatches  = pLocToGloMap->GetNumPatches();
                int nEntriesA = 0; int nEntriesB = 0;
                int nEntriesC = 0; int nEntriesD = 0;

                for(i = 0; i < nPatches; i++)
                {
                    nEntriesA += nBndDofsPerPatch[i]*nBndDofsPerPatch[i];
                    nEntriesB += nBndDofsPerPatch[i]*nIntDofsPerPatch[i];
                    nEntriesC += nIntDofsPerPatch[i]*nBndDofsPerPatch[i];
                    nEntriesD += nIntDofsPerPatch[i]*nIntDofsPerPatch[i];
                }

                // Now create the DNekMatrices and link them to the memory
                // allocated above
                Array<OneD, DNekMatSharedPtr> substructeredMat[4]
                    = {Array<OneD, DNekMatSharedPtr>(nPatches),  //Matrix A
                       Array<OneD, DNekMatSharedPtr>(nPatches),  //Matrix B
                       Array<OneD, DNekMatSharedPtr>(nPatches),  //Matrix C
                       Array<OneD, DNekMatSharedPtr>(nPatches)}; //Matrix D

                // Initialise storage for the matrices. We do this separately
                // for each matrix so the matrices may be independently
                // deallocated when no longer required.
                Array<OneD, NekDouble> storageA(nEntriesA,0.0);
                Array<OneD, NekDouble> storageB(nEntriesB,0.0);
                Array<OneD, NekDouble> storageC(nEntriesC,0.0);
                Array<OneD, NekDouble> storageD(nEntriesD,0.0);

                Array<OneD, NekDouble> tmparray;
                PointerWrapper wType = eWrapper;
                int cntA = 0;
                int cntB = 0;
                int cntC = 0;
                int cntD = 0;

                for(i = 0; i < nPatches; i++)
                {
                    // Matrix A
                    tmparray = storageA+cntA;
                    substructeredMat[0][i] = MemoryManager<DNekMat>
                                    ::AllocateSharedPtr(nBndDofsPerPatch[i],
                                                        nBndDofsPerPatch[i],
                                                        tmparray, wType);
                    // Matrix B
                    tmparray = storageB+cntB;
                    substructeredMat[1][i] = MemoryManager<DNekMat>
                                    ::AllocateSharedPtr(nBndDofsPerPatch[i],
                                                        nIntDofsPerPatch[i],
                                                        tmparray, wType);
                    // Matrix C
                    tmparray = storageC+cntC;
                    substructeredMat[2][i] = MemoryManager<DNekMat>
                                    ::AllocateSharedPtr(nIntDofsPerPatch[i],
                                                        nBndDofsPerPatch[i],
                                                        tmparray, wType);
                    // Matrix D
                    tmparray = storageD+cntD;
                    substructeredMat[3][i] = MemoryManager<DNekMat>
                                    ::AllocateSharedPtr(nIntDofsPerPatch[i],
                                                        nIntDofsPerPatch[i],
                                                        tmparray, wType);

                    cntA += nBndDofsPerPatch[i] * nBndDofsPerPatch[i];
                    cntB += nBndDofsPerPatch[i] * nIntDofsPerPatch[i];
                    cntC += nIntDofsPerPatch[i] * nBndDofsPerPatch[i];
                    cntD += nIntDofsPerPatch[i] * nIntDofsPerPatch[i];
                }

                // Then, project SchurComplement onto
                // the substructured matrices of the next level
                DNekScalBlkMatSharedPtr SchurCompl  = m_schurCompl;
                DNekScalMatSharedPtr schurComplSubMat;
                int       schurComplSubMatnRows;
                int       patchId_i ,patchId_j;
                int       dofId_i   ,dofId_j;
                bool      isBndDof_i,isBndDof_j;
                NekDouble sign_i    ,sign_j;
                NekDouble value;
                int       ABCorD;
                for(n = cnt = 0; n < SchurCompl->GetNumberOfBlockRows(); ++n)
                {
                    schurComplSubMat      = SchurCompl->GetBlock(n,n);
                    schurComplSubMatnRows = schurComplSubMat->GetRows();

                    // Set up  Matrix;
                    for(i = 0; i < schurComplSubMatnRows; ++i)
                    {
                        patchId_i  = pLocToGloMap->GetPatchMapFromPrevLevel(cnt+i)->GetPatchId();
                        dofId_i    = pLocToGloMap->GetPatchMapFromPrevLevel(cnt+i)->GetDofId();
                        isBndDof_i = pLocToGloMap->GetPatchMapFromPrevLevel(cnt+i)->IsBndDof();
                        sign_i     = pLocToGloMap->GetPatchMapFromPrevLevel(cnt+i)->GetSign();

                        for(j = 0; j < schurComplSubMatnRows; ++j)
                        {
                            patchId_j  = pLocToGloMap->GetPatchMapFromPrevLevel(cnt+j)->GetPatchId();
                            dofId_j    = pLocToGloMap->GetPatchMapFromPrevLevel(cnt+j)->GetDofId();
                            isBndDof_j = pLocToGloMap->GetPatchMapFromPrevLevel(cnt+j)->IsBndDof();
                            sign_j     = pLocToGloMap->GetPatchMapFromPrevLevel(cnt+j)->GetSign();

                            ASSERTL0(patchId_i==patchId_j,"These values should be equal");

                            ABCorD = 2*(isBndDof_i?0:1)+(isBndDof_j?0:1);
                            value = substructeredMat[ABCorD][patchId_i]->GetValue(dofId_i,dofId_j) +
                                sign_i*sign_j*(*schurComplSubMat)(i,j);

                            substructeredMat[ABCorD][patchId_i]->SetValue(dofId_i,dofId_j,value);
                        }
                    }
                    cnt += schurComplSubMatnRows;
                }

                // STEP 2: condense the system
                // This can be done elementally (i.e. patch per patch)
                for(i = 0; i < nPatches; i++)
                {
                    if(nIntDofsPerPatch[i])
                    {
                        // 1. D -> InvD
                        substructeredMat[3][i]->Invert();
                        // 2. B -> BInvD
                        (*substructeredMat[1][i]) = (*substructeredMat[1][i])*(*substructeredMat[3][i]);
                        // 3. A -> A - BInvD*C (= schurcomplement)
                        (*substructeredMat[0][i]) = (*substructeredMat[0][i]) -
                            (*substructeredMat[1][i])*(*substructeredMat[2][i]);
                    }
                }

                // STEP 3: fill the blockmatrices
                // however, do note that we first have to convert them to
                // a DNekScalMat in order to be compatible with the first
                // level of static condensation

                const Array<OneD,const unsigned int>& nbdry_size = pLocToGloMap->GetNumLocalBndCoeffsPerPatch();
                const Array<OneD,const unsigned int>& nint_size  = pLocToGloMap->GetNumLocalIntCoeffsPerPatch();
                MatrixStorage blkmatStorage = eDIAGONAL;

                blkMatrices[0] = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nbdry_size, nbdry_size, blkmatStorage);
                blkMatrices[1] = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nbdry_size, nint_size , blkmatStorage);
                blkMatrices[2] = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nint_size , nbdry_size, blkmatStorage);
                blkMatrices[3] = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nint_size , nint_size , blkmatStorage);

                DNekScalMatSharedPtr tmpscalmat;
                for(i = 0; i < nPatches; i++)
                {
                    for(j = 0; j < 4; j++)
                    {
                        tmpscalmat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,substructeredMat[j][i]);
                        blkMatrices[j]->SetBlock(i,i,tmpscalmat);
                    }
                }
            }

            // We've finished with the Schur complement matrix passed to this
            // level, so return the memory to the system.
            // The Schur complement matrix need only be retained at the last
            // level. Save the other matrices at this level though.
            m_schurCompl.reset();

            m_recursiveSchurCompl = MemoryManager<GlobalLinSysIterativeStaticCond>::
                AllocateSharedPtr(m_linSysKey,m_expList,blkMatrices[0],blkMatrices[1],blkMatrices[2],blkMatrices[3],pLocToGloMap);
        }

        /**
         * Diagonal preconditioner computed by evaluating the local matrix
         * acting on each basis vector (0,...,0,1,0,...,0). (deprecated)
         * @param   pLocToGloMap    Local to global mapping.
         */
        void GlobalLinSysIterativeStaticCond::ComputeDiagonalPreconditioner(
                const boost::shared_ptr<LocalToGlobalBaseMap> &pLocToGloMap)
        {
            int nInt = m_gmat->GetRows();
            m_preconditioner = DNekMat(nInt, nInt, eDIAGONAL);

            for (unsigned int i = 0; i < nInt; ++i)
            {
                m_preconditioner.SetValue(i,i,1.0/(*m_gmat)(i,i));
            }
            cout << "Computed preconditioner" << endl;
        }



    }
}