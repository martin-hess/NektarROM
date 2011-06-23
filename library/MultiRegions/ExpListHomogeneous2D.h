///////////////////////////////////////////////////////////////////////////////
//
// File ExpListHomogeneous2D.h
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
// Description: Base class for expansions which are homogeneous in 2
// directions
//
///////////////////////////////////////////////////////////////////////////////

#ifndef EXPLISTHOMO2D_H
#define EXPLISTHOMO2D_H
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <vector>
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList.h>
#include <LibUtilities/FFT/NektarFFT.h>

namespace Nektar
{
    namespace MultiRegions
    {

        enum Homogeneous2DMatType
        {
            eForwardsCoeffSpace2D,
            eBackwardsCoeffSpace2D,
            eForwardsPhysSpace2D,
            eBackwardsPhysSpace2D
        };

        /// A map between homo matrix keys and their associated block
        /// matrices.
        typedef map< Homogeneous2DMatType, DNekBlkMatSharedPtr> Homo2DBlockMatrixMap;
        /// A shared pointer to a BlockMatrixMap.
        typedef boost::shared_ptr<Homo2DBlockMatrixMap> Homo2DBlockMatrixMapShPtr;

        // Forward declaration for typedefs
        class ExpListHomogeneous2D;

        /// Shared pointer to an ExpList3DHomogeneous2D object.
        typedef boost::shared_ptr<ExpListHomogeneous2D>   ExpListHomogeneous2DSharedPtr;
        /// Vector of pointers to ExpList3DHomogeneous2D objects.
        typedef std::vector< ExpListHomogeneous2DSharedPtr > ExpListHomogeneous2DVector;
        /// Iterator for the vector of ExpList3DHomogeneous2D pointers.
        typedef std::vector< ExpListHomogeneous2DSharedPtr >::iterator ExpListHomogeneous2DVectorIter;

        /// Abstraction of a two-dimensional multi-elemental expansion which
        /// is merely a collection of local expansions.
        class ExpListHomogeneous2D: public ExpList
        {
        public:
            /// Default constructor.
            MULTI_REGIONS_EXPORT ExpListHomogeneous2D();

            MULTI_REGIONS_EXPORT ExpListHomogeneous2D(LibUtilities::CommSharedPtr &pComm,
                                                      const LibUtilities::BasisKey &HomoBasis_y,
													  const LibUtilities::BasisKey &HomoBasis_z,
													  const NekDouble ly,
													  const NekDouble lz,
													  bool useFFT);

            /// Copy constructor.
            MULTI_REGIONS_EXPORT ExpListHomogeneous2D(const ExpListHomogeneous2D &In);

            /// Destructor.
            MULTI_REGIONS_EXPORT ~ExpListHomogeneous2D();

            MULTI_REGIONS_EXPORT void Homogeneous2DTrans(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray, bool IsForwards, bool UseContCoeffs = false);

            inline void Homogeneous2DFwdTrans(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray, bool UseContCoeffs = false);

            inline void Homogeneous2DBwdTrans(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray, bool UseContCoeffs = false);

            MULTI_REGIONS_EXPORT void ShuffleIntoHomogeneous2DClosePacked(
                              const Array<OneD, const NekDouble> &inarray,
                              Array<OneD, NekDouble> &outarray,
                              bool UseNumModes = false);

            MULTI_REGIONS_EXPORT void UnshuffleFromHomogeneous2DClosePacked(
                              const Array<OneD, const NekDouble> &inarray,
                              Array<OneD, NekDouble> &outarray,
                              bool UseNumModes = false);
			
			MULTI_REGIONS_EXPORT void Transpose(const Array<OneD, const NekDouble> &inarray,
												Array<OneD, NekDouble> &outarray,
												bool XtoZ);
			
			MULTI_REGIONS_EXPORT void PhysDeriv(const Array<OneD, const NekDouble> &inarray,
								                Array<OneD, NekDouble> &out_d0,
								                Array<OneD, NekDouble> &out_d1, 
								                Array<OneD, NekDouble> &out_d2, bool UseContCoeffs);
			
			MULTI_REGIONS_EXPORT void PhysDeriv(const int dir,
								                const Array<OneD, const NekDouble> &inarray,
								                Array<OneD, NekDouble> &out_d, bool UseContCoeffs);
			
			/// FFT variables
			bool                                    m_useFFT;
			LibUtilities::NektarFFTSharedPtr        m_FFT_y;
			LibUtilities::NektarFFTSharedPtr        m_FFT_z;
			Array<OneD,NekDouble>                   m_tmpIN;
			Array<OneD,NekDouble>                   m_tmpOUT;

        protected:

            /// Definition of the total number of degrees of freedom and
            /// quadrature points. Sets up the storage for \a m_coeff and \a
            ///  m_phys.
            LibUtilities::BasisSharedPtr    m_homogeneousBasis_y;       ///< Base expansion in y direction
			LibUtilities::BasisSharedPtr    m_homogeneousBasis_z;       ///< Base expansion in z direction
            NekDouble                       m_lhom_y;                   ///< Width of homogeneous direction y
			NekDouble                       m_lhom_z;                   ///< Width of homogeneous direction z
            Homo2DBlockMatrixMapShPtr       m_homogeneous2DBlockMat;
            Array<OneD, ExpListSharedPtr>   m_lines;                    ///< Vector of ExpList, will be filled with ExpList1D
			int                             m_ny;                       ///< Number of modes = number of poitns in y direction
			int                             m_nz;                       ///< Number of modes = number of poitns in z direction

            DNekBlkMatSharedPtr GenHomogeneous2DBlockMatrix(Homogeneous2DMatType mattype, bool UseContCoeffs = false) const;

            DNekBlkMatSharedPtr GetHomogeneous2DBlockMatrix(Homogeneous2DMatType mattype, bool UseContCoeffs = false) const;

            //  virtual functions
            virtual int v_GetNumElmts(void)
            {
                return m_lines[0]->GetExpSize();
            }

            virtual void v_FwdTrans(const Array<OneD,const NekDouble> &inarray,
                                    Array<OneD,      NekDouble> &outarray,
                                    bool  UseContCoeffs);

            virtual void v_BwdTrans(const Array<OneD,const NekDouble> &inarray,
                                    Array<OneD,      NekDouble> &outarray,
                                    bool  UseContCoeffs);

            virtual void v_IProductWRTBase(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray, bool UseContCoeffs);

            virtual std::vector<SpatialDomains::FieldDefinitionsSharedPtr> v_GetFieldDefinitions(void);

            virtual void v_GetFieldDefinitions(std::vector<SpatialDomains::FieldDefinitionsSharedPtr> &fielddef);

            virtual void v_AppendFieldData(SpatialDomains::FieldDefinitionsSharedPtr &fielddef, std::vector<NekDouble> &fielddata);
			
			virtual void v_AppendFieldData(SpatialDomains::FieldDefinitionsSharedPtr &fielddef, std::vector<NekDouble> &fielddata, Array<OneD, NekDouble> &coeffs);

            virtual void v_ExtractDataToCoeffs(SpatialDomains::FieldDefinitionsSharedPtr &fielddef, std::vector<NekDouble> &fielddata, std::string &field);

            virtual void v_WriteTecplotHeader(std::ofstream &outfile,std::string var = "v");

            virtual void v_WriteTecplotField(std::ofstream &outfile,int expansion);

            virtual void v_WriteVtkPieceData(std::ofstream &outfile, int expansion, std::string var);
			
			virtual void v_Homogeneous2DFwdTrans(const Array<OneD, const NekDouble> &inarray, 
												 Array<OneD, NekDouble> &outarray, 
												 bool UseContCoeffs = false);
			
			virtual void v_Homogeneous2DBwdTrans(const Array<OneD, const NekDouble> &inarray, 
												 Array<OneD, NekDouble> &outarray, 
												 bool UseContCoeffs = false);
			
			virtual void v_PhysDeriv(const Array<OneD, const NekDouble> &inarray,
									 Array<OneD, NekDouble> &out_d0,
									 Array<OneD, NekDouble> &out_d1, 
									 Array<OneD, NekDouble> &out_d2, bool UseContCoeffs);
			
			virtual void v_PhysDeriv(const int dir,
									 const Array<OneD, const NekDouble> &inarray,
									 Array<OneD, NekDouble> &out_d, bool UseContCoeffs);

        private:
        };

        inline void ExpListHomogeneous2D::Homogeneous2DFwdTrans(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray, bool UseContCoeffs)
        {
			v_Homogeneous2DFwdTrans(inarray,outarray,UseContCoeffs);
        }
		
        inline void ExpListHomogeneous2D::Homogeneous2DBwdTrans(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray, bool UseContCoeffs)
        {
            v_Homogeneous2DBwdTrans(inarray,outarray,UseContCoeffs);
        }

    } //end of namespace
} //end of namespace

#endif//EXPLISTHOMO2D_H

/**
* $Log: v $
*
**/
