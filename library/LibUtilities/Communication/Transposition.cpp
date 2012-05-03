///////////////////////////////////////////////////////////////////////////////
//
// File Transposition.cpp
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
// Description: Data manager for homogeneous transpositions
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Communication/Transposition.h>

namespace Nektar
{
    namespace LibUtilities
    {
		//=================================================================================
		// Constructors
		
		Transposition::Transposition(const LibUtilities::BasisKey &HomoBasis0,
									 LibUtilities::CommSharedPtr hcomm)
		{
			m_hcomm = hcomm;
			m_num_homogeneous_directions = 1;
			
			m_num_points_per_proc    = Array<OneD,int>(m_num_homogeneous_directions);
			m_num_homogeneous_points = Array<OneD,int>(m_num_homogeneous_directions);
			m_num_homogeneous_coeffs = Array<OneD,int>(m_num_homogeneous_directions);
			m_num_points_per_proc    = Array<OneD,int>(m_num_homogeneous_directions);
			m_num_processes          = Array<OneD,int>(m_num_homogeneous_directions);
			
			m_num_homogeneous_points[0]  = HomoBasis0.GetNumPoints();
			m_num_homogeneous_coeffs[0]  = HomoBasis0.GetNumModes();
			m_num_processes[0]           = m_hcomm->GetSize();
			m_num_points_per_proc[0]     = m_num_homogeneous_points[0]/m_num_processes[0];
			m_rank_id                    = m_hcomm->GetRank();
			
			//=================================================================================
			// TODO: Need to be generalised for 1D, 2D and 3D
			
			m_planes_IDs = Array<OneD, unsigned int>(m_num_points_per_proc[0]);
			m_K          = Array<OneD, unsigned int>(m_num_points_per_proc[0]);
			
			for(int i = 0 ; i < m_num_points_per_proc[0] ; i++)
			{
				m_planes_IDs[i] = m_rank_id*m_num_points_per_proc[0] + i;
			}
			
			if(HomoBasis0.GetBasisType() == LibUtilities::eFourier)
			{
				for(int i = 0 ; i < m_num_points_per_proc[0] ; i++)
				{
					m_K[i] = m_planes_IDs[i]/2;
				}
			}
			
			if(HomoBasis0.GetBasisType() == LibUtilities::eFourierSingleMode)
			{
				m_K[0] = 1;
				m_K[1] = 1;
			}
			
			//=================================================================================
        }
		
		Transposition::Transposition(const LibUtilities::BasisKey &HomoBasis0,
									 const LibUtilities::BasisKey &HomoBasis1,
									 LibUtilities::CommSharedPtr hcomm)
		{
			m_hcomm = hcomm;
			m_num_homogeneous_directions = 2;
			
			m_num_points_per_proc    = Array<OneD,int>(m_num_homogeneous_directions);
			m_num_homogeneous_points = Array<OneD,int>(m_num_homogeneous_directions);
			m_num_homogeneous_coeffs = Array<OneD,int>(m_num_homogeneous_directions);
			m_num_points_per_proc    = Array<OneD,int>(m_num_homogeneous_directions);
			m_num_processes          = Array<OneD,int>(m_num_homogeneous_directions);
			
			//=================================================================================
			// TODO: Need set up for 2D
			//=================================================================================
		}
		
		Transposition::Transposition(const LibUtilities::BasisKey &HomoBasis0,
									 const LibUtilities::BasisKey &HomoBasis1,
									 const LibUtilities::BasisKey &HomoBasis2,
									 LibUtilities::CommSharedPtr hcomm)
		{
			m_hcomm = hcomm;
			m_num_homogeneous_directions = 3;
			
			m_num_points_per_proc    = Array<OneD,int>(m_num_homogeneous_directions);
			m_num_homogeneous_points = Array<OneD,int>(m_num_homogeneous_directions);
			m_num_homogeneous_coeffs = Array<OneD,int>(m_num_homogeneous_directions);
			m_num_points_per_proc    = Array<OneD,int>(m_num_homogeneous_directions);
			m_num_processes          = Array<OneD,int>(m_num_homogeneous_directions);
			
			//=================================================================================
			// TODO: Need set up for 3D
			//=================================================================================
		}
		
		
		//=================================================================================
		// Destructors
        Transposition::~Transposition()
        {

        }
		
		//=================================================================================
		// TODO: Need to generalise the following methods for 1D, 2D and 3D
		
		unsigned int Transposition::GetK(int i)
		{
			return m_K[i];
		}
		
		Array<OneD, unsigned int> Transposition::GetKs(void)
		{
			return m_K;
		}
		
		unsigned int Transposition::GetPlaneID(int i)
		{
			return m_planes_IDs[i];
		}
		
		Array<OneD, unsigned int> Transposition::GetPlanesIDs(void)
		{
			return m_planes_IDs;
		}
		
		//=================================================================================
		// Main method: General transposition, the dir parameters define if 1D,2D,3D 
	    // and which transposition is required at the same time 
		
		void Transposition::Transpose(const Array<OneD,const NekDouble> &inarray,Array<OneD, NekDouble> &outarray,bool UseNumMode,TranspositionDir dir)
		{
			switch(dir)
			{
				case eXYtoZ:
				{
					TransposeXYtoZ(inarray,outarray,UseNumMode);
				}
				break;
				case eZtoXY:
				{
					TransposeZtoXY(inarray,outarray,UseNumMode);
				}
				break;
				case eXtoYZ:
				{
					ASSERTL0(false,"Transposition not implemented yet.");
				}
				break;
				case eYZtoX:
				{
					ASSERTL0(false,"Transposition not implemented yet.");
				}
				break;
				case eYZtoZY:
				{
					ASSERTL0(false,"Transposition not implemented yet.");
				}
				break;
				case eZYtoYZ:
				{
					ASSERTL0(false,"Transposition not implemented yet.");
				}
				break;
				case eXtoY:
				{
					ASSERTL0(false,"Transposition not implemented yet.");
				}
				break;
				case eYtoZ:
				{
					ASSERTL0(false,"Transposition not implemented yet.");
				}
				break;
				case eZtoX:
				{
					ASSERTL0(false,"Transposition not implemented yet.");
				}
				break;
				default:
				{
					ASSERTL0(false,"Transposition type does not exist.");
				}
			}
		}
		
		//=================================================================================
		// Homogeneous 1D transposition from SEM to Homogeneous ordering
		
		void Transposition::TransposeXYtoZ(const Array<OneD,const NekDouble> &inarray,Array<OneD, NekDouble> &outarray,bool UseNumMode)
		{
			if(m_num_processes[0] > 1)
			{
				int i,j;
				int num_dofs             = inarray.num_elements();
				int num_points_per_plane = num_dofs/m_num_points_per_proc[0];
			    int num_pencil_per_proc  = num_points_per_plane/m_num_processes[0] + (num_points_per_plane%m_num_processes[0] > 0);
			    int copy_len             = num_pencil_per_proc;
				int packed_len;
			
			    m_SizeMap = Array<OneD,int> (m_num_processes[0],0);
			    m_OffsetMap = Array<OneD,int> (m_num_processes[0],0);
			
			    Array< OneD, NekDouble> tmp_outarray(num_pencil_per_proc*m_num_homogeneous_points[0],0.0);
			
			    for(i = 0; i < m_num_processes[0] ; i++)
			    {
					if(i == m_num_processes[0]-1)
					{
						copy_len = num_pencil_per_proc - ((num_pencil_per_proc*m_num_processes[0]) - num_points_per_plane);
					}
					
					for(j = 0; j < m_num_points_per_proc[0];j++)
					{
						Vmath::Vcopy(copy_len,
								 &(inarray[i*num_pencil_per_proc+j*num_points_per_plane]),1,
								 &(outarray[i*num_pencil_per_proc*m_num_points_per_proc[0]+j*num_pencil_per_proc]),1);
					}
				
					m_SizeMap[i]   = num_pencil_per_proc*m_num_points_per_proc[0];
					m_OffsetMap[i] = i*num_pencil_per_proc*m_num_points_per_proc[0];				
				}
			
				m_hcomm->AlltoAllv(outarray,m_SizeMap,m_OffsetMap,tmp_outarray,m_SizeMap,m_OffsetMap);
				
				if(UseNumMode)
				{
					packed_len = m_num_homogeneous_coeffs[0];
				}
				else
				{
					packed_len = m_num_homogeneous_points[0];
				}
			
				ASSERTL1(&inarray[0] != &outarray[0],"Inarray and outarray cannot be the same");
			
				for(i = 0; i < packed_len; ++i)
				{
					Vmath::Vcopy(num_pencil_per_proc,&(tmp_outarray[i*num_pencil_per_proc]),1,
							 &(outarray[i]),packed_len);
				}
			}
			else 
			{
				int i, pts_per_plane;
				int n = inarray.num_elements();
				int packed_len;
				
				pts_per_plane = n/m_num_points_per_proc[0];
				
				if(UseNumMode)
				{
					packed_len = m_num_homogeneous_coeffs[0];
				}
				else
				{
					packed_len = m_num_homogeneous_points[0];
				}
				
				ASSERTL1(&inarray[0] != &outarray[0],"Inarray and outarray cannot be the same");
				
				for(i = 0; i < packed_len; ++i)
				{
					Vmath::Vcopy(pts_per_plane,&(inarray[i*pts_per_plane]),1,
								 &(outarray[i]),packed_len);
				}
				
			}
		}

		//=================================================================================
		// Homogeneous 1D transposition from Homogeneous to SEM ordering
		
		void Transposition::TransposeZtoXY(const Array<OneD,const NekDouble> &inarray,Array<OneD, NekDouble> &outarray,bool UseNumMode)
		{
			if(m_num_processes[0] > 1)
			{
				int i,j;
				int num_dofs             = inarray.num_elements();
				int num_points_per_plane = num_dofs/m_num_points_per_proc[0];
				int num_pencil_per_proc  = num_points_per_plane/m_num_processes[0] + (num_points_per_plane%m_num_processes[0] > 0);
				int copy_len             = num_pencil_per_proc;
				int packed_len;
			
				m_SizeMap = Array<OneD,int> (m_num_processes[0],0);
				m_OffsetMap = Array<OneD,int> (m_num_processes[0],0);
			
				Array< OneD, NekDouble> tmp_inarray(num_pencil_per_proc*m_num_homogeneous_points[0],0.0);
				Array< OneD, NekDouble> tmp_outarray(num_pencil_per_proc*m_num_homogeneous_points[0],0.0);
			
				if(UseNumMode)
				{
					packed_len = m_num_homogeneous_coeffs[0];
				}
				else
				{
					packed_len = m_num_homogeneous_points[0];
				}
			
				ASSERTL1(&inarray[0] != &outarray[0],"Inarray and outarray cannot be the same");
			
				for(i = 0; i < packed_len; ++i)
				{
					Vmath::Vcopy(num_pencil_per_proc,&(inarray[i]),packed_len,
                             &(tmp_inarray[i*num_pencil_per_proc]),1);
				}
			
				for(i = 0; i < m_num_processes[0] ; i++)
				{				
					m_SizeMap[i]   = num_pencil_per_proc*m_num_points_per_proc[0];
					m_OffsetMap[i] = i*num_pencil_per_proc*m_num_points_per_proc[0];				
				}
			
				m_hcomm->AlltoAllv(tmp_inarray,m_SizeMap,m_OffsetMap,tmp_outarray,m_SizeMap,m_OffsetMap);

				for(i = 0; i < m_num_processes[0];i++)
				{	
					if(i == m_num_processes[0]-1)
					{
						copy_len = num_pencil_per_proc - ((num_pencil_per_proc*m_num_processes[0])-num_points_per_plane);
					}
					for(j = 0; j < m_num_points_per_proc[0];j++)
					{
						Vmath::Vcopy(copy_len,
								 &(tmp_outarray[i*num_pencil_per_proc*m_num_points_per_proc[0]+j*num_pencil_per_proc]),1,
								 &(outarray[i*num_pencil_per_proc+j*num_points_per_plane]),1);
					}
				}
			}
			else 
			{
				int i,pts_per_plane;
				int n = inarray.num_elements();
				int packed_len;
				
				// use length of inarray to determine data storage type
				// (i.e.modal or physical).
				pts_per_plane = n/m_num_points_per_proc[0];
				
				if(UseNumMode)
				{
					packed_len = m_num_homogeneous_coeffs[0];
				}
				else
				{
					packed_len = m_num_homogeneous_points[0];
				}
				
				ASSERTL1(&inarray[0] != &outarray[0],"Inarray and outarray cannot be the same");
				
				
				for(i = 0; i < packed_len; ++i)
				{
					Vmath::Vcopy(pts_per_plane,&(inarray[i]),packed_len,
								 &(outarray[i*pts_per_plane]),1);
				}
			}
		}
		
		//=================================================================================
		// TODO: Impelement 2D and 3D transposition routines 
    }
}
