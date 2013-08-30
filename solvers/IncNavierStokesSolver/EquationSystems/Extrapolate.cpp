///////////////////////////////////////////////////////////////////////////////
//
// File: Extrapolate.cpp
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
// Description: Abstract base class for Extrapolate.
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/EquationSystems/Extrapolate.h>
#include <LibUtilities/Communication/Comm.h>

namespace Nektar
{
    ExtrapolateFactory& GetExtrapolateFactory()
    {
        typedef Loki::SingletonHolder<ExtrapolateFactory,
                                      Loki::CreateUsingNew,
                                      Loki::NoDestroy > Type;
        return Type::Instance();
    }

    Extrapolate::Extrapolate(
        const LibUtilities::SessionReaderSharedPtr pSession,
        Array<OneD, MultiRegions::ExpListSharedPtr> pFields,
        Array<OneD, int> pVel)
        : m_session(pSession),
          m_fields(pFields),
          m_velocity(pVel),
          m_comm(pSession->GetComm())
    {
		// count number of fields to get which one is the pressure
		int numfields = pFields.num_elements();
		
		// generating the HOPBC map and setting dimensionality
		// for following calculation (CurlCurl, tec.)
		GenerateHOPBCMap(pFields[numfields-1]);
    }

    Extrapolate::~Extrapolate()
    {
    }

        
    std::string Extrapolate::def =
        LibUtilities::SessionReader::RegisterDefaultSolverInfo(
            "NoSubStepping", "No SubStepping");
	
	/** 
	 * Function to extrapolate the new pressure boundary condition.
	 * Based on the velocity field and on the advection term.
	 * Acceleration term is also computed.
	 * This routine is a general one for 2d and 3D application and it can be called
	 * directly from velocity correction scheme. Specialisation on dimensionality is
	 * redirected to the CalcPressureBCs method.
	 */
	void Extrapolate::EvaluatePressureBCs(const MultiRegions::ExpListSharedPtr &pField,
										  const Array<OneD, const Array<OneD, NekDouble> > &fields,
										  const Array<OneD, const Array<OneD, NekDouble> >  &N,
										  const int kinvis)
	{		
		
		Array<OneD, NekDouble> tmp;
		Array<OneD, NekDouble> accelerationTerm;
		
		int  n,cnt;
		int  nint    = min(m_pressureCalls++,m_intSteps);
		int  nlevels = m_pressureHBCs.num_elements();
		
		int acc_order = 0;
        
		accelerationTerm = Array<OneD, NekDouble>(m_acceleration[0].num_elements(), 0.0);
		
		// Rotate HOPBCs storage
		RollOver(m_pressureHBCs);
		
		// Rotate acceleration term
		RollOver(m_acceleration);
		
		// Calculate BCs at current level
		CalcPressureBCs(pField,fields,N,kinvis);
		
		// Copy High order values into storage array 
		for(cnt = n = 0; n < PBndConds.num_elements(); ++n)
		{
			// High order boundary condition;
			if(m_PBndConds[n]->GetUserDefined() == SpatialDomains::eHigh)
			{
				int nq = m_PBndExp[n]->GetNcoeffs();
				Vmath::Vcopy(nq,&(m_PBndExp[n]->GetCoeffs()[0]),1,&(m_pressureHBCs[0])[cnt],1);
				cnt += nq;
			}
		}
		
		//Calculate acceleration term at level n based on previous steps
		if (m_pressureCalls > 2)
		{
			acc_order = min(m_pressureCalls-2,m_intSteps);
			Vmath::Smul(cnt, StifflyStable_Gamma0_Coeffs[acc_order-1],
						m_acceleration[0], 1,
						accelerationTerm,  1);
			
			for(int i = 0; i < acc_order; i++)
			{
				Vmath::Svtvp(cnt, -1*StifflyStable_Alpha_Coeffs[acc_order-1][i],
							 m_acceleration[i+1], 1,
							 accelerationTerm,    1,
							 accelerationTerm,    1);
			}
		}
		
		// Adding acceleration term to HOPBCs
		Vmath::Svtvp(cnt, -1.0/m_timestep,
					 accelerationTerm,  1,
					 m_pressureHBCs[0], 1,
					 m_pressureHBCs[0], 1);
		
		// Extrapolate to n+1
		Vmath::Smul(cnt, StifflyStable_Betaq_Coeffs[nint-1][nint-1],
					m_pressureHBCs[nint-1],    1,
					m_pressureHBCs[nlevels-1], 1);
		
		for(n = 0; n < nint-1; ++n)
		{
			Vmath::Svtvp(cnt,StifflyStable_Betaq_Coeffs[nint-1][n],
						 m_pressureHBCs[n],1,m_pressureHBCs[nlevels-1],1,
						 m_pressureHBCs[nlevels-1],1);
		}
		
		// Copy values of [dP/dn]^{n+1} in the pressure bcs storage.
		// m_pressureHBCS[nlevels-1] will be cancelled at next time step
		for(cnt = n = 0; n < m_PBndConds.num_elements(); ++n)
		{
			// High order boundary condition;
			if(m_PBndConds[n]->GetUserDefined() == SpatialDomains::eHigh)
			{
				int nq = m_PBndExp[n]->GetNcoeffs();
				Vmath::Vcopy(nq,&(m_pressureHBCs[nlevels-1])[cnt],1,&(m_PBndExp[n]->UpdateCoeffs()[0]),1);
				cnt += nq;
			}
		}
	}
    
	
	/**
     * Unified routine for calculation high-oder terms
     */
    void Extrapolate::CalcPressureBCs(const MultiRegions::ExpListSharedPtr &pField,
									  const Array<OneD, const Array<OneD, NekDouble> > &fields,
									  const Array<OneD, const Array<OneD, NekDouble> >  &N,
									  const int kinvis)
    {	
        Array<OneD, NekDouble> Pvals;
        Array<OneD, NekDouble> Uvals;
        StdRegions::StdExpansionSharedPtr Pbc;
		
        Array<OneD, <Array<OneD, const NekDouble> > Velocity(m_curl_dim);
        Array<OneD, <Array<OneD, const NekDouble> > Advection(m_bnd_dim);
        
        Array<OneD, <Array<OneD, NekDouble> > BndValues(m_bnd_dim);
        Array<OneD, <Array<OneD, NekDouble> > Q(m_bnd_dim);
		
        for(int i = 0; i < m_bnd_dim; i++)
        {
            BndValues[i] = Array<OneD, NekDouble> (m_pressureBCsMaxPts,0.0);
            Q[i]         = Array<OneD, NekDouble> (m_pressureBCsMaxPts,0.0);
        }
		
        for(int j = 0 ; j < m_HBCdata.num_elements() ; j++)
        {
            /// Casting the boundary expansion to the specific case
            Pbc =  boost::dynamic_pointer_cast<StdRegions::StdExpansion> (m_PBndExp[m_HBCdata[j].m_bndryElmtID]->GetExp(m_HBCdata[j].m_bndElmtOffset));
			
			/// Picking up the element where the HOPBc is located
			m_elmt = pField->GetExp(m_HBCdata[j].m_globalElmtID);
			
			/// Assigning 
			for(int i = 0; i < m_bnd_dim; i++)
			{
				Velocity[i]  = fields[i] + m_HBCdata[j].m_physOffset;
				Advection[i] = N[i]      + m_HBCdata[j].m_physOffset;
			}
			
			// for the 3DH1D case we need to grap the conjugate mode
			if(pField->GetExpType() == MultiRegions::e3DH1D)
			{
				Velocity[2]  = fields[2] + m_HBCdata[j].m_assPhysOffset;
			}
			
			/// Calculating the curl-curl and storing it in Q
			CurlCurl(pField,Velocity,Q,j);
			
			// Mounting advection component into the high-order condition
			for(int i = 0; m_bnd_dim; i++)
			{
				MountHOPBCs(m_HBCdata[j].m_ptsInElmt,kinvis,Q[i],Advection[i]);
			}
			
			Pvals = m_PBndExp[m_HBCdata[j].m_bndryElmtID]->UpdateCoeffs()+m_PBndExp[m_HBCdata[j].m_bndryElmtID]->GetCoeff_Offset(m_HBCdata[j].m_bndElmtOffset);
			Uvals = (m_acceleration[0]) + m_HBCdata[j].m_coeffOffset;
			
			// Getting values on the edge and filling the pressure boundary expansion
			// and the acceleration term. Multiplication by the normal is required
			switch(pField->GetExpType())
			{
				case MultiRegions::e2D:
				case MultiRegions::e3DH1D:
				{
					m_elmt->GetEdgePhysVals(m_HBCdata[j].m_elmtTraceID,Pbc,Q[0],BndValues[0]);
					m_elmt->GetEdgePhysVals(m_HBCdata[j].m_elmtTraceID,Pbc,Q[1],BndValues[1]);
					Pbc->NormVectorIProductWRTBase(BndValues[0],BndValues[1],Pvals);
					
					m_elmt->GetEdgePhysVals(m_HBCdata[j].m_elmtTraceID,Pbc,Velocity[0],BndValues[0]);
					m_elmt->GetEdgePhysVals(m_HBCdata[j].m_elmtTraceID,Pbc,Velocity[1],BndValues[1]);
					Pbc->NormVectorIProductWRTBase(BndValues[0],BndValues[1],Uvals);
				}
					break;
					
				case MultiRegions::e3DH2D:
				{
					if(m_HBCdata[j].m_elmtTraceID == 0)
					{
						(m_PBndExp[m_HBCdata[j].m_bndryElmtID]->UpdateCoeffs()+m_PBndExp[m_HBCdata[j].m_bndryElmtID]->GetCoeff_Offset(m_HBCdata[j].m_bndElmtOffset))[0] = -1.0*Q[0][0];
					}
					else if (m_HBCdata[j].m_elmtTraceID == 1)
					{
						(m_PBndExp[m_HBCdata[j].m_bndryElmtID]->UpdateCoeffs()+m_PBndExp[m_HBCdata[j].m_bndryElmtID]->GetCoeff_Offset(m_HBCdata[j].m_bndElmtOffset))[0] = Q[0][m_HBCdata[j].m_ptsInElmt-1];
					}
					else 
					{
						ASSERTL0(false,"In the 3D homogeneous 2D approach BCs edge ID can be just 0 or 1 ");
					}
					
				}
					break;
					
				case MultiRegions::e3D:
				{
					m_elmt->GetFacePhysVals(m_HBCdata[j].m_elmtTraceID,Pbc,Q[0],BndValues[0]);
					m_elmt->GetFacePhysVals(m_HBCdata[j].m_elmtTraceID,Pbc,Q[1],BndValues[1]);
					m_elmt->GetFacePhysVals(m_HBCdata[j].m_elmtTraceID,Pbc,Q[2],BndValues[2]);
					Pbc->NormVectorIProductWRTBase(BndValues[0],BndValues[1],BndValues[2],Pvals);
					
					m_elmt->GetEdgePhysVals(m_HBCdata[j].m_elmtTraceID,Pbc,Velocity[0],BndValues[0]);
					m_elmt->GetEdgePhysVals(m_HBCdata[j].m_elmtTraceID,Pbc,Velocity[1],BndValues[1]);
					m_elmt->GetEdgePhysVals(m_HBCdata[j].m_elmtTraceID,Pbc,Velocity[2],BndValues[2]);
					Pbc->NormVectorIProductWRTBase(BndValues[0],BndValues[1],BndValues[2],Uvals);
				}
					break;
			}
		}
	}
	
	
	/**
	 * Curl Curl routine - dimension dependent
	 */
	void Extrapolate::CurlCurl(const MultiRegionss::ExpListSharedPtr &pField,
							   const Array<OneD, Array<OneD, NekDouble> > &Vel,
							   Array<OneD, Array<OneD, NekDouble> > &Q,
							   const int j)
	{
		m_elmt = pField->GetExp(m_HBCdata[j].m_globalElmtID);
		
		Array<OneD,NekDouble> Vx(m_pressureBCsMaxPts);
		Array<OneD,NekDouble> Uy(m_pressureBCsMaxPts);
		
		switch(pField->GetExpType())
		{
			case MultiRegions::e2D:
			{
				m_elmt->PhysDeriv(MultiRegionss::DirCartesianMap[0],Vel[1],Vx);
				m_elmt->PhysDeriv(MultiRegionss::DirCartesianMap[1],Vel[2],Uy);  
				
				Vmath::Vsub(m_HBCdata[j].m_ptsInElmt,Vx,1,Uy,1,Dummy,1);
				
				m_elmt->PhysDeriv(Dummy,Q[1],Q[0]);
				
				Vmath::Smul(m_HBCdata[j].m_ptsInElmt,-1.0,Q[1],1,Q[1],1);
			}
				break;
				
			case MultiRegions::e3DH1D:
			{
				Array<OneD,NekDouble> Wz(m_pressureBCsMaxPts);
				Array<OneD,NekDouble> Dummy1(m_pressureBCsMaxPts);
				Array<OneD,NekDouble> DUmmy2(m_pressureBCsMaxPts);
				
				m_elmt->PhysDeriv(MultiRegionss::DirCartesianMap[0],Vel[1],Vx);
				m_elmt->PhysDeriv(MultiRegionss::DirCartesianMap[1],Vel[0],Uy);
				Vmath::Smul(m_HBCdata[j].m_ptsInElmt,m_wavenumber[j],Vel[2],1,Wz,1);
				
				m_elmt->PhysDeriv(MultiRegionss::DirCartesianMap[1],Vx,Dummy1);
				m_elmt->PhysDeriv(MultiRegionss::DirCartesianMap[1],Uy,Dummy2);
				Vmath::Vsub(m_HBCdata[j].m_ptsInElmt,Dummy1,1,Dummy2,1,Q[0],1);
				Vmath::Smul(m_HBCdata[j].m_ptsInElmt,m_negWavenumberSq[j],Vel[0],1,Dummy1,1);
				m_elmt->PhysDeriv(MultiRegionss::DirCartesianMap[0],Wz,Dummy2);
				Vmath::Vsub(m_HBCdata[j].m_ptsInElmt,Q[0],1,Dummy1,1,Q[0],1);
				Vmath::Vadd(m_HBCdata[j].m_ptsInElmt,Q[0],1,Dummy2,1,Q[0],1);
				
				m_elmt->PhysDeriv(MultiRegionss::DirCartesianMap[1],Wz,Dummy1);
				Vmath::Smul(m_HBCdata[j].m_ptsInElmt,m_negWavenumberSq[j],Vel[1],1,Dummy2,1);
				Vmath::Vsub(m_HBCdata[j].m_ptsInElmt,Dummy1,1,Dummy2,1,Q[1],1);
				m_elmt->PhysDeriv(MultiRegionss::DirCartesianMap[0],Vx,Dummy1);
				m_elmt->PhysDeriv(MultiRegionss::DirCartesianMap[0],Uy,Dummy2);
				Vmath::Vsub(m_HBCdata[j].m_ptsInElmt,Q[1],1,Dummy1,1,Q[1],1);
				Vmath::Vadd(m_HBCdata[j].m_ptsInElmt,Q[1],1,Dummy2,1,Q[1],1);			
			}
				break;
				
			case MultiRegions::e3DH2D:
			{
				m_elmt->PhysDeriv(MultiRegionss::DirCartesianMap[0],Vel[2],Wx);
				m_elmt->PhysDeriv(MultiRegionss::DirCartesianMap[0],Vel[1],Vx);
				
				//m_elmt->PhysDeriv(MultiRegionss::DirCartesianMap[1],Vel[0],Uy);
				Vmath::Smul(m_HBCdata[j].m_ptsInElmt,m_negWavenumberSq[j],Vel[0],1,Uy,1);
				
				//m_elmt->PhysDeriv(MultiRegionss::DirCartesianMap[2],Vel[0],Uz);
				Vmath::Smul(m_HBCdata[j].m_ptsInElmt,m_wavenumber[j],Vel[0],1,Uz,1);
				
				Vmath::Vsub(m_HBCdata[j].m_ptsInElmt,Wz,1,Wx,1,qy,1);
				Vmath::Vsub(m_HBCdata[j].m_ptsInElmt,Vx,1,Uy,1,qz,1);
				
				//m_elmt->PhysDeriv(MultiRegionss::DirCartesianMap[1],qz,Uy);
				Vmath::Smul(m_HBCdata[j].m_ptsInElmt,m_negWavenumberSq[j],qz,1,Uy,1);
				
				//m_elmt->PhysDeriv(MultiRegionss::DirCartesianMap[2],qy,Uz);
				Vmath::Smul(m_HBCdata[j].m_ptsInElmt,m_wavenumber[j],qy,1,Uz,1);
				
				Vmath::Vsub(m_HBCdata[j].m_ptsInElmt,Uy,1,Uz,1,Q[0],1);
			}
				break;
				
			case MultiRegions::e3D:
			{
				Array<OneD,NekDouble> Dummy(m_pressureBCsMaxPts);
				Array<OneD,NekDouble> Vz(m_pressureBCsMaxPts);
				Array<OneD,NekDouble> Uz(m_pressureBCsMaxPts);
				Array<OneD,NekDouble> Wx(m_pressureBCsMaxPts);
				Array<OneD,NekDouble> Wy(m_pressureBCsMaxPts);
				
				m_elmt->PhysDeriv(Vel[0],Dummy,Uy,Uz);
				m_elmt->PhysDeriv(Vel[1],Vx,Dummy,Vz);
				m_elmt->PhysDeriv(Vel[2],Wx,Wy,Dummy);
				
				Vmath::Vsub(m_HBCdata[j].m_ptsInElmt,Wy,1,Vz,1,Q[0],1);
				Vmath::Vsub(m_HBCdata[j].m_ptsInElmt,Uz,1,Wx,1,Q[1],1);
				Vmath::Vsub(m_HBCdata[j].m_ptsInElmt,Vx,1,Uy,1,Q[2],1);
				
				m_elmt->PhysDeriv(Q[0],Dummy,Wy,Vx);
				m_elmt->PhysDeriv(Q[1],Wx,Dummy,Uz);
				m_elmt->PhysDeriv(Q[2],Vz,Uy,Dummy);
				
				Vmath::Vsub(m_HBCdata[j].m_ptsInElmt,Uy,1,Uz,1,Q[0],1);
				Vmath::Vsub(m_HBCdata[j].m_ptsInElmt,Vx,1,Vz,1,Q[1],1);
				Vmath::Vsub(m_HBCdata[j].m_ptsInElmt,Wx,1,Wy,1,Q[2],1);
			}
				break;
		}
	}
	
	
	/** 
	 * Function to roll time-level storages to the next step layout.
	 * The stored data associated with the oldest time-level 
	 * (not required anymore) are moved to the top, where they will
	 * be overwritten as the solution process progresses.
	 */
	void Extrapolate::RollOver(Array<OneD, Array<OneD, NekDouble> > &input)
	{
        int  nlevels = input.num_elements();
        
		Array<OneD, NekDouble> tmp;
		
        tmp = input[nlevels-1];
		
        for(int n = nlevels-1; n > 0; --n)
        {
            input[n] = input[n-1];
        }
		
        input[0] = tmp;
    }
	
	
	/**
	 * Map to directly locate HOPBCs position and offsets in all scenarios
	 */
	void Extrapolate::GenerateHOPBCMap(const MultiRegionss::ExpListSharedPtr &pField)
	{
		
		m_PBndConds   = pField->GetBndConditions();
		m_PBndExp     = pField->GetBndCondExpansions();
		
		// Set up mapping from pressure boundary condition to pressure element details.
		pField->GetBoundaryToElmtMap(m_pressureBCtoElmtID,m_pressureBCtoTraceID);
		
		// find the maximum values of points  for pressure BC evaluation
		m_pressureBCsMaxPts = 0; 
		int cnt, n;
		for(cnt = n = 0; n < m_PBndConds.num_elements(); ++n)
		{
			for(int i = 0; i < m_PBndExp[n]->GetExpSize(); ++i)
			{
				m_pressureBCsMaxPts = max(m_pressureBCsMaxPts, PressureField->GetExp(m_pressureBCtoElmtID[cnt++])->GetTotPoints());
			}
		}
		
		// Storage array for high order pressure BCs
		m_pressureHBCs = Array<OneD, Array<OneD, NekDouble> > (m_intSteps);
		m_acceleration = Array<OneD, Array<OneD, NekDouble> > (m_intSteps + 1);
		
		int HBCnumber = 0;
		for(cnt = n = 0; n < m_PBndConds.num_elements(); ++n)
		{
			// High order boundary condition;
			if(m_PBndConds[n]->GetUserDefined() == SpatialDomains::eHigh)
			{
                cnt += m_PBndExp[n]->GetNcoeffs();
                HBCnumber += m_PBndExp[n]->GetExpSize();
			}
		}
		
		ASSERTL0(HBCnumber > 0 ,"At least one high-order pressure boundary condition is required for scheme consistency");
		
		m_acceleration[0] = Array<OneD, NekDouble>(cnt, 0.0);
		for(n = 0; n < m_intSteps; ++n)
		{
			m_pressureHBCs[n]   = Array<OneD, NekDouble>(cnt, 0.0);
			m_acceleration[n+1] = Array<OneD, NekDouble>(cnt, 0.0);
		}
		
		switch(pField->GetExpType())
		{
			case MultiRegions::e2D:
			{
				m_curl_dim = 2;
				m_bnd_dim  = 2;
			}
				break;
			case MultiRegions::e3DH1D:
			{
				m_curl_dim = 3;
				m_bnd_dim  = 2;
			}
				break;
			case MultiRegions::e3DH2D:
			{
				m_curl_dim = 3;
				m_bnd_dim  = 1;
			}
				break;
			case MultiRegions::e3D:
			{
				m_curl_dim = 3;
				m_bnd_dim  = 3;
			}
				break;
		}
		
		
		m_HBCdata = Array<OneD, HBCInfo>(HOPBCnumber);
		
        switch(pField->GetExpType())
        {
			case MultiRegions::e2D:
			case MultiRegions::e3D:
			{
				int coeff_count = 0;
				int exp_size;
				int j=0;
				int cnt = 0;
				for(int n = 0 ; n < m_PBndConds.num_elements(); ++n)
				{
					exp_size = m_PBndExp[n]->GetExpSize();
                    
					if(m_PBndConds[n]->GetUserDefined() == SpatialDomains::eHigh)
					{
						for(int i = 0; i < exp_size; ++i,cnt++)
						{
							m_HBCdata[j].m_globalElmtID = m_pressureBCtoElmtID[cnt];   
							m_elmt      = pField->GetExp(m_HBCdata[j].m_globalElmtID);
							m_HBCdata[j].m_ptsInElmt = m_elmt->GetTotPoints();         
							m_HBCdata[j].m_physOffset = pField->GetPhys_Offset(m_HBCdata[j].m_globalElmtID);
							m_HBCdata[j].m_bndElmtOffset = i;       
							m_HBCdata[j].m_elmtTraceID = m_pressureBCtoTraceID[cnt];      
							m_HBCdata[j].m_bndryElmtID = n;
							m_HBCdata[j].m_coeffOffset = coeff_count;
							coeff_count += m_elmt->GetEdgeNcoeffs(m_HBCdata[j].m_elmtTraceID);
							j = j+1;
						}
					}
					else // setting if just standard BC no High order
					{
						cnt += exp_size;
					}
				}
			}
				break;
				
			case MultiRegions::e3DH1D:
			{
				Array<OneD, unsigned int> planes;
				planes = pField->GetZIDs();
				int num_planes = planes.num_elements();            
				int num_elm_per_plane = (pField->GetExpSize())/num_planes;
				
				m_wavenumber      = Array<OneD, NekDouble>(HOPBCnumber);
				m_negWavenumberSq = Array<OneD, NekDouble>(HOPBCnumber);
				
				int coeff_count = 0;
				int exp_size, exp_size_per_plane;
				int j=0;
				int K;
				NekDouble sign = -1.0;
				int cnt = 0;
				
				for(int k = 0; k < num_planes; k++)
				{
					K = planes[k]/2;
					for(int n = 0 ; n < m_PBndConds.num_elements(); ++n)
					{
						exp_size = m_PBndExp[n]->GetExpSize();
						exp_size_per_plane = exp_size/num_planes;
						
						if(m_PBndConds[n]->GetUserDefined() == SpatialDomains::eHigh)
						{
							for(int i = 0; i < exp_size_per_plane; ++i,cnt++)
							{
								m_HBCdata[j].m_globalElmtID = m_pressureBCtoElmtID[cnt];   
								m_elmt      = pField->GetExp(m_HBCdata[j].m_globalElmtID);
								m_HBCdata[j].m_ptsInElmt = m_elmt->GetTotPoints();         
								m_HBCdata[j].m_physOffset = pField->GetPhys_Offset(m_HBCdata[j].m_globalElmtID);
								m_HBCdata[j].m_bndElmtOffset = i+k*exp_size_per_plane;       
								m_HBCdata[j].m_elmtTraceID = m_pressureBCtoTraceID[cnt];      
								m_HBCdata[j].m_bndryElmtID = n;
								m_HBCdata[j].m_coeffOffset = coeff_count;
								coeff_count += m_elmt->GetEdgeNcoeffs(m_HBCdata[j].m_elmtTraceID);
								
								if(m_SingleMode)
								{
									m_wavenumber[j]      = -2*M_PI/m_LhomZ;       
									m_negWavenumberSq[j] = -1.0*m_wavenumber[j]*m_wavenumber[j];
								}
								else if(m_HalfMode || m_MultipleModes)
								{
									m_wavenumber[j]      = 2*M_PI/m_LhomZ;       
									m_negWavenumberSq[j] = -1.0*m_wavenumber[j]*m_wavenumber[j];
								}
								else
								{
									m_wavenumber[j]     = 2*M_PI*sign*(NekDouble(K))/m_LhomZ; 
									m_negWavenumberSq[j] = -1.0*m_wavenumber[j]*m_wavenumber[j];
								}
								
								int assElmtID;
								
								if(k%2==0)
								{
									if(m_HalfMode)
									{
										assElmtID = m_HBCdata[j].m_globalElmtID;
										
									}
									else
									{
										assElmtID = m_HBCdata[j].m_globalElmtID + num_elm_per_plane;
									}
								}
								else 
								{
									assElmtID = m_HBCdata[j].m_globalElmtID - num_elm_per_plane;
								}
								
								m_HBCdata[j].m_assPhysOffset = pField->GetPhys_Offset(assElmtID);
								
								j = j+1;
							}
						}
						else // setting if just standard BC no High order
						{
							cnt += exp_size_per_plane;
						}
					}
					sign = -1.0*sign;
				}
			}
				break;
				
			case MultiRegions::e3DH2D:
			{
				int cnt = 0;
				int exp_size, exp_size_per_line;
				int j=0;
				
				for(int k1 = 0; k1 < m_npointsZ; k1++)
				{
					for(int k2 = 0; k2 < m_npointsY; k2++)
					{
						for(int n = 0 ; n < m_PBndConds.num_elements(); ++n)
						{
							exp_size = m_PBndExp[n]->GetExpSize();
							
							exp_size_per_line = exp_size/(m_npointsZ*m_npointsY);
							
							if(m_PBndConds[n]->GetUserDefined() == SpatialDomains::eHigh)
							{
								for(int i = 0; i < exp_size_per_line; ++i,cnt++)
								{
									// find element and edge of this expansion. 
									// calculate curl x curl v;
									m_HBCdata[j].m_globalElmtID = m_pressureBCtoElmtID[cnt];
									m_elmt      = pField->GetExp(m_HBCdata[j].m_globalElmtID);
									m_HBCdata[j].m_ptsInElmt = m_elmt->GetTotPoints();
									m_HBCdata[j].m_physOffset = pField->GetPhys_Offset(m_HBCdata[j].m_globalElmtID);
									m_HBCdata[j].m_bndElmtOffset = i+(k1*m_npointsY+k2)*exp_size_per_line;
									m_HBCdata[j].m_elmtTraceID = m_pressureBCtoTraceID[cnt];                
									m_HBCdata[j].m_bndryElmtID = n;
									m_wavenumber[j] = 2*M_PI*sign*(NekDouble(k1))/m_LhomZ;
									m_negWavenumberSq[j] = 2*M_PI*sign*(NekDouble(k2))/m_LhomY;
								}
							}
							else
							{
								cnt += exp_size_per_line;
							}
						}
					}
				}
			}
				break;
		}
    }
}

