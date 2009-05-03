///////////////////////////////////////////////////////////////////////////////
//
// File AdvectionDiffusionReaction.cpp
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
// Description: Advection Diffusion Reaction class definition built on
// ADRBase class
//
///////////////////////////////////////////////////////////////////////////////

#include <AdvectionDiffusionReactionSolver/AdvectionDiffusionReaction.h>
#include <cstdio>
#include <cstdlib>
namespace Nektar
{
    /**
     * Basic construnctor
     */
    AdvectionDiffusionReaction::AdvectionDiffusionReaction(void):
        ADRBase(),
        m_infosteps(100),
        m_explicitAdvection(true),
        m_explicitDiffusion(true),
        m_explicitReaction(true)
    {     
    }
    
    /**
     * Constructor. Creates ...
     *
     * \param 
     * \param
     */
    AdvectionDiffusionReaction::AdvectionDiffusionReaction(string &fileNameString):
        ADRBase(fileNameString,true),
        m_infosteps(10),
        m_explicitDiffusion(true),
        m_explicitReaction(true)
    {

        int i;

        // Set up equation type enum using kEquationTypeStr
        std::string typeStr = m_boundaryConditions->GetSolverInfo("EQTYPE");

        for(i = 0; i < (int) eEquationTypeSize; ++i)
        {
            if(NoCaseStringCompare(kEquationTypeStr[i],typeStr) == 0 )
            {
                m_equationType = (EquationType)i; 
                break;
            }
        }

        ASSERTL0(i != (int) eEquationTypeSize, "Invalid expansion type.");
        
        // Equation specific Setups 
        switch(m_equationType)
        {
        case eHelmholtz: case eLaplace: case ePoisson:
            break;
        case eSteadyAdvection:
            m_velocity = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
        
            for(int i = 0; i < m_spacedim; ++i)
            {
                m_velocity[i] = Array<OneD, NekDouble> (GetNpoints());
            }
            
            EvaluateAdvectionVelocity();
            break;
			
        case eUnsteadyAdvection:
            m_velocity = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
        
            for(int i = 0; i < m_spacedim; ++i)
            {
                m_velocity[i] = Array<OneD, NekDouble> (GetNpoints());
            }
            
            EvaluateAdvectionVelocity();
            m_timeIntMethod = LibUtilities::eClassicalRungeKutta4;		
            goto UnsteadySetup;
            break;
        case eUnsteadyDiffusion:
            // default explicit scheme 
            m_timeIntMethod = LibUtilities::eClassicalRungeKutta4; 
            goto UnsteadySetup;
        case eUnsteadyDiffusionReaction:
            m_timeIntMethod = LibUtilities::eIMEXdirk_3_4_3;	
            goto UnsteadySetup;
            break;

        UnsteadySetup:
            {
                std::string Implicit = "Implicit"; 
                if(m_boundaryConditions->CheckForParameter("IO_InfoSteps") == true)
                {
                    m_infosteps =  m_boundaryConditions->GetParameter("IO_InfoSteps");
                }
                
                // check that any user defined boundary condition is indeed implemented
                for(int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
                {	
                    // Time Dependent Boundary Condition (if no use defined then this is empty)
                    if (m_fields[0]->GetBndConditions()[n]->GetUserDefined().GetEquation() != "")
                    {
                        if (m_fields[0]->GetBndConditions()[n]->GetUserDefined().GetEquation() != "TimeDependent")
                        {
                            ASSERTL0(false,"Unknown USERDEFINEDTYPE boundary condition");
                        }
                    }
                }
                
                // Check for definition of Implicit/Explicit terms in solverinfo
                if(m_boundaryConditions->SolverInfoExists("ADVECTIONADVANCEMENT"))
                {
                    std::string AdvStr = m_boundaryConditions->GetSolverInfo("ADVECTIONADVANCEMENT");
                    
                    if(NoCaseStringCompare(AdvStr,Implicit) == 0)
                    {
                        m_explicitAdvection = false;
                    }
                    else
                    {
                        m_explicitAdvection = true;
                    }
                }
                else
                {
                    m_explicitAdvection = true;
                }
                
                
                if(m_boundaryConditions->SolverInfoExists("DIFFUSIONADVANCEMENT"))
                {
                    std::string AdvStr = m_boundaryConditions->GetSolverInfo("DIFFUSIONADVANCEMENT");
                    
                    if(NoCaseStringCompare(AdvStr,Implicit) == 0 )
                    {
                        m_explicitDiffusion = false;
                        // Reset default for implicit diffusion
                        if(m_equationType == eUnsteadyDiffusion)
                        {
                            m_timeIntMethod = LibUtilities::eDIRKOrder3;		
                        }
                    }
                    else
                    {
                        m_explicitDiffusion = true;
                    }
                }
                else
                {
                    m_explicitDiffusion = true;
                }
                
                if(m_boundaryConditions->SolverInfoExists("REACTIONADVANCEMENT"))
                {
                    std::string AdvStr = m_boundaryConditions->GetSolverInfo("REACTIONADVANCEMENT");
                    
                    if(NoCaseStringCompare(AdvStr,Implicit) == 0)
                    {
                        m_explicitReaction = false;
                    }
                    else
                    {
                        m_explicitReaction = true;
                    }
                }
                else
                {
                    m_explicitReaction = true;
                }

                // check to see if time stepping has been reset
                if(m_boundaryConditions->SolverInfoExists("TIMEINTEGRATIONMETHOD"))
                {
                    std::string TimeIntStr = m_boundaryConditions->GetSolverInfo("TIMEINTEGRATIONMETHOD");
                    int i;
                    for(i = 0; i < (int) LibUtilities::SIZE_TimeIntegrationMethod; ++i)
                    {
                        if(NoCaseStringCompare(LibUtilities::TimeIntegrationMethodMap[i],TimeIntStr) == 0 )
                        {
                            m_timeIntMethod = (LibUtilities::TimeIntegrationMethod)i; 
                            break;
                        }
                    }

                    ASSERTL0(i != (int) LibUtilities::SIZE_TimeIntegrationMethod, "Invalid time integration type.");
                }

                break;
            }
        case eNoEquationType:
        default:
            ASSERTL0(false,"Unknown or undefined equation type");
        }
    }

    void AdvectionDiffusionReaction::EvaluateAdvectionVelocity()
    {
        int nq = m_fields[0]->GetNpoints();
        
        std::string velStr[3] = {"Vx","Vy","Vz"};

        Array<OneD,NekDouble> x0(nq);
        Array<OneD,NekDouble> x1(nq);
        Array<OneD,NekDouble> x2(nq);
      
        // get the coordinates (assuming all fields have the same
        // discretisation)
        m_fields[0]->GetCoords(x0,x1,x2);

        for(int i = 0 ; i < m_velocity.num_elements(); i++)
	{
            SpatialDomains::ConstUserDefinedEqnShPtr ifunc = m_boundaryConditions->GetUserDefinedEqn(velStr[i]);
            
            for(int j = 0; j < nq; j++)
	    {
                m_velocity[i][j] = ifunc->Evaluate(x0[j],x1[j],x2[j]);
	    }
	}
    }
    
    void AdvectionDiffusionReaction::ODErhs(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
                                                  Array<OneD,       Array<OneD, NekDouble> >&outarray, 
                                            const NekDouble time) 
    {
        int i;
        int nvariables = inarray.num_elements();
        int ncoeffs    = inarray[0].num_elements();

	
        switch(m_projectionType)
        {
        case eDiscontinuousGalerkin:
            {	  
                switch(m_equationType)
                {
#ifdef OLDENUM
                case eAdvection:
#else
                case eUnsteadyAdvection:
#endif
                    {
                        SetBoundaryConditions(time);
                        WeakDGAdvection(inarray, outarray);
                        for(i = 0; i < nvariables; ++i)
                        {
                            m_fields[i]->MultiplyByElmtInvMass(outarray[i],outarray[i]);
                            Vmath::Neg(ncoeffs,outarray[i],1);
                        }
                    }
                    break;
                    
#ifdef OLDENUM
                case eDiffusion:
#else
                case eUnsteadyDiffusion:
#endif
                    {
			// BoundaryConditions are imposed weakly at the Diffusion operator
                        WeakDGDiffusion(inarray,outarray);
                        for(i = 0; i < nvariables; ++i)
                        {
                            m_fields[i]->MultiplyByElmtInvMass(outarray[i],outarray[i]);
                        }						
                    }
                    break;
                }
                break;
            }
        case eGalerkin:
            {
	  	SetBoundaryConditions(time);
                Array<OneD, NekDouble> physfield(GetNpoints());
		
                for(i = 0; i < nvariables; ++i)
                {
                    m_fields[i]->MultiplyByInvMassMatrix(inarray[i],  
                                                         outarray[i],
                                                         false);
                    // Calculate -(\phi, V\cdot Grad(u))
                    m_fields[i]->BwdTrans_IterPerExp(outarray[i],physfield);
                    
                    WeakAdvectionNonConservativeForm(m_velocity,
                                                     physfield, outarray[i]);
                    
                    Vmath::Neg(ncoeffs,outarray[i],1);		   		    
                }
            }
            break;
        default:
            ASSERTL0(false,"Unknown projection scheme");
            break;
        }
    }
    
	void AdvectionDiffusionReaction::ODEeReaction(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
												  Array<OneD, Array<OneD, NekDouble> >&outarray, 
											const NekDouble time)
	
	{
		int i,k;
                int nvariables = inarray.num_elements();
                int ncoeffs    = inarray[0].num_elements();
		const NekDouble coeff = 3.14159265*3.14159265;
					
		for (i = 0; i < nvariables; ++i)
		{
			Vmath::Smul(ncoeffs, coeff, inarray[i], 1, outarray[i], 1);
		}
	}
	
	
    void AdvectionDiffusionReaction::ODElhs(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
						  Array<OneD,       Array<OneD, NekDouble> >&outarray, 
                                            const NekDouble time) 
    {
        int nvariables = inarray.num_elements();
        MultiRegions::GlobalMatrixKey key(StdRegions::eMass);

        for(int i = 0; i < nvariables; ++i)
        {
            m_fields[i]->MultiRegions::ExpList::GeneralMatrixOp(key,inarray[i],outarray[i]);
        }
    }
    
    void AdvectionDiffusionReaction::ODElhsSolve(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
                                                       Array<OneD,       Array<OneD, NekDouble> >&outarray, 
                                                 const NekDouble time)   
    {
	SetBoundaryConditions(time);
        int i;
        int nvariables = inarray.num_elements();
	
        switch(m_projectionType)
        {
        case eDiscontinuousGalerkin:

            for(i = 0; i < nvariables; ++i)
            {
                m_fields[i]->MultiplyByElmtInvMass(inarray[i],outarray[i]);
            }
	  break;
        case eGalerkin:
	  {
              for(i = 0; i < nvariables; ++i)
              {
                  m_fields[i]->MultiplyByInvMassMatrix(inarray[i],outarray[i],false);
              }
          }
          break;
        default:
            ASSERTL0(false,"Unknown projection scheme");
            break;
        }
    }


	void AdvectionDiffusionReaction::ODEhelmSolve(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
												  Array<OneD, Array<OneD, NekDouble> >&outarray,
												  NekDouble time, 
												  NekDouble lambda)
    {
        int nvariables = inarray.num_elements();
        int ncoeffs    = inarray[0].num_elements();
							
        // We solve ( \nabla^2 - HHlambda ) Y[i] = rhs [i]
        // inarray = input: \hat{rhs} -> output: \hat{Y}    
        // outarray = output: nabla^2 \hat{Y}       
        // where \hat = modal coeffs
	
        MultiRegions::GlobalMatrixKey key(StdRegions::eMass);
	
        for (int i = 0; i < nvariables; ++i)
        {
            // Multiply by inverse of mass matrix
            m_fields[i]->MultiplyByInvMassMatrix(inarray[i],outarray[i],false);
            
            // Multiply rhs[i] with -1.0/gamma/timestep
            Vmath::Smul(ncoeffs, -1.0/lambda, outarray[i], 1, outarray[i], 1);
            
            // Update coeffs to m_fields
            m_fields[i]->UpdateCoeffs() = outarray[i];
            
            // Backward Transformation to nodal coefficients
            m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(), m_fields[i]->UpdatePhys());
            
            NekDouble kappa = 1.0/lambda;
            
            // Solve a system of equations with Helmholtz solver
            m_fields[i]->HelmSolve(m_fields[i]->GetPhys(),m_fields[i]->UpdateCoeffs(),kappa);
            
            // The solution is Y[i]
            outarray[i] = m_fields[i]->GetCoeffs();	  
            
            // Multiply back by mass matrix
            m_fields[i]->MultiRegions::ExpList::GeneralMatrixOp(key,outarray[i],outarray[i]);
        }
    }
	

    void AdvectionDiffusionReaction::SolveHelmholtz(NekDouble lambda)
    {
        for(int i = 0; i < m_fields.num_elements(); ++i)
        {
            m_fields[i]->HelmSolve(m_fields[i]->GetPhys(),m_fields[i]->UpdateCoeffs(),lambda);
            m_fields[i]->SetPhysState(false);
        }
    }

    // For Continuous Galerkin projections with time-dependent dirichlet boundary conditions,
    // the time integration can be done as follows:
    // The ODE resulting from the PDE can be formulated as:
    // 
    // M du/dt = F(u)  or du/dt = M^(-1) F(u)
    //
    // Now suppose that M does not depend of time, the ODE can than be written as:
    //
    // d(Mu)/dt = F(u)
    //
    // Introducing the variable u* = Mu, this yields
    //
    // du*/dt = F( M^(-1) u* ) = F*(u*)
    //
    // So rather than solving the initial ODE, it is advised to solve this new ODE for u*
    // as this allows for an easier treatment of the dirichlet boundary conditions.
    // However, note that at the end of every time step, the actual solution u can
    // be calculated as:
    // 
    // u = M^(-1) u*;
    //
    // This can be viewed as projecting the solution u* onto the known boundary conditions.
    // Note that this step is also done inside the ODE rhs function F*.
    //
    // In order for all of this to work appropriately, make sure that the operator M^(-1)
    // does include the enforcment of the dirichlet boundary conditionst

    void AdvectionDiffusionReaction::GeneralTimeIntegration(int nsteps, 
	                                                     LibUtilities::TimeIntegrationMethod IntMethod,
							     LibUtilities::TimeIntegrationSchemeOperators ode)
    {
        int i,n,nchk = 0;
        int ncoeffs = m_fields[0]->GetNcoeffs();
        int nvariables = m_fields.num_elements();

        // Set up wrapper to fields data storage. 
        Array<OneD, Array<OneD, NekDouble> >   fields(nvariables);
        Array<OneD, Array<OneD, NekDouble> >   tmp(nvariables);
        
        for(i = 0; i < nvariables; ++i)
        {
            fields[i]  = m_fields[i]->UpdateCoeffs();
        }

        if(m_projectionType==eGalerkin)
        {
            // calculate the variable u* = Mu
            // we are going to TimeIntegrate this new variable u*
            MultiRegions::GlobalMatrixKey key(StdRegions::eMass);
            for(int i = 0; i < nvariables; ++i)
            {
                tmp[i] = Array<OneD, NekDouble>(ncoeffs);
                m_fields[i]->MultiRegions::ExpList::GeneralMatrixOp(key,fields[i],fields[i]);
            }

            for(int i = 0; i < nvariables; ++i)
            {
                m_fields[i]->SetPhysState(false);
            }
        }

        // Declare an array of TimeIntegrationSchemes
        // For multi-stage methods, this array will have just one entry containing
        // the actual multi-stage method...
        // For multi-steps method, this can have multiple entries
        //  - the first scheme will used for the first timestep (this is an initialization scheme)
        //  - the second scheme will used for the first timestep (this is an initialization scheme)
        //  - ...
        //  - the last scheme will be used for all other time-steps (this will be the actual scheme)
        Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr> IntScheme;
        LibUtilities::TimeIntegrationSolutionSharedPtr u;
        int numMultiSteps;

        switch(IntMethod)
        {
	case LibUtilities::eIMEXdirk_3_4_3:
	case LibUtilities::eDIRKOrder3:
        case LibUtilities::eBackwardEuler:      
        case LibUtilities::eForwardEuler:      
        case LibUtilities::eClassicalRungeKutta4:
            {
                numMultiSteps = 1;

                IntScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr>(numMultiSteps);

                LibUtilities::TimeIntegrationSchemeKey IntKey(IntMethod);
                IntScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey];

                u = IntScheme[0]->InitializeScheme(m_timestep,fields,m_time,ode);
            }
            break;
        case LibUtilities::eAdamsBashforthOrder2:
            {
                numMultiSteps = 2;

                IntScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr>(numMultiSteps);

                // Used in the first time step to initalize the scheme
                LibUtilities::TimeIntegrationSchemeKey IntKey0(LibUtilities::eForwardEuler);
				
                // Used for all other time steps 
                LibUtilities::TimeIntegrationSchemeKey IntKey1(IntMethod); 
                IntScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey0];
                IntScheme[1] = LibUtilities::TimeIntegrationSchemeManager()[IntKey1];

                // Initialise the scheme for the actual time integration scheme
                u = IntScheme[1]->InitializeScheme(m_timestep,fields,m_time,ode);
            }
            break;
        default:
            {
                ASSERTL0(false,"populate switch statement for integration scheme");
            }
        }
					          
        for(n = 0; n < nsteps; ++n)
        {
            //----------------------------------------------
            // Perform time step integration
            //----------------------------------------------
            if( n < numMultiSteps-1)
            {
                // Use initialisation schemes
                fields = IntScheme[n]->TimeIntegrate(m_timestep,u,ode);
            }
            else
            {
                fields = IntScheme[numMultiSteps-1]->TimeIntegrate(m_timestep,u,ode);
            }

            m_time += m_timestep;

            if(m_projectionType==eGalerkin)
            {
                // Project the solution u* onto the boundary conditions to
                // obtain the actual solution
                SetBoundaryConditions(m_time);
                for(i = 0; i < nvariables; ++i)
                {
                    m_fields[i]->MultiplyByInvMassMatrix(fields[i],tmp[i],false);
                    fields[i] = tmp[i];	   		    
                }
            }

            //----------------------------------------------
            // Dump analyser information
            //----------------------------------------------
            if(!((n+1)%m_infosteps))
            {
	      cout << "Steps: " << n+1 << "\t Time: " << m_time << endl;
            }
            
            if(n&&(!((n+1)%m_checksteps)))
            {
	      for(i = 0; i < nvariables; ++i)
		{
		  (m_fields[i]->UpdateCoeffs()) = fields[i];
		}
	      Checkpoint_Output(nchk++);
            }
        }
        
        for(i = 0; i < nvariables; ++i)
        {
	  (m_fields[i]->UpdateCoeffs()) = fields[i];
        }
    }
	
    
  //----------------------------------------------------
  void AdvectionDiffusionReaction::SetBoundaryConditions(NekDouble time)
  {
    int nvariables = m_fields.num_elements();
    for (int i = 0; i < nvariables; ++i)
    {
        m_fields[i]->EvaluateBoundaryConditions(time);
    }
  }
  
  // Evaulate flux = m_fields*ivel for i th component of Vu 
  void AdvectionDiffusionReaction::GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield, 
						 Array<OneD, Array<OneD, NekDouble> > &flux)
  {
        ASSERTL1(flux.num_elements() == m_velocity.num_elements(),"Dimension of flux array and velocity array do not match");

        for(int j = 0; j < flux.num_elements(); ++j)
        {
            Vmath::Vmul(GetNpoints(),physfield[i],1,
                        m_velocity[j],1,flux[j],1);
        }
    }
	
	
 // Evaulate flux = m_fields*ivel for i th component of Vu for direction j
  void AdvectionDiffusionReaction::GetFluxVector(const int i, const int j, Array<OneD, Array<OneD, NekDouble> > &physfield, 
						 Array<OneD, Array<OneD, NekDouble> > &flux)
  {
        ASSERTL1(flux.num_elements() == m_velocity.num_elements(),"Dimension of flux array and velocity array do not match");

        for(int k = 0; k < flux.num_elements(); ++k)
        {
		Vmath::Zero(GetNpoints(),flux[k],1);
        }
		
		Vmath::Vcopy(GetNpoints(),physfield[i],1,flux[j],1);
    }
	

    void AdvectionDiffusionReaction::NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, 
						   Array<OneD, Array<OneD, NekDouble> > &numflux)
    {
        int i;

        int nTraceNumPoints = GetTraceNpoints();
        int nvel = m_velocity.num_elements();

        Array<OneD, NekDouble > Fwd(nTraceNumPoints);
        Array<OneD, NekDouble > Bwd(nTraceNumPoints);
        Array<OneD, NekDouble > Vn (nTraceNumPoints,0.0);

	// Get Edge Velocity - Could be stored if time independent
        for(i = 0; i < nvel; ++i)
	  {
            m_fields[0]->ExtractTracePhys(m_velocity[i], Fwd);
            Vmath::Vvtvp(nTraceNumPoints,m_traceNormals[i],1,Fwd,1,Vn,1,Vn,1);
	  }

        for(i = 0; i < numflux.num_elements(); ++i)
        {
            m_fields[i]->GetFwdBwdTracePhys(physfield[i],Fwd,Bwd);
            //evaulate upwinded m_fields[i]
	    m_fields[i]->GetTrace()->Upwind(Vn,Fwd,Bwd,numflux[i]);
	    // calculate m_fields[i]*Vn
            Vmath::Vmul(nTraceNumPoints,numflux[i],1,Vn,1,numflux[i],1);
        }
    }

  void AdvectionDiffusionReaction::NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, 
						 Array<OneD, Array<OneD, NekDouble> > &numfluxX, 
						 Array<OneD, Array<OneD, NekDouble> > &numfluxY)
    {
        int i;

        int nTraceNumPoints = GetTraceNpoints();
        int nvel = m_velocity.num_elements();

        Array<OneD, NekDouble > Fwd(nTraceNumPoints);
        Array<OneD, NekDouble > Bwd(nTraceNumPoints);
	Array<OneD, NekDouble > tmp(nTraceNumPoints,0.0);

	Array<OneD, Array<OneD, NekDouble > > traceVelocity(2);

	traceVelocity[0] = Array<OneD,NekDouble>(nTraceNumPoints,0.0);
	traceVelocity[1] = Array<OneD,NekDouble>(nTraceNumPoints,0.0);

	// Get Edge Velocity - Could be stored if time independent
	m_fields[0]->ExtractTracePhys(m_velocity[0], traceVelocity[0]);
	m_fields[0]->ExtractTracePhys(m_velocity[1], traceVelocity[1]);

	m_fields[0]->GetFwdBwdTracePhys(physfield[0],Fwd,Bwd);
	
	m_fields[0]->GetTrace()->Upwind(traceVelocity,Fwd,Bwd,tmp);
	
	Vmath::Vmul(nTraceNumPoints,tmp,1,traceVelocity[0],1,numfluxX[0],1);
	Vmath::Vmul(nTraceNumPoints,tmp,1,traceVelocity[1],1,numfluxY[0],1);
  
    }
	
	
// Compute the fluxes of q and u vector fields for discontinuous diffusion term
// Input:   ufield : 1 by # of total trace points
// Output:  uflux  : 2 by # of total trace points
	
  void AdvectionDiffusionReaction::NumFluxforDiff(Array<OneD, Array<OneD, NekDouble> > &ufield, 
						   Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &uflux)
    {
        int i,j;

        int nTraceNumPoints = GetTraceNpoints();
        int nvel = m_velocity.num_elements();
		
        Array<OneD, NekDouble > Fwd(nTraceNumPoints);
        Array<OneD, NekDouble > Bwd(nTraceNumPoints);
	Array<OneD, NekDouble > Vn (nTraceNumPoints,0.0);
	Array<OneD, NekDouble > fluxtemp (nTraceNumPoints);
	  		  
     // Get the sign of (v \cdot n), v = an arbitrary vector
        for(i = 0; i < nvel; ++i)
	    {
            m_fields[0]->ExtractTracePhys(m_velocity[i], Fwd);
            Vmath::Vvtvp(nTraceNumPoints,m_traceNormals[i],1,Fwd,1,Vn,1,Vn,1);
	    }

        for(i = 0; i < uflux.num_elements(); ++i)
        {
		
		   //  Evaulate upwind flux of uflux = \hat{u} \phi \cdot u = u^{(+,-)} n

		   //  Fwd = u+,  Bwd = u-
               m_fields[i]->GetFwdBwdTracePhys(ufield[i],Fwd,Bwd);

		   //  uflux[i] = minmod (u+, u-)   
		       m_fields[i]->GetTrace()->Upwind(Vn,Fwd,Bwd,fluxtemp,1);  
			   
		   // Imposing weak boundary condition with flux
			   m_fields[i]->ExtractTracePhys(ufield[i],Fwd);
		       WeakPenaltyBoundary(Fwd,fluxtemp);
			   
		       for(j = 0; j < nvel; ++j)
	            {
                   Vmath::Vmul(nTraceNumPoints,m_traceNormals[j],1,fluxtemp,1,uflux[i][j],1);
	            }
        }
    }

			
// Compute the fluxes of q and u vector fields for discontinuous diffusion term
// Input:   qfield : 2 by # of total trace points
// Output:  qflux  : 2 by # of total trace points
	
  void AdvectionDiffusionReaction::NumFluxforDiff(Array<OneD, Array<OneD, NekDouble> > &ufield,
	                                          Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &qfield,
						   Array<OneD, Array<OneD, NekDouble> >  &qflux)
    {
        int nTraceNumPoints = GetTraceNpoints();
        int nvel = m_velocity.num_elements();
	NekDouble C11 = 1.0;
		
	int numBDEdge = m_fields[0]->GetBndCondExpansions()[0]->GetExpSize();
	int Nfps = m_fields[0]->GetBndCondExpansions()[0]->GetExp(0)->GetNumPoints(0) ;
	int id1, id2;
		
        Array<OneD, NekDouble > Fwd(nTraceNumPoints);
        Array<OneD, NekDouble > Bwd(nTraceNumPoints);
	Array<OneD, NekDouble > Vn (nTraceNumPoints,0.0);	
			
	Array<OneD, NekDouble > qFwd(nTraceNumPoints);
        Array<OneD, NekDouble > qBwd(nTraceNumPoints);
	Array<OneD, NekDouble > qfluxtemp(nTraceNumPoints,0.0);
	Array<OneD, NekDouble > uterm(nTraceNumPoints);
	  		  
     // Get the sign of (v \cdot n), v = an arbitrary vector
        for(int j = 0; j < nvel; ++j)
	    {
            m_fields[0]->ExtractTracePhys(m_velocity[j], qFwd);
            Vmath::Vvtvp(nTraceNumPoints,m_traceNormals[j],1,qFwd,1,Vn,1,Vn,1);
	    }

        for(int i = 0; i < qfield.num_elements(); ++i)
        {
		   // Evaulate upwind flux of qflux = \hat{q} \cdot u = q \cdot n - C_(11)*(u^+ - u^-)		   			
			for(int j = 0; j < nvel; ++j)
			{
			    //  qFwd = q+,  qBwd = q-
			      m_fields[i]->GetFwdBwdTracePhys(qfield[i][j],qFwd,qBwd);

                            //  qfluxtemp = maxmod(q-,q+)
	                      m_fields[i]->GetTrace()->Upwind(Vn,qFwd,qBwd,qfluxtemp,-1);
				
		           //   uterm = - C11 ( u+ - u- )
			      m_fields[i]->GetFwdBwdTracePhys(ufield[i],Fwd,Bwd);
			      Vmath::Vsub(nTraceNumPoints,Fwd,1,Bwd,1,uterm,1);					  
		              Vmath::Smul(nTraceNumPoints,-1.0*C11,uterm,1,uterm,1);
				
		           //   qflux[i][j] = (qhat_x, qhat_y) + temp*(n_x,n_y)  
			      Vmath::Vvtvp(nTraceNumPoints,uterm,1,m_traceNormals[j],1,qfluxtemp,1,qfluxtemp,1);
				
		          //   Imposing weak boundary condition with flux
		              m_fields[i]->ExtractTracePhys(qfield[i][j],qFwd);
			      m_fields[i]->ExtractTracePhys(ufield[i],Fwd);
			      WeakPenaltyBoundary(j,Fwd,qFwd,qfluxtemp,C11);

		          //   Calculate m_fields[i]*Vn
		              Vmath::Vvtvp(nTraceNumPoints,m_traceNormals[j],1,qfluxtemp,1,qflux[i],1,qflux[i],1);	   	
			}			  
        }
   }
	
// Diffusion: Imposing weak boundary condition for u with flux 
//  uflux = g_D  on Dirichlet boundary condition
//  uflux = u_Fwd  on Neumann boundary condition
 void AdvectionDiffusionReaction::WeakPenaltyBoundary(const Array<OneD, const NekDouble> &Fwd, 
									     Array<OneD, NekDouble> &penaltyflux,
									     NekDouble initialtime)
	       {
                     int i, j, e, npoints, id1, id2;
                     int nbnd = m_fields[0]->GetBndCondExpansions().num_elements();
		     int numBDEdge = m_fields[0]->GetBndCondExpansions()[0]->GetExpSize();
	             int Nfps = m_fields[0]->GetBndCondExpansions()[0]->GetExp(0)->GetNumPoints(0) ;
            
                  for(i = 0; i < nbnd; ++i)
                  {                 
			
			  // Evaluate boundary values g_D or g_N from input files
			    SpatialDomains::ConstInitialConditionShPtr ifunc = m_boundaryConditions->GetInitialCondition(i);
                            npoints = m_fields[0]->GetBndCondExpansions()[i]->GetNpoints();

			    Array<OneD,NekDouble> BDphysics(npoints);
			    Array<OneD,NekDouble> x0(npoints,0.0);
                            Array<OneD,NekDouble> x1(npoints,0.0);
                            Array<OneD,NekDouble> x2(npoints,0.0);  
                
                            m_fields[0]->GetBndCondExpansions()[i]->GetCoords(x0,x1,x2);
			   for(j = 0; j < npoints; j++)
                             {
				 BDphysics[j] = ifunc->Evaluate(x0[j],x1[j],x2[j],initialtime);
                             }
			
			       // Weakly impose boundary conditions by modifying flux values
				for (e = 0; e < numBDEdge ; ++e)
				 {
					  id1 = m_fields[i]->GetBndCondExpansions()[0]->GetPhys_Offset(e);
					  id2 = m_fields[i]->GetTrace()->GetPhys_Offset(m_fields[i]->GetTraceMap()->GetBndCondCoeffsToGlobalCoeffsMap(e));
					
					// For Dirichlet boundary condition: uflux = g_D
					if(m_fields[0]->GetBndConditions()[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
					 {
					   Vmath::Vcopy(Nfps,&BDphysics[id1],1,&penaltyflux[id2],1);
					 }
					
					// For Neumann boundary condition: uflux = u+
					else if((m_fields[0]->GetBndConditions()[i])->GetBoundaryConditionType() == SpatialDomains::eNeumann)
					 {
					   Vmath::Vcopy(Nfps,&Fwd[id2],1,&penaltyflux[id2],1);
					 }
				 }
                        }           
	        }
	  
// Diffusion: Imposing weak boundary condition for q with flux 
//  uflux = g_D  on Dirichlet boundary condition
//  uflux = u_Fwd  on Neumann boundary condition
    void AdvectionDiffusionReaction::WeakPenaltyBoundary(const int dir,
				                         const Array<OneD, const NekDouble> &Fwd, 
							 const Array<OneD, const NekDouble> &qFwd, 
							 Array<OneD, NekDouble> &penaltyflux,
							 NekDouble C11,
							 NekDouble initialtime)
	{
		int i, j, e, npoints, id1, id2;
                int nbnd = m_fields[0]->GetBndCondExpansions().num_elements();
		int numBDEdge = m_fields[0]->GetBndCondExpansions()[0]->GetExpSize();
	        int Nfps = m_fields[0]->GetBndCondExpansions()[0]->GetExp(0)->GetNumPoints(0) ;
		int nTraceNumPoints = GetTraceNpoints();
		Array<OneD, NekDouble > uterm(nTraceNumPoints);
            
            for(i = 0; i < nbnd; ++i)
            {                 
			
		// Evaluate boundary values g_D or g_N from input files
		  SpatialDomains::ConstInitialConditionShPtr ifunc = m_boundaryConditions->GetInitialCondition(i);
                  npoints = m_fields[0]->GetBndCondExpansions()[i]->GetNpoints();

		  Array<OneD,NekDouble> BDphysics(npoints);
		  Array<OneD,NekDouble> x0(npoints,0.0);
                  Array<OneD,NekDouble> x1(npoints,0.0);
                  Array<OneD,NekDouble> x2(npoints,0.0);  
                
                  m_fields[0]->GetBndCondExpansions()[i]->GetCoords(x0,x1,x2);
		   for(j = 0; j < npoints; j++)
                    {
		       BDphysics[j] = ifunc->Evaluate(x0[j],x1[j],x2[j],initialtime);
                    }
					
				// Weakly impose boundary conditions by modifying flux values
				for (e = 0; e < numBDEdge ; ++e)
				 {
					  id1 = m_fields[i]->GetBndCondExpansions()[0]->GetPhys_Offset(e);
					  id2 = m_fields[i]->GetTrace()->GetPhys_Offset(m_fields[i]->GetTraceMap()->GetBndCondCoeffsToGlobalCoeffsMap(e));
					
					// For Dirichlet boundary condition: qflux = q+ - C_11 (u+ - g_D) (nx, ny)
					if(m_fields[0]->GetBndConditions()[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
					 {
					    Vmath::Vsub(Nfps,&Fwd[id2],1,&BDphysics[id1],1,&uterm[id2],1);
					    Vmath::Vmul(Nfps,&m_traceNormals[dir][id2],1,&uterm[id2],1,&uterm[id2],1);
					    Vmath::Svtvp(Nfps,-1.0*C11,&uterm[id2],1,&qFwd[id2],1,&penaltyflux[id2],1);
					 }
					
					// For Neumann boundary condition: qflux = g_N
					else if((m_fields[0]->GetBndConditions()[i])->GetBoundaryConditionType() == SpatialDomains::eNeumann)
					{
					  Vmath::Vcopy(Nfps,&BDphysics[id1],1,&penaltyflux[id2],1);
					}
				 }
			  }       
	      }

    void AdvectionDiffusionReaction::Summary(std::ostream &out)
    {   
      cout << "=======================================================================" << endl;
      cout << "\tEquation Type   : "<< kEquationTypeStr[m_equationType] << endl;
      ADRBase::SessionSummary(out);
      switch(m_equationType)
      {
      case eSteadyDiffusion: case eSteadyDiffusionReaction:
      case eHelmholtz: case eLaplace: case ePoisson:
          out << "\tLambda          : " << m_boundaryConditions->GetParameter("Lambda") << endl;
          for(int i = 0; i < m_fields.num_elements(); ++i)
          {
              out << "\tForcing (field " << i << ") : " << m_boundaryConditions->GetForcingFunction(i)->GetEquation() << endl;
          }
          
          break;
      case eUnsteadyAdvection: 
          if(m_explicitAdvection)
          {
              out << "\t\tAdvection Advancement   : Explicit" <<endl;
          }
          else
          {
              out << "\t\tAdvection Advancement   : Implicit" <<endl;
          }
          out << "\t\tTime Integration Method : " << LibUtilities::TimeIntegrationMethodMap[m_timeIntMethod] << endl;
          
          ADRBase::TimeParamSummary(out);
          break;
      case eUnsteadyDiffusion:
          if(m_explicitDiffusion)
          {
              out << "\t\tDiffusion Advancement   : Explicit" <<endl;
          }
          else
          {
              out << "\t\tDiffusion Advancement   : Implicit" <<endl;
          }
          out << "\t\tTime Integration Method : " << LibUtilities::TimeIntegrationMethodMap[m_timeIntMethod] << endl;
          ADRBase::TimeParamSummary(out);
          break;
      case eUnsteadyDiffusionReaction:
          if(m_explicitDiffusion)
          {
            out << "\t\tDiffusion Advancement   : Explicit" <<endl;
          }
          else
          {
              out << "\t\tDiffusion Advancement   : Implicit" <<endl;
          }
          if(m_explicitReaction)
          {
              out << "\t\tReaction Advancement    : Explicit" <<endl;
          }
          else
          {
              out << "\t\tReaction Advancement    : Implicit" <<endl;
          }
          out << "\t\tTime Integration Method : " << LibUtilities::TimeIntegrationMethodMap[m_timeIntMethod] << endl;
          ADRBase::TimeParamSummary(out);
          break;
      }
      cout << "=======================================================================" << endl;

    }
    
} //end of namespace

/**
* $Log: AdvectionDiffusionReaction.cpp,v $
* Revision 1.15  2009/04/27 21:37:14  sherwin
* Updated to dump .fld and .chk file in compressed coefficient format
*
* Revision 1.14  2009/03/06 12:00:10  sehunchun
* Some minor changes on nomenclatures and tabbing errors
*
* Revision 1.13  2009/03/05 14:02:38  pvos
* Fixed bug
*
* Revision 1.12  2009/03/05 11:50:32  sehunchun
* Implicit scheme and IMEX scheme are now implemented
*
* Revision 1.11  2009/03/04 14:17:38  pvos
* Removed all methods that take and Expansion as argument
*
* Revision 1.10  2009/03/03 16:11:26  pvos
* New version of TimeIntegrator classes
*
* Revision 1.9  2009/02/28 21:59:09  sehunchun
* Explicit Diffusion solver is added
*
* Revision 1.8  2009/02/16 16:07:03  pvos
* Update of TimeIntegration classes
*
* Revision 1.7  2009/02/10 16:39:35  sherwin
* Added new format of SolverInfo reader to identify EQTYPE
*
* Revision 1.6  2009/02/08 09:13:08  sherwin
* Updates to go with Multiple matrix/variable solve
*
* Revision 1.6  2009/01/06 21:10:34  sherwin
* Updates for virtual calls to IProductWRTBase and introduced reader to handle SOLVERINFO section to specify different solvers
*
* Revision 1.5  2008/11/19 10:53:51  pvos
* Made 2D CG version working
*
* Revision 1.4  2008/11/17 08:20:14  claes
* Temporary fix for CG schemes. 1D CG working (but not for userdefined BC). 1D DG not working
*
* Revision 1.3  2008/11/12 12:12:26  pvos
* Time Integration update
*
* Revision 1.2  2008/11/02 22:38:51  sherwin
* Updated parameter naming convention
*
* Revision 1.1  2008/10/31 10:50:10  pvos
* Restructured directory and CMakeFiles
*
* Revision 1.3  2008/10/29 22:51:07  sherwin
* Updates for const correctness and ODEforcing
*
* Revision 1.2  2008/10/19 15:59:20  sherwin
* Added Summary method
*
* Revision 1.1  2008/10/16 15:25:45  sherwin
* Working verion of restructured AdvectionDiffusionReactionSolver
*
* Revision 1.1  2008/08/22 09:48:23  pvos
* Added Claes' AdvectionDiffusionReaction, ShallowWater and Euler solver
*
**/