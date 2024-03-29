///////////////////////////////////////////////////////////////////////////////
//
// File VelocityCorrectionScheme_ROM.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: Velocity Correction Scheme header 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_VelocityCorrectionScheme_ROM_H
#define NEKTAR_SOLVERS_VelocityCorrectionScheme_ROM_H

#include <IncNavierStokesSolver/EquationSystems/IncNavierStokes.h>
#include "ThirdParty/Eigen/Dense"

namespace Nektar
{
    class VelocityCorrectionScheme_ROM: public IncNavierStokes
    {
    public:

        /// Creates an instance of this class
        static SolverUtils::EquationSystemSharedPtr create(
            const LibUtilities::SessionReaderSharedPtr& pSession,
            const SpatialDomains::MeshGraphSharedPtr &pGraph)
        {
            SolverUtils::EquationSystemSharedPtr p =
                MemoryManager<VelocityCorrectionScheme_ROM>::AllocateSharedPtr(
                    pSession, pGraph);
            p->InitObject();
            return p;
        }

        /// Name of class
        static std::string className;

        /// Constructor.
        VelocityCorrectionScheme_ROM(
            const LibUtilities::SessionReaderSharedPtr& pSession,
            const SpatialDomains::MeshGraphSharedPtr &pGraph);

        virtual ~VelocityCorrectionScheme_ROM();

        virtual void v_InitObject();
        
        void SetUpPressureForcing(
                    const Array<OneD, const Array<OneD, NekDouble> > &fields,
                    Array<OneD, Array<OneD, NekDouble> > &Forcing,
                    const NekDouble aii_Dt)
        {
            v_SetUpPressureForcing( fields, Forcing, aii_Dt);
        }

        void SetUpViscousForcing(
                    const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                    Array<OneD, Array<OneD, NekDouble> > &Forcing,
                    const NekDouble aii_Dt)
        {
            v_SetUpViscousForcing( inarray, Forcing, aii_Dt);
        }
        
        void SolvePressure( const Array<OneD, NekDouble>  &Forcing)
        {
            v_SolvePressure( Forcing);
        }
        
        void SolveViscous(
                    const Array<OneD, const Array<OneD, NekDouble> > &Forcing,
                    Array<OneD, Array<OneD, NekDouble> > &outarray,
                    const NekDouble aii_Dt)
        {
            v_SolveViscous( Forcing, outarray, aii_Dt);
        }

        void SolveUnsteadyStokesSystem(
                    const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                    Array<OneD, Array<OneD, NekDouble> > &outarray,
                    const NekDouble time,
                    const NekDouble a_iixDt);

        void EvaluateAdvection_SetPressureBCs(
                    const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                    Array<OneD, Array<OneD, NekDouble> > &outarray,
                    const NekDouble time)
        {
            v_EvaluateAdvection_SetPressureBCs( inarray, outarray, time);
        }

    double Linfnorm_ITHACA( Array< OneD, NekDouble > component_x, Array< OneD, NekDouble > component_y );
    double L2norm_ITHACA( Array< OneD, NekDouble > component_x, Array< OneD, NekDouble > component_y );
    double Linfnorm_abs_error_ITHACA( Array< OneD, NekDouble > component1_x, Array< OneD, NekDouble > component1_y, Array< OneD, NekDouble > component2_x, Array< OneD, NekDouble > component2_y );
    double L2norm_abs_error_ITHACA( Array< OneD, NekDouble > component1_x, Array< OneD, NekDouble > component1_y, Array< OneD, NekDouble > component2_x, Array< OneD, NekDouble > component2_y );


    protected:
        /// bool to identify if spectral vanishing viscosity is active.
        bool m_useHomo1DSpecVanVisc;
        /// bool to identify if spectral vanishing viscosity is active.
        bool m_useSpecVanVisc;
        /// cutt off ratio from which to start decayhing modes
        NekDouble m_sVVCutoffRatio;
        /// Diffusion coefficient of SVV modes
        NekDouble m_sVVDiffCoeff;
        NekDouble m_sVVCutoffRatioHomo1D;
        /// Diffusion coefficient of SVV modes in homogeneous 1D Direction
        NekDouble m_sVVDiffCoeffHomo1D;
        /// Array of coefficient if power kernel is used in SVV
        Array<OneD, NekDouble> m_svvVarDiffCoeff;
        /// Identifier for Power Kernel otherwise DG kernel
        bool m_IsSVVPowerKernel;
        /// Diffusion coefficients (will be kinvis for velocities)
        Array<OneD, NekDouble> m_diffCoeff;

        /// Variable Coefficient map for the Laplacian which can be activated as part of SVV or otherwise
        StdRegions::VarCoeffMap m_varCoeffLap;

        /// Desired volumetric flowrate
        NekDouble m_flowrate;
        /// Area of the boundary through which we are measuring the flowrate
        NekDouble m_flowrateArea;
        // Bool to identify 3D1HD with forcing explicitly defined
        bool m_homd1DFlowinPlane;
        /// Flux of the Stokes function solution
        NekDouble m_greenFlux;
        /// Current flowrate correction
        NekDouble m_alpha;
        /// Boundary ID of the flowrate reference surface
        int m_flowrateBndID;
        /// Plane ID for cases with homogeneous expansion
        int m_planeID;
        /// Flowrate reference surface
        MultiRegions::ExpListSharedPtr m_flowrateBnd;
        /// Stokes solution used to impose flowrate
        Array<OneD, Array<OneD, NekDouble> > m_flowrateStokes;
        /// Output stream to record flowrate
        std::ofstream m_flowrateStream;
        /// Interval at which to record flowrate data
        int m_flowrateSteps;
        /// Value of aii_dt used to compute Stokes flowrate solution.
        NekDouble m_flowrateAiidt;

        void SetupFlowrate(NekDouble aii_dt);
        NekDouble MeasureFlowrate(
            const Array<OneD, Array<OneD, NekDouble> > &inarray);

        // Virtual functions
        virtual bool v_PostIntegrate(int step);

        virtual void v_GenerateSummary(SolverUtils::SummaryList& s);

        virtual void v_TransCoeffToPhys(void);

        virtual void v_TransPhysToCoeff(void);

        virtual void v_DoInitialise(void);

        virtual Array<OneD, bool> v_GetSystemSingularChecks();

        virtual int v_GetForceDimension();
        
        virtual void v_SetUpPressureForcing(
                    const Array<OneD, const Array<OneD, NekDouble> > &fields,
                    Array<OneD, Array<OneD, NekDouble> > &Forcing,
                    const NekDouble aii_Dt);
        
        virtual void v_SetUpViscousForcing(
                    const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                    Array<OneD, Array<OneD, NekDouble> > &Forcing,
                    const NekDouble aii_Dt);
        
        virtual void v_SolvePressure( const Array<OneD, NekDouble>  &Forcing);
        
        virtual void v_SolveViscous(
                    const Array<OneD, const Array<OneD, NekDouble> > &Forcing,
                    Array<OneD, Array<OneD, NekDouble> > &outarray,
                    const NekDouble aii_Dt);
        
        virtual void v_EvaluateAdvection_SetPressureBCs(
                    const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                    Array<OneD, Array<OneD, NekDouble> > &outarray,
                    const NekDouble time);

        virtual bool v_RequireFwdTrans()
        {
            return false;
        }

        virtual std::string v_GetExtrapolateStr(void)
        {
            return "Standard";
        }
        
        virtual std::string v_GetSubSteppingExtrapolateStr(
                                               const std::string &instr)
        {
            return instr;
        }
        
        Array<OneD, Array< OneD, NekDouble> > m_F;

        void SetUpSVV(void);
        void SetUpExtrapolation(void);
        
        void SVVVarDiffCoeff(const NekDouble velmag, 
                             Array<OneD, NekDouble> &diffcoeff,
                             const Array<OneD, Array<OneD, NekDouble> >
                             &vel = NullNekDoubleArrayOfArray);
        void AppendSVVFactors(
                              StdRegions::ConstFactorMap &factors,
                              MultiRegions::VarFactorsMap &varFactorsMap);
    private:
    
    Array<OneD, Array<OneD, Array<OneD, NekDouble> > > fields_time_trajectory;  
    Array<OneD, Array<OneD, Array<OneD, NekDouble> > > global_fields_time_trajectory;  
    Array<OneD, Array<OneD, Array<OneD, NekDouble> > > ROM_fields_time_trajectory;  
    Array<OneD, Array<OneD, Array<OneD, NekDouble> > > interp_traj;
    Array<OneD, Array<OneD, NekDouble> > last_added_field; 
	Eigen::MatrixXd POD_modes_x;
	Eigen::MatrixXd POD_modes_y;
	Eigen::MatrixXd interp_traj_x;
	Eigen::MatrixXd interp_traj_y;
    int globalNcoeff;
    int step;
    int no_of_added_ones;
    bool ROM_started;
    int ROM_stage;
    int ROM_size_x;
    int ROM_size_y;
	int interp_traj_start;
	int interp_traj_length;
	int interp_traj_mode_number;
	bool overwrite_with_interp_field;
	bool overwrite_with_interp_field_and_diff;
	int POD_collect_start;
	int max_time_samples;
	Array<OneD, NekDouble> collect_interp_relative_L2_error;  
    Array<OneD, NekDouble> collect_interp_relative_Linf_error;  		

        
    };

    typedef std::shared_ptr<VelocityCorrectionScheme_ROM>
                VelocityCorrectionScheme_ROMSharedPtr;

} //end of namespace


#endif //VELOCITY_CORRECTION_SCHEME_H
