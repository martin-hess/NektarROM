#ifndef NEKTAR_LIBRARY_MF_BWDTRANS_H
#define NEKTAR_LIBRARY_MF_BWDTRANS_H

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <LibUtilities/Foundations/Basis.h>

#include "Operator.hpp"
#include "BwdTransKernels.hpp"

namespace Nektar
{
namespace MatrixFree
{

struct BwdTransQuad : public BwdTrans, public Helper<2>
{
    BwdTransQuad(std::vector<LibUtilities::BasisSharedPtr> basis,
                    int nElmt)
        : BwdTrans(basis, nElmt),
          Helper<2>(basis, nElmt),
          m_nmTot(LibUtilities::StdQuadData::getNumberOfCoefficients(
                      this->m_nm[0], this->m_nm[1]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<BwdTransQuad>(basis, nElmt);
    }

    static NekDouble FlopsPerElement(
        const int nm,
        const int nq0,
        const int nq1)
    {
        int loop_i = nq0*nm*nm*2;
        int loop_j = nq1*nq0*nm*2;
        return loop_i + loop_j;
    }

    NekDouble GFlops() final
    {

        const int nm = m_basis[0]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();

        int flops = m_nElmt * BwdTransQuad::FlopsPerElement(nm, nq0, nq1);
        return flops * 1e-9;
    }

    NekDouble NLoads() final
    {
        const int nm0 = m_basis[0]->GetNumModes();
        const int nm1 = m_basis[1]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();

        int load_iqp = nq0 * nm1 * nm0 * 2;
        int load_jiq = nq1 * nq0 * nm1 * 2;
        int load_expected = load_iqp + load_jiq;

        return m_nElmt * load_expected;

    }

    NekDouble NStores() final
    {
        // const int nm0 = m_basis[0]->GetNumModes();
        const int nm1 = m_basis[1]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();

        int store_iq = nq0 * nm1;
        int store_ij = nq0 * nq1;
        int store_expected = store_iq + store_ij;

        return m_nElmt * store_expected;
    }

    NekDouble Ndof() final
    {
        return m_nmTot * this->m_nElmt;
    }

    void operator()(const Array<OneD, const NekDouble> &in,
        Array<OneD,       NekDouble> &out) final
    {
        // Check preconditions
        ASSERTL0(m_basis[0]->GetNumModes() == m_basis[1]->GetNumModes() &&
            m_basis[0]->GetNumPoints() == m_basis[1]->GetNumPoints(),
            "MatrixFree requires homogenous modes/points");

        switch(m_basis[0]->GetNumModes())
        {
            case 2:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 2: BwdTransQuadImpl<2 ,2 ,2 ,2 >(in, out); break;
                    case 3: BwdTransQuadImpl<2 ,2 ,3 ,3 >(in, out); break;
                    case 4: BwdTransQuadImpl<2 ,2 ,4 ,4 >(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransQuad: # of modes / points combo not implemented.");
                } break;
            case 3:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 3: BwdTransQuadImpl<3 ,3 ,3 ,3 >(in, out); break;
                    case 4: BwdTransQuadImpl<3 ,3 ,4 ,4 >(in, out); break;
                    case 5: BwdTransQuadImpl<3 ,3 ,5 ,5 >(in, out); break;
                    case 6: BwdTransQuadImpl<3 ,3 ,6 ,6 >(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransQuad: # of modes / points combo not implemented.");
                } break;
            case 4:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 4: BwdTransQuadImpl<4 ,4 ,4 ,4 >(in, out); break;
                    case 5: BwdTransQuadImpl<4 ,4 ,5 ,5 >(in, out); break;
                    case 6: BwdTransQuadImpl<4 ,4 ,6 ,6 >(in, out); break;
                    case 7: BwdTransQuadImpl<4 ,4 ,7 ,7 >(in, out); break;
                    case 8: BwdTransQuadImpl<4 ,4 ,8 ,8 >(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransQuad: # of modes / points combo not implemented.");
                } break;
            case 5:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 5: BwdTransQuadImpl<5 ,5 ,5 ,5 >(in, out); break;
                    case 6: BwdTransQuadImpl<5 ,5 ,6 ,6 >(in, out); break;
                    case 7: BwdTransQuadImpl<5 ,5 ,7 ,7 >(in, out); break;
                    case 8: BwdTransQuadImpl<5 ,5 ,8 ,8 >(in, out); break;
                    case 9: BwdTransQuadImpl<5 ,5 ,9 ,9 >(in, out); break;
                    case 10: BwdTransQuadImpl<5 ,5 ,10 ,10 >(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransQuad: # of modes / points combo not implemented.");
                } break;
            case 6:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 6: BwdTransQuadImpl<6 ,6 ,6 ,6 >(in, out); break;
                    case 7: BwdTransQuadImpl<6 ,6 ,7 ,7 >(in, out); break;
                    case 8: BwdTransQuadImpl<6 ,6 ,8 ,8 >(in, out); break;
                    case 9: BwdTransQuadImpl<6 ,6 ,9 ,9 >(in, out); break;
                    case 10: BwdTransQuadImpl<6 ,6 ,10 ,10 >(in, out); break;
                    case 11: BwdTransQuadImpl<6 ,6 ,11 ,11 >(in, out); break;
                    case 12: BwdTransQuadImpl<6 ,6 ,12 ,12 >(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransQuad: # of modes / points combo not implemented.");
                } break;
            case 7:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 7: BwdTransQuadImpl<7 ,7 ,7 ,7 >(in, out); break;
                    case 8: BwdTransQuadImpl<7 ,7 ,8 ,8 >(in, out); break;
                    case 9: BwdTransQuadImpl<7 ,7 ,9 ,9 >(in, out); break;
                    case 10: BwdTransQuadImpl<7 ,7 ,10 ,10 >(in, out); break;
                    case 11: BwdTransQuadImpl<7 ,7 ,11 ,11 >(in, out); break;
                    case 12: BwdTransQuadImpl<7 ,7 ,12 ,12 >(in, out); break;
                    case 13: BwdTransQuadImpl<7 ,7 ,13 ,13 >(in, out); break;
                    case 14: BwdTransQuadImpl<7 ,7 ,14 ,14 >(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransQuad: # of modes / points combo not implemented.");
                } break;
            case 8:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 8: BwdTransQuadImpl<8 ,8 ,8 ,8 >(in, out); break;
                    case 9: BwdTransQuadImpl<8 ,8 ,9 ,9 >(in, out); break;
                    case 10: BwdTransQuadImpl<8 ,8 ,10 ,10 >(in, out); break;
                    case 11: BwdTransQuadImpl<8 ,8 ,11 ,11 >(in, out); break;
                    case 12: BwdTransQuadImpl<8 ,8 ,12 ,12 >(in, out); break;
                    case 13: BwdTransQuadImpl<8 ,8 ,13 ,13 >(in, out); break;
                    case 14: BwdTransQuadImpl<8 ,8 ,14 ,14 >(in, out); break;
                    case 15: BwdTransQuadImpl<8 ,8 ,15 ,15 >(in, out); break;
                    case 16: BwdTransQuadImpl<8 ,8 ,16 ,16 >(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransQuad: # of modes / points combo not implemented.");
                } break;;
            default: NEKERROR(ErrorUtil::efatal,
                "BwdTransQuad: # of modes / points combo not implemented.");
        }
    }

    template<int NM0, int NM1, int NQ0, int NQ1>
    void BwdTransQuadImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &output)
    {
        using namespace tinysimd;
        using vec_t = simd<NekDouble>;

        auto *inptr = &input[0];
        auto *outptr = &output[0];

        constexpr auto nqTot = NQ0 * NQ1;
        constexpr auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        vec_t p_sums[NM0 * NQ0]; //Sums over q for each quadpt p_i
        std::vector<vec_t, allocator<vec_t>> tmpIn(m_nmTot), tmpOut(nqTot);

        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            // Load and transpose data
            load_interleave(inptr, m_nmTot, tmpIn);

            BwdTransQuadKernel<NM0, NM1, NQ0, NQ1>(
                tmpIn, this->m_bdata[0], this->m_bdata[1],
                p_sums, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, nqTot, outptr);

            inptr += nmBlocks;
            outptr += nqBlocks;
        }
    }

private:
    int m_nmTot;
};

struct BwdTransTri : public BwdTrans, public Helper<2>
{
    BwdTransTri(std::vector<LibUtilities::BasisSharedPtr> basis, int nElmt)
        : BwdTrans(basis, nElmt),
          Helper<2>(basis, nElmt),
          m_nmTot(LibUtilities::StdTriData::getNumberOfCoefficients(
                      this->m_nm[0], this->m_nm[1]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<BwdTransTri>(basis, nElmt);
    }

    static NekDouble FlopsPerElement(
        const int nm,
        const int nq0,
        const int nq1)
    {
        const int pq_loop = nm * (nm+1) / 2;
        const int inner = 2 * nm + 3;

        return nq1 * (2*pq_loop + nq0*inner);
    }

    NekDouble GFlops() final
    {
        const int nm = m_basis[0]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();

        int flops = m_nElmt * BwdTransTri::FlopsPerElement(nm, nq0, nq1);

        return flops * 1e-9;
    }

    NekDouble NLoads() final
    {
        const int nm0 = m_basis[0]->GetNumModes();
        // const int nm1 = m_basis[1]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();

        int load_jpq = nq1 * (nm0 * (nm0 + 1) / 2) * 2;
        int load_jip =  nq1 * nq0 * nm0 * 2;
        int load_corr = nq1 * nq0 * 3;
        int load_expected = load_jpq + load_jip + load_corr;

        return m_nElmt * load_expected;

    }

    NekDouble NStores() final
    {
        const int nm0 = m_basis[0]->GetNumModes();
        // const int nm1 = m_basis[1]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();

        int store_jp = nq1 * nm0;
        int store_ji = nq1 * nq0;
        int store_expected = store_jp + store_ji;

        return m_nElmt * store_expected;
    }

    NekDouble Ndof() final
    {
        return m_nmTot * this->m_nElmt;
    }

    void operator()(const Array<OneD, const NekDouble> &in,
        Array<OneD,       NekDouble> &out) final
    {
        // Check preconditions
        ASSERTL0(m_basis[0]->GetNumModes() == m_basis[1]->GetNumModes() &&
            m_basis[0]->GetNumPoints() == m_basis[1]->GetNumPoints()+1,
            "MatrixFree requires homogenous modes/points");

        if (m_basis[0]->GetBasisType() == LibUtilities::eModified_A)
        {
            switch(m_basis[0]->GetNumModes())
            {
                case 2: switch(m_basis[0]->GetNumPoints())
                {
                    case 3: BwdTransTriImpl<2, 2, 3, 2, true>(in, out); break;
                    case 4: BwdTransTriImpl<2, 2, 4, 3, true>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransTri: # of modes / points combo not implemented.");
                } break;
            case 3:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 4: BwdTransTriImpl<3, 3, 4, 3, true>(in, out); break;
                    case 5: BwdTransTriImpl<3, 3, 5, 4, true>(in, out); break;
                    case 6: BwdTransTriImpl<3, 3, 6, 5, true>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransTri: # of modes / points combo not implemented.");
                } break;
            case 4:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 5: BwdTransTriImpl<4, 4, 5, 4, true>(in, out); break;
                    case 6: BwdTransTriImpl<4, 4, 6, 5, true>(in, out); break;
                    case 7: BwdTransTriImpl<4, 4, 7, 6, true>(in, out); break;
                    case 8: BwdTransTriImpl<4, 4, 8, 7, true>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransTri: # of modes / points combo not implemented.");
                } break;
            case 5:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 6: BwdTransTriImpl<5, 5, 6, 5, true>(in, out); break;
                    case 7: BwdTransTriImpl<5, 5, 7, 6, true>(in, out); break;
                    case 8: BwdTransTriImpl<5, 5, 8, 7, true>(in, out); break;
                    case 9: BwdTransTriImpl<5, 5, 9, 8, true>(in, out); break;
                    case 10: BwdTransTriImpl<5, 5, 10, 9, true>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransTri: # of modes / points combo not implemented.");
                } break;
            case 6:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 7: BwdTransTriImpl<6, 6, 7, 6, true>(in, out); break;
                    case 8: BwdTransTriImpl<6, 6, 8, 7, true>(in, out); break;
                    case 9: BwdTransTriImpl<6, 6, 9, 8, true>(in, out); break;
                    case 10: BwdTransTriImpl<6, 6, 10, 9, true>(in, out); break;
                    case 11: BwdTransTriImpl<6, 6, 11, 10, true>(in, out); break;
                    case 12: BwdTransTriImpl<6, 6, 12, 11, true>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransTri: # of modes / points combo not implemented.");
                } break;
            case 7:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 8: BwdTransTriImpl<7, 7, 8, 7, true>(in, out); break;
                    case 9: BwdTransTriImpl<7, 7, 9, 8, true>(in, out); break;
                    case 10: BwdTransTriImpl<7, 7, 10, 9, true>(in, out); break;
                    case 11: BwdTransTriImpl<7, 7, 11, 10, true>(in, out); break;
                    case 12: BwdTransTriImpl<7, 7, 12, 11, true>(in, out); break;
                    case 13: BwdTransTriImpl<7, 7, 13, 12, true>(in, out); break;
                    case 14: BwdTransTriImpl<7, 7, 14, 13, true>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransTri: # of modes / points combo not implemented.");
                } break;
            case 8:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 9: BwdTransTriImpl<8, 8, 9, 8, true>(in, out); break;
                    case 10: BwdTransTriImpl<8, 8, 10, 9, true>(in, out); break;
                    case 11: BwdTransTriImpl<8, 8, 11, 10, true>(in, out); break;
                    case 12: BwdTransTriImpl<8, 8, 12, 11, true>(in, out); break;
                    case 13: BwdTransTriImpl<8, 8, 13, 12, true>(in, out); break;
                    case 14: BwdTransTriImpl<8, 8, 14, 13, true>(in, out); break;
                    case 15: BwdTransTriImpl<8, 8, 15, 14, true>(in, out); break;
                    case 16: BwdTransTriImpl<8, 8, 16, 15, true>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransTri: # of modes / points combo not implemented.");
                } break;
            default: NEKERROR(ErrorUtil::efatal,
                "BwdTransTri: # of modes / points combo not implemented.");

            }
        }
        else
        {
            switch(m_basis[0]->GetNumModes())
            {
                case 2: switch(m_basis[0]->GetNumPoints())
                {
                    case 3: BwdTransTriImpl<2 ,2 ,3 ,2, false >(in, out); break;
                    case 4: BwdTransTriImpl<2 ,2 ,4 ,3, false >(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransTri: # of modes / points combo not implemented.");
                } break;
            case 3:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 4: BwdTransTriImpl<3 ,3 ,4 ,3, false >(in, out); break;
                    case 5: BwdTransTriImpl<3 ,3 ,5 ,4, false >(in, out); break;
                    case 6: BwdTransTriImpl<3 ,3 ,6 ,5, false >(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransTri: # of modes / points combo not implemented.");
                } break;
            case 4:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 5: BwdTransTriImpl<4, 4, 5, 4, false>(in, out); break;
                    case 6: BwdTransTriImpl<4, 4, 6, 5, false>(in, out); break;
                    case 7: BwdTransTriImpl<4, 4, 7, 6, false>(in, out); break;
                    case 8: BwdTransTriImpl<4, 4, 8, 7, false>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransTri: # of modes / points combo not implemented.");
                } break;
            case 5:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 6: BwdTransTriImpl<5, 5, 6, 5, false>(in, out); break;
                    case 7: BwdTransTriImpl<5, 5, 7, 6, false>(in, out); break;
                    case 8: BwdTransTriImpl<5, 5, 8, 7, false>(in, out); break;
                    case 9: BwdTransTriImpl<5, 5, 9, 8, false>(in, out); break;
                    case 10: BwdTransTriImpl<5, 5, 10, 9, false>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransTri: # of modes / points combo not implemented.");
                } break;
            case 6:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 7: BwdTransTriImpl<6, 6, 7, 6, false>(in, out); break;
                    case 8: BwdTransTriImpl<6, 6, 8, 7, false>(in, out); break;
                    case 9: BwdTransTriImpl<6, 6, 9, 8, false>(in, out); break;
                    case 10: BwdTransTriImpl<6, 6, 10, 9, false>(in, out); break;
                    case 11: BwdTransTriImpl<6, 6, 11, 10, false>(in, out); break;
                    case 12: BwdTransTriImpl<6, 6, 12, 11, false>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransTri: # of modes / points combo not implemented.");
                } break;
            case 7:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 8: BwdTransTriImpl<7, 7, 8, 7, false>(in, out); break;
                    case 9: BwdTransTriImpl<7, 7, 9, 8, false>(in, out); break;
                    case 10: BwdTransTriImpl<7, 7, 10, 9, false>(in, out); break;
                    case 11: BwdTransTriImpl<7, 7, 11, 10, false>(in, out); break;
                    case 12: BwdTransTriImpl<7, 7, 12, 11, false>(in, out); break;
                    case 13: BwdTransTriImpl<7, 7, 13, 12, false>(in, out); break;
                    case 14: BwdTransTriImpl<7, 7, 14, 13, false>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransTri: # of modes / points combo not implemented.");
                } break;
            case 8:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 9: BwdTransTriImpl<8, 8, 9, 8, false>(in, out); break;
                    case 10: BwdTransTriImpl<8, 8, 10, 9, false>(in, out); break;
                    case 11: BwdTransTriImpl<8, 8, 11, 10, false>(in, out); break;
                    case 12: BwdTransTriImpl<8, 8, 12, 11, false>(in, out); break;
                    case 13: BwdTransTriImpl<8, 8, 13, 12, false>(in, out); break;
                    case 14: BwdTransTriImpl<8, 8, 14, 13, false>(in, out); break;
                    case 15: BwdTransTriImpl<8, 8, 15, 14, false>(in, out); break;
                    case 16: BwdTransTriImpl<8, 8, 16, 15, false>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransTri: # of modes / points combo not implemented.");
                } break;
            default: NEKERROR(ErrorUtil::efatal,
                "BwdTransTri: # of modes / points combo not implemented.");

            }
        }
    }

    template<int NM0, int NM1, int NQ0, int NQ1, bool CORRECT>
    void BwdTransTriImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &output)
    {
        using namespace tinysimd;
        using vec_t = simd<NekDouble>;

        auto* inptr = &input[0];
        auto* outptr = &output[0];

        constexpr auto nqTot = NQ0 * NQ1;
        constexpr auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        vec_t q_sums[NM0]; //Sums over q for each p
        std::vector<vec_t, allocator<vec_t>> tmpIn(m_nmTot), tmpOut(nqTot);

        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            // Load and transpose data
            load_interleave(inptr, m_nmTot, tmpIn);

            BwdTransTriKernel<NM0, NM1, NQ0, NQ1, CORRECT>(
                tmpIn,
                this->m_bdata[0], this->m_bdata[1],
                q_sums,
                tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, nqTot, outptr);

            inptr += nmBlocks;
            outptr += nqBlocks;
        }
    }

private:
    int m_nmTot;
};

struct BwdTransHex : public BwdTrans, public Helper<3>
{
    BwdTransHex(std::vector<LibUtilities::BasisSharedPtr> basis,
                   int nElmt )
    : BwdTrans(basis, nElmt),
        Helper<3>(basis, nElmt),
        m_nmTot(LibUtilities::StdHexData::getNumberOfCoefficients(
                    this->m_nm[0], this->m_nm[1], this->m_nm[2]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<BwdTransHex>(basis, nElmt);
    }

    static NekDouble FlopsPerElement(
        const int nm,
        const int nq0,
        const int nq1,
        const int nq2)
    {
        int loop_i = nq0 * nm *  nm  * nm  * 2;
        int loop_j = nq1 * nq0 * nm  * nm  * 2;
        int loop_k = nq2 * nq1 * nq0 * nm  * 2;

        return (loop_i + loop_j + loop_k);
    }

    NekDouble GFlops() final
    {
        const int nm = m_basis[0]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        int flops = this->m_nElmt * BwdTransHex::FlopsPerElement(nm,nq0,nq1,nq2);
        return flops * 1e-9;
    }

    NekDouble Ndof()
    {
        return m_nmTot * this->m_nElmt;
    }

    NekDouble NLoads() final
    {
        const int nm0 = m_basis[0]->GetNumModes();
        const int nm1 = m_basis[1]->GetNumModes();
        const int nm2 = m_basis[2]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        int load_irqp = nq0 * nm2 * nm1 * nm0 * 2;
        int load_jirq = nq1 * nq0 * nm2 * nm1 * 2;
        int load_kjir = nq2 * nq1 * nq0 * nm2 * 2;
        int load_expected = load_irqp + load_jirq + load_kjir;

        return load_expected * m_nElmt;
    }

    NekDouble NStores() final
    {
        // const int nm0 = m_basis[0]->GetNumModes();
        const int nm1 = m_basis[1]->GetNumModes();
        const int nm2 = m_basis[2]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        int store_irq = nq0 * nm2 * nm1;
        int store_jir = nq1 * nq0 * nm2;
        int store_kji = nq0 * nq1 * nq2;
        int store_expected = store_irq + store_jir + store_kji;

        return store_expected * m_nElmt;
    }

    void operator()(const Array<OneD, const NekDouble> &in,
                                  Array<OneD,       NekDouble> &out) final
    {
        switch(m_basis[0]->GetNumModes())
        {
            case 2:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 2: BwdTransHexImpl<2, 2, 2, 2, 2, 2>(in, out); break;
                    case 3: BwdTransHexImpl<2, 2, 2, 3, 3, 3>(in, out); break;
                    case 4: BwdTransHexImpl<2, 2, 2, 4, 4, 4>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransHex: # of modes / points combo not implemented.");
                } break;
            case 3:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 3: BwdTransHexImpl<3, 3, 3, 3, 3, 3>(in, out); break;
                    case 4: BwdTransHexImpl<3, 3, 3, 4, 4, 4>(in, out); break;
                    case 5: BwdTransHexImpl<3, 3, 3, 5, 5, 5>(in, out); break;
                    case 6: BwdTransHexImpl<3, 3, 3, 6, 6, 6>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransHex: # of modes / points combo not implemented.");
                } break;
            case 4:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 4: BwdTransHexImpl<4, 4, 4, 4, 4, 4>(in, out); break;
                    case 5: BwdTransHexImpl<4, 4, 4, 5, 5, 5>(in, out); break;
                    case 6: BwdTransHexImpl<4, 4, 4, 6, 6, 6>(in, out); break;
                    case 7: BwdTransHexImpl<4, 4, 4, 7, 7, 7>(in, out); break;
                    case 8: BwdTransHexImpl<4, 4, 4, 8, 8, 8>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransHex: # of modes / points combo not implemented.");
                } break;
            case 5:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 5: BwdTransHexImpl<5, 5, 5, 5, 5, 5>(in, out); break;
                    case 6: BwdTransHexImpl<5, 5, 5, 6, 6, 6>(in, out); break;
                    case 7: BwdTransHexImpl<5, 5, 5, 7, 7, 7>(in, out); break;
                    case 8: BwdTransHexImpl<5, 5, 5, 8, 8, 8>(in, out); break;
                    case 9: BwdTransHexImpl<5, 5, 5, 9, 9, 9>(in, out); break;
                    case 10: BwdTransHexImpl<5, 5, 5, 10, 10, 10>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransHex: # of modes / points combo not implemented.");
                } break;
            case 6:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 6: BwdTransHexImpl<6, 6, 6, 6, 6, 6>(in, out); break;
                    case 7: BwdTransHexImpl<6, 6, 6, 7, 7, 7>(in, out); break;
                    case 8: BwdTransHexImpl<6, 6, 6, 8, 8, 8>(in, out); break;
                    case 9: BwdTransHexImpl<6, 6, 6, 9, 9, 9>(in, out); break;
                    case 10: BwdTransHexImpl<6, 6, 6, 10, 10, 10>(in, out); break;
                    case 11: BwdTransHexImpl<6, 6, 6, 11, 11, 11>(in, out); break;
                    case 12: BwdTransHexImpl<6, 6, 6, 12, 12, 12>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransHex: # of modes / points combo not implemented.");
                } break;
            case 7:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 7: BwdTransHexImpl<7, 7, 7, 7, 7, 7>(in, out); break;
                    case 8: BwdTransHexImpl<7, 7, 7, 8, 8, 8>(in, out); break;
                    case 9: BwdTransHexImpl<7, 7, 7, 9, 9, 9>(in, out); break;
                    case 10: BwdTransHexImpl<7, 7, 7, 10, 10, 10>(in, out); break;
                    case 11: BwdTransHexImpl<7, 7, 7, 11, 11, 11>(in, out); break;
                    case 12: BwdTransHexImpl<7, 7, 7, 12, 12, 12>(in, out); break;
                    case 13: BwdTransHexImpl<7, 7, 7, 13, 13, 13>(in, out); break;
                    case 14: BwdTransHexImpl<7, 7, 7, 14, 14, 14>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransHex: # of modes / points combo not implemented.");
                } break;
            case 8:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 8: BwdTransHexImpl<8, 8, 8, 8, 8, 8>(in, out); break;
                    case 9: BwdTransHexImpl<8, 8, 8, 9, 9, 9>(in, out); break;
                    case 10: BwdTransHexImpl<8, 8, 8, 10, 10, 10>(in, out); break;
                    case 11: BwdTransHexImpl<8, 8, 8, 11, 11, 11>(in, out); break;
                    case 12: BwdTransHexImpl<8, 8, 8, 12, 12, 12>(in, out); break;
                    case 13: BwdTransHexImpl<8, 8, 8, 13, 13, 13>(in, out); break;
                    case 14: BwdTransHexImpl<8, 8, 8, 14, 14, 14>(in, out); break;
                    case 15: BwdTransHexImpl<8, 8, 8, 15, 15, 15>(in, out); break;
                    case 16: BwdTransHexImpl<8, 8, 8, 16, 16, 16>(in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransHex: # of modes / points combo not implemented.");
                } break;
            default: NEKERROR(ErrorUtil::efatal,
                "BwdTransHex: # of modes / points combo not implemented.");
        }
    }

    template<int NM0, int NM1, int NM2, int NQ0, int NQ1, int NQ2>
    void BwdTransHexImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &output)
    {
        using namespace tinysimd;
        using vec_t = simd<NekDouble>;

        auto* inptr = &input[0];
        auto* outptr = &output[0];

        constexpr auto nqTot = NQ0 * NQ1 * NQ2;
        constexpr auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        vec_t sum_irq[nqTot], sum_jir[nqTot];
        std::vector<vec_t, allocator<vec_t>> tmpIn(m_nmTot), tmpOut(nqTot);


        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            // Load and transpose data
            load_interleave(inptr, m_nmTot, tmpIn);

            BwdTransHexKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2>(
                tmpIn,
                this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                sum_irq, sum_jir,
                tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, nqTot, outptr);

            inptr += nmBlocks;
            outptr += nqBlocks;
        }
    }

private:
    int m_nmTot;
};

struct BwdTransTet : public BwdTrans, public Helper<3>
{
    BwdTransTet(std::vector<LibUtilities::BasisSharedPtr> basis,
                   int nElmt)
    : BwdTrans(basis, nElmt),
        Helper<3>(basis, nElmt),
        m_nmTot(LibUtilities::StdTetData::getNumberOfCoefficients(
                    this->m_nm[0], this->m_nm[1], this->m_nm[2]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<BwdTransTet>(basis, nElmt);
    }

    static NekDouble FlopsPerElement(
        const int nm,
        const int nq0,
        const int nq1,
        const int nq2)
    {
        int pqr_loop = 2.0 * nm*(nm+1)*(nm+2) / 6;
        int pq_loop = 2.0 * nm*(nm+1) / 2;
        int i_loop = 2.0 * nm + 8 + 4*(nm - 1);

        return ( nq2 * (pqr_loop + nq1*(pq_loop + nq0*(i_loop) )) );
    }

    NekDouble GFlops() final
    {
        const int nm = m_basis[0]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        int flop_estimate = m_nElmt * BwdTransTet::FlopsPerElement(nm, nq0, nq1, nq2);
        return flop_estimate * 1e-9; //Convert to gigaflops.
    }

    NekDouble Ndof() final
    {
        return m_nmTot * this->m_nElmt;
    }

    NekDouble NLoads() final
    {
        const int nm0 = m_basis[0]->GetNumModes();
        const int nm1 = m_basis[1]->GetNumModes();
        // const int nm2 = m_basis[2]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        int load_kpqr = nq2 * (nm0*(nm0+1)/2) * nm1 * 2;
        int load_kjpq = nq2 * nq1 * nm0 * nm1 * 2;
        int load_kjip = nq2 * nq1 * nq0 * nm0 * 2;
        int corr = nq2 * nq1 * nq0 * (2 + nm1*2);
        int load_expected = load_kpqr + load_kjpq + load_kjip + corr;

        return load_expected * m_nElmt;
    }

    NekDouble NStores() final
    {
        const int nm0 = m_basis[0]->GetNumModes();
        const int nm1 = m_basis[1]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        int store_kpq = nq2 * nm0 * nm1;
        int store_kjp = nq2 * nq1 * nm0;
        int store_kji = nq2 * nq1 * nq0;
        int store_expected = store_kpq + store_kjp + store_kji;

        return store_expected * m_nElmt;
    }

    void operator()(const Array<OneD, const NekDouble> &in,
        Array<OneD, NekDouble> &out) final
    {
        // Check preconditions
        ASSERTL0(m_basis[0]->GetNumModes() == m_basis[1]->GetNumModes() &&
            m_basis[0]->GetNumModes() == m_basis[2]->GetNumModes() &&
            m_basis[0]->GetNumPoints() == m_basis[1]->GetNumPoints()+1 &&
            m_basis[0]->GetNumPoints() == m_basis[2]->GetNumPoints()+1,
            "MatrixFree requires homogenous modes/points");

        if (m_basis[0]->GetBasisType() == LibUtilities::eModified_A)
        {
            switch(m_basis[0]->GetNumModes())
            {
                case 2: switch(m_basis[0]->GetNumPoints())
                {
                    case 3: BwdTransTetImpl<2, 2, 2, 3, 2, 2, true>
                        (in, out); break;
                    case 4: BwdTransTetImpl<2, 2, 2, 4, 3, 3, true>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransTet: # of modes / points combo not implemented.");
                } break;
            case 3:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 4: BwdTransTetImpl<3, 3, 3, 4, 3, 3, true>
                        (in, out); break;
                    case 5: BwdTransTetImpl<3, 3, 3, 5, 4, 4, true>
                        (in, out); break;
                    case 6: BwdTransTetImpl<3, 3, 3, 6, 5, 5, true>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransTet: # of modes / points combo not implemented.");
                } break;
            case 4:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 5: BwdTransTetImpl<4, 4, 4, 5, 4, 4, true>
                        (in, out); break;
                    case 6: BwdTransTetImpl<4, 4, 4, 6, 5, 5, true>
                        (in, out); break;
                    case 7: BwdTransTetImpl<4, 4, 4, 7, 6, 6, true>
                        (in, out); break;
                    case 8: BwdTransTetImpl<4, 4, 4, 8, 7, 7, true>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransTet: # of modes / points combo not implemented.");
                } break;
            case 5:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 6: BwdTransTetImpl<5, 5, 5, 6, 5, 5, true>
                        (in, out); break;
                    case 7: BwdTransTetImpl<5, 5, 5, 7, 6, 6, true>
                        (in, out); break;
                    case 8: BwdTransTetImpl<5, 5, 5, 8, 7, 7, true>
                        (in, out); break;
                    case 9: BwdTransTetImpl<5, 5, 5, 9, 8, 8, true>
                        (in, out); break;
                    case 10: BwdTransTetImpl<5, 5, 5, 10, 9, 9, true>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransTet: # of modes / points combo not implemented.");
                } break;
            case 6:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 7: BwdTransTetImpl<6, 6, 6, 7, 6, 6, true>
                        (in, out); break;
                    case 8: BwdTransTetImpl<6, 6, 6, 8, 7, 7, true>
                        (in, out); break;
                    case 9: BwdTransTetImpl<6, 6, 6, 9, 8, 8, true>
                        (in, out); break;
                    case 10: BwdTransTetImpl<6, 6, 6, 10, 9, 9, true>
                        (in, out); break;
                    case 11: BwdTransTetImpl<6, 6, 6, 11, 10, 10, true>
                        (in, out); break;
                    case 12: BwdTransTetImpl<6, 6, 6, 12, 11, 11, true>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransTet: # of modes / points combo not implemented.");
                } break;
            case 7:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 8: BwdTransTetImpl<7, 7, 7, 8, 7, 7, true>
                        (in, out); break;
                    case 9: BwdTransTetImpl<7, 7, 7, 9, 8, 8, true>
                        (in, out); break;
                    case 10: BwdTransTetImpl<7, 7, 7, 10, 9, 9, true>
                        (in, out); break;
                    case 11: BwdTransTetImpl<7, 7, 7, 11, 10, 10, true>
                        (in, out); break;
                    case 12: BwdTransTetImpl<7, 7, 7, 12, 11, 11, true>
                        (in, out); break;
                    case 13: BwdTransTetImpl<7, 7, 7, 13, 12, 12, true>
                        (in, out); break;
                    case 14: BwdTransTetImpl<7, 7, 7, 14, 13, 13, true>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransTet: # of modes / points combo not implemented.");
                } break;
            case 8:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 9: BwdTransTetImpl<8, 8, 8, 9, 8, 8, true>
                        (in, out); break;
                    case 10: BwdTransTetImpl<8, 8, 8, 10, 9, 9, true>
                        (in, out); break;
                    case 11: BwdTransTetImpl<8, 8, 8, 11, 10, 10, true>
                        (in, out); break;
                    case 12: BwdTransTetImpl<8, 8, 8, 12, 11, 11, true>
                        (in, out); break;
                    case 13: BwdTransTetImpl<8, 8, 8, 13, 12, 12, true>
                        (in, out); break;
                    case 14: BwdTransTetImpl<8, 8, 8, 14, 13, 13, true>
                        (in, out); break;
                    case 15: BwdTransTetImpl<8, 8, 8, 15, 14, 14, true>
                        (in, out); break;
                    case 16: BwdTransTetImpl<8, 8, 8, 16, 15, 15, true>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransTet: # of modes / points combo not implemented.");
                } break;
            default: NEKERROR(ErrorUtil::efatal,
                "BwdTransTet: # of modes / points combo not implemented.");

            }
        }
        else
        {
            switch(m_basis[0]->GetNumModes())
            {
                case 2: switch(m_basis[0]->GetNumPoints())
                {
                    case 3: BwdTransTetImpl<2, 2, 2, 3, 2, 2, false>
                        (in, out); break;
                    case 4: BwdTransTetImpl<2, 2, 2, 4, 3, 3, false>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransTet: # of modes / points combo not implemented.");
                } break;
            case 3:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 4: BwdTransTetImpl<3, 3, 3, 4, 3, 3, false>
                        (in, out); break;
                    case 5: BwdTransTetImpl<3, 3, 3, 5, 4, 4, false>
                        (in, out); break;
                    case 6: BwdTransTetImpl<3, 3, 3, 6, 5, 5, false>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransTet: # of modes / points combo not implemented.");
                } break;
            case 4:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 5: BwdTransTetImpl<4, 4, 4, 5, 4, 4, false>
                        (in, out); break;
                    case 6: BwdTransTetImpl<4, 4, 4, 6, 5, 5, false>
                        (in, out); break;
                    case 7: BwdTransTetImpl<4, 4, 4, 7, 6, 6, false>
                        (in, out); break;
                    case 8: BwdTransTetImpl<4, 4, 4, 8, 7, 7, false>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransTet: # of modes / points combo not implemented.");
                } break;
            case 5:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 6: BwdTransTetImpl<5, 5, 5, 6, 5, 5, false>
                        (in, out); break;
                    case 7: BwdTransTetImpl<5, 5, 5, 7, 6, 6, false>
                        (in, out); break;
                    case 8: BwdTransTetImpl<5, 5, 5, 8, 7, 7, false>
                        (in, out); break;
                    case 9: BwdTransTetImpl<5, 5, 5, 9, 8, 8, false>
                        (in, out); break;
                    case 10: BwdTransTetImpl<5, 5, 5, 10, 9, 9, false>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransTet: # of modes / points combo not implemented.");
                } break;
            case 6:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 7: BwdTransTetImpl<6, 6, 6, 7, 6, 6, false>
                        (in, out); break;
                    case 8: BwdTransTetImpl<6, 6, 6, 8, 7, 7, false>
                        (in, out); break;
                    case 9: BwdTransTetImpl<6, 6, 6, 9, 8, 8, false>
                        (in, out); break;
                    case 10: BwdTransTetImpl<6, 6, 6, 10, 9, 9, false>
                        (in, out); break;
                    case 11: BwdTransTetImpl<6, 6, 6, 11, 10, 10, false>
                        (in, out); break;
                    case 12: BwdTransTetImpl<6, 6, 6, 12, 11, 11, false>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransTet: # of modes / points combo not implemented.");
                } break;
            case 7:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 8: BwdTransTetImpl<7, 7, 7, 8, 7, 7, false>
                        (in, out); break;
                    case 9: BwdTransTetImpl<7, 7, 7, 9, 8, 8, false>
                        (in, out); break;
                    case 10: BwdTransTetImpl<7, 7, 7, 10, 9, 9, false>
                        (in, out); break;
                    case 11: BwdTransTetImpl<7, 7, 7, 11, 10, 10, false>
                        (in, out); break;
                    case 12: BwdTransTetImpl<7, 7, 7, 12, 11, 11, false>
                        (in, out); break;
                    case 13: BwdTransTetImpl<7, 7, 7, 13, 12, 12, false>
                        (in, out); break;
                    case 14: BwdTransTetImpl<7, 7, 7, 14, 13, 13, false>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransTet: # of modes / points combo not implemented.");
                } break;
            case 8:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 9: BwdTransTetImpl<8, 8, 8, 9, 8, 8, false>
                        (in, out); break;
                    case 10: BwdTransTetImpl<8, 8, 8, 10, 9, 9, false>
                        (in, out); break;
                    case 11: BwdTransTetImpl<8, 8, 8, 11, 10, 10, false>
                        (in, out); break;
                    case 12: BwdTransTetImpl<8, 8, 8, 12, 11, 11, false>
                        (in, out); break;
                    case 13: BwdTransTetImpl<8, 8, 8, 13, 12, 12, false>
                        (in, out); break;
                    case 14: BwdTransTetImpl<8, 8, 8, 14, 13, 13, false>
                        (in, out); break;
                    case 15: BwdTransTetImpl<8, 8, 8, 15, 14, 14, false>
                        (in, out); break;
                    case 16: BwdTransTetImpl<8, 8, 8, 16, 15, 15, false>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransTet: # of modes / points combo not implemented.");
                } break;
            default: NEKERROR(ErrorUtil::efatal,
                "BwdTransTet: # of modes / points combo not implemented.");

            }

        }
    }

    template<int NM0, int NM1, int NM2, int NQ0, int NQ1, int NQ2, bool CORRECT>
    void BwdTransTetImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &output)
    {
        using namespace tinysimd;
        using vec_t = simd<NekDouble>;

        auto *inptr = input.data();
        auto *outptr = output.data();

        constexpr auto nqTot = NQ0 * NQ1 * NQ2;
        constexpr auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        vec_t fpq[NM0 * NM1], fp[NM0];
        std::vector<vec_t, allocator<vec_t>> tmpIn(m_nmTot), tmpOut(nqTot);

        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            // Load and transpose data
            load_interleave(inptr, m_nmTot, tmpIn);

            BwdTransTetKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, CORRECT>(
                tmpIn,
                this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                fpq, fp,
                tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, nqTot, outptr);

            inptr += nmBlocks;
            outptr += nqBlocks;
        }
    }
private:
    int m_nmTot;
};

struct BwdTransPrism : public BwdTrans, public Helper<3>
{
    BwdTransPrism(std::vector<LibUtilities::BasisSharedPtr> basis, int nElmt)
    : BwdTrans(basis, nElmt),
        Helper<3>(basis, nElmt),
        m_nmTot(LibUtilities::StdPrismData::getNumberOfCoefficients(
                    this->m_nm[0], this->m_nm[1], this->m_nm[2]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<BwdTransPrism>(basis, nElmt);
    }

    static NekDouble FlopsPerElement(
        const int nm,
        const int nq0,
        const int nq1,
        const int nq2)
    {
        int loop_pqr = 2 * nm * nm * (nm + 1) / 2.0;
        int loop_pq = nm * nm * 2.0;
        int loop_p = nm * 2.0;
        int q_correction = nm * 6.0;

        return nq2 * (loop_pqr + nq1 *( loop_pq + nq0*(loop_p + q_correction)));
    }

    double GFlops() final
    {
        const int nm = m_basis[0]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        int flops = m_nElmt * BwdTransPrism::FlopsPerElement(nm, nq0, nq1, nq2);

        return flops * 1e-9;

    }

    NekDouble Ndof() final
    {
        return m_nmTot * this->m_nElmt;
    }

    NekDouble NLoads() final
    {
        const int nm0 = m_basis[0]->GetNumModes();
        const int nm2 = m_basis[2]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        int load_kpqr = nq2 * (nm0*nm0*(nm0+1) / 2) * 2;
        int load_kjpq = nq2 * nq1 * nm0 * nm0 * 2;
        int load_kjip = nq2 * nq1 * nq0 * nm0 * 2;
        int corr = nq2 * nq1 * nq0 * (8 + 4 *(nm2 - 1));
        int load_expected = load_kpqr + load_kjpq + load_kjip + corr;

        return load_expected * m_nElmt;
    }

    NekDouble NStores() final
    {
        const int nm0 = m_basis[0]->GetNumModes();
        // const int nm1 = m_basis[1]->GetNumModes();
        // const int nm2 = m_basis[2]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        int store_kpr = nq2 * nm0 * nm0;
        int store_kjp = nq2 * nq1 * nm0;
        int store_kji = nq2 * nq1 * nq0;
        int store_expected = store_kpr + store_kjp + store_kji;

        return store_expected * m_nElmt;
    }

    void operator()(const Array<OneD, const NekDouble> &in,
                                    Array<OneD,       NekDouble> &out)
    {
        // Check preconditions
        ASSERTL0(m_basis[0]->GetNumModes() == m_basis[1]->GetNumModes() &&
            m_basis[0]->GetNumModes() == m_basis[2]->GetNumModes() &&
            m_basis[0]->GetNumPoints() == m_basis[1]->GetNumPoints() &&
            m_basis[0]->GetNumPoints() == m_basis[2]->GetNumPoints()+1,
            "MatrixFree requires homogenous modes/points");

        if (m_basis[0]->GetBasisType() == LibUtilities::eModified_A)
        {
            switch(m_basis[0]->GetNumModes())
            {
                case 2: switch(m_basis[0]->GetNumPoints())
                {
                    case 3: BwdTransPrismImpl<2, 2, 2, 3, 3, 2, true>
                        (in, out); break;
                    case 4: BwdTransPrismImpl<2, 2, 2, 4, 4, 3, true>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransPrism: # of modes / points combo not implemented.");
                } break;
            case 3:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 4: BwdTransPrismImpl<3, 3, 3, 4, 4, 3, true>
                        (in, out); break;
                    case 5: BwdTransPrismImpl<3, 3, 3, 5, 5, 4, true>
                        (in, out); break;
                    case 6: BwdTransPrismImpl<3, 3, 3, 6, 6, 5, true>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransPrism: # of modes / points combo not implemented.");
                } break;
            case 4:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 5: BwdTransPrismImpl<4, 4, 4, 5, 5, 4, true>
                        (in, out); break;
                    case 6: BwdTransPrismImpl<4, 4, 4, 6, 6, 5, true>
                        (in, out); break;
                    case 7: BwdTransPrismImpl<4, 4, 4, 7, 7, 6, true>
                        (in, out); break;
                    case 8: BwdTransPrismImpl<4, 4, 4, 8, 8, 7, true>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransPrism: # of modes / points combo not implemented.");
                } break;
            case 5:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 6: BwdTransPrismImpl<5, 5, 5, 6, 6, 5, true>
                        (in, out); break;
                    case 7: BwdTransPrismImpl<5, 5, 5, 7, 7, 6, true>
                        (in, out); break;
                    case 8: BwdTransPrismImpl<5, 5, 5, 8, 8, 7, true>
                        (in, out); break;
                    case 9: BwdTransPrismImpl<5, 5, 5, 9, 9, 8, true>
                        (in, out); break;
                    case 10: BwdTransPrismImpl<5, 5, 5, 10, 10, 9, true>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransPrism: # of modes / points combo not implemented.");
                } break;
            case 6:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 7: BwdTransPrismImpl<6, 6, 6, 7, 7, 6, true>
                        (in, out); break;
                    case 8: BwdTransPrismImpl<6, 6, 6, 8, 8, 7, true>
                        (in, out); break;
                    case 9: BwdTransPrismImpl<6, 6, 6, 9, 9, 8, true>
                        (in, out); break;
                    case 10: BwdTransPrismImpl<6, 6, 6, 10, 10, 9, true>
                        (in, out); break;
                    case 11: BwdTransPrismImpl<6, 6, 6, 11, 11, 10, true>
                        (in, out); break;
                    case 12: BwdTransPrismImpl<6, 6, 6, 12, 12, 11, true>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransPrism: # of modes / points combo not implemented.");
                } break;
            case 7:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 8: BwdTransPrismImpl<7, 7, 7, 8, 8, 7, true>
                        (in, out); break;
                    case 9: BwdTransPrismImpl<7, 7, 7, 9, 9, 8, true>
                        (in, out); break;
                    case 10: BwdTransPrismImpl<7, 7, 7, 10, 10, 9, true>
                        (in, out); break;
                    case 11: BwdTransPrismImpl<7, 7, 7, 11, 11, 10, true>
                        (in, out); break;
                    case 12: BwdTransPrismImpl<7, 7, 7, 12, 12, 11, true>
                        (in, out); break;
                    case 13: BwdTransPrismImpl<7, 7, 7, 13, 13, 12, true>
                        (in, out); break;
                    case 14: BwdTransPrismImpl<7, 7, 7, 14, 14, 13, true>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransPrism: # of modes / points combo not implemented.");
                } break;
            case 8:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 9: BwdTransPrismImpl<8, 8, 8, 9, 9, 8, true>
                        (in, out); break;
                    case 10: BwdTransPrismImpl<8, 8, 8, 10, 10, 9, true>
                        (in, out); break;
                    case 11: BwdTransPrismImpl<8, 8, 8, 11, 11, 10, true>
                        (in, out); break;
                    case 12: BwdTransPrismImpl<8, 8, 8, 12, 12, 11, true>
                        (in, out); break;
                    case 13: BwdTransPrismImpl<8, 8, 8, 13, 13, 12, true>
                        (in, out); break;
                    case 14: BwdTransPrismImpl<8, 8, 8, 14, 14, 13, true>
                        (in, out); break;
                    case 15: BwdTransPrismImpl<8, 8, 8, 15, 15, 14, true>
                        (in, out); break;
                    case 16: BwdTransPrismImpl<8, 8, 8, 16, 16, 15, true>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransPrism: # of modes / points combo not implemented.");
                } break;
            default: NEKERROR(ErrorUtil::efatal,
                "BwdTransPrism: # of modes / points combo not implemented.");

            }
        }
        else
        {
            switch(m_basis[0]->GetNumModes())
            {
                case 2: switch(m_basis[0]->GetNumPoints())
                {
                    case 3: BwdTransPrismImpl<2, 2, 2, 3, 3, 2, false>
                        (in, out); break;
                    case 4: BwdTransPrismImpl<2, 2, 2, 4, 4, 3, false>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransPrism: # of modes / points combo not implemented.");
                } break;
            case 3:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 4: BwdTransPrismImpl<3, 3, 3, 4, 4, 3, false>
                        (in, out); break;
                    case 5: BwdTransPrismImpl<3, 3, 3, 5, 5, 4, false>
                        (in, out); break;
                    case 6: BwdTransPrismImpl<3, 3, 3, 6, 6, 5, false>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransPrism: # of modes / points combo not implemented.");
                } break;
            case 4:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 5: BwdTransPrismImpl<4, 4, 4, 5, 5, 4, false>
                        (in, out); break;
                    case 6: BwdTransPrismImpl<4, 4, 4, 6, 6, 5, false>
                        (in, out); break;
                    case 7: BwdTransPrismImpl<4, 4, 4, 7, 7, 6, false>
                        (in, out); break;
                    case 8: BwdTransPrismImpl<4, 4, 4, 8, 8, 7, false>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransPrism: # of modes / points combo not implemented.");
                } break;
            case 5:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 6: BwdTransPrismImpl<5, 5, 5, 6, 6, 5, false>
                        (in, out); break;
                    case 7: BwdTransPrismImpl<5, 5, 5, 7, 7, 6, false>
                        (in, out); break;
                    case 8: BwdTransPrismImpl<5, 5, 5, 8, 8, 7, false>
                        (in, out); break;
                    case 9: BwdTransPrismImpl<5, 5, 5, 9, 9, 8, false>
                        (in, out); break;
                    case 10: BwdTransPrismImpl<5, 5, 5, 10, 10, 9, false>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransPrism: # of modes / points combo not implemented.");
                } break;
            case 6:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 7: BwdTransPrismImpl<6, 6, 6, 7, 7, 6, false>
                        (in, out); break;
                    case 8: BwdTransPrismImpl<6, 6, 6, 8, 8, 7, false>
                        (in, out); break;
                    case 9: BwdTransPrismImpl<6, 6, 6, 9, 9, 8, false>
                        (in, out); break;
                    case 10: BwdTransPrismImpl<6, 6, 6, 10, 10, 9, false>
                        (in, out); break;
                    case 11: BwdTransPrismImpl<6, 6, 6, 11, 11, 10, false>
                        (in, out); break;
                    case 12: BwdTransPrismImpl<6, 6, 6, 12, 12, 11, false>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransPrism: # of modes / points combo not implemented.");
                } break;
            case 7:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 8: BwdTransPrismImpl<7, 7, 7, 8, 8, 7, false>
                        (in, out); break;
                    case 9: BwdTransPrismImpl<7, 7, 7, 9, 9, 8, false>
                        (in, out); break;
                    case 10: BwdTransPrismImpl<7, 7, 7, 10, 10, 9, false>
                        (in, out); break;
                    case 11: BwdTransPrismImpl<7, 7, 7, 11, 11, 10, false>
                        (in, out); break;
                    case 12: BwdTransPrismImpl<7, 7, 7, 12, 12, 11, false>
                        (in, out); break;
                    case 13: BwdTransPrismImpl<7, 7, 7, 13, 13, 12, false>
                        (in, out); break;
                    case 14: BwdTransPrismImpl<7, 7, 7, 14, 14, 13, false>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransPrism: # of modes / points combo not implemented.");
                } break;
            case 8:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 9: BwdTransPrismImpl<8, 8, 8, 9, 9, 8, false>
                        (in, out); break;
                    case 10: BwdTransPrismImpl<8, 8, 8, 10, 10, 9, false>
                        (in, out); break;
                    case 11: BwdTransPrismImpl<8, 8, 8, 11, 11, 10, false>
                        (in, out); break;
                    case 12: BwdTransPrismImpl<8, 8, 8, 12, 12, 11, false>
                        (in, out); break;
                    case 13: BwdTransPrismImpl<8, 8, 8, 13, 13, 12, false>
                        (in, out); break;
                    case 14: BwdTransPrismImpl<8, 8, 8, 14, 14, 13, false>
                        (in, out); break;
                    case 15: BwdTransPrismImpl<8, 8, 8, 15, 15, 14, false>
                        (in, out); break;
                    case 16: BwdTransPrismImpl<8, 8, 8, 16, 16, 15, false>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransPrism: # of modes / points combo not implemented.");
                } break;
            default: NEKERROR(ErrorUtil::efatal,
                "BwdTransPrism: # of modes / points combo not implemented.");

            }

        }
    }

    template<int NM0, int NM1, int NM2, int NQ0, int NQ1, int NQ2, bool CORRECT>
    void BwdTransPrismImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &output)
    {
        using namespace tinysimd;
        using vec_t = simd<NekDouble>;

        auto *inptr = input.data();
        auto *outptr = output.data();

        constexpr auto nqTot = NQ0 * NQ1 * NQ2;
        constexpr auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        vec_t fpq[NM0 * NM1], fp[NM0];
        std::vector<vec_t, allocator<vec_t>> tmpIn(m_nmTot), tmpOut(nqTot);

        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            // Load and transpose data
            load_interleave(inptr, m_nmTot, tmpIn);

            BwdTransPrismKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, CORRECT>(
                tmpIn,
                this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                fpq, fp,
                tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, nqTot, outptr);

            inptr += nmBlocks;
            outptr += nqBlocks;
        }
    }


private:
    int m_nmTot;
};


struct BwdTransPyr : public BwdTrans, public Helper<3>
{
    BwdTransPyr(std::vector<LibUtilities::BasisSharedPtr> basis, int nElmt)
    : BwdTrans(basis, nElmt),
        Helper<3>(basis, nElmt),
        m_nmTot(LibUtilities::StdPyrData::getNumberOfCoefficients(
                    this->m_nm[0], this->m_nm[1], this->m_nm[2]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<BwdTransPyr>(basis, nElmt);
    }

    static NekDouble FlopsPerElement(
        const int nm,
        const int nq0,
        const int nq1,
        const int nq2)
    {
        int loop_pqr = 2 * nm * (nm + 1) *(2*nm +1)/ 6.0;
        int loop_pq = nm * nm * 2.0;
        int loop_p = nm * 2.0;
        int q_correction = nm * 6.0;

        return nq2 * (loop_pqr + nq1 *( loop_pq + nq0*(loop_p + q_correction)));
    }

    double GFlops() final
    {
        const int nm  = m_basis[0]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        int flops = m_nElmt * BwdTransPyr::FlopsPerElement(nm, nq0, nq1, nq2);

        return flops * 1e-9;

    }

    NekDouble Ndof() final
    {
        return m_nmTot * this->m_nElmt;
    }
    NekDouble NLoads() final
    {
        const int nm0 = m_basis[0]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        int load_kpqr = nq2 * (nm0*(nm0+1)*(2*nm0+1) / 6) * 2;
        int load_kjpq = nq2 * nq1 * nm0 * nm0 * 2;
        int load_kjip = nq2 * nq1 * nq0 * nm0 * 2;
        int corr = nq2 * nq1 * nq0 * 8;
        int load_expected = load_kpqr + load_kjpq + load_kjip + corr;

        return load_expected * m_nElmt;
    }

    NekDouble NStores() final
    {
        const int nm0 = m_basis[0]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        int store_kpr = nq2 * nm0 * nm0;
        int store_kjp = nq2 * nq1 * nm0;
        int store_kji = nq2 * nq1 * nq0;
        int store_expected = store_kpr + store_kjp + store_kji;

        return store_expected * m_nElmt;
    }

    void operator()(const Array<OneD, const NekDouble> &in,
                          Array<OneD,       NekDouble> &out)
    {
        // Check preconditions
        ASSERTL0(m_basis[0]->GetNumModes() == m_basis[1]->GetNumModes() &&
                 m_basis[0]->GetNumModes() == m_basis[2]->GetNumModes() &&
                 m_basis[0]->GetNumPoints() == m_basis[1]->GetNumPoints() &&
                 m_basis[0]->GetNumPoints() == m_basis[2]->GetNumPoints()+1,
                 "MatrixFree requires homogenous modes/points");

        if (m_basis[0]->GetBasisType() == LibUtilities::eModified_A)
        {
            switch(m_basis[0]->GetNumModes())
            {
                case 2: switch(m_basis[0]->GetNumPoints())
                {
                    case 3: BwdTransPyrImpl<2, 2, 2, 3, 3, 2, true>
                        (in, out); break;
                    case 4: BwdTransPyrImpl<2, 2, 2, 4, 4, 3, true>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransPyr: # of modes / points combo not implemented.");
                } break;
            case 3:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 4: BwdTransPyrImpl<3, 3, 3, 4, 4, 3, true>
                        (in, out); break;
                    case 5: BwdTransPyrImpl<3, 3, 3, 5, 5, 4, true>
                        (in, out); break;
                    case 6: BwdTransPyrImpl<3, 3, 3, 6, 6, 5, true>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransPyr: # of modes / points combo not implemented.");
                } break;
            case 4:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 5: BwdTransPyrImpl<4, 4, 4, 5, 5, 4, true>
                        (in, out); break;
                    case 6: BwdTransPyrImpl<4, 4, 4, 6, 6, 5, true>
                        (in, out); break;
                    case 7: BwdTransPyrImpl<4, 4, 4, 7, 7, 6, true>
                        (in, out); break;
                    case 8: BwdTransPyrImpl<4, 4, 4, 8, 8, 7, true>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransPyr: # of modes / points combo not implemented.");
                } break;
            case 5:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 6: BwdTransPyrImpl<5, 5, 5, 6, 6, 5, true>
                        (in, out); break;
                    case 7: BwdTransPyrImpl<5, 5, 5, 7, 7, 6, true>
                        (in, out); break;
                    case 8: BwdTransPyrImpl<5, 5, 5, 8, 8, 7, true>
                        (in, out); break;
                    case 9: BwdTransPyrImpl<5, 5, 5, 9, 9, 8, true>
                        (in, out); break;
                    case 10: BwdTransPyrImpl<5, 5, 5, 10, 10, 9, true>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransPyr: # of modes / points combo not implemented.");
                } break;
            case 6:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 7: BwdTransPyrImpl<6, 6, 6, 7, 7, 6, true>
                        (in, out); break;
                    case 8: BwdTransPyrImpl<6, 6, 6, 8, 8, 7, true>
                        (in, out); break;
                    case 9: BwdTransPyrImpl<6, 6, 6, 9, 9, 8, true>
                        (in, out); break;
                    case 10: BwdTransPyrImpl<6, 6, 6, 10, 10, 9, true>
                        (in, out); break;
                    case 11: BwdTransPyrImpl<6, 6, 6, 11, 11, 10, true>
                        (in, out); break;
                    case 12: BwdTransPyrImpl<6, 6, 6, 12, 12, 11, true>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransPyr: # of modes / points combo not implemented.");
                } break;
            case 7:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 8: BwdTransPyrImpl<7, 7, 7, 8, 8, 7, true>
                        (in, out); break;
                    case 9: BwdTransPyrImpl<7, 7, 7, 9, 9, 8, true>
                        (in, out); break;
                    case 10: BwdTransPyrImpl<7, 7, 7, 10, 10, 9, true>
                        (in, out); break;
                    case 11: BwdTransPyrImpl<7, 7, 7, 11, 11, 10, true>
                        (in, out); break;
                    case 12: BwdTransPyrImpl<7, 7, 7, 12, 12, 11, true>
                        (in, out); break;
                    case 13: BwdTransPyrImpl<7, 7, 7, 13, 13, 12, true>
                        (in, out); break;
                    case 14: BwdTransPyrImpl<7, 7, 7, 14, 14, 13, true>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransPyr: # of modes / points combo not implemented.");
                } break;
            case 8:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 9: BwdTransPyrImpl<8, 8, 8, 9, 9, 8, true>
                        (in, out); break;
                    case 10: BwdTransPyrImpl<8, 8, 8, 10, 10, 9, true>
                        (in, out); break;
                    case 11: BwdTransPyrImpl<8, 8, 8, 11, 11, 10, true>
                        (in, out); break;
                    case 12: BwdTransPyrImpl<8, 8, 8, 12, 12, 11, true>
                        (in, out); break;
                    case 13: BwdTransPyrImpl<8, 8, 8, 13, 13, 12, true>
                        (in, out); break;
                    case 14: BwdTransPyrImpl<8, 8, 8, 14, 14, 13, true>
                        (in, out); break;
                    case 15: BwdTransPyrImpl<8, 8, 8, 15, 15, 14, true>
                        (in, out); break;
                    case 16: BwdTransPyrImpl<8, 8, 8, 16, 16, 15, true>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransPyr: # of modes / points combo not implemented.");
                } break;
            default: NEKERROR(ErrorUtil::efatal,
                "BwdTransPyr: # of modes / points combo not implemented.");

            }
        }
        else
        {
            switch(m_basis[0]->GetNumModes())
            {
                case 2: switch(m_basis[0]->GetNumPoints())
                {
                    case 3: BwdTransPyrImpl<2, 2, 2, 3, 3, 2, false>
                        (in, out); break;
                    case 4: BwdTransPyrImpl<2, 2, 2, 4, 4, 3, false>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransPyr: # of modes / points combo not implemented.");
                } break;
            case 3:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 4: BwdTransPyrImpl<3, 3, 3, 4, 4, 3, false>
                        (in, out); break;
                    case 5: BwdTransPyrImpl<3, 3, 3, 5, 5, 4, false>
                        (in, out); break;
                    case 6: BwdTransPyrImpl<3, 3, 3, 6, 6, 5, false>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransPyr: # of modes / points combo not implemented.");
                } break;
            case 4:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 5: BwdTransPyrImpl<4, 4, 4, 5, 5, 4, false>
                        (in, out); break;
                    case 6: BwdTransPyrImpl<4, 4, 4, 6, 6, 5, false>
                        (in, out); break;
                    case 7: BwdTransPyrImpl<4, 4, 4, 7, 7, 6, false>
                        (in, out); break;
                    case 8: BwdTransPyrImpl<4, 4, 4, 8, 8, 7, false>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransPyr: # of modes / points combo not implemented.");
                } break;
            case 5:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 6: BwdTransPyrImpl<5, 5, 5, 6, 6, 5, false>
                        (in, out); break;
                    case 7: BwdTransPyrImpl<5, 5, 5, 7, 7, 6, false>
                        (in, out); break;
                    case 8: BwdTransPyrImpl<5, 5, 5, 8, 8, 7, false>
                        (in, out); break;
                    case 9: BwdTransPyrImpl<5, 5, 5, 9, 9, 8, false>
                        (in, out); break;
                    case 10: BwdTransPyrImpl<5, 5, 5, 10, 10, 9, false>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransPyr: # of modes / points combo not implemented.");
                } break;
            case 6:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 7: BwdTransPyrImpl<6, 6, 6, 7, 7, 6, false>
                        (in, out); break;
                    case 8: BwdTransPyrImpl<6, 6, 6, 8, 8, 7, false>
                        (in, out); break;
                    case 9: BwdTransPyrImpl<6, 6, 6, 9, 9, 8, false>
                        (in, out); break;
                    case 10: BwdTransPyrImpl<6, 6, 6, 10, 10, 9, false>
                        (in, out); break;
                    case 11: BwdTransPyrImpl<6, 6, 6, 11, 11, 10, false>
                        (in, out); break;
                    case 12: BwdTransPyrImpl<6, 6, 6, 12, 12, 11, false>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransPyr: # of modes / points combo not implemented.");
                } break;
            case 7:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 8: BwdTransPyrImpl<7, 7, 7, 8, 8, 7, false>
                        (in, out); break;
                    case 9: BwdTransPyrImpl<7, 7, 7, 9, 9, 8, false>
                        (in, out); break;
                    case 10: BwdTransPyrImpl<7, 7, 7, 10, 10, 9, false>
                        (in, out); break;
                    case 11: BwdTransPyrImpl<7, 7, 7, 11, 11, 10, false>
                        (in, out); break;
                    case 12: BwdTransPyrImpl<7, 7, 7, 12, 12, 11, false>
                        (in, out); break;
                    case 13: BwdTransPyrImpl<7, 7, 7, 13, 13, 12, false>
                        (in, out); break;
                    case 14: BwdTransPyrImpl<7, 7, 7, 14, 14, 13, false>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransPyr: # of modes / points combo not implemented.");
                } break;
            case 8:
                switch(m_basis[0]->GetNumPoints())
                {
                    case 9: BwdTransPyrImpl<8, 8, 8, 9, 9, 8, false>
                        (in, out); break;
                    case 10: BwdTransPyrImpl<8, 8, 8, 10, 10, 9, false>
                        (in, out); break;
                    case 11: BwdTransPyrImpl<8, 8, 8, 11, 11, 10, false>
                        (in, out); break;
                    case 12: BwdTransPyrImpl<8, 8, 8, 12, 12, 11, false>
                        (in, out); break;
                    case 13: BwdTransPyrImpl<8, 8, 8, 13, 13, 12, false>
                        (in, out); break;
                    case 14: BwdTransPyrImpl<8, 8, 8, 14, 14, 13, false>
                        (in, out); break;
                    case 15: BwdTransPyrImpl<8, 8, 8, 15, 15, 14, false>
                        (in, out); break;
                    case 16: BwdTransPyrImpl<8, 8, 8, 16, 16, 15, false>
                        (in, out); break;
                    default: NEKERROR(ErrorUtil::efatal,
                "BwdTransPyr: # of modes / points combo not implemented.");
                } break;
            default: NEKERROR(ErrorUtil::efatal,
                "BwdTransPyr: # of modes / points combo not implemented.");

            }
        }
    }

    template<int NM0, int NM1, int NM2, int NQ0, int NQ1, int NQ2, bool CORRECT>
    void BwdTransPyrImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &output)
    {
        using namespace tinysimd;
        using vec_t = simd<NekDouble>;

        auto *inptr = input.data();
        auto *outptr = output.data();

        constexpr auto nqTot = NQ0 * NQ1 * NQ2;
        constexpr auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        vec_t fpq[NM0 * NM1], fp[NM0];
        std::vector<vec_t, allocator<vec_t>> tmpIn(m_nmTot), tmpOut(nqTot);

        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            // Load and transpose data
            load_interleave(inptr, m_nmTot, tmpIn);

            BwdTransPyrKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, CORRECT>(
                tmpIn, this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                fpq, fp, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, nqTot, outptr);

            inptr += nmBlocks;
            outptr += nqBlocks;
        }
    }

private:
    int m_nmTot;
};

} // namespace MatrixFree
} // namespace Nektar

#endif
