///////////////////////////////////////////////////////////////////////////////
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
/// The above copyright notice and this permission notice shall be included
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
// Description:
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/SimdLib/tinysimd.hpp>

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <array>
#include <cmath>
#include <iostream>

namespace Nektar
{
namespace SimdLibTests
{
    using namespace tinysimd;
    using vec_t = simd<double>;

    BOOST_AUTO_TEST_CASE(SimdLib_type)
    {
        // specific checks
        #if defined(__AVX2__)
        auto width = vec_t::width;
        auto alignment = vec_t::alignment;
        BOOST_CHECK_EQUAL(width, 4);
        BOOST_CHECK_EQUAL(alignment, 32);
        #endif

    }

    BOOST_AUTO_TEST_CASE(SimdLib_mem_size)
    {
        BOOST_CHECK_EQUAL(sizeof(vec_t), sizeof(double)*vec_t::width);
    }

    BOOST_AUTO_TEST_CASE(SimdLib_load)
    {
        alignas(vec_t::alignment) std::array<double, vec_t::width>
            ascalararr = {1.0, 2.0, 3.0, 4.0};
        vec_t avec;
        avec.load(ascalararr.data());
    }

    BOOST_AUTO_TEST_CASE(SimdLib_load_implicit)
    {
        alignas(vec_t::alignment) std::array<double, vec_t::width>
            ascalararr = {1.0, 2.0, 3.0, 4.0};
        vec_t avec;
        avec = *(reinterpret_cast<vec_t*>(ascalararr.data()));
    }

    BOOST_AUTO_TEST_CASE(SimdLib_load_aligned)
    {
        alignas(vec_t::alignment) std::array<double, vec_t::width>
            ascalararr = {1.0, 2.0, 3.0, 4.0};
        vec_t avec;
        avec.load(ascalararr.data(), is_aligned);
    }

    BOOST_AUTO_TEST_CASE(SimdLib_load_unaligned)
    {
        std::array<double, vec_t::width> ascalararr = {1.0, 2.0, 3.0, 4.0};
        vec_t avec;
        avec.load(ascalararr.data(), is_not_aligned);
    }

    BOOST_AUTO_TEST_CASE(SimdLib_store)
    {
        double val = 4.0;
        vec_t avec(val);
        alignas(vec_t::alignment) std::array<double, vec_t::width> ascalararr{};
        avec.store(ascalararr.data());

        for (size_t i = 0; i < vec_t::width; ++i)
        {
            BOOST_CHECK_EQUAL(ascalararr[i], val);
        }
    }

    BOOST_AUTO_TEST_CASE(SimdLib_store_aligned)
    {
        double val = 4.0;
        vec_t avec(val);
        alignas(vec_t::alignment) std::array<double, vec_t::width> ascalararr{};
        avec.store(ascalararr.data(), is_aligned);

        for (size_t i = 0; i < vec_t::width; ++i)
        {
            BOOST_CHECK_EQUAL(ascalararr[i], val);
        }
    }

    BOOST_AUTO_TEST_CASE(SimdLib_store_unaligned)
    {
        double val = 4.0;
        vec_t avec(val);
        std::array<double, vec_t::width> ascalararr{};
        avec.store(ascalararr.data(), is_not_aligned);

        for (size_t i = 0; i < vec_t::width; ++i)
        {
            BOOST_CHECK_EQUAL(ascalararr[i], val);
        }
    }

    BOOST_AUTO_TEST_CASE(SimdLib_store_non_temporal)
    {
        double val = 4.0;
        vec_t avec(val);
        alignas(vec_t::alignment) std::array<double, vec_t::width> ascalararr{};
        avec.store(ascalararr.data(), is_not_reused);

        for (size_t i = 0; i < vec_t::width; ++i)
        {
            BOOST_CHECK_EQUAL(ascalararr[i], val);
        }
    }

    BOOST_AUTO_TEST_CASE(SimdLib_subscript_assign_read)
    {
        vec_t avec;
        std::array<double, vec_t::width> ascalararr{1, 2, 3, 4};

        for (size_t i = 0; i < vec_t::width; ++i)
        {
            avec[i] = ascalararr[i];
        }

        for (size_t i = 0; i < vec_t::width; ++i)
        {
            BOOST_CHECK_EQUAL(ascalararr[i], avec[i]);
        }
    }

    #if defined(__AVX2__)
    BOOST_AUTO_TEST_CASE(SimdLib_load_any)
    {
        vec_t avec;
        double* p0, * p1, * p2, * p3;
        std::array<double, 16> ascalararr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};

        p0 =  ascalararr.data();
        p1 =  ascalararr.data() + 3;
        p2 =  ascalararr.data() + 5;
        p3 =  ascalararr.data() + 6;

        avec.load(p0, p1, p2, p3);

        BOOST_CHECK_EQUAL(ascalararr[0], avec[0]);
        BOOST_CHECK_EQUAL(ascalararr[3], avec[1]);
        BOOST_CHECK_EQUAL(ascalararr[5], avec[2]);
        BOOST_CHECK_EQUAL(ascalararr[6], avec[3]);

    }

    BOOST_AUTO_TEST_CASE(SimdLib_gather64)
    {
        using index_t = simd<size_t>;
        vec_t avec;
        index_t aindexvec;
        aindexvec[0] = 0;
        aindexvec[1] = 3;
        aindexvec[2] = 5;
        aindexvec[3] = 6;
        std::array<double, 16> ascalararr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};

        avec.gather(ascalararr.data(), aindexvec);

        BOOST_CHECK_EQUAL(ascalararr[0], avec[0]);
        BOOST_CHECK_EQUAL(ascalararr[3], avec[1]);
        BOOST_CHECK_EQUAL(ascalararr[5], avec[2]);
        BOOST_CHECK_EQUAL(ascalararr[6], avec[3]);

    }

    BOOST_AUTO_TEST_CASE(SimdLib_scatter64)
    {
        using index_t = simd<size_t>;
        vec_t avec;
        index_t aindexvec;
        aindexvec[0] = 0;
        aindexvec[1] = 3;
        aindexvec[2] = 5;
        aindexvec[3] = 6;
        std::array<double, 16> ascalararr{};

        avec[0] = 10;
        avec[1] =  9;
        avec[2] =  8;
        avec[3] =  7;

        avec.scatter(ascalararr.data(), aindexvec);

        BOOST_CHECK_EQUAL(avec[0], ascalararr[0]);
        BOOST_CHECK_EQUAL(avec[1], ascalararr[3]);
        BOOST_CHECK_EQUAL(avec[2], ascalararr[5]);
        BOOST_CHECK_EQUAL(avec[3], ascalararr[6]);

    }
    #endif


    BOOST_AUTO_TEST_CASE(SimdLib_add_mul)
    {
        double val1 = -4.0;
        double val2 =  2.0;
        double val3 =  2.0;
        vec_t avec1(val1);
        vec_t avec2(val2);
        vec_t avec3(val3);
        vec_t res = avec1 + avec2 * avec3;
        alignas(vec_t::alignment) std::array<double, vec_t::width>
            ascalararr{};
        res.store(ascalararr.data());

        for (size_t i = 0; i < vec_t::width; ++i)
        {
            BOOST_CHECK_EQUAL(ascalararr[i], val1 + val2 + val3);
        }

    }

    BOOST_AUTO_TEST_CASE(SimdLib_sqrt)
    {
        double val = 4.0;
        vec_t avec(val);
        vec_t asqrt = sqrt(avec);
        alignas(vec_t::alignment) std::array<double, vec_t::width>
            ascalararr{};
        asqrt.store(ascalararr.data());

        for (size_t i = 0; i < vec_t::width; ++i)
        {
            BOOST_CHECK_EQUAL(ascalararr[i], std::sqrt(val));
        }
    }


    BOOST_AUTO_TEST_CASE(SimdLib_abs)
    {
        double val = -4.0;
        vec_t avec(val);
        vec_t aabs = abs(avec);
        alignas(vec_t::alignment) std::array<double, vec_t::width>
            ascalararr{};
        aabs.store(ascalararr.data());

        for (size_t i = 0; i < vec_t::width; ++i)
        {
            BOOST_CHECK_EQUAL(ascalararr[i], std::abs(val));
        }

    }

    BOOST_AUTO_TEST_CASE(SimdLib_load_interleave_unload)
    {
        constexpr size_t nDof{5};
        // no padding in load_interleave deinterleave_store
        constexpr size_t nEle{vec_t::width * 4};
        constexpr size_t nDofBlock = nDof * vec_t::width;

        constexpr size_t size{nDof*nEle};
        std::array<double,size> dofScalarArr{};
        for (size_t i = 0; i < size; ++i)
        {
            dofScalarArr[i] = i;
        }

        // number of blocks
        size_t nBlock = nEle / vec_t::width;

        // aligned vector
        std::vector<vec_t, allocator<vec_t>> dofVectorArr(nDof);

        double* dataPtr = dofScalarArr.data();
        // loop over blocks vec_t::width elements at the time
        for (size_t b = 0; b < nBlock; ++b)
        {
            // load
            load_interleave(dataPtr, nDof, dofVectorArr);

            // manipulate each block
            for (size_t j = 0; j < nDof; ++j)
            {
                dofVectorArr[j] = dofVectorArr[j] + j;
            }

            // store
            deinterleave_store(dofVectorArr, nDof, dataPtr);
            dataPtr += nDofBlock;
        }

        // check
        // for (size_t i = 0; i < size; ++i)
        // {
        //     BOOST_CHECK_EQUAL(dofScalarArr[i], i + 1.0);
        // }

        for (size_t b = 0, i = 0; b < nBlock; ++b)
        {
            for (size_t j = 0; j < nDof; ++j, ++i)
            {
                BOOST_CHECK_EQUAL(dofScalarArr[i], i + j);
            }
        }

    }


}
}