////////////////////////////////////////////////////////////////////////////////
//
//  File: TetGenInterface.h
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: class for interacting with tetgen
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_MESHUTILS_EXTLIBINTERFACE_TETGENINTERFACE_H
#define NEKTAR_MESHUTILS_EXTLIBINTERFACE_TETGENINTERFACE_H

#include <boost/shared_ptr.hpp>

#define TETLIBRARY
#include <tetgen.h>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

#include <NekMeshUtils/MeshElements/MeshElements.h>

#include <NekMeshUtils/NekMeshUtilsDeclspec.h>

namespace Nektar
{
namespace NekMeshUtils
{

/**
 * @brief class for interacting with the external library tetgen
 */
class TetGenInterface
{
    public:
        friend class MemoryManager<TetGenInterface>;

        /**
         * @brief default constructor
         */
        NEKMESHUTILS_EXPORT TetGenInterface()
        {
        };

        /**
         * @brief assign parameters for meshing
         */
        NEKMESHUTILS_EXPORT void InitialMesh(std::map<int, NodeSharedPtr> tgidton,
                         std::vector<Array<OneD, int> > tri);

        /**
         * @brief gets the locations of the stiener points added by tetgen
         */
        NEKMESHUTILS_EXPORT void GetNewPoints(int num, std::vector<Array<OneD, NekDouble> > &newp);

        /**
         * @brief refines a previously made tetmesh with node delta information from the Octree
         */
        NEKMESHUTILS_EXPORT void RefineMesh(std::map<int, NekDouble> delta);

        /**
         * @brief get the list of connectivites of the nodes
         */
        NEKMESHUTILS_EXPORT std::vector<Array<OneD, int> > Extract();

        /**
         * @brief clear previous mesh
         */
        NEKMESHUTILS_EXPORT void freetet();

    private:

        ///tetgen objects
        tetgenio surface, output, input, additional;
};

typedef boost::shared_ptr<TetGenInterface> TetGenInterfaceSharedPtr;

}
}
#endif
