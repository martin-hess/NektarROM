////////////////////////////////////////////////////////////////////////////////
//
//  File:  Geometry.cpp
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
//  Description:  This file contains the base class implementation for the
//                Geometry class.
//
//
////////////////////////////////////////////////////////////////////////////////
#include "pchSpatialDomains.h"

#include <SpatialDomains/Geometry.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        Geometry::Geometry():
            m_coordim(0),
            m_state(eNotFilled),
            m_geomShapeType(eNoGeomShapeType),
            m_globalID(-1)
        {
        }

        Geometry::Geometry(const int coordim):
            m_coordim(coordim),
            m_state(eNotFilled),
            m_geomShapeType(eNoGeomShapeType),
            m_globalID(-1)
        {
        }

        Geometry::~Geometry()
        {
        }

        GeomFactorsVector Geometry::m_regGeomFactorsManager;
        GeomFactorsSharedPtr Geometry::ValidateRegGeomFactor(GeomFactorsSharedPtr geomFactor)
        {
            GeomFactorsSharedPtr returnval = geomFactor;
#if 0
            bool found = false;
            if (geomFactor->GetGtype() == eRegular)
            {
                for (GeomFactorsVectorIter iter = m_regGeomFactorsManager.begin();
                    iter != m_regGeomFactorsManager.end();
                    ++iter)
                {
                    if (**iter == *geomFactor)
                    {
                        returnval = *iter;
                        found = true;
                        break;
                    }
                }

                if (!found)
                {
                    m_regGeomFactorsManager.push_back(geomFactor);
                    returnval = geomFactor;
                }
            }
#endif
            return returnval;
        }

        bool SortByGlobalId(const boost::shared_ptr<Geometry>& lhs, 
            const boost::shared_ptr<Geometry>& rhs)
        {
            return lhs->GetGlobalID() < rhs->GetGlobalID();
        }

        bool GlobalIdEquality(const boost::shared_ptr<Geometry>& lhs, 
            const boost::shared_ptr<Geometry>& rhs)
        {
            return lhs->GetGlobalID() == rhs->GetGlobalID();
        }
    }; //end of namespace
}; //end of namespace

//
// $Log: Geometry.cpp,v $
// Revision 1.11  2008/05/29 19:00:55  delisi
// Constructors initialize m_globalID to -1.
//
// Revision 1.10  2008/05/28 21:52:27  jfrazier
// Added GeomShapeType initialization for the different shapes.
//
// Revision 1.9  2007/07/26 19:38:47  jfrazier
// Added general equation evaluation.
//
// Revision 1.8  2007/07/26 18:02:42  jfrazier
// Manage the storage of geofactors.
//
// Revision 1.7  2007/07/20 02:15:08  bnelson
// Replaced boost::shared_ptr with Nektar::ptr
//
// Revision 1.6  2007/07/10 22:21:00  jfrazier
// Revision of geo fac manager to test for equality.
//
// Revision 1.5  2007/07/10 17:06:31  jfrazier
// Added method and underlying structure to manage geomfactors.
//
// Revision 1.4  2007/03/14 21:24:08  sherwin
// Update for working version of MultiRegions up to ExpList1D
//
// Revision 1.3  2006/08/24 18:50:00  jfrazier
// Completed error checking on permissable composite item combinations.
//
// Revision 1.2  2006/06/01 14:15:30  sherwin
// Added typdef of boost wrappers and made GeoFac a boost shared pointer.
//
// Revision 1.1  2006/05/04 18:59:00  kirby
// *** empty log message ***
//
// Revision 1.14  2006/04/09 02:08:34  jfrazier
// Added precompiled header.
//
// Revision 1.13  2006/03/13 11:17:03  sherwin
//
// First compiing version of Demos in SpatialDomains and LocalRegions. However they do not currently seem to execute properly
//
// Revision 1.12  2006/03/12 14:20:42  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.11  2006/03/12 07:42:02  sherwin
//
// Updated member names and StdRegions call. Still has not been compiled
//
// Revision 1.10  2006/02/19 01:37:33  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//
