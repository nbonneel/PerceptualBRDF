/* Copyright (C) 2021 Cyril Soler <cyril.soler@inria.fr>
 *
 *                    BrdfIO.h
 *
 * This code is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This code is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this code; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#pragma once

#include <map>
#include <string>
#include <stdint.h>

#include "IsotropicMERLBRDF.h"

class BrdfIO
{
public:
    /*!
     * \brief loadIsotropicBrdf
     * 			Loads an isotropic MERL Brdf that corresponds to a given material name. The file is searched into a collection
     * 			of directories that is listed in the cpp file.
     *
     * \param material_name			should be just the name of the material, e.g. "gold-paint", or "brass"
     * \return
     */
    static const IsotropicMERLBRDF *loadIsotropicBrdf(const std::string& material_name) ;

    static std::map<std::string,const IsotropicMERLBRDF *> mBrdfCache ;

    static void saveToTitopo(const IsotropicMERLBRDF& b,int inc,const std::string& brdfFileName);
    static void convertToTitopo(const IsotropicMERLBRDF& b,int inc,float *data);

protected:
    static IsotropicMERLBRDF *loadFromTitopoh(const std::string& brdfFileName,int inc);
    static IsotropicMERLBRDF *loadAllFormats(const std::string& brdfFileName);

#ifdef USE_BRDF_CACHE
	static threads::Mutex mBrdfCacheMutex ;
#endif
};
