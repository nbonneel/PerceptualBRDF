/* Copyright (C) 2021 Cyril Soler <cyril.soler@inria.fr>
 *
 *                    main.cpp
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

#include <iostream>
#include <sstream>
#include <fenv.h>
#define _USE_MATH_DEFINES // for C++  
#include <cmath>  

//Modif JPF -  FPE handling on osX is handled by SSE instructions
#ifdef __APPLE__
#include <xmmintrin.h>
#endif

#include "argstream.h"
#include "IsotropicMERLBRDF.h"
#include "BrdfIO.h"

bool has_extension(const std::string& fname,const std::string& extension)
{
    return fname.length() > extension.length() && !strcmp( fname.c_str() + fname.length() - extension.length(),extension.c_str());
}

int main(int argc, char *argv[])
{
    try
    {
        argstream as(argc,argv) ;

        std::string input_brdf,output_brdf;

        as >> parameter('i',"input" ,input_brdf ,"Input brdf (MERL .binary / .titopo)",false)
                        >> parameter('o',"output",output_brdf,"Output brdf (MERL .titopo / .binary)",false)
                        >> help('h',"help","display this help");

        as.defaultErrorHandling() ;

        const IsotropicMERLBRDF *pB = BrdfIO::loadIsotropicBrdf(input_brdf);
        const IsotropicMERLBRDF& b(*pB);

        if(has_extension(output_brdf,std::string(".titopo")))
            BrdfIO::saveToTitopo(b,1,output_brdf);
        else if(has_extension(output_brdf,std::string(".binary")))
            b.save(output_brdf) ;
        else
            throw std::runtime_error("Only .binary or .titopo is accepted as output extension");

#ifdef DEBUG
        std::cerr << "Modified BRDF saved to file " << output_brdf << std::endl;
#endif

        return 0;
    }
    catch(std::exception& e)
    {
        std::cerr << "Exception never handled: " << e.what() << std::endl;
        return 1;
    }
}






