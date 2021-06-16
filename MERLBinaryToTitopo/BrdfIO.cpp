/* Copyright (C) 2021 Cyril Soler <cyril.soler@inria.fr>
 *
 *                    BrdfIO.cpp
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

#include <math.h>
#include <stdlib.h>
#include "BrdfIO.h"

#ifdef USE_BRDF_CACHE
threads::Mutex BrdfIO::mBrdfCacheMutex("BRDF Cache mutex") ;
std::map<std::string,const IsotropicMERLBRDF *> BrdfIO::mBrdfCache ;
#endif

const IsotropicMERLBRDF *BrdfIO::loadIsotropicBrdf(const std::string& material_name)
{
    // Put your own directories there.

    static const uint32_t nb_brdf_places = 2 ;
    static const char *brdf_places[nb_brdf_places] = {
        ".",
        "/media/disc/MERL"
    };

    // now search for the file.

    try
    {
        return loadAllFormats(material_name);
    }
    catch(std::exception&)
    {
    }

    for(uint32_t i=0;i<nb_brdf_places;++i)
    {
        try
        {
            std::string fname = std::string(brdf_places[i]) + "/" + material_name + ".binary";

#if defined(__APPLE__)|defined(WIN32)
			// canonicalise_file_name not available on macOS...
			std::string real_name(fname.c_str());
#else
			char *res = canonicalize_file_name(fname.c_str());

			std::string real_name(res);
			free(res);
#endif

#ifdef USE_BRDF_CACHE
            STACK_MUTEX(mBrdfCacheMutex) ;


            if(mBrdfCache.find(real_name) != mBrdfCache.end())
            {
                std::cerr << "(II) brdf found in cache for name " << real_name << std::endl;
                return mBrdfCache[real_name] ;
            }
            else
            {
                std::cerr << "(II) new brdf added to cache: " << real_name << std::endl;
                return mBrdfCache[real_name] = loadAllFormats(real_name) ;
            }
#else
            return loadAllFormats(real_name);
#endif
        }
        catch(std::exception& e)
        {
        }
    }
    std::cerr << "(EE) Cannot find BRDF binary file for material " << material_name << " in any of the supplied directories. Please add your own dir in the brdf_places tab in GplvmXmlIO.cpp" << std::endl;
    return NULL ;
}

struct TitopoData
{
    int inc;
    void *data;
};

static arzh::Spectrum titopoh_proxy_function(double theta_in,double theta_out,double phi_out,void *data)
{
    TitopoData& dat(*static_cast<TitopoData*>(data));

    int inc = dat.inc;
    int n_theta_in  = 90  / inc ;
    int n_theta_out = 90  / inc ;
    int n_phi_out   = 360 / inc ;

    if(phi_out < 0)
        phi_out += 2.0*M_PI ;

    float *titopoh_data = static_cast<float*>(dat.data);

    int theta_in_index  = floor(theta_in /(M_PI/2.0) * n_theta_in );
    int theta_out_index = floor(theta_out/(M_PI/2.0) * n_theta_out);
    int phi_out_index   = floor(phi_out  /(M_PI*2.0) * n_phi_out  );

    float delta_theta_in  = theta_in /(M_PI/2.0) * n_theta_in  - theta_in_index ;
    float delta_theta_out = theta_out/(M_PI/2.0) * n_theta_out - theta_out_index;
    float delta_phi_out   = phi_out  /(M_PI*2.0) * n_phi_out   - phi_out_index  ;

    arzh::Spectrum res(0.0,0.0,0.0);

    for(int i=0;i<2;++i)
    {
       float d1 = i?(1-delta_theta_in):delta_theta_in ;

       assert(d1 >= 0 && d1 <= 1.0);

       for(int j=0;j<2;++j)
       {
          float d2 = j?(1-delta_theta_out):delta_theta_out ;

          assert(d2 >= 0 && d2 <= 1.0);

          for(int k=0;k<2;++k)
          {
             float d3 = k?(1-delta_phi_out):delta_phi_out ;

             assert(d3 >= 0 && d3 <= 1.0);

             int index = std::min(n_phi_out-1,k+phi_out_index) +  n_phi_out * (std::min(n_theta_out-1,j+theta_out_index) + n_theta_out*std::min(n_theta_in-1,i+theta_in_index) ) ;

             res += d1*d2*d3*arzh::Spectrum( titopoh_data[3*index + 0],titopoh_data[3*index + 1],titopoh_data[3*index + 2]);
          }
       }
    }

    return res ;
}

IsotropicMERLBRDF *BrdfIO::loadAllFormats(const std::string& brdfFileName)
{
    if(!strcmp(brdfFileName.c_str() + brdfFileName.length() - std::string(".titopoh").length(), ".titopoh"))
        return loadFromTitopoh(brdfFileName,2);

    if(!strcmp(brdfFileName.c_str() + brdfFileName.length() - std::string(".titopo").length(), ".titopo"))
        return loadFromTitopoh(brdfFileName,1);

    if(!strcmp(brdfFileName.c_str() + brdfFileName.length() - std::string(".binary").length(), ".binary"))
        return new IsotropicMERLBRDF(brdfFileName);

    if(!strcmp(brdfFileName.c_str() + brdfFileName.length() - std::string(".gz").length(), ".gz"))
        return new IsotropicMERLBRDF(brdfFileName);

    throw std::runtime_error("Cannot find proper BRDF format from extension for file \"" + brdfFileName +"\"");
    return NULL;
}

void BrdfIO::convertToTitopo(const IsotropicMERLBRDF& b,int inc,float *data)
{
    uint32_t n_theta_in  = 90 ;
    uint32_t n_theta_out = 90 ;
    uint32_t n_phi_out   = 360 ;

    uint32_t n=0;

    for(uint32_t i=0;i<n_theta_in;i+=inc)
    {
        float theta_in = M_PI/2.0 * i/(float)n_theta_in;

        for(uint32_t j=0;j<n_theta_out;j+=inc)
        {
            float theta_out = M_PI/2.0 * j/(float)n_theta_out;

            for(uint32_t k=0;k<n_phi_out;k+=inc)
            {
                float phi_out = k/(float)n_phi_out * 2*M_PI ;

                arzh::Spectrum s = b(theta_in,0.0,theta_out,phi_out,true) ;	// we allow negative to keep the data unchanged

                data[n++] = s[0] ;
                data[n++] = s[1] ;
                data[n++] = s[2] ;
            }
        }
    }
}

void BrdfIO::saveToTitopo(const IsotropicMERLBRDF& b,int inc,const std::string& brdfFileName)
{
    uint32_t n_theta_in  = 90 ;
    uint32_t n_theta_out = 90 ;
    uint32_t n_phi_out   = 360 ;

    float *data = new float[3*n_theta_in*n_theta_out*n_phi_out/(inc*inc*inc)] ;
    uint32_t n=0;

    for(uint32_t i=0;i<n_theta_in;i+=inc)
    {
        float theta_in = M_PI/2.0 * i/(float)n_theta_in;

        for(uint32_t j=0;j<n_theta_out;j+=inc)
        {
            float theta_out = M_PI/2.0 * j/(float)n_theta_out;

            for(uint32_t k=0;k<n_phi_out;k+=inc)
            {
                float phi_out = k/(float)n_phi_out * 2*M_PI ;

                arzh::Spectrum s = b(theta_in,0.0,theta_out,phi_out,true) ;	// we allow negative to keep the data unchanged

                data[n++] = s[0] ;
                data[n++] = s[1] ;
                data[n++] = s[2] ;
            }
        }
    }

    FILE *f = fopen(brdfFileName.c_str(),"wb") ;

    if(!f)
        throw std::runtime_error("Cannot open file "+brdfFileName+" for writing.") ;

    if(fwrite(data,sizeof(float),3*n_theta_in*n_theta_out*n_phi_out/(inc*inc*inc),f) != 3*n_theta_in*n_theta_out*n_phi_out/(inc*inc*inc))
        throw std::runtime_error("Write error for file "+brdfFileName+". Check disk space and permissions.") ;

    fclose(f);
    delete[] data ;
}


IsotropicMERLBRDF *BrdfIO::loadFromTitopoh(const std::string& brdfFileName,int inc)
{
    uint32_t n_theta_in  = 90 ;
    uint32_t n_theta_out = 90 ;
    uint32_t n_phi_out   = 360 ;

    const int Tin  = n_theta_in  / inc ;
    const int Tout = n_theta_out / inc ;
    const int Pout = n_phi_out   / inc ;

    int brdf_data_size = 3*Tin*Tout*Pout;

    float *data = new float[brdf_data_size];
    uint32_t n=0;

    if(inc == 2 && strcmp(brdfFileName.c_str() + brdfFileName.length() - std::string(".titopoh").length(), ".titopoh"))
        throw std::runtime_error("Input file "+brdfFileName+" should end with .titopoh.") ;
    else if(inc == 1 && strcmp(brdfFileName.c_str() + brdfFileName.length() - std::string(".titopo").length(), ".titopo"))
        throw std::runtime_error("Input file "+brdfFileName+" should end with .titopo.") ;

    FILE *f = fopen(brdfFileName.c_str(),"rb") ;

    if(!f)
        throw std::runtime_error("Cannot open file "+brdfFileName+" for reading.") ;

    if(fread(data, sizeof(float),brdf_data_size,f) != brdf_data_size)
        throw std::runtime_error("Read error for file "+brdfFileName + ".") ;

    fclose(f);

    // std::cerr << __PRETTY_FUNCTION__ << ": Data[223291] = " << data[223291] << std::endl;

    TitopoData dat;
    dat.inc = inc;
    dat.data = data;

    IsotropicMERLBRDF *b = new IsotropicMERLBRDF(titopoh_proxy_function,&dat);

    delete[] data ;

    return b;
}
