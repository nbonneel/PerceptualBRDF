/* Copyright (C) 2021 Cyril Soler <cyril.soler@inria.fr>
 *
 *                    IsotropicMERLBRDF.cpp
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

#include <stdexcept>
#include "MERLBRDFRead.h"
#include "IsotropicMERLBRDF.h"

using arzh::Spectrum ;

IsotropicMERLBRDF::IsotropicMERLBRDF(const std::string& filename)
    : _filename(filename), _albedo(0,0,0)
{
	data = NULL ;
	
    if(!MERL::read_brdf_SB(filename.c_str(),data,data_size,_assume_reciprocity))
        throw std::runtime_error("Missing file "+filename) ;
}

struct tmpData
{
   void *proxy_param ;
   arzh::Spectrum (*proxy_func)(double,double,double,void *) ;
};

std::string IsotropicMERLBRDF::material_name() const
{
    return extractMaterialName(_filename);
}

static bool brdf_vals_func(double theta_in,double theta_out,double phi_out,double& red,double& green,double& blue,void *param)
{
	const tmpData& dat(*static_cast<tmpData*>(param)) ;

	arzh::Spectrum res = (*dat.proxy_func)(theta_in,theta_out,phi_out,dat.proxy_param) ;

	red = res[0] ;
	green = res[1] ;
	blue = res[2] ;

	return true ;
}

IsotropicMERLBRDF::IsotropicMERLBRDF(arzh::Spectrum (*func)(double,double,double,void *),void *param)
{
	data = NULL ;
	tmpData proxdata ;
        _assume_reciprocity = false ;

	proxdata.proxy_func = func ;
	proxdata.proxy_param = param ;

        MERL::allocate_brdf_SB(data,data_size,_assume_reciprocity) ;
    MERL::fill_brdf_val_SB(data,brdf_vals_func,_assume_reciprocity,&proxdata) ;
}

void IsotropicMERLBRDF::save(const std::string& fname) const
{
    if(!MERL::write_brdf_SB(fname.c_str(),(const double *)data,data_size))
        throw std::runtime_error("Could not write file "+fname) ;
}


IsotropicMERLBRDF::~IsotropicMERLBRDF()
{
	MERL::free_brdf_SB(data) ;
}

arzh::Spectrum IsotropicMERLBRDF::getIntegratedReflectance(float theta_in, float) const
{
	if(_precomp_int.empty())
		precomputeSums() ;

	int n_theta_in = std::min(N_PRECOMP_THETA-1,(int)floor(theta_in * N_PRECOMP_THETA / (M_PI/2.0))) ;

	if(n_theta_in < 0)
		return Spectrum(0,0,0) ;

	return _precomp_int[n_theta_in] ;
}
arzh::Spectrum IsotropicMERLBRDF::getAlbedo() const
{
    if(_precomp_int.empty())
        precomputeSums();

    return _albedo;
}

std::string IsotropicMERLBRDF::getBRDFName() const
{
    std::string BRDFName = _filename.substr(_filename.find_last_of('/', _filename.length())+1, _filename.length()) ;
    BRDFName = BRDFName.substr(0, BRDFName.find_last_of('.', _filename.length()));
    return BRDFName;
}

void IsotropicMERLBRDF::precomputeSums() const
{
#ifdef DEBUG
	std::cout << "Precomputing integrals for ntheta=" << N_PRECOMP_THETA << "." << std::endl ;
#endif

	_precomp_int.resize(N_PRECOMP_THETA,Spectrum(0,0,0)) ;
    _albedo = arzh::Spectrum(0,0,0);

	for(int T=0;T<N_PRECOMP_THETA;++T)
	{
		double theta_in = T*M_PI/2.0 / N_PRECOMP_THETA ;
		arzh::Spectrum total(0,0,0);

		for(int i=0;i<50;++i)
		{
			double theta_out = i*M_PI/2.0/50 ;

			double sin_theta_out = sin(theta_out) ;

			for(int j=0;j<200;++j)
			{
				double   phi_out = j*M_PI*2.0/200 ;

				// operator() includes cos_theta
				total += sin_theta_out * operator()(theta_in,0.0,theta_out,phi_out);
			}
		}

		_precomp_int[T] = total / (50*200) * M_PI * M_PI ;
        _albedo += _precomp_int[T] * cos(theta_in) / (float)N_PRECOMP_THETA * M_PI ;
	}
}

Spectrum IsotropicMERLBRDF::operator()(float theta_in,float phi_in,float theta_out,float phi_out) const
{
	if(theta_in > M_PI/2.0 || theta_out > M_PI/2.0)
		return Spectrum(0,0,0) ;

	double r,g,b ;
    MERL::lookup_brdf_val_SB(data, theta_in, phi_in, theta_out, phi_out, r,g,b,_assume_reciprocity) ;

    return Spectrum(r,g,b) ;// theta_out is the light
}
Spectrum IsotropicMERLBRDF::operator()(float theta_in,float phi_in,float theta_out,float phi_out,bool allow_negative) const
{
    if(theta_in > M_PI/2.0 || theta_out > M_PI/2.0)
        return Spectrum(0,0,0) ;

    double r,g,b ;
    MERL::lookup_brdf_val_SB(data, theta_in, phi_in, theta_out, phi_out, r,g,b,_assume_reciprocity,true,allow_negative) ;

    return Spectrum(r,g,b) ;// theta_out is the light
}

IsotropicMERLBRDF::IsotropicMERLBRDF(const IsotropicMERLBRDF& brdf)
        :_filename ( brdf._filename )
{
    _assume_reciprocity = brdf._assume_reciprocity;
    MERL::allocate_brdf_SB(data,data_size,_assume_reciprocity) ;
	memcpy(data,brdf.data,brdf.data_size * sizeof(double)) ;

	_precomp_int = brdf._precomp_int ;
}

IsotropicMERLBRDF& IsotropicMERLBRDF::operator=(const IsotropicMERLBRDF& brdf)
{ 
    if(data != NULL)
        MERL::free_brdf_SB(data) ;
    
    _assume_reciprocity = brdf._assume_reciprocity;

    MERL::allocate_brdf_SB(data,data_size,_assume_reciprocity) ;
    memcpy(data,brdf.data,brdf.data_size*sizeof(double)) ;

    _filename = brdf._filename ;
    _precomp_int = brdf._precomp_int ;

    return *this ;
}
