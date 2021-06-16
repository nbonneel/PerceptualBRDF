/* Copyright (C) 2021 Cyril Soler <cyril.soler@inria.fr>
 *
 *                    IsotropicMERLBRDF.h
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

#include <vector>
#include "BRDF.h"

class IsotropicMERLBRDF: public BRDF
{
public:
    IsotropicMERLBRDF(const std::string& filename) ;
    IsotropicMERLBRDF(arzh::Spectrum (*func)(double,double,double,void*),void*) ;	// func(theta_in,theta_out,phi_out)

    virtual ~IsotropicMERLBRDF() ;

    virtual arzh::Spectrum operator()(float theta_in,float phi_in,float theta_out,float phi_out) const ;
    virtual arzh::Spectrum operator()(float theta_in,float phi_in,float theta_out,float phi_out,bool allow_negative) const ;

#ifdef TO_REMOVE
    virtual H_BRDF *convertToSH(int LMax) const ;
#endif
    virtual std::string getFilename() const {return _filename;}

    virtual double* getBRDFData() const { return data ; }
    virtual int getBRDFDataSize() const { return data_size ; }
    const std::string& filename() const { return _filename ; }
    std::string material_name() const ;
    bool assumeReciprocal() const { return _assume_reciprocity ; }

    virtual arzh::Spectrum getIntegratedReflectance(float theta_in,float phi_in) const;
    virtual arzh::Spectrum getAlbedo() const;

    std::string getBRDFName() const ;

    /*!
     * \brief IsotropicMERLBRDF
     * \param brdf
     * 		Saves to MERL .binary format. If reciprocity is assumed, the output will be identical to the
     * 		original MERL data.
     */
    void save(const std::string& filename) const ;

    IsotropicMERLBRDF(const IsotropicMERLBRDF& brdf); 	// deep copy operators
    IsotropicMERLBRDF& operator=(const IsotropicMERLBRDF& brdf);
protected:
    static const int N_PRECOMP_THETA = 50 ;

    void precomputeSums() const ;

    double *data ;	// MERL data.
    uint32_t data_size ;	// number of doubles
    std::string _filename;
    bool _assume_reciprocity ;

    mutable std::vector<arzh::Spectrum> _precomp_int ;
    mutable arzh::Spectrum _albedo;

};


