/* Copyright (C) 2021 Cyril Soler <cyril.soler@inria.fr>
 *
 *                    BRDF.h
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

#include <stdexcept>
#include <Spectrum.h>

class BRDF
{
public:
    BRDF() {}
    BRDF(const std::string& filename) : _filename(filename) {}
    virtual ~BRDF() {}

    // Returns the value of the BRDF (without the cos_theta) for the given incoming/outgoing directions.
    // Vectors are expressed in the local coordinate system of course.
    //
    virtual arzh::Spectrum operator()(float theta_in,float phi_in,float theta_out,float phi_out) const=0;

    virtual std::string getFilename() const {return _filename;}

    virtual double* getBRDFData() const { return NULL ;}

    virtual arzh::Spectrum getIntegratedReflectance(float theta_in,float phi_in) const =0;
    virtual arzh::Spectrum getAlbedo() const
    {
        throw std::runtime_error("Albedo calculation not implemented in your base class?" );
        return arzh::Spectrum(0,0,0);
    }

    static std::string extractMaterialName(const std::string& filename)
    {
        if(filename.length() < strlen(".binary") || strcmp(&filename[filename.length() - strlen(".binary")],".binary"))
            return std::string( "unknown" );

        int i;

        for(i=filename.length();i > 0;--i)
            if(filename[i] == '/' || filename[i] == '\\')
            {
                i++;
                break ;
            }

        return filename.substr(i,filename.length()-i-strlen(".binary")) ;
    }
protected:
    std::string _filename;
};


