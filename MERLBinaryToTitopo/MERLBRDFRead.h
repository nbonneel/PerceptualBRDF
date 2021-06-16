// Copyright 2005 Mitsubishi Electric Research Laboratories All Rights Reserved.

// Permission to use, copy and modify this software and its documentation without
// fee for educational, research and non-profit purposes, is hereby granted, provided
// that the above copyright notice and the following three paragraphs appear in all copies.

// To request permission to incorporate this software into commercial products contact:
// Vice President of Marketing and Business Development;
// Mitsubishi Electric Research Laboratories (MERL), 201 Broadway, Cambridge, MA 02139 or 
// <license@merl.com>.

// IN NO EVENT SHALL MERL BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL,
// OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND
// ITS DOCUMENTATION, EVEN IF MERL HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.

// MERL SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE PROVIDED
// HEREUNDER IS ON AN "AS IS" BASIS, AND MERL HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT,
// UPDATES, ENHANCEMENTS OR MODIFICATIONS.
#pragma once

#include <algorithm>
#include <stdexcept>
#include <assert.h>
#include <iostream>
#include "stdlib.h"
#include "stdint.h"
#include "math.h"
#include <wchar.h>
#include <stdio.h>

#ifdef _MSC_VER // Visual Studio
#define __PRETTY_FUNCTION__ __FUNCTION__ 
#endif

#define BRDF_SAMPLING_RES_THETA_H       90
#define BRDF_SAMPLING_RES_THETA_D       90
#define BRDF_SAMPLING_RES_PHI_D         360

#define RED_SCALE (1.0/1500.0)
#define GREEN_SCALE (1.15/1500.0)
#define BLUE_SCALE (1.66/1500.0)
#ifndef M_PI
#	define M_PI	3.1415926535897932384626433832795
#endif

namespace MERL
{
// cross product of two vectors
inline void cross_product_SB (double* v1, double* v2, double* out)
{
	out[0] = v1[1]*v2[2] - v1[2]*v2[1];
	out[1] = v1[2]*v2[0] - v1[0]*v2[2];
	out[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

// normalize_SB vector
inline void normalize_SB(double* v)
{
	// normalize_SB
	double len = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	v[0] = v[0] / len;
	v[1] = v[1] / len;
	v[2] = v[2] / len;
}

// rotate vector along one axis
inline void rotate_vector_SB(double* vector, double* axis, double angle, double* out)
{
	double temp;
	double cross[3];
	double cos_ang = cos(angle);
	double sin_ang = sin(angle);

	out[0] = vector[0] * cos_ang;
	out[1] = vector[1] * cos_ang;
	out[2] = vector[2] * cos_ang;

	temp = axis[0]*vector[0]+axis[1]*vector[1]+axis[2]*vector[2];
	temp = temp*(1.0-cos_ang);

	out[0] += axis[0] * temp;
	out[1] += axis[1] * temp;
	out[2] += axis[2] * temp;

	cross_product_SB(axis,vector,cross);
	
	out[0] += cross[0] * sin_ang;
	out[1] += cross[1] * sin_ang;
	out[2] += cross[2] * sin_ang;
}

// convert standard coordinates to half vector/difference vector coordinates
inline void std_coords_to_half_diff_coords_SB(double theta_in, double fi_in, double theta_out, double fi_out,
								double& theta_half,double& fi_half,double& theta_diff,double& fi_diff )
{

	// compute in vector
	double in_vec_z = cos(theta_in);
	double proj_in_vec = sin(theta_in);
	double in_vec_x = proj_in_vec*cos(fi_in);
	double in_vec_y = proj_in_vec*sin(fi_in);
	double in[3]= {in_vec_x,in_vec_y,in_vec_z};
	normalize_SB(in);


	// compute out vector
	double out_vec_z = cos(theta_out);
	double proj_out_vec = sin(theta_out);
	double out_vec_x = proj_out_vec*cos(fi_out);
	double out_vec_y = proj_out_vec*sin(fi_out);
	double out[3]= {out_vec_x,out_vec_y,out_vec_z};
	normalize_SB(out);


	// compute halfway vector
	double half_x = (in_vec_x + out_vec_x)/2.0f;
	double half_y = (in_vec_y + out_vec_y)/2.0f;
	double half_z = (in_vec_z + out_vec_z)/2.0f;
	double half[3] = {half_x,half_y,half_z};
	normalize_SB(half);

	// compute  theta_half, fi_half
	theta_half = acos(std::min(1.0,std::max(-1.0,half[2])));
	fi_half = atan2(half[1], half[0]);


	double bi_normal[3] = {0.0, 1.0, 0.0};
	double normal[3] = { 0.0, 0.0, 1.0 };
	double temp[3];
	double diff[3];

	// compute diff vector
	rotate_vector_SB(in, normal , -fi_half, temp);
	rotate_vector_SB(temp, bi_normal, -theta_half, diff);
	
	// compute  theta_diff, fi_diff	
	theta_diff = acos(std::min(1.0,std::max(-1.0,diff[2])));
	fi_diff = atan2(diff[1], diff[0]);

}

// convert standard coordinates to half vector/difference vector coordinates
inline void half_diff_coords_to_std_coords_SB( double theta_half,double fi_half,double theta_diff,double fi_diff ,
                                               double& theta_in, double& fi_in, double& theta_out, double& fi_out)
{
    // compute half-vector

    double h[3] = { cos(fi_half)*sin(theta_half), sin(fi_half)*sin(theta_half), cos(theta_half) } ;

    // compute diff direction in H CS. 

    double diff[3] = { cos(fi_diff)*sin(theta_diff), sin(fi_diff)*sin(theta_diff), cos(theta_diff) } ;

    // then rotate it to put in in global CS

    double bi_normal[3] = {0.0, 1.0, 0.0};
    double normal[3] = { 0.0, 0.0, 1.0 };

    double temp[3];
    double in[3];

    rotate_vector_SB(diff, bi_normal, theta_half, temp);
    rotate_vector_SB(temp, normal , fi_half, in);

    // finally reflect in around H to get out vector

    double dot = h[0]*in[0] + h[1]*in[1] + h[2]*in[2] ;
    double out[3] = { 2*dot*h[0] - in[0], 2*dot*h[1] - in[1], 2*dot*h[2] - in[2]  }; // out = 2*(h*in)*h - in;

    theta_out = acos(std::min(1.0,std::max(-1.0,out[2])));
    fi_out = atan2(out[1],out[0]) ;

    theta_in = acos(std::min(1.0,std::max(-1.0,in[2])));
    fi_in = atan2(in[1],in[0]) ;
}

// Lookup theta_half index
// This is a non-linear mapping!
// In:  [0 .. pi/2]
// Out: [0 .. 89]
inline int theta_half_index_SB(double theta_half,float *delta=NULL)
{
	if(delta)
		*delta=0.0;

	if (theta_half <= 0.0)
		return 0;
	double theta_half_deg = ((theta_half / (M_PI/2.0))*BRDF_SAMPLING_RES_THETA_H);
	double temp = theta_half_deg*BRDF_SAMPLING_RES_THETA_H;
	temp = sqrt(temp);
	int ret_val = (int)floor(temp);

	if(delta)
		*delta = temp - ret_val ;

	if (ret_val < 0) ret_val = 0;
	if (ret_val >= BRDF_SAMPLING_RES_THETA_H)
		ret_val = BRDF_SAMPLING_RES_THETA_H-1;
	return ret_val;
}

// This function inverses the previous function

inline double theta_half_from_index(int index)
{
	double temp = index*index ;
    	double theta_half_deg = temp / BRDF_SAMPLING_RES_THETA_H ;
        double theta_half = theta_half_deg * M_PI/2.0 / BRDF_SAMPLING_RES_THETA_H + 0.0001;
        
        return theta_half ;
}

// Lookup theta_diff index
// In:  [0 .. pi/2]
// Out: [0 .. 89]
inline int theta_diff_index_SB(double theta_diff,float *delta=NULL)
{
	int tmp = int(floor(theta_diff / (M_PI * 0.5) * BRDF_SAMPLING_RES_THETA_D));

	if(delta)
		*delta = theta_diff / (M_PI * 0.5) * BRDF_SAMPLING_RES_THETA_D - tmp;

	if (tmp < 0)
		return 0;
	else if (tmp < BRDF_SAMPLING_RES_THETA_D - 1)
		return tmp;
	else
		return BRDF_SAMPLING_RES_THETA_D - 1;
}

inline double theta_diff_from_index(int index)
{
    return index * (M_PI*0.5) / BRDF_SAMPLING_RES_THETA_D + 0.0001;
}
// Lookup phi_diff index
inline int phi_diff_index_SB(double phi_diff,bool reciprocity,float *delta=NULL)
{
	// Because of reciprocity, the BRDF is unchanged under
	// phi_diff -> phi_diff + M_PI

	int fx = reciprocity?1:2 ;

    if(phi_diff < 0.0)
        phi_diff += fx*M_PI;

    assert(phi_diff <= fx*M_PI) ;

	// In: phi_diff in [0 .. pi]
	// Out: tmp in [0 .. 179]
    int tmp = int(floor(phi_diff / (fx*M_PI) * BRDF_SAMPLING_RES_PHI_D / 2 * fx));

	if(delta)
        *delta = phi_diff / (fx*M_PI) * BRDF_SAMPLING_RES_PHI_D / 2 * fx - tmp ;

	if (tmp < 0)	
		return 0;
	else if (tmp < BRDF_SAMPLING_RES_PHI_D / 2 * fx - 1)
		return tmp;
	else
		return BRDF_SAMPLING_RES_PHI_D / 2 * fx - 1;
}
inline double phi_diff_from_index(int index,bool reciprocity)
{
	int fx=reciprocity?1:2 ;
    return index * fx * M_PI / (BRDF_SAMPLING_RES_PHI_D/2*fx) + 0.0001;	// this offset ensures that n==phi_diff_index_SB(phi_diff_from_index(n)) for all n
}

// Given a pair of incoming/outgoing angles, look up the BRDF.
// Modif JPF - default parameter is needed for reciprocity 
inline void lookup_brdf_val_SB(double* brdf, double theta_in, double fi_in, double theta_out, double fi_out, double& red_val,double& green_val,double& blue_val,bool reciprocity=true,bool interpolate=true,bool allow_negative=false)
{
	// Convert to half-angle / difference angle coordinates
	double theta_half, fi_half, theta_diff, fi_diff;
	
	std_coords_to_half_diff_coords_SB(theta_in, fi_in, theta_out, fi_out, theta_half, fi_half, theta_diff, fi_diff);

	int fx = reciprocity?1:2 ;

	if(interpolate)
	{
		float dfd, dtd, dth ;

        int fi_diff_index    = phi_diff_index_SB(fi_diff,reciprocity,&dfd) ;
		int theta_diff_index = theta_diff_index_SB(theta_diff,&dtd) ;
		int theta_half_index = theta_half_index_SB(theta_half,&dth) ;

		red_val = 0.0 ;
		green_val = 0.0 ;
		blue_val = 0.0 ;

		for(int i=0;i<2;++i)
			for(int j=0;j<2;++j)
				for(int k=0;k<2;++k)
				{
					int ind_ijk = std::min(BRDF_SAMPLING_RES_PHI_D/2*fx-1,fi_diff_index + i) 
									+ std::min(BRDF_SAMPLING_RES_THETA_D-1,theta_diff_index+j) * BRDF_SAMPLING_RES_PHI_D / 2 * fx 
									+ std::min(BRDF_SAMPLING_RES_THETA_H-1,theta_half_index+k) * BRDF_SAMPLING_RES_PHI_D / 2 * fx * BRDF_SAMPLING_RES_THETA_D;

					float delta = (i?dfd:(1-dfd)) * (j?dtd:(1-dtd)) * (k?dth:(1-dth));

					red_val   += delta * brdf[ind_ijk                                                                                ] ;
					green_val += delta * brdf[ind_ijk + BRDF_SAMPLING_RES_THETA_H*BRDF_SAMPLING_RES_THETA_D*BRDF_SAMPLING_RES_PHI_D/2 * fx] ;
					blue_val  += delta * brdf[ind_ijk + BRDF_SAMPLING_RES_THETA_H*BRDF_SAMPLING_RES_THETA_D*BRDF_SAMPLING_RES_PHI_D   * fx] ;
				}
	}
	else
	{
		// Find index.
		// Note that phi_half is ignored, since isotropic BRDFs are assumed
        int ind = phi_diff_index_SB(fi_diff,reciprocity) 
        + theta_diff_index_SB(theta_diff) * BRDF_SAMPLING_RES_PHI_D / 2 * fx 
        + theta_half_index_SB(theta_half) * BRDF_SAMPLING_RES_PHI_D / 2 * fx * BRDF_SAMPLING_RES_THETA_D;

		red_val   = brdf[ind                                                                                      ] ;
		green_val = brdf[ind + BRDF_SAMPLING_RES_THETA_H*BRDF_SAMPLING_RES_THETA_D*BRDF_SAMPLING_RES_PHI_D/2 * fx ] ;
		blue_val  = brdf[ind + BRDF_SAMPLING_RES_THETA_H*BRDF_SAMPLING_RES_THETA_D*BRDF_SAMPLING_RES_PHI_D   * fx ] ;
	}

	red_val   *=   RED_SCALE;
	green_val *= GREEN_SCALE;	
	blue_val  *=  BLUE_SCALE; 	

    if(!allow_negative)
    {
       if(  red_val < 0.0)   red_val = 0 ;
       if(green_val < 0.0) green_val = 0 ;
       if( blue_val < 0.0)  blue_val = 0 ;
    }
}

inline void free_brdf_SB(double *& brdf)
{
    free(brdf) ;
    brdf = NULL ;
}
inline void allocate_brdf_SB(double *& brdf,uint32_t& size,bool reciprocity)
{
	int n = BRDF_SAMPLING_RES_THETA_H * BRDF_SAMPLING_RES_THETA_D * BRDF_SAMPLING_RES_PHI_D / 2 ;

	if(!reciprocity)
		n *= 2 ;
    
	brdf = (double*) malloc (sizeof(double)*3*n);
    	size = 3*n ;
}

inline void fill_brdf_val_SB(double *brdf, bool (*brdf_vals_func)(double,double,double,double&,double&,double&,void *),bool reciprocity,void *param)
{
	int fx = reciprocity?1:2 ;
#pragma omp parallel for schedule(static)
	for(int i=0;i<BRDF_SAMPLING_RES_THETA_H;++i)
		for(int j=0;j<BRDF_SAMPLING_RES_THETA_D;++j)
            for(int k=0;k<BRDF_SAMPLING_RES_PHI_D/2*fx;++k)
			{
				int ind = k + j * BRDF_SAMPLING_RES_PHI_D / 2 * fx + i * BRDF_SAMPLING_RES_PHI_D / 2 * fx * BRDF_SAMPLING_RES_THETA_D;

				double theta_in, fi_in, theta_out, fi_out ;
				double fi_half = 0.0 ;

				double red_val, green_val, blue_val ;

                double    fi_diff =   phi_diff_from_index(k,reciprocity) ;
				double theta_diff = theta_diff_from_index(j) ;
				double theta_half = theta_half_from_index(i) ;

				half_diff_coords_to_std_coords_SB(theta_half, fi_half, theta_diff, fi_diff,theta_in, fi_in, theta_out, fi_out);
                
				fi_out -= fi_in ; // this is because on isotropic BRDFs, fi_in is implicitly 0.

				if(!brdf_vals_func(theta_in,theta_out,fi_out,red_val,green_val,blue_val,param))
				{
					std::cerr << "Lookup error when querying BRDF values for theta_in=" << theta_in << ", theta_out=" << theta_out << ", phi_out=" << fi_out << std::endl;
					continue ;
				}

                brdf[ind                                                                                     ] =   red_val / RED_SCALE;
                brdf[ind + BRDF_SAMPLING_RES_THETA_H*BRDF_SAMPLING_RES_THETA_D*BRDF_SAMPLING_RES_PHI_D/2 * fx] = green_val / GREEN_SCALE;
                brdf[ind + BRDF_SAMPLING_RES_THETA_H*BRDF_SAMPLING_RES_THETA_D*BRDF_SAMPLING_RES_PHI_D   * fx] =  blue_val / BLUE_SCALE;
			}
}

// Given a pair of incoming/outgoing angles, look up the BRDF.
inline void lookup_brdf_val_SB(double* brdf, int index,bool reciprocity, double& red_val,double& green_val,double& blue_val)
{
	if(reciprocity)
	{
		red_val   = brdf[index                                                                                ] *   RED_SCALE;
		green_val = brdf[index + BRDF_SAMPLING_RES_THETA_H*BRDF_SAMPLING_RES_THETA_D*BRDF_SAMPLING_RES_PHI_D/2] * GREEN_SCALE;
		blue_val  = brdf[index + BRDF_SAMPLING_RES_THETA_H*BRDF_SAMPLING_RES_THETA_D*BRDF_SAMPLING_RES_PHI_D  ] *  BLUE_SCALE;
	}
	else
	{
		red_val   = brdf[index                                                                                ] *   RED_SCALE;
		green_val = brdf[index + BRDF_SAMPLING_RES_THETA_H*BRDF_SAMPLING_RES_THETA_D*BRDF_SAMPLING_RES_PHI_D  ] * GREEN_SCALE;
		blue_val  = brdf[index + BRDF_SAMPLING_RES_THETA_H*BRDF_SAMPLING_RES_THETA_D*BRDF_SAMPLING_RES_PHI_D*2] *  BLUE_SCALE;
	}
	if (red_val < 0.0 || green_val < 0.0 || blue_val < 0.0)
		fprintf(stderr, "Below horizon.\n");
}
// write BRDF data

inline bool write_brdf_SB(const char *filename, const double* brdf,uint32_t data_size)
{
    uint32_t dims[3] = { BRDF_SAMPLING_RES_THETA_H, BRDF_SAMPLING_RES_THETA_D, BRDF_SAMPLING_RES_PHI_D/2 } ;

	if(data_size == 3*BRDF_SAMPLING_RES_THETA_H * BRDF_SAMPLING_RES_THETA_D * BRDF_SAMPLING_RES_PHI_D/2)
		dims[2] = BRDF_SAMPLING_RES_PHI_D/2 ;
	else if(data_size == 3*BRDF_SAMPLING_RES_THETA_H * BRDF_SAMPLING_RES_THETA_D * BRDF_SAMPLING_RES_PHI_D)
		dims[2] = BRDF_SAMPLING_RES_PHI_D ;
	else
		throw std::runtime_error("Dimension errors!") ;

    if(data_size != dims[0]*dims[1]*dims[2] * 3)
	{
		std::cerr << __PRETTY_FUNCTION__ << ": wrong sizes. data_size=" << data_size << ", dims[0]*dims[1]*dims[2]=" <<  dims[0]*dims[1]*dims[2]<< std::endl;
		return false;
	}

	FILE *f = fopen(filename, "wb");

	if (!f)
	{
		std::cerr << __PRETTY_FUNCTION__ << ": cannot open file " << filename << " in writing mode." << std::endl;
		return false;
	}

	if(3 != fwrite(dims, sizeof(uint32_t), 3, f))
	{
		std::cerr << __PRETTY_FUNCTION__ << ": writing a file dimensions failed." << std::endl;

		fclose(f) ;
		return false ;
	}

	uint32_t n = dims[0]*dims[1]*dims[2];

	if(3*n != fwrite(brdf, sizeof(double), 3*n, f))
	{
		std::cerr << __PRETTY_FUNCTION__ << ": writing a file content failed." << std::endl;
		fclose(f) ;
		return false ;
	}

	fclose(f);
	return true;
}
// Read BRDF data
inline bool read_brdf_SB(const char *filename, double* &brdf,uint32_t& data_size,bool& reciprocity)
{
	FILE *f = fopen(filename, "rb");
	if (!f)
		return false;

	int dims[3];
	fread(dims, sizeof(int), 3, f);
	int n = dims[0] * dims[1] * dims[2];

	if(dims[2] == BRDF_SAMPLING_RES_PHI_D / 2)
		reciprocity = true ;
	else if(dims[2] == BRDF_SAMPLING_RES_PHI_D)
		reciprocity = false ;
	else
	{
        fprintf(stderr, "Dimensions don't match on Phi. Dims[2]=%d, expected is %d or %d\n",dims[2],BRDF_SAMPLING_RES_PHI_D, BRDF_SAMPLING_RES_PHI_D/2 );
		fclose(f);
		return false; 
	}

    if(dims[0] != BRDF_SAMPLING_RES_THETA_H)
	{
        fprintf(stderr, "Dimensions don't match on ThetaH. dims[0]=%d, expected value=%d\n",dims[0],BRDF_SAMPLING_RES_THETA_H);
        fclose(f);
		return false;
	}
        if(dims[1] != BRDF_SAMPLING_RES_THETA_D)
    {
        fprintf(stderr, "Dimensions don't match on ThetaD. dims[1]=%d, expected value=%d\n",dims[1],BRDF_SAMPLING_RES_THETA_D);
        fclose(f);
        return false ;
    }

    allocate_brdf_SB(brdf,data_size,reciprocity) ;

	fread(brdf, sizeof(double), 3*n, f);

	fclose(f);
	return true;
}
/*
int main_SB(int argc, char *argv[])
{
	const char *filename = argv[1];
	double* brdf;

	// read brdf
	if (!read_brdf_SB(filename, brdf)) 
	{
		fprintf(stderr, "Error reading %s\n", filename);
		exit(1);
	}

	// print out a 16x64x16x64 table table of BRDF values
	const int n = 16;
	for (int i = 0; i < n; i++) 
	{
	    double theta_in = i * 0.5 * M_PI / n;
	    for (int j = 0; j < 4*n; j++) 
		{
			double phi_in = j * 2.0 * M_PI / (4*n);
			for (int k = 0; k < n; k++) 
			{
				double theta_out = k * 0.5 * M_PI / n;
				for (int l = 0; l < 4*n; l++) 
				{
					double phi_out = l * 2.0 * M_PI / (4*n);
					double red,green,blue;
					lookup_brdf_val_SB(brdf, theta_in, phi_in, theta_out, phi_out, red, green, blue);
					printf("%f %f %f\n", (float)red, (float)green, (float)blue);
				}
			}
	    }
	}
	return 0;
}
*/
}
