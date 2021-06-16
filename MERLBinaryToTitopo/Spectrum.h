/****************************************************************************

 Copyright (C) 2002-2010 Cyril Soler. All rights reserved.

 This file is part of the HQR (High Quality Rendering) plateform.

 http://artis.imag.fr/~Cyril.Soler/HQR

 This file may be used under the terms of the GNU General Public License 
 versions 2.0 or 3.0 as published by the Free Software Foundation and
 appearing in the LICENSE file included in the packaging of this file.

 This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

*****************************************************************************/


#pragma once

#if defined(__GNUC__) && defined(__SSE__)
//#define USE_SSE_IN_SPECTRUM
#endif

#include <iostream>
#include <string.h>
#include <assert.h>

#ifdef DEBUG
#include <stdexcept>
#endif

#ifdef USE_SSE_IN_SPECTRUM
#include <xmmintrin.h>
//# warning Compiling using sse/sse2 instructions
#endif

// -------------------------------------------------------------------------- //
//                                                                            //
//  INTERFACE DECLARATION                                                     //
//                                                                            //
// -------------------------------------------------------------------------- //

namespace arzh
{
	template<class T>
		class t_Spectrum
		{
			public:
				t_Spectrum()
#ifdef USE_SSE_IN_SPECTRUM
				{ 
					this->m128 = _mm_set1_ps(0.0) ;
				}
#else
				{
					memset(_values,0,3*sizeof(T));
				}
#endif
#ifdef USE_SSE_IN_SPECTRUM
				inline t_Spectrum(const __m128& v)
				{
					m128 = v ;
				}
#endif
				inline t_Spectrum(T r, T g, T b)
				{
					_values[0] = r ;
					_values[1] = g ;
					_values[2] = b ;
#ifdef DEBUG
					//		CheckSpectrum(_values) ;
#endif
				}

				inline explicit t_Spectrum(T intensity) { _values[0]=_values[1]=_values[2]=intensity ; }

				void setRGB(T r,T g,T b)
				{
					_values[0] = r ;
					_values[1] = g ;
					_values[2] = b ;
				}

				inline T infNorm() const { return std::max(fabs(_values[0]),std::max(fabs(_values[1]),fabs(_values[2]))) ; }
				inline T square_L2_norm() const { return _values[0]*_values[0] + _values[1]*_values[1] + _values[2]*_values[2] ; }
				inline T getR() const { return _values[0] ; }
				inline T getG() const { return _values[1] ; }
				inline T getB() const { return _values[2] ; }
#ifdef A_VIRER
				inline void setR(T r) { _r = r ; }
				inline void setG(T g) { _g = g ; }
				inline void setB(T b) { _b = b ; }
#endif

				inline T operator[](unsigned int c) const { return _values[c] ; }
				inline T& operator[](unsigned int c) { return _values[c] ; }

				inline t_Spectrum  operator*(T a) const
				{
#ifdef USE_SSE_IN_SPECTRUM
					__m128 f = _mm_set1_ps(a) ;
					return t_Spectrum(_mm_mul_ps(f,m128)) ;
#else
					return t_Spectrum(a*_values[0],a*_values[1],a*_values[2]);
#endif
				}

				inline t_Spectrum  operator+(const t_Spectrum &s) const
				{
#ifdef USE_SSE_IN_SPECTRUM
					return t_Spectrum(_mm_add_ps(m128,s.m128)) ;
#else
					return t_Spectrum(_values[0]+s._values[0], _values[1]+s._values[1], _values[2]+s._values[2]);
#endif
				}

				inline t_Spectrum  operator*(const t_Spectrum &s) const
				{
#ifdef USE_SSE_IN_SPECTRUM
					return t_Spectrum(_mm_mul_ps(m128,s.m128)) ;
#else
					return t_Spectrum(_values[0]*s._values[0], _values[1]*s._values[1], _values[2]*s._values[2]);
#endif
				}

				inline t_Spectrum  operator-(const t_Spectrum &s) const
				{
#ifdef USE_SSE_IN_SPECTRUM
					return t_Spectrum(_mm_sub_ps(m128,s.m128)) ;
#else
					return t_Spectrum(_values[0]-s._values[0], _values[1]-s._values[1], _values[2]-s._values[2]);
#endif
				}

				inline t_Spectrum& operator+=(const t_Spectrum &s)
				{
#ifdef USE_SSE_IN_SPECTRUM
					m128 = _mm_add_ps(m128,s.m128) ;
#else
					_values[0] += s._values[0];
					_values[1] += s._values[1];
					_values[2] += s._values[2];
#endif
					return *this;
				}

				inline t_Spectrum& operator*=(const t_Spectrum &s)
				{
#ifdef USE_SSE_IN_SPECTRUM
					m128 = _mm_mul_ps(m128,s.m128) ;
#else
					_values[0] *= s._values[0];
					_values[1] *= s._values[1];
					_values[2] *= s._values[2];
#endif
					return *this;
				}

				inline t_Spectrum& operator*=(T a)
				{
#ifdef USE_SSE_IN_SPECTRUM
					__m128 f = _mm_set1_ps(a) ;
					m128 = _mm_mul_ps(m128,f) ;
#else
					_values[0] *= a ;
					_values[1] *= a ;
					_values[2] *= a ;
#endif
					return *this;
				}

				inline t_Spectrum& operator-=(const t_Spectrum& a)
				{
#ifdef USE_SSE_IN_SPECTRUM
					m128 = _mm_sub_ps(m128,s.m128) ;
#else
					_values[0] -= a._values[0];
					_values[1] -= a._values[1];
					_values[2] -= a._values[2];
#endif
					return *this;
				}

				inline t_Spectrum& operator/=(T a)
				{
#ifdef DEBUG
					if(a == 0.0)
						throw std::runtime_error("Division by zero") ;
#endif
#ifdef USE_SSE_IN_SPECTRUM
					__m128 f = _mm_set1_ps(a) ;
					m128 = _mm_div_ps(m128,f) ;
#else
					_values[0] /= a;
					_values[1] /= a;
					_values[2] /= a;
#endif
					return *this;
				}

				inline friend t_Spectrum operator*(T a, const t_Spectrum &s) { return s*a ; }
				inline friend t_Spectrum operator/(const t_Spectrum &s, T a)
				{ 
#ifdef DEBUG
					assert(a>0.0) ;
#endif
#ifdef USE_SSE_IN_SPECTRUM
					return t_Spectrum(_mm_div_ps(s.m128,_mm_set1_ps(a))) ;
#else
					return t_Spectrum(s._values[0]/a,s._values[1]/a,s._values[2]/a) ;
#endif
				}

				inline T intensity() const
				{
#ifdef DEBUG
					//		CheckSpectrum(_r,_g,_b) ;
#endif
#ifdef USE_SSE_IN_SPECTRUM
					__m128 fact = _mm_set_ps(0.30,0.59,0.11,0.0) ;
					__m128 res = _mm_mul_ps(fact,m128) ;
					template<class T>
						return ((T*)&res)[0]+((T*)&res)[1]+((T*)&res)[2] ;
#else
					return 0.30*_values[0] + 0.59*_values[1] + 0.11*_values[2];
#endif
				}

				friend std::ostream &operator<< (std::ostream& out,const t_Spectrum &s)
				{
					return out << s._values[0] << " " << s._values[1] << " " << s._values[2] ;
				}

			private:
#ifdef USE_SSE_IN_SPECTRUM
				union 
				{
					__m128  m128;
					T   _values[4];
				};
#else
				T _values[3] ;
#endif
		}; // interface of Spectrum

#ifdef DEBUG
#define CheckSpectrum(v) { { assert(v[0]>=0.0) ; assert(v[1]>=0.0) ; assert(v[2]>=0.0) ; } }
#endif

	typedef t_Spectrum<double> Spectrum ;
}

