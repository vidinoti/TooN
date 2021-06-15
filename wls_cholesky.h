/*                       
	Copyright (C) 2005 Tom Drummond

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	This library is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
	Lesser General Public License for more details.

	You should have received a copy of the GNU Lesser General Public
	License along with this library; if not, write to the Free Software
	Foundation, Inc.
     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
//-*- c++ -*-
// A WLS class using Cholesky decomposition and sparse JtJ
// Also stores the sum squared error and can compute the residual

#ifndef CVD_WLS_CHOLESKY_H
#define CVD_WLS_CHOLESKY_H

#include "VDARMatrix.h"

/// Performs weighted least squares using Cholesky decomposition and sparse JtJ.
/// Much faster (but less robust) than the standard WLS.
/// Also stores the sum squares error and can compute the residual.
/// @param The number of parameters in the system
/// @ingroup gEquations
template <int Size, class TYPE = float>
class WLSCholesky {
public:
  /// Default constructor
  WLSCholesky(){clear();}
  /// Construct using a given regularisation prior
  WLSCholesky(TYPE prior){clear(prior);}

  /// Clear all the measurements and apply a constant regularisation term. 
  /// Equates to a prior that says all the parameters are zero with \f$\sigma^2 = \frac{1}{\text{val}}\f$.
  /// @param prior The strength of the prior
    void clear(TYPE prior){
        for(int i = 0; i < Size; ++i){
            my_C_inv(i,i) = prior;
            for(int j = i + 1; j < Size; ++j)
                my_C_inv(i,j) = TYPE(0);
        }
        
        
		my_vector.clear();
        //my_err = TYPE(0);
    }
	
	void clear(){
        
        my_C_inv.clear();
	        
		my_vector.clear();
        //my_err = TYPE(0);
    }
	
	/// Add a single measurement 
	/// @param m The value of the measurement
	/// @param J The Jacobian for the measurement \f$\frac{\partial\text{m}}{\partial\text{param}_i}\f$
	/// @param weight The inverse variance of the measurement (default = 1)
    inline void add_df(const TYPE m, const VDARVector<Size, TYPE> & J) {
        //VDARVector<Size> Jw = J;
        for(unsigned int i=0; i<(Size); i++){
            for(unsigned int j=i; j<(Size); j++){
	            my_C_inv(i,j)+=J(i)*J(j);
            }
            my_vector(i)+=J(i)*m;
        }
        //my_err+=m*weight*m;
    }
  
    void compute(){
        // Homegrown Cholesky
        VDARMatrix<Size, Size, TYPE> L;
        VDARVector<Size, TYPE> invLdiag;
        for(unsigned int i=0;i<Size;i++) {
			
            TYPE a=my_C_inv(i,i);
			
            for(unsigned int k=0;k<i;k++) 
				a-=L(k,i)*L(k,i);
			
            invLdiag(i) = TYPE(1)/sqrtf(a); 
			
            L(i,i)= a * invLdiag(i);
			
            for(unsigned int j=i;j<Size;j++) {
				
	            a=my_C_inv(i,j);
				
	            for(unsigned int k=0;k<i;k++) a-=L(k,j)*L(k,i);
	                L(i,j)=a * invLdiag(i);
            }
        }
		
        VDARVector<Size, TYPE> y;
		
        for(unsigned int i=0;i<Size;i++) {
            TYPE a=my_vector(i);
            for(unsigned int j=0;j<i;j++) a-=L(j,i)*y(j);
            y(i)=a*invLdiag(i);
        }
        for(int i=Size-1;i>-1;i--) {
            TYPE a=y(i);
            for(int j=i+1;j<Size;j++) a-=L(i,j)*my_mu(j);
            my_mu(i)=a*invLdiag(i);
        } 

    }


    /// Returns the inverse covariance matrix
    VDARMatrix<Size,Size,TYPE>& get_C_inv() {return my_C_inv;}
  
    /// Returns the inverse covariance matrix
    const VDARMatrix<Size,Size,TYPE>& get_C_inv() const {return my_C_inv;}
  
    VDARVector<Size, TYPE> & get_mu(){return my_mu;}
    const VDARVector<Size, TYPE>& get_mu() const {return my_mu;}

private:
  VDARVector<Size, TYPE> my_mu;
  VDARMatrix<Size,Size, TYPE> my_C_inv;
  VDARVector<Size, TYPE> my_vector;
  //TYPE my_err;    // error before optimisation
};

#endif
