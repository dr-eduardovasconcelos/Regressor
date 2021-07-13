/*
	This file is part of Regressor.

    Regressor is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Regressor is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Regressor.  If not, see <https://www.gnu.org/licenses/>
*/

#include "completepolynomialregression128.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

long double *qrbetascalculator(long double *T, double* Y, int rows, int columnsT);
long double hypotlong(long double a, long double b);
long double pow2(long double a, int e);
long double sqrt2(long double x);


/*
* This function is responsible for create the matrix M, used to perform the regression.
*/

long double* performRegression(double *Y, double *X, int degree, int rows, int XColumns)
{
	
	/*
	* this variable is basically the number of columns of columns of X, 
	* that represents the indent 
	*/
	int variables = XColumns;
	
	/*
	* Adjust the polynomial degree, since the number columns of M cannot be greater
	* than the number of rows
	*/
	
    while (pow(degree + 1, variables) > rows)
    {
        degree--;
    }

	/*
	* The number of columns of M is (degree+1)^variables
	*/
    int columns = pow(degree + 1, variables);
    
    #ifndef NOVERBOSE
    	printf("Creating Regression Matrix \n");
	#endif
	
	/*
	* rMatrix is the matrix used for regression. It´s a transposed vertion of M.
	* We have used a transposed metrix to increase the algorithm performance,
	* since, Jama matrix implementation makes rows iterations to solve the QR
	* factorization. 
	*/
    long double* rMatrix = (long double*) malloc (columns * rows * sizeof(long double));

    if(rMatrix == NULL){
        printf("rMatrix null \n");

        exit(1);
    }

    for (int i = 0; i < rows; i++)
    {
    	#ifndef NOVERBOSE
			printf("%d from %d lines \n",i,rows);
		#endif
		
        *(rMatrix + i) = 1.0;

		/*
		* This code calculate each term of the polynomial.
		* The term structure is: beta_i x prod{X_i,j}
		*/
        for (int j = 1; j < columns; j++)
        {
            long double product = 1.0;
            
            /*
            * By using limit variable, it's possible to decrease the complexity of building
			* M. As the value of M_i,j = prod{X_i,k} com k < ln(j)/ln(degree+1). 
            */
            int limit = (int)(log(j) / log(degree+1)) + 1;

            for (int k = 0; k < limit; k++)
            {
                product *= pow2((long double) *(X + i*variables + k) , (((int) ((j) / (pow(degree + 1, k)))) % (degree + 1)));
            }

			/*
			* rMatrix is a transposed of M, so M_i,j = rMatrix_j,i with M beeing a rows x columns matrix and
			* rMatrix a columns x rows.
			*/
            *(rMatrix + j * rows + i) = product;
            
        }

    }

    return qrbetascalculator(rMatrix, Y, rows, columns);
    
}

/*
* This function is a modification of Jama System QR solving.
*
* In this function, matrix T is a transposed version of M, so, the modification
* consists in making iterations on T by accessing its columns instead of its rows.
*
* Access https://math.nist.gov/javanumerics/jama/ for more information about Jama project.
*/
long double *qrbetascalculator(long double *T, double* Y, int rows, int columnsT){

    long double* diagonal = (long double*)malloc(columnsT * sizeof(long double));
    
    #ifndef NOVERBOSE
    	printf("Making QR factorization of M \n");
    #endif

    for(int k = 0; k < columnsT; k++){
    	
    	#ifndef NOVERBOSE
    		printf("%d from %d columns \n",k,columnsT);
        #endif
        
        long double norm = 0.0;

        for(int i = k; i < rows; i++){
            norm = hypotlong(norm, *(T + k*rows +i));
            
        }

        if(norm != 0){
            if(*(T + k*rows + k) < 0){
                norm = -norm;
            }

            for(int i = k; i < rows; i++){
                *(T + k*rows +i) = *(T + k*rows +i)/norm;
            }

            *(T + k*rows + k) = *(T + k*rows + k) + 1.0; 

            for(int j = k+1; j < columnsT; j++){
                long double s = 0.0;

                for(int i = k; i < rows; i++){
                    s += *(T + k*rows +i) * *(T + j*rows +i);
                }

                s = -s/(*(T + k*rows + k));

                for(int i = k; i < rows; i++){
                    *(T + j*rows +i) = *(T + j*rows +i) + (s * *(T + k*rows +i));
                }
            }
        }

        *(diagonal + k) = -norm;
    }

    long double* B = (long double*) malloc(rows * sizeof(long double));

    for(int i = 0; i < rows; i++){
        *(B + i) = (long double)*(Y + i);
    }

	#ifndef NOVERBOSE
    	printf("Solving Linear System \n");
    #endif

    for(int k = 0; k < columnsT; k++){
    	
        long double s = 0;

        for(int i = k; i < rows; i++){
            s += *(T + k*rows +i) * *(B + i);
        }

        s = -s/(*(T + k*rows + k));

        for(int i = k; i < rows; i++){
            *(B + i) = *(B + i) + (s * *(T + k*rows +i));
        }
    }

    for(int k = columnsT - 1; k >= 0; k--){
        *(B + k) = *(B + k) / (*(diagonal + k));

        for(int i = 0; i < k; i++){
            *(B + i) = *(B + i) - (*(B + k) * *(T + k*rows +i));
        }
    }
    
    #ifndef NOVERBOSE
    	printf("Creating Coefficients \n");
    #endif
    
    long double* betas = (long double*)malloc(columnsT * sizeof(long double));
    
    for(int i = 0; i < columnsT; i++){    	
    	
    	*(betas + i) = *(B + i);
	}

    return betas;

}

/*
* This funcion represents the long double ajustment for the 
* hipotenuse funcion.
*/
long double hypotlong(long double a, long double b){
	return (long double) sqrt((double)(a*a + b*b));
}

/*
* This funcion is a simple implementation of power funcion
* with long double  precision.
*/
long double pow2(long double a, int e){
	long double product = 1;
	for(int i = 0; i < e; i++){
		product *= a;
	}
	
	return product;
}
