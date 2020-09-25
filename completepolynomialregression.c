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

#include "completepolynomialregression.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double *qrbetascalculator(double *T, double* Y, int rows, int columnsT);

/*
* This function is responsible for create the matrix M, used to perform the regression.
*/

double* performRegression(double *Y, double *X, int degree, int rows, int XColumns)
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
    double* rMatrix = (double*) malloc (columns * rows * sizeof(double));

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
            double product = 1.0;
            
            /*
            * By using limit variable, it's possible to decrease the complexity of building
			* M. As the value of M_i,j = prod{X_i,k} com k < ln(j)/ln(degree+1). 
            */
            int limit = (int)(log(j) / log(degree+1)) + 1;

            for (int k = 0; k < limit; k++)
            {
                product *= pow(*(X + i*variables + k) , (((int) ((j) / (pow(degree + 1, k)))) % (degree + 1)));
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
double *qrbetascalculator(double *T, double* Y, int rows, int columnsT){

    double* diagonal = (double*)malloc(columnsT * sizeof(double));
    
    #ifndef NOVERBOSE
    	printf("Making QR factorization of M \n");
    #endif

    for(int k = 0; k < columnsT; k++){
    	
    	#ifndef NOVERBOSE
    		printf("%d from %d columns \n",k,columnsT);
        #endif
        
        double norm = 0.0;

        for(int i = k; i < rows; i++){
            norm = hypot(norm, *(T + k*rows +i));
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
                double s = 0.0;

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

    double* B = (double*) malloc(rows * sizeof(double));

    for(int i = 0; i < rows; i++){
        *(B + i) = *(Y + i);
    }

	#ifndef NOVERBOSE
    	printf("Solving Linear System \n");
    #endif

    for(int k = 0; k < columnsT; k++){
    	
        double s = 0;

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
    
    double* betas = (double*)malloc(columnsT * sizeof(double));
    
    for(int i = 0; i < columnsT; i++){    	
    	
    	*(betas + i) = *(B + i);
	}

    return betas;

}
