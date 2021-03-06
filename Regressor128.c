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

#include "csvReader.h"
#include "completepolynomialregression128.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/*
* This program receives arguments, and the first one is obrigatory.
*  
* This program generates two files, the first is a .txt with the RSquared value
* and the regression formula, to be used in a spreadsheet tool,
* and the second is a .csv file with the residues obtined from the model.
*/
int main(int argc, char* argv[]){
	
	if(argc < 2){
		
		printf("Required arguments are missing! \n");
		printf("The correct form to start Regressor is: Regressor csvpath degree buffer_size \n");
		printf("You have to pass at last the csvpath argument \n");
		printf("if you don't define the degree argument, Regressor will consider the value 1 \n");
		printf("if you don't define the buffer_size, Regressor will consider the value 20'\n");
		
		exit(0);
	}
	
	int degree = 1;
	int buffer = 20;
	
	if(argc > 2){
	
		sscanf(argv[2], "%d", &degree);
	
		if(argc > 3)
			sscanf(argv[3], "%d", &buffer);
	}
	
	char* path = argv[1];
	
	double* X;
	double* Y;
	
	int rows;
	int variables;
	
	/*
	* Retrieving data from CSV file. X, y, rows and columns are passed as reference
	* parameters.
	*/
	getDataFromCSV(path, &X, &Y, &rows, &variables, buffer);
	
	/*
	* retrieving the beta coeficients
	*/
	long double* B = performRegression(Y, X, degree, rows, variables);
    
    while (pow(degree + 1, variables) > rows)
    {
        degree--;
    }
    
    int columns = pow(degree + 1, variables);

	/*
	* The following procedure is used to compute predicted array.
	* This array is used to compute the RSquared value and to genetate the 
	* residues.csv file.
	*/
    long double* predicted = (long double*) malloc(rows * sizeof(long double));

	/*
	* This variable is used to calculate the total variation of array Y.
	*/
    double averageY = 0.0; 

    for (int i = 0; i < rows; i++)
    {
        averageY += *(Y + i);

		/*
		* Adding beta_0 to the polynomial
		*/
        *(predicted + i) = *B;

		/*
		* This code calculate each term of the polynomial.
		* The term structure is: beta_i x prod{X_i,j}
		*/
        for(int j = 1; j < columns; j++){
            long double product = 1.0;
            
             /*
            * By using limit variable, it's possible to reduce the complexity of building
			* M. As the value of M_i,j = prod{X_i,k} com k < ln(j)/ln(degree+1). 
            */
            int limit = (int)(log(j)/log(degree+1))+1;

            for (int k = 0; k < limit; k++)
            {
                product *= (long double)pow( *(X + i*variables+ k) , (((int)((j) / (pow(degree + 1, k)))) % (degree + 1)));
            }

            *(predicted + i) += (product * *(B + j));
        }
    }

    averageY /= rows;

    double sstot = 0.0;
    double ssres = 0.0;
    double sqe = 0.0;

    for (int i = 0; i < rows; i++)
    {
        sstot += pow((*(Y + i)) - averageY,2);
        ssres += pow((*(Y + i)) - (double) *(predicted + i),2);
    }
    
    double RSquared = 1.0 - ssres/sstot;
    
    /*
    * Creating the spreadsheet formated formula and inserting it into a txt file.
    */
    printf("RSquared value is: %lf \n", RSquared);
    
    FILE * file;

    file = fopen("spreadsheet_equation.txt","w");

    fprintf(file, "this equation has a RSquared value of %lf \n", RSquared);
    fprintf(file, "put the following equation on cel A%i \n", variables + 1);

    fprintf(file,"= %G + ",(double) (*B));

    for(int i = 1; i < columns; i++){

        int j;
        
        int limit = (int)(log(i)/log(degree+1))+1;
        
        fprintf(file,"(%G)*",(double) *(B + i));
        
        for(j = 0; j < limit; j++) {

            int exp = ( ((int)((i) / (pow(degree + 1, j)))) % (degree + 1));

            if(exp == 0) continue;

            if(exp == 1)
                fprintf(file, "(A%d)", j);
            else
                fprintf(file, "(A%d^%d)", j, exp);

            if(j < limit - 1)
                fprintf(file, "*");

        }
        if(i < columns - 1)
            fprintf(file, " + ");
        
    }

    fclose(file);
    
    /*
    * Writes a binary file with the coeficients. 
    * This file can be useful when the equation is large.
    * In some cases, with large equations, some coeficientes written into file, spreadsheet_equation.txt,
    * can suffer from loss of precision, and these equations can return imprecise outcomes.
    */
    
    file = fopen("binary_coeficients.bin","wb");
    
    fwrite(&variables, sizeof(int), 1, file);
    fwrite(&degree, sizeof(int), 1, file);
    
    double* aux =  (double*) malloc(sizeof(double) * columns);
    
    for(int i = 0; i<columns; i++){
    	
    	*(aux + i) = (double) *(B + i);
    	
	}
	fwrite(&(*aux), sizeof(double) * columns, 1, file);
    
    fclose(file);
    
    /*
    * Writing a binary file with 128 bits precision.
    * In some cases, coeficients are large enouth to fit into a
    * 64 bits variable. 
    * This pack contains a program entitled read128BinaryModel,
    * used to execute models generated by this program.
    */
    
    file = fopen("binary_coeficients128.bin","wb");
    
    fwrite(&variables, sizeof(int), 1, file);
    fwrite(&degree, sizeof(int), 1, file);
    
   
	fwrite(&(*B), sizeof(long double) * columns, 1, file);
    
    fclose(file);
    
    /*
    * Creating the residues file.
    */
    file = fopen("residues.csv","w");
    
    fprintf(file, "Y,predicted,residuo\n");
    
    for(int i = 0; i < rows; i++){ 
    	fprintf(file, "%lf,%lf,%lf\n", *(Y + i), (double)*(predicted + i), *(Y + i) - (double) *(predicted + i));
	}
	
	fclose(file);
	
}
