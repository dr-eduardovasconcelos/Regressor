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
#include <stdlib.h>
#include <stdio.h>

void cleanBuff(char* buff, int buffersize);

/*
* This is an auxuliary procedure, that reads a specific format of csv value.
* This procedure accept only file with columns separated with "," and separate
*
* X and Y are regular arrays passed as reference, rows and columns are regular 
* integers passed as reference
*/
void getDataFromCSV(char* path, double** X, double** Y, int* rows, int* columns, int buffersize){
	
	FILE* file = fopen(path, "r");
	
	/*
	* The following procedure is used to count the number of rows and columns.
	*/
	
	char curr;
	(*columns) = 0;
	(*rows) = 1;
	
	do{
		curr = getc(file);
		
		if(curr == ',')
			++(*columns);
	}while(curr!='\n');
	
	do{
		
		curr = getc(file);
		
		if(curr == '\n')
			++(*rows);
		
	}while(curr!=EOF);
	
	fclose(file);
	
	/*
	* The following procedure is used for fill matrices Y and X.
	*/
	
	char buff[buffersize];
	
	curr = ' ';
	
	file = fopen(path, "r");
	
	*Y = (double*) malloc((*rows)*sizeof(double));
	int yIndex = 0;
	*X = (double*) malloc((*columns)*(*rows)*sizeof(double));
	int xRowIndex = 0;
	int xColumnIndex = 0;
	
	for(xRowIndex = 0; xRowIndex < (*rows); xRowIndex++){
		
		int  i = 0;
		
		curr = getc(file);
		
		cleanBuff(buff, buffersize);
		
		/*
		* This loop runs through the string by adding caracteres into the buffer variable.
		*/
		while(curr != ','){
			
			if(curr == '\n'){
				printf("problem with the CSV file");
				exit(1);
			}
			
			buff[i++] = curr;
			curr = getc(file);
			
		}
		
		double value;
		
		sscanf(buff, "%lf", &value);
		
		*(*Y + yIndex++) = value;
		
		while(curr != '\n'){
			
			i = 0;
			
			cleanBuff(buff, buffersize);
			
			curr = getc(file);
			
			while(curr != ',' && curr != '\n'){
				buff[i++] = curr;
				curr = getc(file);
			}
		
			sscanf(buff, "%lf", &value);
			
			*(*X + xRowIndex*(*columns) + xColumnIndex++) = value;
			
		}
		
		xColumnIndex = 0;
		
	};
	
}

/*
* This procedure is used to remove residues of buffer.
*/
void cleanBuff(char* buff, int buffersize){
	for(int i = 0; i < buffersize; i++){
		
		buff[i] = 0;
		
	}
}

