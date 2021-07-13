# Regressor

How do I execute Regressor and Regressor Polynomial?

- Open the command line (CMD or Terminal) and access the folder of the program.
- type: "Regressor.exe csvpath degree buffersize"
-- exemple: Regressor.exe C:\covid_9_var.csv 2 20
- do the same for RegressorPolynomial.exe
- argument csvpath is obrigatory.
- arguments degree and buffersize are both optional. If not specified, these argments will be defined as 1 and 20, respectivelly
- if your CSV file has large numbers, just increase the buffersize.

How do I Compile Regressor and Regressor Polynomial?

- first, compile the libraries: gcc.exe -c -g completepolynomialregression.c -std=c11

  gcc.exe -c -g csvReader.c -std=c11

- after, compile the program: gcc.exe Regressor.c -o Regressor csvReader.o completepolynomialregression.o -std=c11

  do the same for Regressor Polynomial. Note that, for Polynomial Regressor, the correct library is polynomialregression.c

  if you are using linux, use -lm to compile the programs.

  example: gcc -c -g completepolynomialregression.c -std=c11 -lm
  
- Do the same to use and compile Regressor128


for more information about combinatorial regression, please, access the manuscript: https://arxiv.org/abs/2009.12386

--- Notes about v1.1

This version contains some new files that will help you to execute your models.

First: Regressor now generate a binary file containing the beta coeficients of your model. There is also a html file able to read this file and generate responses for these models.

Second: We release a 128 bit program to generate models from bit data. Besides Regressor128, we have developed another program called read128BinaryModel, since, javascript hasn´t a native structure to represent 128 bit floats.

to use program read128BinaryModel.exe, place the binary file generated by Regressor128.exe in the same folder and create a file called inputs.txt. Inside inputs.txt, insert the predictor values separated by comma.
