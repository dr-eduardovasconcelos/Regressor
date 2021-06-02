# Regressor

How do I execute Regressor and Regressor Polynomial?

- Open the command line (CMD or Terminal) and access the folder of the program.
- type: "Regressor.exe csvpath degree buffersize"
-- exemple: Regressor.exe C:\covid_9_var.csv 2 20
--- do the same for RegressorPolynomial.exe
- argument csvpath is obrigatory.
- arguments degree and buffersize are both optional. If not specified, these argments will be defined as 1 and 20, respectivelly
- if your CSV file has large numbers, just increase the buffersize.

How do I Compile Regressor and Regressor Polynomial?

- first, compile the libraries: gcc.exe -c -g completepolynomialregression.c -std=c11
-- gcc.exe -c -g csvReader.c -std=c11
- after, compile the program: gcc.exe Regressor.c -o Regressor csvReader.o completepolynomialregression.o -std=c11
-- do the same for Regressor Polynomial. Note that, for Polynomial Regressor, the correct library is polynomialregression.c
--- if you are using linux, use -lm to compile the programs.
---- example: gcc -c -g completepolynomialregression.c -std=c11 -lm

for more information about combinatorial regression, please, access the manuscript: https://arxiv.org/abs/2009.12386
