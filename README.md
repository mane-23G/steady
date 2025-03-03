# steady
This is a parallel algorithm that computes the temperature in a steady state of the point specifed by the user. 

### Build with
mpicc -Wall -o steady steady.c

### Usage 
mpirun \[MPI options\] steady \<plate\> \<row\> \<column\>

The user must supply a pathname of a text file that consists of six numbers in the following order, separated by
white space:
    – an integer representing the number of rows in the grid
    – an integer representing the number of columns in the grid
    – a double precision float representing the temperature in Celsius on the north edge
    – a double precision float representing the temperature in Celsius on the east edge
    – a double precision float representing the temperature in Celsius on the south edge
    – a double precision float representing the temperature in Celsius on the west edge
The second and third are a pair of non-negative intgers representing the row and column of the grid.

### Output and Features
The program outputs the temperature at the specefied point. The steady state occurs when the maximum change is less than the convergence threshold.

### Defects/Shortcomings
-The programs does not produce the correct temperature all the time, the reason has yet to be found

Thank you to Professor Weiss for the assigment!
