prog: main.cpp initial_data_and_functions.cpp integrals.cpp 
	gcc -o prog main.cpp initial_data_and_functions.cpp integrals.cpp -I. -lm -Wno-write-strings
