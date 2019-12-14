echo "Build ex_3_1."
mpic++ -g -Wall -std=c++11 src/ex_3_1.cpp -o build/ex_3_1

echo "Build ex_3_1_omp."
mpic++ -g -Wall -std=c++11 src/ex_3_1_omp.cpp -o build/ex_3_1_omp -fopenmp

echo "Build ex_3_2."
mpicc src/ex_3_2.c -o build/ex_3_2
