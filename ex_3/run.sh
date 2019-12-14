# ex_3_1
if [ $1 == "1" ]
then
	echo "ex_3_1 running..."
	mpirun -np 3 build/ex_3_1
# ex_3_1_omp
elif [ $1 == "1omp" ]
then
	echo "ex_3_1_omp running..."
	mpirun -np 3 build/ex_3_1_omp
# ex_3_2
elif [ $1 == "2" ]
then
	echo "ex_3_2 running..."
	mpirun -np 2 build/ex_3_2
# error
else
	exit 2
fi
