# ex_1_a_1
if [ $1 == "1" ]
then
	echo "ex_3_1 running..."
	mpirun -np 3 build/ex_3_1
# ex_1_a_2
elif [ $1 == "a2" ]
then
	echo "ex_1_a_2 running..."
	mpirun -np 34 build/ex_1_a_2
# ex_1_b
elif [ $1 == "b" ]
then
	echo "ex_1_b running..."
	mpirun -np 9 build/ex_1_b
# ex_1_c
elif [ $1 == "c" ]
then
	echo "ex_1_c running..."
	mpirun -np 50 build/ex_1_c
# error
else
	exit 2
fi
