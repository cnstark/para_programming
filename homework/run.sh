# homework1
if [ $1 == "1" ]
then
	echo "homework1 running..."
	mpirun -np 5 build/homework1
# homework2
elif [ $1 == "2" ]
then
	echo "homework2 running..."
	mpirun -np 5 build/homework2
# homework3
elif [ $1 == "3" ]
then
    echo "Copy build/homework3 to other nodes."
    scp build/homework3 pp11@node3:/home/pp11/SA19011084/homework/build/
    scp build/homework3 pp11@node4:/home/pp11/SA19011084/homework/build/

	echo "homework3 running..."
	mpirun -np 12 -f mpi_config build/homework3
# homework4
elif [ $1 == "4" ]
then
    echo "homework4 running..."
	mpirun -np 4 build/homework4
# error
else
	exit 2
fi
