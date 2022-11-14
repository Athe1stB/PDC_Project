# PDC_Project

##MPI

**Install MPI**

	sudo apt-get mpich

**Compile**
	
	mpicc mpi.c -lm
	mpicc mpi_pthread.c -lm -lpthread
	mpicc mpi_openmp.c -lm -fopenmp

**Run**

	mpirun ./a.out

