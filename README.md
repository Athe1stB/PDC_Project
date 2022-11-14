# PDC_Project

## Sequential

**Compile**
	
	g++ sequential.cpp

**Run**
	
	./a.out

## Pthread
	
**Compile**
	
	g++ pthread_algo1.cpp -lpthread

**Run**

	./a.out

## OpenMP

**Compile**
	g++ openmp_algo1.cpp -fopenmp

**Run**
	./a.out


## MPI

**Install MPI**

	sudo apt-get mpich

**Compile**
	
	mpicc mpi.c -lm
	mpicc mpi_pthread.c -lm -lpthread
	mpicc mpi_openmp.c -lm -fopenmp

**Run**

	mpirun ./a.out

