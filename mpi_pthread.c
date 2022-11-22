#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>

// defining the dimensions of the matrix.
int n, n_total, m, m_total, valCount;

// defining the number of threads
int NUM_THREADS = 16;

// arrays to store the sparse matrix in CSR format [value, column and rowOffset]
float *value;
int *column, *rowOffset;

// store dense vector b and answer matrix.
float *b, *ans, *ansSequential;

// local variable for child processes
int index2;
float value2[150], b2[150];
int column2[150], rowOffset2[150];

// data structure to store all necessary data in a single entity
struct ds{
    long long thID;
    int index, n, n_total;
    float *value, *b, *ans;
    int *column, *rowOffset; 
};

// define lock
pthread_mutex_t lockk;

// constructor for the data structure to initialize with values
struct ds* construct(long long a, int b, int c, int d, float *x, float *y, int *xx, int *yy, float *zz)
{
    struct ds *s = malloc(sizeof (struct ds));
    s->thID = a;
    s->index = b;
    s->n = c;
    s->n_total = d;
    s->value = x;
    s->b = y;
    s->column = xx;
    s->rowOffset = yy;
    s->ans = zz;
    return s;
}

// this function will sum up the answer for lth row to (r-1)th element of matrix 1.
void *calculateRowsRange(void *arg)
{
    struct ds *d = ((struct ds*)arg);
    
    int valSize = d->n;
    long long thID = d->thID;
    
    int l = (thID*valSize)/NUM_THREADS;
    int r = ((thID+1)*valSize)/NUM_THREADS;
    
    // for each element of matrix 1 in range l to r
    for(int i=l; i<r; i++)
    {
        int c = d->column[i];
        int low = 0, high = d->n_total+1;

        // here we find the index of the row in which the element lies
        while(high>low+1)
        {
            int mid = low + (high-low)/2;

            if(d->rowOffset[mid]<=i+1)
              low = mid;
            else
              high = mid;
        }
        
        int r = low-1;
        
        pthread_mutex_lock(&lockk);
        d->ans[d->index + r]+= d->value[i] * d->b[c];
        pthread_mutex_unlock(&lockk);
    }
}

// function to calculate sequential answer to check for correctness
void sequential()
{
    ansSequential = (float*)calloc(n, sizeof(float));
    for(int i=0; i<n; i++)
    {
        ansSequential[i] = 0;
        for(int k=rowOffset[i]; k<rowOffset[i+1]; k++)
        {
            int index = column[k];
            ansSequential[i]+= (value[k]*b[index]);
        }
    }
}

// function to calculate the partial answer in parallel
void calculatePartialAnswerParallel(int index, int elements_per_process, int n, float *value, float *b, int *column, int *rowOffset, float *ans)
{
    pthread_t th[NUM_THREADS];
    
    // creating threads
    for(int i=0; i<NUM_THREADS; i++)
    {
        struct ds* d = construct(i, index, elements_per_process, n, value, b, column, rowOffset, ans);
        pthread_create(&th[i], NULL, calculateRowsRange, (void *)d);
    }
        
    // joining threads
    for(int i=0; i<NUM_THREADS; i++)
        pthread_join(th[i], NULL);
}

int main()
{
    int pid, np, elements_per_process, n_elements_recieved;
    // np -> no. of processes
    // pid -> process id

    MPI_Status status;

    // Creation of parallel processes
    MPI_Init(NULL, NULL);
    
    //inititalize lock
    pthread_mutex_init(&lockk, NULL);

    // find out process ID, and how many processes were started
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    // master process
    if(pid == 0) 
    {
        /*** read the sparse matrix from file ***/
        FILE *inputMatrix;
        inputMatrix = fopen("inputfile.mtx", "r");
        
        int num_lines;
        fscanf(inputMatrix, "%d", &n);
        fscanf(inputMatrix, "%d", &m);
        fscanf(inputMatrix, "%d", &num_lines);
        valCount = num_lines;
        
        // initialize all the matrices with the input dimension
        ans = (float*)calloc(n, sizeof(float));
        value = (float*)calloc(valCount, sizeof(float));
        column = (int*)calloc(valCount, sizeof(int));
        rowOffset = (int*)calloc((n+3), sizeof(int));
        b = (float*)calloc(m, sizeof(float));
        
        for(int l = 0; l < valCount; l++)
        {
            float val;
            int r, c;
            
            fscanf(inputMatrix, "%d", &r);
            fscanf(inputMatrix, "%d", &c);
            fscanf(inputMatrix, "%f", &val);
            
            //storing in csr format
            value[l] = val;
            column[l] = (c-1);
            rowOffset[r]++;
        }
        
        printf("The number of rows are: %d\nThe number of Columns are: %d\n\n", n, m);
        
        fclose(inputMatrix);
        
        // make rowOffset by prefix sum
        for(int i=2; i<n+1; i++)
            rowOffset[i]+=rowOffset[i-1];
        
        rowOffset[n+1] = rowOffset[n+2] = n+10;
        
        /*** read the vector from file ***/
        FILE *inputVector;
        inputVector = fopen("vector.txt", "r");

        if(inputVector == NULL)
        {
            printf("Error Reading input vector\n");
            exit (0);
        }

        for(int i=0; i<m; i++)
            fscanf(inputVector, "%f, ", &b[i]);
        
        fclose(inputVector);
        
        int index, i;
        elements_per_process = n / np;
        
        double startTime = MPI_Wtime();
        // check if more than 1 processes are running
        if(np > 1)
        {
            // distributes the portion of array to child processes to calculate their partial sums
            for(i = 1; i<np-1; i++) 
            {
                // starting index of the array sending to child process
                index = i * elements_per_process;

                // send all the required data to child processes to calculate partial answer 
                MPI_Send(&elements_per_process, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&index, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&n, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&m, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&rowOffset[0], n+3, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&column[index], elements_per_process, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&value[index], elements_per_process, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&b[0], m, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
            }

            // last process adds remaining elements
            index = i * elements_per_process;
            int elements_left = valCount - index;

            MPI_Send(&elements_left, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&index, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&n, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&m, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&rowOffset[0], n+3, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&column[index], elements_left, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&value[index], elements_left, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&b[0], m, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
        }
            
        // master process operates on its own partial answer in parallel
        calculatePartialAnswerParallel(0, elements_per_process, n, value, b, column, rowOffset, ans);

        // collects partial answers from other processes
        float tmp[n];
        
        for(i = 1; i < np; i++)
        {
            MPI_Recv(&tmp, n, MPI_FLOAT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
            
            // add partial answer to final ans
            for(int j=0; j<n; j++)
                ans[j]+=tmp[j];
                
            int sender = status.MPI_SOURCE;
        }
        
        double endTime = MPI_Wtime();
        
        // output the dense vector
        printf("The input dense vector is: \n");
        for(int i=0; i<m; i++)printf("%f ",b[i]); printf("\n\n");
        
        // output final answer using distributed computing
        printf("The final answer is: \n");
        for(int i=0; i<n; i++)printf("%f ", ans[i]); printf("\n\n");
        
        // calculate answer using sequential algorithm and print it
        sequential();
        printf("The sequential answer is: \n");
        for(int i=0; i<n; i++)printf("%f ", ansSequential[i]); printf("\n\n");
        
        // check for correctness using sequential answer
        bool correct = true;
        for(int i=0; i<n; i++)
        {
            if(ans[i]!=ansSequential[i])
                correct = false;
        }
        
        printf((correct?"***Correct***\n\n":"***InCorrect***\n\n"));
        printf("Time taken: %f seconds\n\n", endTime-startTime);
    }
    // slave processes
    else
    {
        // recieves all the required data in local arrays and variables to calculate partial answer
        MPI_Recv(&n_elements_recieved, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&index2, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&n_total, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&m_total, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&rowOffset2, n_total+3, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&column2, n_elements_recieved, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&value2, n_elements_recieved, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&b2, m_total, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status);

        // calculates partial answer for the current child process
        float* partial_ans = (float*)calloc(n_total, sizeof(float));
        
        // calculate partial answer in parallel using pthread
        calculatePartialAnswerParallel(index2, n_elements_recieved, n_total, value2, b2, column2, rowOffset2, partial_ans);
        
        // sends the partial answe to the root process
        MPI_Send(&partial_ans[0], n_total, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
    }

    // cleans up all MPI state before exit of process
    MPI_Finalize();
    
    return 0;
}