#include<bits/stdc++.h>
#include<chrono>
using namespace std;

#define nline "\n"

// defining the number of threads
int NUM_THREADS = 1;

// defining the dimensions of the matrix.
int n, m, valCount;

// arrays to store the sparse matrix in CSR format [value, column and rowOffset]
vector<double> value;
vector<int> column, row, rowOffset;

// store matrices a, b and answer matrix.
vector<double> b, ans, ansSequential;

// define lock
pthread_mutex_t lockk, barrierLock;

// define barrier
pthread_barrier_t barrier; 

void sequential()
{
    ansSequential.assign(n,0);
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

struct data{
    int thID;
    int pow;
};

// this function will calculate prefix sum using hillis algorithm
void *hillis(void *arg)
{
    struct data *d = ((struct data*)arg);
    int tid = d->thID, pow = d->pow;
    
    rowOffset[tid+pow]+= row[tid];
}

// this function will sum up the answer for lth row to (r-1)th element of matrix 1.
void *calculateRowsRange(void *threadID)
{
    long long thID = (long long)threadID;
    
    int valSize = value.size();
    
    int l = (thID*valSize)/NUM_THREADS;
    int r = ((thID+1)*valSize)/NUM_THREADS;
    
    // for each element of matrix 1 in range l to r
    for(int i=l; i<r; i++)
    {
        int c = column[i];
        auto rp = lower_bound(rowOffset.begin(), rowOffset.end(), i+1) - rowOffset.begin();
        int r = rp-1;
        
        // critical section
        pthread_mutex_lock(&lockk);
        ans[r]+= value[i] * b[c];
        pthread_mutex_unlock(&lockk);
    }
}

void solve(int threads)
{    
    NUM_THREADS = threads;
    
    //initialize all variables
    value.clear();
    rowOffset.clear();
    column.clear();
    b.clear();
    ans.clear();
    
    //inititalize lock
    pthread_mutex_init(&lockk, NULL);
    pthread_mutex_init(&barrierLock, NULL);
    
    // read the sparse matrix from file
    ifstream inputMatrix("inputfile.mtx");
    int num_row, num_col, num_lines;

    while (inputMatrix.peek() == '%') inputMatrix.ignore(2048, '\n');
    inputMatrix >> num_row>> num_col >> num_lines;    
    n = num_row;
    m = num_col;
    valCount = num_lines;
    
    // initialize all arrays
    ans.assign(n,0.0);
    value.assign(num_lines,0);
    column.assign(num_lines,0);
    row.assign(n+1,0);    
    
    for(int l = 0; l < num_lines; l++)
    {
        double val;
        int r, c;
        inputMatrix >> r >> c >> val;
        
        //storing in csr format
        value[l] = val;
        column[l] = (c-1);
        row[r]++;
    }
    
    cout<<"The number of rows are: "<<num_row<<"\nThe number of columns are: "<<num_col<<"\nThe number of threads are: "<<NUM_THREADS<<nline;

    inputMatrix.close();
    
    int pow = 1;
    
    rowOffset = row;
    
    // for each stride we are working on it
    while(pow<n+1)
    {
        int n_threads = n+1-pow;
        
        // make rowOffset by prefix sum using Hillis Algorithm        
        pthread_t ht[n_threads];
        
        for(int i=0; i<n_threads; i++)
        {
            struct data *s = (struct data*)malloc(sizeof (struct data));
            s->thID = i;
            s->pow = pow;
            pthread_create(&ht[i], NULL, hillis, (void *)s);
        }
        
        for(int i=0; i<n_threads; i++)
            pthread_join(ht[i], NULL);
        
        row = rowOffset;
        
        pow*=2;
    }
    
    // read the vector from file
    string s, temp="";
    ifstream inputVector("vector.txt");
    while(inputVector>>s)
    {
        temp = "";
        if(s.back()==',')
            s.pop_back();
        int number = stoi(s);
        double xx = number*1.0;
        b.push_back(xx);
    }
    
    inputVector.close();
    
    auto start = chrono::high_resolution_clock::now();
    pthread_t th[NUM_THREADS];
    
    // creating threads
    for(int i=0; i<NUM_THREADS; i++)
        pthread_create(&th[i], NULL, calculateRowsRange, (void *)i);
        
    // joining threads
    for(int i=0; i<NUM_THREADS; i++)
        pthread_join(th[i], NULL);       
             
    auto end = chrono::high_resolution_clock::now();
          
    // Output the dense vector
    cout<<nline<<"The input dense vector is: "<<nline;
    for(auto it:b)
        cout<<it<<" "; cout<<nline<<nline;
    
    // Output the final answer using pthread
    cout<<"The final answer is: "<<nline;
    for(int i=0; i<n; i++)cout<<ans[i]<<" "; cout<<nline<<nline;
        
    // calculate sequential answer and print it
    sequential();
    cout<<"The sequential answer is: "<<nline;
    for(int i=0; i<n; i++)cout<<ans[i]<<" "; cout<<nline<<nline;
        
    // check for correctness
    cout<<((ans==ansSequential)?"Correct":"InCorrect")<<nline;
    
    // Calculating total time taken by the program.
    double time_taken = 
      chrono::duration_cast<chrono::nanoseconds>(end - start).count();
  
    time_taken *= 1e-9;
  
    cout << "Time taken by program using "<< NUM_THREADS <<" threads is " << fixed 
         << time_taken << setprecision(9);
    cout << " sec" << nline;
}

int main()
{
    ios_base::sync_with_stdio(false);    cin.tie(NULL);    cout.tie(NULL);
    
    // run for threads
    solve(4);
    
    return 0;
}
