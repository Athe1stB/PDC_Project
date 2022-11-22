#include<bits/stdc++.h>
#include<chrono>
using namespace std;

#define nline "\n"

// defining the dimensions of the matrix.
int n, m, valCount;

// arrays to store the sparse matrix in CSR format [value, column and rowOffset]
vector<double> value;
vector<int> column, rowOffset;

// store matrices a, b and answer matrix.
vector<double> b, ansSequential;

void sequential()
{
    ansSequential.assign(n,0);
    // calculate answer at ith index of answer array
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

void solve()
{
    
    // read the sparse matrix from file
    ifstream inputMatrix("inputfile.mtx");
    int num_row, num_col, num_lines;

    while (inputMatrix.peek() == '%') inputMatrix.ignore(2048, '\n');
    inputMatrix >> num_row>> num_col >> num_lines;    
    n = num_row;
    m = num_col;
    valCount = num_lines;
    
    // initialize all arrays
    value.assign(valCount,0);
    column.assign(valCount,0);
    rowOffset.assign(n+1,0);
    
    for(int l = 0; l < num_lines; l++)
    {
        double val;
        int r, c;
        inputMatrix >> r >> c >> val;
        
        //storing in csr format
        value[l] = val;
        column[l] = (c-1);
        rowOffset[r]++;
    }
    
    cout<<"The number of rows are: "<<num_row<<"\nThe number of columns are: "<<num_col<<nline;

    inputMatrix.close();
    
    // make rowOffset by prefix sum
    for(int i=2; i<n+1; i++)
        rowOffset[i]+=rowOffset[i-1];
    
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
    sequential();       
             
    auto end = chrono::high_resolution_clock::now();
          
    // Output the dense vector
    cout<<nline<<"The input dense vector is: "<<nline;
    for(auto it:b)
        cout<<it<<" "; cout<<nline<<nline;
    
    // Output the sequential answer
    cout<<"The sequential answer is: "<<nline;
    for(int i=0; i<n; i++)cout<<ansSequential[i]<<" "; cout<<nline<<nline;
    
    // Calculating total time taken by the program.
    double time_taken = 
      chrono::duration_cast<chrono::nanoseconds>(end - start).count();
  
    time_taken *= 1e-9;
  
    cout << "Time taken by program is " << fixed 
         << time_taken << setprecision(9);
    cout << " sec" << nline;
}

int main()
{
    ios_base::sync_with_stdio(false);    cin.tie(NULL);    cout.tie(NULL);
    
    // run for threads
    solve();
    
    return 0;
}
