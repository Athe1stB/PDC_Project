#include<bits/stdc++.h>
using namespace std;
using namespace std::chrono;

#define MAX_THREADS 2
#define MAX_BLOCKS 8

float** h_A;
float* h_csr_vals;
int* h_csr_rows;
int* h_csr_col;
float* h_B;
float* h_C1;
float* h_C2;
int n, m;

//sequential multiplication
void sequentialMulti(){
	for(int i=0; i<n; i++)
    {
        h_C1[i] = 0;
		for(int k = h_csr_rows[i]; k < h_csr_rows[i+1]; k++)
        {
            int index = h_csr_col[k];
            h_C1[i]+= (h_csr_vals[k]*h_B[index]);
        }
    }
}

__global__ void parallelMulti(int d_n, int* d_csr_rows, float* d_C2, float* d_csr_vals, float* d_B, int* d_csr_col){

	long long int thread_id = blockIdx.x * MAX_THREADS + threadIdx.x;
	int m = d_n;
	int lx = (thread_id*(d_n+m+1))/MAX_THREADS;
	int rx = (thread_id+1)*(d_n+m+1); rx/=MAX_THREADS;
	
	int l = 0, r = lx;
    int r_el_f_A = 0;
    
    while(l<=r){
        int mid = l+(r-l)/2;
        int x = mid -1;
        int y = lx - mid - 1;
        
        if(x>d_n)
            r = mid-1;
        else if(x<0 || y<0 || d_csr_rows[x]<=y)
            r_el_f_A = mid,
            l = mid + 1;
        else
            r = mid -1;
    }
	
	int ind_el_A = lx - r_el_f_A;
	int r_el_f_B = d_n+1;
	if(thread_id!=MAX_THREADS-1)
 	{
		int l = 0, r = rx;
		r_el_f_B = 0;
		
		while(l<=r){
			int mid = l+(r-l)/2;
			int x = mid -1;
			int y = rx - mid - 1;
			
			if(x>d_n)
				r = mid-1;
			else if(x<0 || y<0 || d_csr_rows[x]<=y)
				r_el_f_B = mid,
				l = mid + 1;
			else
				r = mid -1;
		}
	 }

	int ind_el_B = rx - r_el_f_B;
	
	int i = r_el_f_A; int j = ind_el_A;
	int N = r_el_f_B; int M = ind_el_B;
	
	while(i<N && j<M){
	    if(d_csr_rows[i]<=j) //no entries in the row (moving down).
	        i++;
	    else{
	       //moving right and adding row_element*corresponding_element_of_vector to answer vector.
	       d_C2[i-1]+= d_csr_vals[j] * d_B[d_csr_col[j]];
	       j++;
	   }
	}
	while(i<N){
	    i++;
	}
	while(j<M){
	    d_C2[i-1]+= d_csr_vals[j] * d_B[d_csr_col[j]];
	    j++;
	}
}

int main(){
	auto start = high_resolution_clock::now();
	
	vector<double> vals;
	vector<int> rows;
	vector<int> col;

    
	//reading inputfile.mtx
	ifstream file("inputfile.mtx");
	
	//removing header comments
	while(file.peek() == '%')
	file.ignore(1000, '\n');

	int x, y;
	double z;
	int r, c, nz;
	
	//reading number of rows, cols, and non-zero values.
	file>>r>>c>>nz;
	n=r;
	m=c;
	//A.resize(r, vector<double>(c, 0.0));

	//storing data in matrix
	for(int i = 0; i < nz; i++){
		file>>x>>y>>z;
		rows.push_back(x);
		col.push_back(y);
		vals.push_back(z);
		//A[x-1][y-1] = z;
	}
	file.close();
	
	//converting coo format to csr with 3 arrays: values, columns and number of non zero values till a particular row

	h_csr_vals = (float*)malloc(nz*sizeof(float));
	h_csr_col = (int*)malloc(c*sizeof(int));
	h_csr_rows = (int*)malloc((r+1)*sizeof(int));
	h_B = (float*)malloc(r*sizeof(float));
	for(int i = 0; i < nz; i++){
		h_csr_vals[i] = vals[i];
		h_csr_col[i] = col[i]-1;
		h_csr_rows[rows[i]]++;
	}
	for(int i = 0; i < r; i++)
	h_csr_rows[i+1] += h_csr_rows[i];
	
	// reading vector.txt
	ifstream file1("vector.txt");
	int val_B; int p = 0;
	while(file1>>val_B){
		h_B[p] = val_B;
		file1.ignore(1,' ');
		p++;
	}
	file1.close();


	h_C1 = (float*)malloc(r*sizeof(float));
	h_C2 = (float*)malloc(r*sizeof(float));

	
	//A*B sequentially
	sequentialMulti();
	cout<<endl;


	//Multi-threading (A*B in parallel using CUDA)
	//Memory allocation in GPU
	int* d_csr_rows;
	cudaMalloc(&d_csr_rows, (n+1)*sizeof(int));
	float* d_C2;
	cudaMalloc(&d_C2, n*sizeof(float));
	float* d_csr_vals;
	cudaMalloc(&d_csr_vals, n*sizeof(float));
	float* d_B;
	cudaMalloc(&d_B, n*sizeof(float));
	int* d_csr_col;
	cudaMalloc(&d_csr_col, n*sizeof(int));

	cudaMemcpy(d_csr_rows, h_csr_rows, (n+1)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_csr_col, h_csr_col, (n)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, h_B, (n)*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_csr_vals, h_csr_vals, (n)*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_C2, h_C2, (n)*sizeof(float), cudaMemcpyHostToDevice);
    
	//SpMV using CUDA
	parallelMulti<<<MAX_BLOCKS, MAX_THREADS>>>(n, d_csr_rows, d_C2, d_csr_vals, d_B, d_csr_col);
	
	cudaMemcpy(h_C2, d_C2, n*sizeof(float), cudaMemcpyDeviceToHost);
	
	//Free the memory allocated in GPU
	cudaFree(d_csr_rows);
	cudaFree(d_csr_vals);
	cudaFree(d_csr_col);
	cudaFree(d_B);
	cudaFree(d_C2);

	//printing sequential answer
	cout<<"Sequential Answer: \n";
	for(int i = 0; i < r; i++)
	cout<<h_C1[i]<<" ";
	cout<<endl;
	cout<<endl;

	//printing final answer using CUDA
	cout<<"Final Answer using CUDA: \n";
	for(int i = 0; i < r; i++)
	cout<<h_C2[i]<<" ";
	cout<<endl;

	
	//for execution time calculation
	cout<<endl;
	auto stop = high_resolution_clock::now();
	float duration = duration_cast<nanoseconds>(stop - start).count();
	duration*= 1e-9;
	cout << "Time taken using "<< MAX_BLOCKS<<" blocks and "<<MAX_THREADS <<" threads is " << fixed 
         << duration << setprecision(9);
    cout << " sec" << endl;
}
