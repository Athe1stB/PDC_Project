#include<bits/stdc++.h>
using namespace std;
using namespace std::chrono;

int MAX_THREADS;

vector<vector<double>> A;
vector<double> csr_vals;
vector<int> csr_rows;
vector<int> csr_col;
vector<int> B;
vector<double> ans;
vector<double> ansSequential;
int n;


//sequential multiplication
void sequentialMulti(){
	for(int i = 0; i < n; i++){
		int x = csr_rows[i];
		int y = csr_rows[i+1];
		for(int j = x; j < y; j++){
			ansSequential[i] += (csr_vals[j] * B[csr_col[j]]);
		}
	}
}

// this function will calculate the answer for lth row to (r-1)th row for each thread
void* parallelMulti(void *var){
	long long thID = (long long)var;
    
    int l = (thID*n)/MAX_THREADS;
    int r = ((thID+1)*n)/MAX_THREADS;
    
	for(int i = l; i < r; i++){
		int x = csr_rows[i];
		int y = csr_rows[i+1];
		for(int j = x; j < y; j++){
			ans[i] += (csr_vals[j] * B[csr_col[j]]);
		}
	}
}

int main(){
	cout<<"Enter number of threads: ";
	cin>>MAX_THREADS;
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
	n = r;
	A.resize(r, vector<double>(c, 0.0));

	//storing data in matrix
	for(int i = 0; i < nz; i++){
		file>>x>>y>>z;
		rows.push_back(x);
		col.push_back(y);
		vals.push_back(z);
		A[x-1][y-1] = z;
	}
	
	file.close();

	//converting coo format to csr with 3 arrays: values, columns and number of non zero values till a particular row
	csr_vals.resize(nz, 0);
	csr_col.resize(c, 0);
	csr_rows.resize(r+1, 0);
	for(int i = 0; i < nz; i++){
		csr_vals[i] = vals[i];
		csr_col[i] = col[i]-1;
		csr_rows[rows[i]]++;
	}
	for(int i = 0; i < r; i++)
	csr_rows[i+1] += csr_rows[i];
	
	//reading vector.txt
	ifstream file1("vector.txt");
	int val_B;
	while(file1>>val_B){
		B.push_back(val_B);
		file1.ignore(1,' ');
	}
	file1.close();

	//printing the dense vector
	cout<<"Dense Vector:\n";
	for(int i = 0; i < B.size(); i++){
		cout<<B[i]<<" ";	
	}
	cout<<"\n";
	cout<<"\n";
	
	ansSequential.resize(r, 0.0);

	//calculating A*B sequentially
	sequentialMulti();


	ans.resize(r, 0.0);

	auto start = high_resolution_clock::now();
	//Multi-threading (A*B in parallel using pthreads)
	pthread_t th_id[MAX_THREADS];
	for(int i = 0; i < MAX_THREADS; i++){
        pthread_create(&th_id[i], NULL, parallelMulti, (void *)i);
	}
	for(int i = 0; i < MAX_THREADS; i++){
        pthread_join(th_id[i], NULL);
	}

	auto stop = high_resolution_clock::now();

	//printing sequential answer
	cout<<"Sequential answer: \n";
	for(int i = 0; i < ansSequential.size(); i++)
	cout<<ansSequential[i]<<" ";
	cout<<"\n";
	cout<<"\n";

	//printing final answer from parallel code.
	cout<<"Final answer: \n";
	for(int i = 0; i < ans.size(); i++)
	cout<<ans[i]<<" ";
	cout<<"\n";
	cout<<"\n";

	//comparing final answer with sequential answer
	if(ans == ansSequential)
	cout<<"Correct";
	else
	cout<<"Incorrect";
	cout<<"\n";
	cout<<"\n";

	//for time calculation
	double duration = duration_cast<nanoseconds>(stop - start).count();
	duration *= 1e-9;
	cout<<"Time Taken using "<<MAX_THREADS<<" threads: "<<fixed<<duration<<setprecision(9)<<" sec"<<endl;
}
