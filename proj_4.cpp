#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>
#include <fstream>
#include <vector>
#include <sstream>
#include <math.h>
#include <chrono>

#define MASTER  0
#define TAG     0
#define MSGSIZE 100
#define MAX 25

using namespace std;

class sample_stats
{
	public:
		double t_stat;
		int non_empty_test;
		int non_empty_ctrl;
		double converted[60];
};

sample_stats make_ttest(string sample_string[]){
	double converted[60];
	int non_empty_test = 0;
	int non_empty_ctrl = 0;
	double sum_test = 0;
	double sum_ctrl = 0;
	double mean_test, sd_sq_test, mean_ctrl, sd_sq_ctrl, t_stat;
	double temp = 0;
	sample_stats key_stats;
	
	for(int i = 1; i <= 60; i++){
		converted[i - 1] = atof(sample_string[i].c_str());
		// add them to the returned class 
		key_stats.converted[i - 1] = converted[i - 1];
	}
	
	// get the mean for all non 0 items
	for(int i = 0; i < 8; i++){
		if(converted[i] != 0){
			sum_test = sum_test + converted[i];
			non_empty_test++;
		}
	}
	mean_test = sum_test / non_empty_test;
	
	//calculate sd now that we have the mean
	for(int i = 0; i < 8; i++){
		if(converted[i] != 0){
			temp = temp + pow(converted[i] - mean_test, 2);
		}
	}
	
	sd_sq_test = temp / (non_empty_test - 1);
	
	temp = 0;
	for(int i = 8; i < 60; i++){
		if(converted[i] != 0){
			sum_ctrl = sum_ctrl + converted[i];
			non_empty_ctrl++;
		}
	}
	mean_ctrl = sum_ctrl / non_empty_ctrl;
	
	//calculate sd now that we have the mean
	for(int i = 8; i < 60; i++){
		if(converted[i] != 0){
			temp = temp + pow(converted[i] - mean_ctrl, 2);
		}
	}
	sd_sq_ctrl = temp / (non_empty_ctrl - 1);
	
	t_stat = (mean_test - mean_ctrl)/(sqrt((sd_sq_test/non_empty_test) + (sd_sq_ctrl/non_empty_ctrl)));
	
	key_stats.t_stat = t_stat;
	key_stats.non_empty_ctrl = non_empty_ctrl;
	key_stats.non_empty_test = non_empty_test;
		
	return key_stats;
}

double bootstrapper(sample_stats key_stats){
	double bootstrapped[60];
	string bootstrapped_str[61];
	int temp;
	int test_count = 0;
	int ctrl_count = 0;
	bool index_not_empty = false;
	sample_stats bootstrapped_stats;
	
	// load bootstrapped with the bootstrapped index
	for(int i = 0; i < 60; i++){
		index_not_empty = false;
		// get a random index - check to see if it is empty. If so, get a new number
		while(!index_not_empty){
			temp = rand() % 60;
			if(key_stats.converted[temp] != 0){
				index_not_empty = true;
			} 
		}
		// if we're in the test zone (first 8 spots) AND we have not reached
		// the number of non empty test entries, place the value.
		// If not zero out and let make_ttest() handle later.
		if(i < 8 & test_count <= key_stats.non_empty_test){
			bootstrapped[i] = key_stats.converted[temp];
			test_count++;
		} else if(i < 8 & i > key_stats.non_empty_test) {
			bootstrapped[i] = 0;
		} else if(i < 60 & ctrl_count <= key_stats.non_empty_ctrl){
			bootstrapped[i] = key_stats.converted[temp];
			ctrl_count++;
		} else {
			bootstrapped[i] = 0;
		}
		
		bootstrapped_str[i + 1] = to_string(bootstrapped[i]);
	}
	
	bootstrapped_stats = make_ttest(bootstrapped_str);
	return bootstrapped_stats.t_stat;
}

double calc_D_score(double bootstraps[], double t_stat, int bootstrap_size){
	double sum = 0;
	double bootstrap_mean, bootstrap_sd, temp, d_score;
	
	for(int i = 0; i < bootstrap_size; i++){
		sum += bootstraps[i];
	}
	
	bootstrap_mean = sum / bootstrap_size;
	
	//calculate sd now that we have the mean
	for(int i = 0; i < bootstrap_size; i++){
		temp += pow(bootstraps[i] - bootstrap_mean, 2);
	}
	bootstrap_sd = temp / bootstrap_size;
	
	d_score = (t_stat - bootstrap_mean) / bootstrap_sd;
	
	return d_score;
}

int main(int argc, char* argv[])
{
    int my_rank, source, num_nodes, groups, remaining_rows, end_interval, a, b;
	double start_time, end_time, start_interval;
	string sample[61];
    char my_host[MAX];
	sample_stats key_stats; 
	double bootstrapped_distro[500];
	ofstream logFile;
	char filename[64];
	double d_scores[4549];
	
	
	srand(time(NULL)); 
	
	// Read in file
	fstream file_in;
	file_in.open("NCI-60.csv", ios::in);
	vector<string> row; 
    string line, word;
	
	while(getline(file_in, line)){ 
		// used for breaking words 
		stringstream s(line); 

		// read every column data of a row and 
		// store it in a string variable, 'word' 
		while (getline(s, word, ',')) { 
			// add all the column data 
			// of a row to a vector 
			row.push_back(word); 	
		}
	}
	file_in.close();
	
	auto start = chrono::steady_clock::now();
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_nodes);
	
	start_interval = 0.0;
    end_interval = 4549;
	groups = end_interval / num_nodes;
	remaining_rows = end_interval % num_nodes;
	
    a = start_interval + my_rank*groups;
    b = a + groups - 1;
	
	// tack on remaining rows to the last node
	if(my_rank == num_nodes - 1) b = b + remaining_rows;
	
    if (my_rank != MASTER) {
        gethostname (my_host, MAX);
		double local_dscores[b-a];
		
		for(int i = a; i <= b; i++){
			// skip header
			if(i == 0) {
				continue;
			}
			for(int j = 0; j < 61; j++){
				sample[j] = row[i*61 + j];
			}

			key_stats = make_ttest(sample);
			for(int j = 0; j < 500; j++){
				bootstrapped_distro[j] = bootstrapper(key_stats);
			}
			local_dscores[i - a] = calc_D_score(bootstrapped_distro, key_stats.t_stat, 500);
		}
		
		cout << my_host << ": " << local_dscores[0] << endl;
//		MPI_Gather();
    }
    else {		
        gethostname (my_host, MAX);
		
		for(int i = a; i <= b; i++){
			// skip header
			if(i == 0) {
				continue;
			}
			for(int j = 0; j < 61; j++){
				sample[j] = row[i*61 + j];
			}

			key_stats = make_ttest(sample);
			for(int j = 0; j < 500; j++){
				bootstrapped_distro[j] = bootstrapper(key_stats);
			}
			// perhaps gather these
			d_scores[i] = calc_D_score(bootstrapped_distro, key_stats.t_stat, 500);
		}
    }

	logFile.open("d_scores.txt", ios_base::app);
	for(int i = 1; i <= 4549; i++){
		logFile << d_scores[i] << endl;
	}
	logFile.close();
    MPI_Finalize();
	
	auto end = chrono::steady_clock::now();
	
	logFile.open("logging.txt", ios_base::app);
	logFile << num_nodes << "," 
			<< chrono::duration_cast<chrono::milliseconds>(end - start).count() << "," 
			<< "ms" << endl;
	logFile.close();
	


    return 0;
}

