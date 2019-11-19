#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>
#include <fstream>
#include <vector>
#include <sstream>
#include <math.h>

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

sample_stats make_ttest(string sample_string[], int row){
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
	cout << "tstat: " << key_stats.t_stat << endl;
	cout << "first item of sample: " << key_stats.converted[0] << endl;
	cout << "non_zero_test: " << key_stats.non_empty_test << endl;
	cout << "non_zero_ctrl: " << key_stats.non_empty_ctrl << endl;
	
	return 0.0;
}
int main(int argc, char* argv[])
{
    int my_rank, source, num_nodes, groups, remaining_rows, end_interval;
	double start_time, end_time, start_interval, a, b;
	string sample[61];
    char my_host[MAX];
	sample_stats key_stats; 
	
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
	
//	cout << row[0] << endl;
//	cout << row[61*4549 + 60] << endl;
	
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_nodes);
	
	start_interval = 0.0;
    end_interval = 4549;
	groups = end_interval / num_nodes;
	remaining_rows = end_interval % num_nodes;
	
//	cout << "Nodes: " << num_nodes << endl;
//	cout << "Groups: " << groups << endl;
//	cout << "Remainder Rows: " << remaining_rows << endl;
	
    a = start_interval + my_rank*groups;
    b = a + groups - 1;
	
	// tack on remaining rows to the last node
	if(my_rank == num_nodes - 1) b = b + remaining_rows;
	
    if (my_rank != MASTER) {
        gethostname (my_host, MAX);
//		start_send = MPI_Wtime();
//		MPI_Send(an_array, array_size, MPI_INT, MASTER, 0, MPI_COMM_WORLD);
//		end_send = MPI_Wtime();
							

		
//        MPI_Send(message_host, strlen(message_host) + 1, MPI_CHAR, MASTER, 1, MPI_COMM_WORLD);
//		MPI_Send(array_metrics, 2, MPI_DOUBLE , MASTER, 2, MPI_COMM_WORLD);*/
    }
    else {		
        gethostname (my_host, MAX);
        cout << my_host << " - a: " << a << ", b: " << b << endl;
		
		for(int i = a; i <= a + 1; i++){
			// skip header
			if(i == 0) {
				continue;
			}
			for(int j = 0; j < 61; j++){
				sample[j] = row[i*61 + j];
			}
			key_stats = make_ttest(sample, i);
			bootstrapper(key_stats);
		}
       /*   for (source = 1; source < num_nodes; source++) {
			MPI_Recv(an_array, array_size, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(message_host, MSGSIZE, MPI_CHAR, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(array_metrics, 2, MPI_DOUBLE, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			// logging
			ofstream logFile;
			logFile.open("log_seq.txt", ios_base::app);
			logFile << argv[2] << "," << num_nodes << "," << array_size << "," << message_host << "," 
					<< array_metrics[0] << "," << array_metrics[1] << endl;
			logFile.close();
        }*/
    }

    MPI_Finalize();

    return 0;
}

