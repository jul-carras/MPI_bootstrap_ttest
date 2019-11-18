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

int main(int argc, char* argv[])
{
    int my_rank, source, num_nodes, groups, remaining_rows, end_interval;
	double start_time, end_time, start_interval, a, b;
    char my_host[MAX];
    char message_host[MSGSIZE];
	
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
	
    a = start_interval + my_rank*groups + 61;
    b = a + (my_rank + 1)*groups*61;
	
    if (my_rank != MASTER) {
        gethostname (my_host, MAX);
//		start_send = MPI_Wtime();
//		MPI_Send(an_array, array_size, MPI_INT, MASTER, 0, MPI_COMM_WORLD);
//		end_send = MPI_Wtime();
							
		cout << my_host << " - a: " << a << ", b: " << b << endl;
		
//        MPI_Send(message_host, strlen(message_host) + 1, MPI_CHAR, MASTER, 1, MPI_COMM_WORLD);
//		MPI_Send(array_metrics, 2, MPI_DOUBLE , MASTER, 2, MPI_COMM_WORLD);*/
    }
    else {		
        gethostname (my_host, MAX);
        cout << my_host << " - a: " << a << ", b: " << b << endl;
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