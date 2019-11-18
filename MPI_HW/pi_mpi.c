// MPI-based computation of PI
// demonstrates how to use rank to determine workload
// gw

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define MASTER  0
#define TAG     0

// integrates area under positive quadrant of unit circle
float integrate (double a, double b, long int steps)
{
    double x, y, width, psum = 0.0;
    long int i;

    width = (b - a) / steps;
    for (i=0; i < steps; i++) {
        x = a + i * width;
        y = sqrt(1 - x * x);
        psum = psum + y * width;
    }
    return (psum);
}

int main(int argc, char* argv[])
{
    long int num_steps; 
    int my_rank, num_nodes, source;
    double start, end, a, b, psum, sum;

    if (argc != 2) {
        fprintf (stderr, "usage: pi_MPI num_steps\n");
        exit(-1);
    }
    // number of rectangles to use per process
    num_steps = atol (argv[1]);
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_nodes);

    start = 0.0;
    end = 1.0;
    a = start + my_rank*(end-start) / num_nodes;
    b = a + (end-start) / num_nodes;
    psum = integrate(a, b, num_steps);
            
    if (my_rank != MASTER) {
        MPI_Send(&psum, 1, MPI_DOUBLE, MASTER, TAG, MPI_COMM_WORLD);
    }
    else {
        sum = psum;
        for (source = 1; source < num_nodes; source++) {
            MPI_Recv(&psum, 1, MPI_DOUBLE, source, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            sum = sum + psum;
        } 
        // PI = sum of 4 quadrants of unit circle
        printf("Pi: %2.6f\n", 4 * sum);
    }

    MPI_Finalize();

    return 0;
}
