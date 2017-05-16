#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

void adjacencyMatrixGen(int mat[], int n)
{
    int i, j, tmp;
    for (i = 0; i < n; i++) {
        mat[i * n + i] = 0;
        for (j = 0; j < i; j++) {
            tmp = ((i + j) % 10 + 1);
            mat[i * n + j] = tmp;  
            mat[j * n + i] = tmp;  
	    }
    }
}

void adjacencyMatrixPrint(int mat[], int n)
{
    int i, j;
    for (i = 0; i < n; i++) {
        printf("|");
        for (j = 0; j < n; j++) {
            printf("%4d|", mat[i * n + j]);  
	    }
        printf("\n");
    }
}

// shortest path calculation
void calculate_shortest_path(int p, int n, int local_mat[], int my_rank, MPI_Comm comm)
{
    int* tmp; /*k-th row array*/
    int local_index;
    int root; //process controlling row to be bcast
    int i, j, k;
    int rsz;

    tmp = malloc(n * sizeof(int));

    //double *temp;
    //temp= malloc(n*n*sizeof(double)); 
    rsz = n / p; //we expect p is a divisor of n!

    for (k = 0; k < n; k++) { 
        //determine which process owns the kth row
        root = k / rsz; 
        if (my_rank == root) { //copy k-th row to special array
            local_index = k % rsz;
            for(i = 0; i< n; i++) {
                tmp[i] = local_mat[local_index * n + i];
            }
        }
        //broadcast k-th row
        MPI_Bcast(tmp, n, MPI_INT, root, MPI_COMM_WORLD);
        for (i = 0; i < rsz; i++) {
		    for (j = 0; j < n; j++) {
                local_mat[i* n + j] = MIN(local_mat[i * n + j], local_mat[i * n + k] + tmp[j]);
            }
        }
	}

    free(tmp);
    return;
}

int main(int argc, char* argv[])
{
    int p, my_rank, n;
    double start, end;
    n = 16; //n is number of nodes
    int* mat;
    int lm_sz;
    int* local_mat;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    /*Local rowwise block of matrix memory allocation*/
    /*p expected to be a divisor of n*/
    lm_sz = n * (n / p);
	local_mat = malloc(lm_sz * sizeof(int));

    /*Global matrix memory allocation & generation in rank=0*/
    mat = NULL;
    if(my_rank == 0) {
        mat = malloc(n * n * sizeof(int));
        adjacencyMatrixGen(mat, n);
    }
    /*Delivering matrix elements to processes*/
    MPI_Scatter((void*)mat, lm_sz, MPI_INT, 
                (void*)local_mat, lm_sz, MPI_INT,
                0, MPI_COMM_WORLD);
		
    start = MPI_Wtime();//To find the time the path calulation take

	calculate_shortest_path(p, n, local_mat, my_rank, MPI_COMM_WORLD);

	end = MPI_Wtime();

    /*Gathering back local matrix elements to main matrix*/
    MPI_Gather((void*)local_mat, lm_sz, MPI_INT,
                (void*)mat, lm_sz, MPI_INT,
                0, MPI_COMM_WORLD);

	if (my_rank == 0) {
        printf("Path Calculation Runtime: %f \nShortest-Paths matrix\n", end-start);
        printf("Matrix values:\n");
        adjacencyMatrixPrint(mat, n);
	}

    if(mat) {
        free(mat);
    } 
    free(local_mat);
     
    MPI_Finalize();
    return 0;
} 

/*
 * vim: set expandtab sw=4 ts=4 sts=4:
 */
