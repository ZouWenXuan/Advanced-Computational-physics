#include "mpi.h"
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>

// Master-slaves Mode
void segment(int ndim, int nprc, int* vcbg, int* vced, int* vcln);
double fRand(double fMin, double fMax);

int main(int argc, char *argv[])
{
	// counting number;
	int i, j, k; 
	double sum;

	// Computing parameters;
	int n, m, ni;
	int fMin = 0;
	int fMax = 1;
	double norm, normi=0;

	// MPI parameter;
	int myid, namelen, nprc, masternode = 0;
	double startwtime = 0.0, endtime;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	
	// tag for MPI Send and Recv;
	int tag_vec1 = 10;
	int tag_vcln = 30;
	int tag_mat = 50;
	int tag_vec2 = 70;

	// Intermediate variables;
	double* vec1 = NULL;
	double** mat = NULL;
	double* vec2 = NULL;
	double** mati = NULL;
	double* vecj = NULL;

	int* vcbg = NULL;
	int* vced = NULL;
	int* vcln = NULL;

	// Main process;
	MPI_Init(&argc, &argv);
	MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &nprc);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Get_processor_name(processor_name, &namelen);
	

	if (myid == masternode)
	{
		printf("Master-Slave Mode to compute a random matrix by a random vector!\n");
		printf("\nProcess (kernal) involved as follows:\n");
		
	}

	MPI_Barrier(MPI_COMM_WORLD);
	fprintf(stdout, "	Process %d of %d on %s involved!\n", myid, nprc, processor_name);
	fflush(stdout); 
	
	MPI_Barrier(MPI_COMM_WORLD);
	if (myid == masternode) 
	{
		//input the rows and cols
		printf("\nPlease give the rows of the matrix, N = ");
		fflush(stdout);
		scanf_s("%d", &n);
		printf("Please give the cols of the matrix, M = ");
		fflush(stdout);
		scanf_s("%d", &m);
	}

	// Broadcast rows and cols to all process;
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	

	if (myid == masternode) // masternode perform initialization and integration;
	{
		// Goals: Compute Vec2 = Mat*Vec1;
		// Initialized the Vec1;
		vec1 = (double*)malloc(m * sizeof(double));
		for (i = 0; i < m; i++)
		{
			vec1[i] = fRand(fMin, fMax);
		}
		// Initialized the Mat;
		mat = (double**)malloc(n * sizeof(double*));
		mat[0] = (double*)malloc(n * m * sizeof(double));
		for (i = 1; i < n; i++)
		{
			mat[i] = mat[i - 1] + m;
		}
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < m; j++)
			{
				mat[i][j] = fRand(fMin, fMax);
			}
		}
		startwtime = MPI_Wtime();
		
		// Send Vec1 to slave process;
		printf("\nSend data to slave processes.\n ");
		for (i = 1; i < nprc; i++)
		{
			printf("	Send Vec1 to slave processes.\n");
			MPI_Send(&vec1[0], m, MPI_DOUBLE, i, tag_vec1+i, MPI_COMM_WORLD);
		}

		// Matrix segmentation: obtain the task index for different process;
		vcbg = (int*)malloc(nprc * sizeof(int));
		vced = (int*)malloc(nprc * sizeof(int));
		vcln = (int*)malloc(nprc * sizeof(int));
		segment(n, nprc-1, vcbg, vced, vcln);

		// Send numbers of rows, and corresponding data;
		
		for (i = 1; i < nprc; i++)
		{
			printf("	Send segmented data to Process %d, who computes %d rows.\n", i, vcln[i - 1]);
			MPI_Send(&vcln[i-1], 1, MPI_INT, i, tag_vcln + i, MPI_COMM_WORLD);
			MPI_Send(&mat[vcbg[i-1]][0], vcln[i-1]*m, MPI_DOUBLE, i, tag_mat + i, MPI_COMM_WORLD);
		}

		// Receive the results from slave process;
		printf("\nReceive data from slave processes.\n ");
		vec2 = (double*)malloc(n * sizeof(double));

		for (i = 1; i < nprc; i++)
		{
			printf("	Receive data from Process %d, \n", i);
			MPI_Recv(&vec2[vcbg[i - 1]], vcln[i - 1], MPI_DOUBLE, i, tag_vec2 + i, MPI_COMM_WORLD, &status);
		}

	}

	else //Slave mode perform computation.
	{	
		// Receive Vec1;
		vec1 = (double*) malloc (m * sizeof(double));
		MPI_Recv(&vec1[0], m, MPI_DOUBLE, 0, tag_vec1+myid, MPI_COMM_WORLD, &status);
		MPI_Recv(&ni, 1, MPI_INT, 0, tag_vcln + myid, MPI_COMM_WORLD, &status);

		// Receive rows data from Mat;
		mati = (double**)malloc(ni * sizeof(double*));
		mati[0] = (double*)malloc(ni * m * sizeof(double));
		for (i = 1; i < ni; i++)
		{
			mati[i] = mati[i - 1] + m;
		}	
		MPI_Recv(&mati[0][0], ni*m, MPI_DOUBLE, 0, tag_mat + myid, MPI_COMM_WORLD, &status);

		// Compute the multiplitation and send the results.
		vecj = (double*)malloc(ni * sizeof(double));
		for (i = 0; i < ni; i++)
		{
			sum = 0;
			for (j = 0; j < m; j++)
			{
				sum += mati[i][j] * vec1[j];
			}
			vecj[i] = sum;
		}
		MPI_Send(&vecj[0], ni, MPI_DOUBLE, 0, tag_vec2 + myid, MPI_COMM_WORLD);
		
		// Normalization
		for (i = 0; i < ni; i++)
		{
			normi += vecj[i] * vecj[i];
		}
		
	}
	MPI_Reduce(&normi, &norm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	

	if (myid == masternode)
	{
		for (i = 0; i < n; i++)
		{
			vec2[i] = vec2[i] / sqrt(norm);
		}	 
		endtime = MPI_Wtime();
		printf("\nTotal time cost = %f seconds\n", endtime - startwtime);
		fflush(stdout);
		free(vec1);
		free(mat[0]);
		free(vcbg);
		free(vced);
		free(vcln);
		free(vec2);
	}
	else
	{
		free(vec1);
		free(mati[0]);
		free(mati);
		free(vecj);
	}
	MPI_Finalize();

}


void segment(int ndim, int nprc, int* vcbg, int* vced, int* vcln)
{
	int i, j, k, m;

	k = ndim / nprc;
	j = ndim % nprc;
	m = 0;
	for (i = 0; i < nprc; i++)
	{
		vcbg[i] = m;
		m = m + k;
		if (j > i)
		{
			m = m + 1;
		}
		vcln[i] = m - vcbg[i];
		vced[i] = vcbg[i] + vcln[i] - 1;
	}
}



double fRand(double fMin, double fMax)
{
	double f = (double)rand() / ((double)RAND_MAX + 1.0);
	return fMin + f * (fMax - fMin);
}
