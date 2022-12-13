// #include "mpi.h"
// #include <stdio.h>
// #include <math.h>
// #include <limits.h>
// #include <stdlib.h>

// // SPMD Mode
// void segment(int ndim, int nprc, int* vcbg, int* vced, int* vcln);
// double fRand(double fMin, double fMax);


// int main(int argc, char *argv[])
// {
// 	// counting number;
// 	int i, j;
// 	double sum;

// 	// Computing parameters;
// 	int n, m;
// 	int fMin = 0;
// 	int fMax = 1;
// 	double norm, normi = 0;

// 	// MPI parameter;
// 	int myid, namelen, nprc, masternode = 0;
// 	double startwtime = 0.0, endtime;
// 	char processor_name[MPI_MAX_PROCESSOR_NAME];

// 	// Intermediate variables;
// 	double* vec1 = NULL;
// 	double** mat = NULL;
// 	double* vec2 = NULL;
// 	double* vecj = NULL;
// 	double** mati = NULL;

// 	int* vcbg = NULL;
// 	int* vced = NULL;
// 	int* vcln = NULL;
// 	int* ndata = NULL;
// 	int* ndispl = NULL;

// 	// Main process:
// 	MPI_Init(&argc, &argv);
// 	MPI_Comm_size(MPI_COMM_WORLD, &nprc);
// 	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
// 	MPI_Get_processor_name(processor_name, &namelen);
	

// 	if (myid == masternode)
// 	{
// 		printf("SPMD Mode to compute a random matrix by a random vector!\n");
// 		printf("\nProcess (kernal) involved as follows:\n");
// 	}

// 	MPI_Barrier(MPI_COMM_WORLD);
// 	fprintf(stdout, "Process %d of %d on %s\n", myid, nprc, processor_name);
// 	fflush(stdout); 


// 	MPI_Barrier(MPI_COMM_WORLD);
// 	if (myid == masternode) 
// 	{
// 		//input the rows and cols
// 		printf("\nPlease give the rows of the matrix, N = ");
// 		fflush(stdout);
// 		scanf_s("%d", &n);
// 		printf("Please give the cols of the matrix, M = ");
// 		fflush(stdout);
// 		scanf_s("%d", &m);
// 	}

// 	// Broadcast rows and cols to all process;
// 	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
// 	MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
// 	MPI_Barrier(MPI_COMM_WORLD);
	

// 	// Goals: Compute Vec2 = Mat*Vec1;
// 	// Apply Vec1 in all process;
// 	vec1 = (double*)malloc(m * sizeof(double));
// 	if (myid == masternode) 
// 	{
		
// 		// Initialzed the Vec1;
// 		for (i = 0; i < m; i++)
// 		{
// 			vec1[i] = fRand(fMin, fMax);
// 		}
// 		// Apply Mat in master process;
// 		mat = (double**)malloc(n * sizeof(double*));
// 		mat[0] = (double*)malloc(n * m * sizeof(double));
// 		for (i = 1; i < n; i++)
// 		{
// 			mat[i] = mat[i - 1] + m;
// 		}
// 		// Initialzed the Mat;
// 		for (i = 0; i < n; i++)
// 		{
// 			for (j = 0; j < m; j++)
// 			{
// 				mat[i][j] = fRand(fMin, fMax);
// 			}
// 		}
// 		startwtime = MPI_Wtime();
// 		printf("\nSegment data to all the process involved.\n");
// 		printf("	Broadcast Vec1 to all processes.\n");
// 	}

// 	// Broadcast Vec1;
// 	MPI_Bcast(&vec1[0], m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
// 	MPI_Barrier(MPI_COMM_WORLD);

// 	// Matrix segmentation
// 	vcbg = (int* )malloc(nprc * sizeof(int));
// 	vced = (int* )malloc(nprc * sizeof(int));
// 	vcln = (int* )malloc(nprc * sizeof(int));
// 	segment(n, nprc, vcbg,  vced, vcln);

// 	// MPI_Scatterv parameter
// 	ndata = (int* )malloc(nprc* sizeof(int));
// 	for (i = 0; i < nprc; i++)
// 	{
// 		ndata[i] = vcln[i] * n;

// 	}
// 	ndispl = (int*)malloc(nprc * sizeof(int));
// 	ndispl[0] = 0;
// 	for (i = 1; i < nprc; i++)
// 	{
// 		ndispl[i] = ndispl[i - 1] + vcln[i - 1] * n;
// 	}

// 	// Scatter different rows of the matrix
// 	mati = (double**)malloc(vcln[myid] * sizeof(double*));
// 	mati[0] = (double*)malloc(vcln[myid] * m * sizeof(double));
// 	for (i = 1; i < vcln[myid]; i++)
// 	{
// 		mati[i] = mati[i - 1] + m;
// 	}
	
// 	if (myid == masternode)
// 	{
// 		MPI_Scatterv(&mat[0][0], ndata, ndispl, MPI_DOUBLE, &mati[0][0], ndata[myid], MPI_DOUBLE, 0, MPI_COMM_WORLD);
// 		printf("	Scatterv Mat to all processes.\n");
// 	}
// 	else 
// 	{
// 		MPI_Scatterv(NULL, ndata, ndispl, MPI_DOUBLE, &mati[0][0], ndata[myid], MPI_DOUBLE, 0, MPI_COMM_WORLD);
// 	}
// 	MPI_Barrier(MPI_COMM_WORLD);


// 	// Computation
// 	vecj = (double*)malloc(vcln[myid] * sizeof(double));
// 	for (i = 0; i < vcln[myid]; i++)
// 	{
// 		sum = 0;
// 		for (j = 0; j < m; j++)
// 		{
// 			sum += mati[i][j] * vec1[j];
// 		}
// 		vecj[i] = sum;
// 	}

// 	// Normalization
// 	for (i = 0; i < vcln[myid]; i++)
// 	{
// 		normi += vecj[i] * vecj[i];
// 	}

// 	// MPI_Gatherv parameter
// 	for (i = 0; i < nprc; i++)
// 	{
// 		ndata[i] = vcln[i];

// 	}
// 	ndispl[0] = 0;
// 	for (i = 1; i < nprc; i++)
// 	{
// 		ndispl[i] = ndispl[i - 1] + vcln[i - 1];
// 	}

// 	// Gatherv the vecj to Vec2;
// 	if (myid == masternode)
// 	{
// 		printf("Gatherv vecj to Vec2.\n");
// 	}
// 	vec2 = (double*)malloc(m * sizeof(double));
// 	MPI_Gatherv(&vecj[0], vcln[myid], MPI_DOUBLE, &vec2[0], ndata, ndispl, MPI_DOUBLE, 0, MPI_COMM_WORLD);
// 	MPI_Barrier(MPI_COMM_WORLD);
// 	MPI_Reduce(&normi, &norm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


// 	if (myid == masternode)
// 	{
// 		for (i = 0; i < m; i++)
// 		{
// 			vec2[i] = vec2[i] / sqrt(norm);
// 		}
// 		endtime = MPI_Wtime();
// 		printf("\nFinish. Total time cost = %f seconds\n", endtime - startwtime);
// 		fflush(stdout);
// 	}

// 	free(vec1);
// 	if (myid == 0)
// 	{
// 		free(mat[0]);
// 		free(mat);
// 	}
// 	free(vcbg);
// 	free(vced);
// 	free(vcln);
// 	free(ndata);
// 	free(ndispl);
// 	free(mati[0]); 
// 	free(mati);
// 	free(vecj);
// 	free(vec2);

// 	MPI_Finalize();
	
// }


// void segment(int ndim, int nprc, int* vcbg, int* vced, int* vcln)
// {
// 	int i, j, k, m;

// 	k = ndim / nprc;
// 	j = ndim % nprc;
// 	m = 0;
// 	for (i = 0; i < nprc; i++)
// 	{
// 		vcbg[i] = m;
// 		m = m + k;
// 		if (j > i)
// 		{
// 			m = m + 1;
// 		}
// 		vcln[i] = m - vcbg[i];
// 		vced[i] = vcbg[i] + vcln[i] - 1;
// 	}

// }

// double fRand(double fMin, double fMax)
// {
// 	double f = (double)rand() / ((double)RAND_MAX + 1.0);
// 	return fMin + f * (fMax - fMin);
// }

