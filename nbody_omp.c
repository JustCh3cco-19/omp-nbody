#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include "student.h"

int main(int argc, char **argv)
{
	char fn[FILENAME_MAX];
	FILE *fin, *fout;
	double (*pos)[3];
	double ene, cut2 = 1000.0;
	unsigned i, j, nbodies;

	if (argc != 2)
	{
		fprintf(stderr, "Usage: %s <nbodies>\n", argv[0]);
		exit(-1);
	}

	nbodies = atoi(argv[1]);

	if (nbodies == 0 || (nbodies != 20000 && nbodies != 55000))
	{
		fprintf(stderr, "Invalid number of bodies: %s\n", argv[1]);
		fprintf(stderr, "Valid options are 20000 or 55000\n");
		exit(-1);
	}

	sprintf(fn, "positions.xyz.%u", nbodies);
	fin = fopen(fn, "r");

	if (fin == NULL)
	{
		perror(fn);
		exit(-1);
	}

	pos = calloc(nbodies, sizeof(*pos));
	if (pos == NULL)
	{
		fprintf(stderr, "Not enough memory!\n");
		exit(-2);
	}

	for (i = 0; i < nbodies; ++i)
		fscanf(fin, "%lf%lf%lf", pos[i] + 0, pos[i] + 1, pos[i] + 2);

	fclose(fin);

	double start_time = omp_get_wtime();

	double (*forces)[3] = calloc(nbodies, sizeof(*forces));
	if (forces == NULL) {
		fprintf(stderr, "Not enough memory!\n");
		exit(-2);
	}
	ene = student_nbody(nbodies, pos, cut2, forces);

	double end_time = omp_get_wtime();
	printf("Computation completed in %lf seconds.\n", end_time - start_time);

	// saving results to file
	fout = fopen("results", "w");
	if (fout == NULL)
	{
		perror("results");
		exit(-1);
	}

	double *linear_forces = (double *)forces;

	fprintf(fout, "%20.10lE\n", ene);
	for (i = 0; i < nbodies; ++i)
	{
		fprintf(fout, "%5d ", i);
		for (j = 0; j < 3; ++j)
			fprintf(fout, "%20.10lE", linear_forces[i * 3 + j]);
		fprintf(fout, "\n");
	}

	fclose(fout);
	free(forces);
	free(pos);
	return 0;
}