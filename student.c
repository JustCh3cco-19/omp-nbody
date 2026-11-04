#include "student.h"
#include <math.h>
#include <omp.h>
#include <stdlib.h>

/* Parallelize this implementation. You may add includes and helper functions. */
double student_nbody(unsigned nbodies, double (*pos)[3], double cut2,
                     double (*forces)[3])
{
    int nthreads = omp_get_max_threads();
    size_t force_count = (size_t)nthreads * nbodies * 3;
    double *private_forces = calloc(force_count, sizeof(*private_forces));
    double *private_energy = calloc((size_t)nthreads, sizeof(*private_energy));

    /* The allocations are normally small compared with the input.  Keep a
       correct serial fallback if memory is unavailable. */
    if (private_forces == NULL || private_energy == NULL) {
        free(private_forces);
        free(private_energy);
        for (unsigned i = 0; i < nbodies; ++i) {
            for (unsigned j = i + 1; j < nbodies; ++j) {
                double rij[3], d2 = 0.0;
                for (unsigned k = 0; k < 3; ++k) {
                    rij[k] = pos[i][k] - pos[j][k];
                    d2 += rij[k] * rij[k];
                }
                if (d2 <= cut2) {
                    double d = sqrt(d2), d3 = d * d2;
                    for (unsigned k = 0; k < 3; ++k) {
                        double f = -rij[k] / d3;
                        forces[i][k] += f;
                        forces[j][k] -= f;
                    }
                    private_energy = NULL;
                }
            }
        }
        /* Recompute the energy in the fallback without changing the force
           accumulation above. */
        double ene = 0.0;
        for (unsigned i = 0; i < nbodies; ++i)
            for (unsigned j = i + 1; j < nbodies; ++j) {
                double dx = pos[i][0] - pos[j][0];
                double dy = pos[i][1] - pos[j][1];
                double dz = pos[i][2] - pos[j][2];
                double d2 = dx * dx + dy * dy + dz * dz;
                if (d2 <= cut2) ene -= 1.0 / sqrt(d2);
            }
        return ene;
    }

#pragma omp parallel
    {
        int thread = omp_get_thread_num();
        double (*local)[3] = (double (*)[3])
            (private_forces + (size_t)thread * nbodies * 3);

#pragma omp for schedule(guided, 64)
        for (unsigned i = 0; i < nbodies; ++i) {
            for (unsigned j = i + 1; j < nbodies; ++j) {
                double dx = pos[i][0] - pos[j][0];
                double dy = pos[i][1] - pos[j][1];
                double dz = pos[i][2] - pos[j][2];
                double d2 = dx * dx + dy * dy + dz * dz;

                if (d2 > 0.0 && d2 <= cut2) {
                    double d = sqrt(d2);
                    double factor = -1.0 / (d * d2);

                    local[i][0] += dx * factor;
                    local[i][1] += dy * factor;
                    local[i][2] += dz * factor;
                    local[j][0] -= dx * factor;
                    local[j][1] -= dy * factor;
                    local[j][2] -= dz * factor;
                    private_energy[thread] -= 1.0 / d;
                }
            }
        }
    }

    double ene = 0.0;
    for (int thread = 0; thread < nthreads; ++thread) {
        ene += private_energy[thread];
        double (*local)[3] = (double (*)[3])
            (private_forces + (size_t)thread * nbodies * 3);
        for (unsigned i = 0; i < nbodies; ++i)
            for (unsigned k = 0; k < 3; ++k)
                forces[i][k] += local[i][k];
    }
    free(private_forces);
    free(private_energy);
    return ene;
}
