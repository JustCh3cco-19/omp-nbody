#include "student.h"
#include <math.h>
#include <omp.h>

/* Parallelize this implementation. You may add includes and helper functions. */
double student_nbody(unsigned nbodies, double (*pos)[3], double cut2,
                     double (*forces)[3])
{
    double ene = 0.0;
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
                ene += -1.0 / d;
            }
        }
    }
    return ene;
}
