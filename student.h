#ifndef STUDENT_H
#define STUDENT_H

/* Fill the initially zero forces[nbodies][3] and return potential energy.
 * pos contains coordinates; cut2 is the squared interaction cutoff. */
double student_nbody(unsigned nbodies, double (*pos)[3], double cut2,
                     double (*forces)[3]);

#endif
