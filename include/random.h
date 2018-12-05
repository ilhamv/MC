#ifndef _RANDOM_H_
#define _RANDOM_H_

// call to return a uniform random number
double Urand(void);

// call at the beginning of each history, nps = index for the history
void RN_init_particle(unsigned long long* nps); 

#endif // RANDOM_H
