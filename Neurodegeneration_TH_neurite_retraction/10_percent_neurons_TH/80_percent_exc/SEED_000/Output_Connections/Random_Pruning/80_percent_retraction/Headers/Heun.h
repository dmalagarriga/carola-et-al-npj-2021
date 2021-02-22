#ifndef _HEUN_H_
#define _HEUN_H_
#include <iostream>
#include "../Headers/Parameters.h"
#include "../Headers/rhs.h"
#include <math.h>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

extern gsl_rng * r;

void Heun (Parameters *);

#endif
