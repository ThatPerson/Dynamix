//
// Created by ben on 19/01/2022.
//

#ifndef DYNAMIX_PREDICT_RATE_H
#define DYNAMIX_PREDICT_RATE_H

#include "datatypes.h"
#include "read_data.h"

#include "global_gaf.h"
#include "runners.h"
#include "chisq.h"

int print_prate(struct Model *m);
int read_relaxation_set(struct Model *m, char *fn);
void clean_relaxation(struct Model *m);
int read_fit_data(struct Model *m);

#endif //DYNAMIX_PREDICT_RATE_H
