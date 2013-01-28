#ifndef BanmiUtil_h
#define BanmiUtil_h

#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>

typedef struct {
    int n;
    int size;

    int *dim;
    int *step;
    double *dat;
} tab_t;

tab_t* alloc_tab(int, const int*);
void free_tab(tab_t*);
void tab_array_index(int*, const tab_t*, int);
int tab_flat_index(int*, const tab_t*, const int*);
double tab_get(const tab_t*, const int*);
void tab_set(tab_t*, const int*, double);
tab_t* new_marginal_tab(const tab_t*, const int*, int);
tab_t* new_conditional_tab(const tab_t*, const int*, const int*, int);

void order(int*, const int*, int);
void order_d(double*, const int*, int);
void order_blocks(int*, const int*, int, int);
void order_blocks_d(double*, const int*, int, int);

int sample(gsl_rng*, int, const int*);
int sample_d(gsl_rng*, int, const double*);

#endif
