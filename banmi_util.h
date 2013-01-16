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
    int *dat;
} tab_t;

tab_t* alloc_tab(int n, const int *dim);
void free_tab(tab_t *t);
void tab_array_index(int *arr_ix, const tab_t *t, int flat_ix);
int tab_flat_index(int *flat_index, const tab_t *t, const int *arr_ix);
int tab_get(const tab_t *t, const int *arr_ix);
void tab_set(tab_t *t, const int *arr_ix, int val);
tab_t* new_marginal_tab(const tab_t *t, const int *margins, int n_margins);
tab_t* new_conditional_tab(const tab_t *t, const int *free_margins, const int *conditioned_on, int n_free);

void order(int *values, const int *by, int n);
void order_d(double *values, const int *by, int n);
void order_blocks(int *values, const int *by, int n_blocks, int block_size);
void order_rows_int(gsl_matrix_int *mat, const int *by, int n_rows);
void order_rows(gsl_matrix *mat, const int *by, int n_rows);

int sample(gsl_rng *rng, int n, const int *weight);
int sample_d(gsl_rng *rng, int n, const double *weight);
