#ifndef Banmi_h
#define Banmi_h

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_vector.h>

#include "banmi_util.h"

#define BoundaryCorrection

typedef struct {
    // Data: original and with imputed values
    gsl_matrix_int *disc, *disc_imp;
    gsl_matrix *cont, *cont_imp;

    // Latent variables
    gsl_matrix_int *x;
    gsl_matrix *mu;
    gsl_vector *sigma, *lambda;

    // hyperparameters
    double dp_weight;
    double init_crosstab;
    tab_t *crosstab;
    gsl_vector *mu_a, *mu_b, *sigma_a, *sigma_b;
    double lambda_a, lambda_b;

    // data properties
    int n_rows, n_complete, n_missing, n_disc, n_cont, *mis_pat;
    int mask_missing_cont_data, mask_missing_disc_data;
    gsl_vector_int *bds_disc;

    // misc properties
    int is_initialized;

} banmi_model_t;

banmi_model_t* 
new_banmi_model(int, gsl_vector_int*, int, double, double, double, double);

void
banmi_add_row(banmi_model_t*, gsl_vector_int*, gsl_vector*);

void 
banmi_data_augmentation(gsl_rng*, banmi_model_t*, int);

int 
banmi_to_ordered_value(double, int);

double 
banmi_from_ordered_value(int, int);

void
banmi_print_shape_variables(banmi_model_t*);

int
banmi_count_unique_modes(banmi_model_t*);

#endif
