#ifndef Mico_h
#define Mico_h

#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics.h>

#include "banmi_util.h"

typedef struct {
    // original data and imputed data
    gsl_matrix *y, *yi;

    // latent variables
    gsl_matrix *xi;
    double sigma;
    
    // margin parameters
    gsl_vector *mu;
    gsl_vector *theta;

    // hyperparameters
    int n_comp;
    double sigma_a, sigma_b;
    gsl_vector *mu_a, *mu_b, *theta_a, *theta_b;

    // data properties
    int n_rows, n_cols, n_complete, n_missing, *mis_pat, is_initialized;

} mico_model_t;

mico_model_t*
new_mico_model(int, int, int, double, double);

void
mico_add_row(mico_model_t*, gsl_vector*);

#endif
