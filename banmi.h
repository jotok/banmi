#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>

#include "banmi_util.h"

typedef struct {
    // Data: multivariate discrete variable and univariate continuous variable,
    // including storage to contain both the original observations (with missing
    // fields) and the original data with imputed fields.
    tab_t *disc;      
    tab_t *disc_imp;
    double *cont;     
    double *cont_imp;

    // Latent variables
    tab_t *x;
    double *mu, sigma, lambda;

    // hyperparameters
    tab_t *crosstab;
    double dp_weight, sigma_a, sigma_b, lambda_a, lambda_b, mu_a, mu_b;

    // data properties
    int n_rows, n_complete, n_missing, n_disc, *mis_pat;
    int mask_missing_cont_data;
    gsl_vector *bds_disc;
} banmi_model_t;

banmi_model_t* new_banmi_model(int max_rows, gsl_vector *bds_disc,
                               double dp_weight, double sigma_a, double sigma_b,
                               double lambda_a, double lambda_b);
void banmi_data_augmentation(gsl_rng *rng, banmi_model_t *model, int n_iter);
