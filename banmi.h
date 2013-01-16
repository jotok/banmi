#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_vector.h>

#include "banmi_util.h"

typedef struct {
    // Data: multivariate discrete variable and univariate continuous variable,
    // including storage to contain both the original observations (with missing
    // fields) and the original data with imputed fields.
    tab_t *disc, *disc_imp;      
    gsl_matrix *cont, *cont_imp;

    // Latent variables
    tab_t *x;
    gsl_matrix *mu;
    gsl_vector *sigma, *lambda;

    // hyperparameters
    double dp_weight;
    tab_t *crosstab;
    gsl_vector *mu_a, *mu_b, *sigma_a, *sigma_b;
    double lambda_a, lambda_b;

    // data properties
    int n_rows, n_complete, n_missing, n_disc, n_cont, *mis_pat;
    int mask_missing_cont_data, mask_missing_disc_data;
    gsl_vector_int *bds_disc;

} banmi_model_t;

banmi_model_t* new_banmi_model(int max_rows, gsl_vector_int *bds_disc, int n_cont, 
                               double dp_weight, double lambda_a, double lambda_b);
void banmi_data_augmentation(gsl_rng *rng, banmi_model_t *model, int n_iter);

int banmi_to_ordered_value(double x, int max_value);
double banmi_from_ordered_value(int o, int max_value);
