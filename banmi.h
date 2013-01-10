#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "banmi_util.h"

#define MaxRows 300 // maximum number of rows of data
#define NDiscrete 5 // number of discrete variables
#define NIter 50

#define DPWeight 50        // weight paramter to the DP
#define SigmaA 17.5        // alpha parameter to sigma prior
#define SigmaB 1.0/129838  // beta parameter to sigma prior
#define LambdaA 3.0        // alpha parameter to lambda prior
#define LambdaB 2.0        // beta parameter to lambda prior

#define MissingContinuousData 32  // bit field mask to indicate missing continuous data

const int DiscDim[5] = {2, 2, 2, 2, 3};

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
    double mu_a, mu_b, sigma_a, sigma_b, lambda_a, lambda_b, dp_weight;

    // data properties
    int n_rows, n_complete, n_missing, *mis_pat;
} banmi_model_t;

