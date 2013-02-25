#include "mico.h"

#include <assert.h>
#include <math.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>

/**
 * Allocate a new mico model
 */
mico_model_t*
new_mico_model(int max_rows, int n_cols, int n_comp, double sigma_a, double sigma_b,
               double x_nu, double mcmc_proposal_nu) 
{
    mico_model_t *model = malloc(sizeof(mico_model_t));

    model->y = gsl_matrix_alloc(max_rows, n_cols);
    model->yi = gsl_matrix_alloc(max_rows, n_cols);

    model->xi = gsl_matrix_alloc(max_rows, n_cols);
    model->z = gsl_matrix_alloc(n_comp, n_cols);
    model->zind = gsl_vector_int_alloc(max_rows);

    model->mu = gsl_vector_alloc(n_cols);
    model->theta = gsl_vector_alloc(n_cols);

    model->n_comp = n_comp;
    model->sigma_a = sigma_a;
    model->sigma_b = sigma_b;
    model->x_nu = x_nu;
    model->mcmc_proposal_nu = mcmc_proposal_nu;
    model->mu_a = gsl_vector_alloc(n_cols);
    model->mu_b = gsl_vector_alloc(n_cols);
    model->theta_a = gsl_vector_alloc(n_cols);
    model->theta_b = gsl_vector_alloc(n_cols);

    model->n_rows = 0;
    model->n_cols = n_cols;
    model->mis_pat = malloc(max_rows * sizeof(int));
    model->is_initialized = 0;

    return model;
}

void
mico_add_row(mico_model_t *model, gsl_vector *row) {
        gsl_matrix_set_row(model->y, model->n_rows, row);
        gsl_matrix_set_row(model->yi, model->n_rows, row);
        model->n_rows++;
}

static int
missingness_pattern(const mico_model_t *model, int row) {
    int j, pattern = 0, mask = 1;

    for (j = 0; j < model->n_cols; j++) {
        if (isnan(gsl_matrix_get(model->y, row, j)))
            pattern |= mask;

        mask = mask << 1;
    }

    return pattern;
}

/**
 * Sort the rows of data by missingness pattern and set the n_complete and
 * n_missing fields. Returns the number of complete rows.
 */
static int
sort_data_by_missingness_pattern(mico_model_t *model) {
    int i, n_complete = 0;

    for (i = 0; i < model->n_rows; i++) {
        model->mis_pat[i] = missingness_pattern(model, i);

        if (model->mis_pat[i] == 0)
            n_complete++;
    }

    banmi_order_blocks(model->y->data, model->mis_pat, model->n_rows, model->y->tda);
    banmi_order_blocks(model->yi->data, model->mis_pat, model->n_rows, model->yi->tda);
    banmi_order_int(model->mis_pat, model->mis_pat, model->n_rows);

    model->n_complete = n_complete;
    model->n_missing = model->n_rows - n_complete;

    return n_complete;
}

static void
init_parameters(gsl_rng *rng, mico_model_t *model) {
    // place a prior on the marginal variance so that the expected value is
    // the variance of the observed values, and the weight is the number of
    // complete values.
    int i, j, insert_index = 0;
    double this_y, mean, var, alpha, beta, theta, temp[model->n_rows];

    for (j = 0; j < model->n_cols; j++) {
        insert_index = 0;
        for (i = 0; i < model->n_rows; i++) {
            if (!isnan(this_y = gsl_matrix_get(model->y, i, j)))
                temp[insert_index++] = this_y;
        }

        mean = gsl_stats_mean(temp, 1, insert_index);
        var = gsl_stats_variance(temp, 1, insert_index);

        // the rule of thumb is to take alpha = n/2. since we will be reusing
        // data to compute the posterior, take n/20 instead.
        alpha = (insert_index+0.0) / 20.0;
        beta = alpha * var;
        gsl_vector_set(model->theta_a, j, alpha);
        gsl_vector_set(model->theta_b, j, beta);

        theta = 1.0 / gsl_ran_gamma(rng, alpha, 1.0 / beta);
        gsl_vector_set(model->theta, j, theta);

        alpha = mean;
        beta = (insert_index+0.0) / 10.0;
        gsl_vector_set(model->mu_a, j, alpha);
        gsl_vector_set(model->mu_b, j, beta);

        gsl_vector_set(model->mu, j,
                       alpha + gsl_ran_gaussian(rng, theta / sqrt(beta)));
    }
}

/**
 * Set initial values in the imputed data set. 
 */
static void
init_missing_values(gsl_rng *rng, mico_model_t *model) {
    // draw from marginal models assuming independence
    int i, j;
    double this_y;

    for (i = model->n_complete; i < model->n_rows; i++) {
        for (j = 0; j < model->n_cols; j++) {
            if (isnan(gsl_matrix_get(model->y, i, j))) {
                this_y = gsl_vector_get(model->mu, j) + 
                         gsl_ran_gaussian(rng, gsl_vector_get(model->theta, j));
                gsl_matrix_set(model->yi, i, j, this_y);
            }
        }
    }
}

static double
x_cdf(double x, const gsl_matrix *z, double sigma, int j) {
    int i;
    double result = 0.0;
    for (i = 0; i < z->size1; i++) {
        result += gsl_cdf_gaussian_P(x - gsl_matrix_get(z, i, j), sigma);
        
    }

    return result / z->size1;
}

static double
x_pdf(double x, const gsl_matrix *z, double sigma, int j) {
    int i;
    double result = 0.0;
    for (i = 0; i < z->size1; i++) {
        result += gsl_ran_gaussian_pdf(x - gsl_matrix_get(z, i, j), sigma);
        
    }

    return result / z->size1;
}

static double
x_pdf_multivariate(const gsl_vector *x, const gsl_matrix *z, double sigma) {
    int i, j;
    double result = 0.0;

    for (i = 0; i < z->size1; i++) {
        double component_result = 1.0;
        for (j = 0; j < x->size; j++) {
            component_result *= gsl_ran_gaussian_pdf(gsl_vector_get(x, j) - 
                                                     gsl_matrix_get(z, i, j), 
                                                     sigma);
        }
        result += component_result;        
    }

    return result / z->size1;
}

#define QuantileTolerance 1e-6

static double
u_quantile(double u, const gsl_matrix *z, double sigma, int j) {
    assert(0 < u && u < 1);

    double x = 0.0, try_x, step;
    double u_of_x = x_cdf(x, z, sigma, j);

    double sign;
    if (u_of_x <= u)
        sign = 1.0;
    else
        sign = -1.0;

    step = sign * 1.0;
    try_x = x + step;
    while (sign * (x_cdf(try_x, z, sigma, j) - u) < QuantileTolerance) {
        x = try_x;
        step *= 2.0;
        try_x = x + step;
    }
    
    u_of_x = x_cdf(x, z, sigma, j);
    step /= 2.0;

    while (sign * (u - u_of_x) >= QuantileTolerance) {
        try_x = x + step;
        while (sign * (x_cdf(try_x, z, sigma, j) - u) < QuantileTolerance) {
            x = try_x;
            try_x += step;
        }
        u_of_x = x_cdf(x, z, sigma, j);
        step /= 2.0;
    }

    return x;
}

static double
x_copula_pdf(const gsl_vector *x, const gsl_matrix *z, double sigma) {
    double result = x_pdf_multivariate(x, z, sigma);

    int j;
    for (j = 0; j < x->size; j++) {
        result /= x_pdf(gsl_vector_get(x, j), z, sigma, j);
    }

    return result;
}

static void
init_latent_variables(gsl_rng *rng, mico_model_t *model) {
    model->sigma = gsl_ran_gamma(rng, model->sigma_a, 1.0 / model->sigma_b);

    // initially draw mixture means from multivariate t distribution
    int i, j;
    double u, x;
    for (i = 0; i < model->n_comp; i++) {
        u = gsl_ran_chisq(rng, model->x_nu);

        for (j = 0; j < model->n_cols; j++) {
            x = gsl_ran_ugaussian(rng) * sqrt(model->x_nu / u);
            gsl_matrix_set(model->z, i, j, x);
        }

        gsl_vector_int_set(model->zind, i, (int)(model->n_comp * gsl_rng_uniform(rng)));
    }

    // map y's to x's
    for (i = 0; i < model->n_rows; i++) {
        for (j = 0; j < model->n_cols; j++) {
            u = gsl_cdf_gaussian_P(gsl_matrix_get(model->yi, i, j) -
                                   gsl_vector_get(model->mu, j),
                                   gsl_vector_get(model->theta, j));

            gsl_matrix_set(model->xi, i, j,
                           u_quantile(u, model->z, model->sigma, j));
        }
    }
}

static void
draw_new_missing_values(gsl_rng *rng, mico_model_t *model) {
    int i, j, k;
    double yc, fc, yp, up, xp, fp; // current and proposed values
    double ratio, uniform_draw;

    gsl_vector *xvecc = gsl_vector_alloc(model->n_cols);
    gsl_vector *xvecp = gsl_vector_alloc(model->n_cols);

    for (i = model->n_complete; i < model->n_rows; i++) {
        for (j = 0; j < model->n_cols; j++) {
            if (isnan(gsl_matrix_get(model->y, i, j))) {
                yc = gsl_matrix_get(model->yi, i, j);
                yp = yc + gsl_ran_tdist(rng, model->mcmc_proposal_nu);
                up = gsl_cdf_gaussian_P(yp - gsl_vector_get(model->mu, j),
                                        gsl_vector_get(model->theta, j));
                xp = u_quantile(up, model->z, model->sigma, j);

                // now compute the density at each point

                for (k = 0; k < model->n_cols; k++) {
                    if (k == j) {
                        gsl_vector_set(xvecp, k, xp);
                    } else {
                        gsl_vector_set(xvecp, k, gsl_matrix_get(model->xi, i, k));
                    }
                }

                gsl_matrix_get_row(xvecc, model->xi, i);

                fc = x_copula_pdf(xvecc, model->z, model->sigma);
                fp = x_copula_pdf(xvecp, model->z, model->sigma);

                fc *= gsl_ran_gaussian_pdf(yc - gsl_vector_get(model->mu, j),
                                           gsl_vector_get(model->theta, j));
                fp *= gsl_ran_gaussian_pdf(yp - gsl_vector_get(model->mu, j),
                                           gsl_vector_get(model->theta, j));

                ratio = fp / fc;
                if (ratio < 1) {
                    uniform_draw = gsl_rng_uniform(rng);
                    if (uniform_draw > ratio) {
                        gsl_matrix_set(model->yi, i, j, yp);
                        gsl_matrix_set(model->xi, i, j, xp);
                    }
                }
            }
        }
    }

    gsl_vector_free(xvecc);
    gsl_vector_free(xvecp);
}

static void
draw_new_margin_parameters(gsl_rng *rng, mico_model_t *model) {
    // draw marginal means / standard deviations from a normal-inverse-gamma
    // distribution.    

    double this_mu, next_mu, this_theta, next_theta;
    double theta_a, theta_b, mu_a, mu_b;
    double this_f, next_f, next_u, ratio, uniform_draw;
    gsl_vector *this_x = gsl_vector_alloc(model->n_cols);
    gsl_vector *next_x = gsl_vector_alloc(model->n_cols);
    int i, j;

    for (j = 0; j < model->n_cols; j++) {
        this_mu = gsl_vector_get(model->mu, j);
        this_theta = gsl_vector_get(model->theta, j);

        // what should be the proposal density for (mu, theta)?
        // for now, just drawing values from prior
        theta_a = gsl_vector_get(model->theta_a, j);
        theta_b = gsl_vector_get(model->theta_b, j);
        next_theta = 1.0 / gsl_ran_gamma(rng, theta_a, 1.0 / theta_b);

        mu_a = gsl_vector_get(model->mu_a, j);
        mu_b = gsl_vector_get(model->mu_b, j);
        next_mu = mu_a + gsl_ran_gaussian(rng, next_theta / sqrt(mu_b));

        //the RHS should be proportional to the normal-inverse-gamma density
        this_f = exp(-(2 * theta_b + mu_b * pow(this_mu - mu_a, 2.0)) /
                     (2 * pow(this_theta, 2.0)));
        next_f = exp(-(2 * theta_b + mu_b * pow(next_mu - mu_a, 2.0)) /
                     (2 * pow(next_theta, 2.0)));

        for (i = 0; i < model->n_rows; i++) {
            this_f *= gsl_ran_gaussian_pdf(gsl_matrix_get(model->yi, i, j) - this_mu,
                                           this_theta);
            next_f *= gsl_ran_gaussian_pdf(gsl_matrix_get(model->yi, i, j) - next_mu,
                                           next_theta);

            gsl_matrix_get_row(this_x, model->xi, i);
            gsl_matrix_get_row(next_x, model->xi, i);

            next_u = gsl_cdf_gaussian_P(gsl_matrix_get(model->yi, i, j) - this_mu,
                                        next_theta);
            gsl_vector_set(next_x, j, u_quantile(next_u, model->z, model->sigma, j));

            this_f *= x_copula_pdf(this_x, model->z, model->sigma);
            next_f *= x_copula_pdf(next_x, model->z, model->sigma);
        }

        ratio = next_f / this_f;
        if (ratio < 1) {
            uniform_draw = gsl_rng_uniform(rng);
            if (uniform_draw > ratio) {
                gsl_vector_set(model->mu, j, next_mu);
                gsl_vector_set(model->theta, j, next_theta);
            }
        }
    }

    gsl_vector_free(this_x);
    gsl_vector_free(next_x);

    // draw new mixture means
}

static void
draw_new_mixture_means(gsl_rng *rng, mico_model_t *model) {
    int i, j, k;
    double weights[model->n_comp];

    // Assign x's to z's
    for (i = 0; i < model->z->size1; i++) {
        for (k = 0; k < model->n_comp; k++) {
            // compute the probability that x[i] was drawn from z[k]
            weights[j] = 1;
            for (j = 0; j < model->n_cols; j++) {
                weights[k] *= gsl_ran_gaussian_pdf(gsl_matrix_get(model->xi, i, j) -
                                                   gsl_matrix_get(model->z, k, j),
                                                   model->sigma);
            }
        }

        gsl_vector_int_set(model->zind, i, banmi_sample(rng, model->n_comp, weights));
    }

    // propose new centers

    double u, x, 
           zp[model->n_comp][model->n_cols]; // proposed z

    for (k = 0; k < model->n_comp; k++) {
        u = gsl_ran_chisq(rng, model->x_nu);
        for (j = 0; j < model->n_cols; j++) {
            x = gsl_ran_gaussian(rng, 1.0) * sqrt(model->x_nu / u);
            zp[k][j] = x + gsl_matrix_get(model->z, k, j);
        }
    }

    // Take one M-H step using a t-distribution proposal
    double wc[model->n_comp], // weight at current z
           wp[model->n_comp]; // weight at proposed z

    for (k = 0; k < model->n_comp; k++) {
        wc[k] = wp[k] = 1.0;
        for (j = 0; j < model->n_cols; j++) {
            wc[k] += pow(gsl_matrix_get(model->z, k, j), 2.0) / model->x_nu;
            wp[k] += pow(zp[k][j], 2.0) / model->x_nu;
        }
        wc[k] = pow(wc[k], -(model->x_nu + model->n_cols) / 2.0);
        wp[k] = pow(wp[k], -(model->x_nu + model->n_cols) / 2.0);
    }

    for (i = 0; i < model->n_rows; i++) {
        k = gsl_vector_int_get(model->zind, i);
        for (j = 0; j < model->n_cols; j++) {
            wc[k] *= gsl_ran_gaussian_pdf(gsl_matrix_get(model->xi, i, j) -
                                          gsl_matrix_get(model->z, k, j),
                                          model->sigma);
            wp[k] *= gsl_ran_gaussian_pdf(gsl_matrix_get(model->xi, i, j) -
                                          zp[k][j],
                                          model->sigma);
        }
    }

    double ratio, uniform_draw;
    for (k = 0; k < model->n_comp; k++) {
        ratio = wp[k] / wc[k];
        if (ratio < 1) {
            uniform_draw = gsl_rng_uniform(rng);
            if (uniform_draw > ratio) {
                for (j = 0; j < model->n_cols; j++) {
                    gsl_matrix_set(model->z, k, j, zp[k][j]);
                }
            }
        }
    }

    // draw a new sigma

    double sigma_a_p = model->sigma_a + model->n_rows / 2.0;
    double sigma_b_p = model->sigma_b;

    for (i = 0; i < model->n_rows; i++) {
        k = gsl_vector_int_get(model->zind, i);
        for (j = 0; j < model->n_cols; j++) {
            sigma_b_p += pow(gsl_matrix_get(model->xi, i, j) - 
                             gsl_matrix_get(model->z, k, j), 2.0) / 2.0;
        }
    }

    model->sigma = gsl_ran_gamma(rng, sigma_a_p, 1.0 / sigma_b_p);
}

static void
draw_new_parameters(gsl_rng *rng, mico_model_t *model) {
    draw_new_margin_parameters(rng, model);
    draw_new_mixture_means(rng, model);
}


void
mico_data_augmentation(gsl_rng *rng, mico_model_t *model, int n_iter) {

    if (!model->is_initialized) {
        sort_data_by_missingness_pattern(model);
        init_parameters(rng, model);
        init_missing_values(rng, model);
        init_latent_variables(rng, model);
    }

    int t;
    for (t = 0; t < n_iter; t++) {
        draw_new_missing_values(rng, model);
        draw_new_parameters(rng, model);
    }
}

int
main(void) {
    gsl_matrix *z = gsl_matrix_alloc(2, 2);
    double sigma = 1.0;

    gsl_matrix_set(z, 0, 0, 0);
    gsl_matrix_set(z, 0, 1, 0);
    gsl_matrix_set(z, 1, 0, 1);
    gsl_matrix_set(z, 1, 1, 1);

    double u, diff;
    for (u = 0.001; u < 1; u += 0.001) {
        printf("%f ", u);
        diff = u - x_cdf(u_quantile(u, z, sigma, 0), z, sigma, 0);
        assert(abs(diff) < QuantileTolerance);
    }
    printf("\n");

    gsl_matrix_free(z);

    return 1;
}
