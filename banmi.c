#include "banmi.h"

int
banmi_to_ordered_value(double x, int max_value) {
    return (int) ((max_value + 1) * x);
}

double
banmi_from_ordered_value(int o, int max_value) {
    return (o + 0.5) / (max_value + 1.0);
}

// Allocate storage for a banmi_model struct, set explicit parameters.
//
banmi_model_t*
new_banmi_model(int max_rows, gsl_vector_int *bds_disc, int n_cont, 
                double dp_weight, double init_crosstab, 
                double lambda_a, double lambda_b) 
{

    banmi_model_t *model = malloc(sizeof(banmi_model_t));

    model->disc = gsl_matrix_int_alloc(max_rows, bds_disc->size);
    model->disc_imp = gsl_matrix_int_alloc(max_rows, bds_disc->size);
    model->cont = gsl_matrix_alloc(max_rows, n_cont);
    model->cont_imp = gsl_matrix_alloc(max_rows, n_cont);

    model->x = gsl_matrix_int_alloc(max_rows, bds_disc->size);
    model->mu = gsl_matrix_alloc(max_rows, n_cont);
    model->lambda = gsl_vector_alloc(bds_disc->size);
    model->sigma = gsl_vector_alloc(n_cont);

    int disc_dim[bds_disc->size];
    int i;
    for (i = 0; i < bds_disc->size; i++)
        disc_dim[i] = gsl_vector_int_get(bds_disc, i);

    model->crosstab = alloc_tab(bds_disc->size, disc_dim);

    model->dp_weight = dp_weight;
    model->init_crosstab = init_crosstab;
    model->mu_a = gsl_vector_alloc(n_cont);
    model->mu_b = gsl_vector_alloc(n_cont);
    model->sigma_a = gsl_vector_alloc(n_cont);
    model->sigma_b = gsl_vector_alloc(n_cont);
    model->lambda_a = lambda_a;
    model->lambda_b = lambda_b;

    model->n_rows = 0;
    model->bds_disc = bds_disc;
    model->n_disc = bds_disc->size;
    model->n_cont = n_cont;

    // this needs to be in sync with missingness_pattern
    int mask = 1;
    model->mask_missing_disc_data = 0;
    model->mask_missing_cont_data = 0;
    for (i = 0; i < model->n_disc; i++) {
        model->mask_missing_disc_data |= mask;
        mask = mask << 1;
    }
    for (i = 0; i < model->n_cont; i++) {
        model->mask_missing_cont_data |= mask;
        mask = mask << 1;
    }
    model->mis_pat = malloc(sizeof(int) * max_rows);

    model->is_initialized = 0;

    return model;
}

// Copy values from from foo to foo_imp
//
void
init_imputed_data(banmi_model_t *model) {
    int i, j;
    
    for (i = 0; i < model->n_rows; i++) {
        for (j = 0; j < model->cont->size2; j++)
            gsl_matrix_set(model->cont_imp, i, j, gsl_matrix_get(model->cont, i, j));

        for (j = 0; j < model->disc->size2; j++)
            gsl_matrix_int_set(model->disc_imp, i, j, 
                               gsl_matrix_int_get(model->disc, i, j));
    }
}

// Compute the missingness pattern for the given row of data.
//
int 
missingness_pattern(const banmi_model_t *model, int row) {
    int i, pattern = 0, mask = 1;

    for (i = 0; i < model->n_disc; i++) {
        if (gsl_matrix_int_get(model->disc, row, i) < 0)
            pattern |= mask;

        mask = mask << 1;
    }

    for (i = 0; i < model->n_cont; i++) {
        if (gsl_matrix_get(model->cont, row, i) < 0)
            pattern |= mask;

        mask = mask << 1;
    }

    return pattern;
}

// convert a missingness pattern to an array of margin indices.
// return the number of margins with missing data.
// Example:
//   [0 -1  0 -1  1]  => [1  3]
int 
missingness_pattern_to_index(int *ix, int n_disc, int pattern) {
    int i, j = 0, mask = 1;
    for (i = 0; i < n_disc; i++) {
        if ((pattern & mask) > 0) 
            ix[j++] = i;

        mask = mask << 1;
    }

    return j;
}

// Sort the rows of data by missingness pattern. Additionally, set the
// n_complete and n_missing fields. Returns the number of complete rows.
//
int
sort_data_by_missingness_pattern(banmi_model_t *model) {
    int i, n_complete = 0;

    for (i = 0; i < model->n_rows; i++) {
        model->mis_pat[i] = missingness_pattern(model, i);

        if (model->mis_pat[i] == 0)
            n_complete++;
    }

    order_blocks(model->disc->data, model->mis_pat, model->n_rows, model->disc->tda);
    order_blocks(model->disc_imp->data, model->mis_pat, model->n_rows, model->disc_imp->tda);
    order_blocks_d(model->cont->data, model->mis_pat, model->n_rows, model->cont->tda);
    order_blocks_d(model->cont_imp->data, model->mis_pat, model->n_rows, model->cont_imp->tda);
    order(model->mis_pat, model->mis_pat, model->n_rows);

    model->n_complete = n_complete;
    model->n_missing = model->n_rows - n_complete;

    return n_complete;
}

// Compute crosstab, mu_a, mu_b, sigma_a, sigma_b from the data. Recall that
// lambda_a, lambda_b, and dp_weight are set explicitly by the user in
// new_banmi_model.
//
// This function should be called after sort_data_by_missingness_pattern, since
// that function sets model->n_complete.
//
void init_hyperparameters(banmi_model_t *model) { 
    int i, j, crosstab_ix[model->n_disc];

    // compute the crosstab of the complete rows

    for (i = 0; i < model->crosstab->size; i++)
        model->crosstab->dat[i] = model->init_crosstab;

    for (i = 0; i < model->n_complete; i++) {

        for (j = 0; j < model->n_disc; j++)
            crosstab_ix[j] = gsl_matrix_int_get(model->disc, i, j);

        tab_set(model->crosstab, crosstab_ix, 
                tab_get(model->crosstab, crosstab_ix) + 1);
    }

    // take the distribution of means for the continuous variables to have a
    // beta distribution.  Note the mean computation uses all known values, not
    // just the ones in complete data rows.
    //

    double mean, mean_shape, var, u, temp[model->n_rows];
    int insert_index;
    for (j = 0; j < model->n_cont; j++) {
        insert_index = 0;
        for (i = 0; i < model->n_rows; i++) {
            if ((u = gsl_matrix_get(model->cont, i, j)) >= 0)
                temp[insert_index++] = u;
        }

        mean = gsl_stats_mean(temp, 1, insert_index);
        var = gsl_stats_variance(temp, 1, insert_index);
        
        gsl_vector_set(model->mu_a, j,
                       mean * (mean * (1 - mean) / var - 1));
        gsl_vector_set(model->mu_b, j,
                       (1 - mean) * (mean * (1 - mean) / var - 1));

        // set sigma_a and sigma_b to have a mean given by Silverman's rule
        // and a weight equal to the number of complete observations
        mean_shape = 1.06 * pow(var, 0.5) * pow(insert_index, -0.2);
        gsl_vector_set(model->sigma_a, j, (insert_index+0.0) / 2.0);
        gsl_vector_set(model->sigma_b, j, 2.0 * mean_shape / insert_index);
    }
}

// Set initial values for missing values in the imputed data set. This function
// should be called after init_hyperparams.
//
void
init_missing_values(gsl_rng *rng, banmi_model_t *model) {
    int i, j, last_pattern, this_pattern, n_missing, flat_ix,
        missing_ix[model->n_disc], marg_ix[model->n_disc];
    tab_t *margtab = NULL;

    // fill in missing discrete values

    last_pattern = 0;
    for (i = model->n_complete; i < model->n_rows; i++) {

        this_pattern = model->mis_pat[i];

        if (this_pattern & model->mask_missing_disc_data) {
            if (this_pattern != last_pattern) {
                // rebuild the marginal count table
                last_pattern = this_pattern;
                n_missing = missingness_pattern_to_index(missing_ix, 
                                                         model->n_disc, 
                                                         this_pattern);

                if (margtab != NULL)
                    free_tab(margtab);

                margtab = new_marginal_tab(model->crosstab, missing_ix, n_missing);
            }

            // sample from the marginal distribution and copy result to disc_imp
            flat_ix = sample_d(rng, margtab->size, margtab->dat);
            tab_array_index(marg_ix, margtab, flat_ix);
            for (j = 0; j < margtab->n; j++) {
                gsl_matrix_int_set(model->disc_imp, i, missing_ix[j], marg_ix[j]);
            }
        }

        if (this_pattern & model->mask_missing_cont_data) {
            for (j = 0; j < model->n_cont; j++) {
                if (gsl_matrix_get(model->cont, i, j) < 0)
                    gsl_matrix_set(model->cont_imp, i, j, 
                                   gsl_ran_beta(rng, 
                                                gsl_vector_get(model->mu_a, j), 
                                                gsl_vector_get(model->mu_b, j)));
            }
        }
    }

}

double 
kernel(const gsl_vector_int *bds_disc,
       const gsl_vector *u, const gsl_vector_int *w, 
       const gsl_vector *mu, const gsl_vector_int *x, 
       const gsl_vector *sigma, const gsl_vector *lambda) 
{
    int i;
    double result = 1.0;

    for (i = 0; i < bds_disc->size; i++) {
        if (gsl_vector_int_get(w, i) == gsl_vector_int_get(x, i))
            result *= 1 - gsl_vector_get(lambda, i);
        else
            result *= gsl_vector_get(lambda, i) / (gsl_vector_int_get(bds_disc, i) - 1);
    }

    double ui, y, mean, sd;
    for (i = 0; i < u->size; i++) {
        ui = gsl_vector_get(u, i);
        mean = gsl_vector_get(mu, i);
        sd = gsl_vector_get(sigma, i);
        y = gsl_ran_gaussian_pdf(ui - mean, sd);
#ifdef BoundaryCorrection
        // for now, assume reflecting once in each direction is sufficient
        y += gsl_ran_gaussian_pdf(2 - ui - mean, sd) + 
             gsl_ran_gaussian_pdf(-ui - mean, sd);
#endif
        result *= y;
    }

    return result;
}

// Initialize the means/modes of the mixture model. This function should be 
// called after init_hyperparams.
//
void
init_latent_variables(gsl_rng *rng, banmi_model_t *model) {
    int i, j, flat_ix, crosstab_ix[model->n_disc];

    for (i = 0; i < model->n_rows; i++) {
        for (j = 0; j < model->n_cont; j++) {
            gsl_matrix_set(model->mu, i, j,
                           gsl_ran_beta(rng, gsl_vector_get(model->mu_a, j), 
                                             gsl_vector_get(model->mu_b, j)));
        }

        flat_ix = sample_d(rng, model->crosstab->size, model->crosstab->dat);
        tab_array_index(crosstab_ix, model->crosstab, flat_ix);
        for (j = 0; j < model->n_disc; j++) {
            gsl_matrix_int_set(model->x, i, j, crosstab_ix[j]);
        }
    }

    double s;
    for (i = 0; i < model->n_cont; i++) {
        s = sqrt(1.0 / gsl_ran_gamma(rng, gsl_vector_get(model->sigma_a, i), 
                                          gsl_vector_get(model->sigma_b, i)));
        gsl_vector_set(model->sigma, i, s);
    }

    for (i = 0; i < model->n_disc; i++) {
        s = gsl_ran_beta(rng, model->lambda_a, model->lambda_b);
        gsl_vector_set(model->lambda, i, s);
    }
}

void
draw_new_latent_variables(gsl_rng *rng, banmi_model_t *model) {

    int i, j, k, flat_ix, choice, crosstab_ix[model->n_disc];
    double weight[model->n_rows];

    gsl_vector_int *w_row = gsl_vector_int_alloc(model->n_disc);
    gsl_vector_int *x_row = gsl_vector_int_alloc(model->n_disc);
    gsl_vector_int *i_row = gsl_vector_int_alloc(model->n_disc);
    gsl_vector *u_row = gsl_vector_alloc(model->n_cont);
    gsl_vector *mu_row = gsl_vector_alloc(model->n_cont);
    gsl_vector *a_row = gsl_vector_alloc(model->n_cont);

    // draw new means/modes

    for (i = 0; i < model->n_rows; i++) {
        gsl_matrix_int_get_row(w_row, model->disc_imp, i);
        gsl_matrix_int_get_row(x_row, model->x, i);
        gsl_matrix_get_row(u_row, model->cont_imp, i);
        gsl_matrix_get_row(mu_row, model->mu, i);

        weight[0] = model->dp_weight * kernel(model->bds_disc,
                                              u_row, w_row,
                                              mu_row, x_row,
                                              model->sigma,
                                              model->lambda);

        k = 1; // insertion index
        for (j = 0; j < model->n_rows; j++) {
            if (j == i)
                continue;

            gsl_matrix_int_get_row(x_row, model->x, j);
            gsl_matrix_get_row(mu_row, model->mu, j);

            weight[k] = kernel(model->bds_disc,
                               u_row, w_row,
                               mu_row, x_row,
                               model->sigma,
                               model->lambda);

            k++;
        }

        choice = sample_d(rng, model->n_rows, weight);

        if (choice == 0) {
            // sample from G_0
            for (j = 0; j < model->n_cont; j++) 
                gsl_matrix_set(model->mu, i, j,
                               gsl_ran_beta(rng, gsl_vector_get(model->mu_a, j), 
                                                 gsl_vector_get(model->mu_b, j)));

            flat_ix = sample_d(rng, model->crosstab->size, model->crosstab->dat);
            tab_array_index(crosstab_ix, model->crosstab, flat_ix);
            for (j = 0; j < model->n_disc; j++) {
                gsl_matrix_int_set(model->x, i, j, crosstab_ix[j]);
            }
        } else {
            // set (x[i], mu[i]) to (x[choice], mu[choice])
            if (choice <= i)
                choice--;

            gsl_matrix_int_get_row(i_row, model->x, choice);
            gsl_matrix_int_set_row(model->x, i, i_row);

            gsl_matrix_get_row(a_row, model->mu, choice);
            gsl_matrix_set_row(model->mu, i, a_row);
        }

    }

    gsl_vector_int_free(w_row);
    gsl_vector_int_free(x_row);
    gsl_vector_int_free(i_row);
    gsl_vector_free(u_row);
    gsl_vector_free(mu_row);
    gsl_vector_free(a_row);

    // draw new sigma and lambda

    double sigma_a_post, sigma_b_post;
    double diff;

    for (j = 0; j < model->n_cont; j++) {
        sigma_a_post = gsl_vector_get(model->sigma_a, j) + model->n_rows / 2.0;
        sigma_b_post = gsl_vector_get(model->sigma_b, j);

        for (i = 0; i < model->n_rows; i++) {
            diff = gsl_matrix_get(model->cont_imp, i, j) - 
                   gsl_matrix_get(model->mu, i, j);
            sigma_b_post += pow(diff, 2.0) / 2.0;
        }
        sigma_b_post = 1.0 / sigma_b_post;

        gsl_vector_set(model->sigma, j,
                       sqrt(1.0 / gsl_ran_gamma(rng, sigma_a_post, sigma_b_post)));
    }

    double lambda_a_post, lambda_b_post;
    for (j = 1; j < model->n_disc; j++) {
        lambda_a_post = model->lambda_a;
        lambda_b_post = model->lambda_b;

        for (i = 0; i < model->n_rows; i++) {
            if (gsl_matrix_int_get(model->disc_imp, i, j) != gsl_matrix_int_get(model->x, i, j))
                lambda_a_post++;
            else
                lambda_b_post++;
        }

        gsl_vector_set(model->lambda, j,
                       gsl_ran_beta(rng, lambda_a_post, lambda_b_post));
    }
}

void
draw_new_missing_values(gsl_rng *rng, banmi_model_t *model) {
    int i, j, choice;
    double u, mu, sigma, cont_val;

    for (i = model->n_complete; i < model->n_rows; i++) {

        if (model->mis_pat[i] & model->mask_missing_cont_data) {
            // there are missing discrete values
            for (j = 0; j < model->n_disc; j++) {
                if (gsl_matrix_int_get(model->disc, i, j) < 0) {
                    u = gsl_rng_uniform(rng);
                    if (u > 1 - gsl_vector_get(model->lambda, j)) {
                        gsl_matrix_int_set(model->disc_imp, i, j,
                                           gsl_matrix_int_get(model->x, i, j));
                    } else {
                        choice = (int) (gsl_rng_uniform(rng) * (gsl_vector_int_get(model->bds_disc, j) - 1));

                        if (choice >= gsl_matrix_int_get(model->x, i, j))
                            choice++;

                        gsl_matrix_int_set(model->disc_imp, i, j, choice);
                    }
                }
            }
        }

        if (model->mis_pat[i] & model->mask_missing_cont_data) {
            // there are missing continuous values
            for (j = 0; j < model->n_cont; j++) {
                sigma = gsl_vector_get(model->sigma, j);
                if (gsl_matrix_get(model->cont, i, j) < 0) {
                    mu = gsl_matrix_get(model->mu, i, j);
                    cont_val = mu + gsl_ran_gaussian(rng, sigma);
#ifdef BoundaryCorrection
                    cont_val = fmod(cont_val, 2.0);
                    if (cont_val > 1)
                        cont_val = 2 - cont_val;
                    else if (-1 <= cont_val && cont_val < 0)
                        cont_val = -cont_val;
                    else if (cont_val < -1)
                        cont_val = cont_val - 2;
#endif
                    gsl_matrix_set(model->cont_imp, i, j, cont_val);
                }
            }
        }
    }
}

void 
banmi_data_augmentation(gsl_rng *rng, banmi_model_t *model, int n_iter) {

    if (!model->is_initialized) {
        // initialize the model before first augmentation
        init_imputed_data(model);
        sort_data_by_missingness_pattern(model);
        init_hyperparameters(model);
        init_missing_values(rng, model);
        init_latent_variables(rng, model);
        model->is_initialized = 1;
    }

    // perform the data augmentation

    int t;
    for (t = 0; t < n_iter; t++) {
        draw_new_latent_variables(rng, model);
        draw_new_missing_values(rng, model);
    }
}

void
banmi_print_shape_variables(banmi_model_t *model) {
    int i;
    printf("lambda ");

    for (i = 0; i < model->n_disc; i++)
        printf("%4.2f ", gsl_vector_get(model->lambda, i));

    printf("\nsigma ");
    for (i = 0; i < model->n_cont; i++)
        printf("%4.2f ", gsl_vector_get(model->sigma, i));

    printf("\n");
}

int
banmi_count_unique_modes(banmi_model_t *model) {
    gsl_matrix *mu = gsl_matrix_alloc(model->n_rows, model->n_cont);
    gsl_matrix_int *x = gsl_matrix_int_alloc(model->n_rows, model->n_disc);

    int i, j, k;

    // copy the first row of mu, x
    for (j = 0; j < model->n_cont; j++)
        gsl_matrix_set(mu, 0, j, gsl_matrix_get(model->mu, 0, j));

    for (j = 0; j < model->n_disc; j++)
        gsl_matrix_int_set(x, 0, j, gsl_matrix_int_get(model->x, 0, j));

    int match_found, insert_index = 1;

    for (i = 1; i < model->n_rows; i++) {

        for (k = 0; k < insert_index; k++) {
            match_found = 1;

            for (j = 0; j < model->n_cont; j++) {
                if (gsl_matrix_get(model->mu, i, j) != gsl_matrix_get(mu, k, j)) {
                    match_found = 0;
                    break;
                }
            }

            for (j = 0; j < model->n_disc; j++) {
                if (gsl_matrix_int_get(model->x, i, j) != gsl_matrix_int_get(x, k, j)) {
                    match_found = 0;
                    break;
                }
            }

            if (match_found)
                break;
        }

        if (!match_found) {
            for (j = 0; j < model->n_cont; j++)
                gsl_matrix_set(mu, insert_index, j, 
                               gsl_matrix_get(model->mu, i, j));

            for (j = 0; j < model->n_disc; j++)
                gsl_matrix_int_set(x, insert_index, j, 
                                   gsl_matrix_int_get(model->x, i, j));

            insert_index++;
        }
    }

    gsl_matrix_free(mu);
    gsl_matrix_int_free(x);

    return insert_index;
}
