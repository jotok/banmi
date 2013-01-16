#include "banmi.h"

// Allocate storage for a banmi_model struct, set explicit parameters.
//
banmi_model_t*
new_banmi_model(int max_rows, gsl_vector_int *bds_disc, gsl_vector_int *bds_orde,
                int n_cont, double dp_weight, 
                double lambda_a, double lambda_b) 
{

    banmi_model_t *model = malloc(sizeof(banmi_model_t));

    int data_dim[2];
    data_dim[0] = max_rows; data_dim[1] = bds_disc->size;

    model->disc = alloc_tab(2, data_dim);
    model->disc_imp = alloc_tab(2, data_dim);
    model->orde = gsl_matrix_int_alloc(max_rows, bds_orde->size);
    model->orde_imp = gsl_matrix_int_alloc(max_rows, bds_orde->size);
    model->cont = gsl_matrix_alloc(max_rows, n_cont);
    model->cont_imp = gsl_matrix_alloc(max_rows, n_cont);

    model->x = alloc_tab(2, data_dim);
    model->xo = gsl_matrix_int_alloc(max_rows, bds_orde->size);
    model->mu = gsl_matrix_alloc(max_rows, n_cont);
    model->lambda = gsl_vector_alloc(bds_disc->size);
    model->kappa = gsl_vector_alloc(bds_orde->size);
    model->sigma = gsl_vector_alloc(n_cont);

    int disc_dim[bds_disc->size];
    int i;
    for (i = 0; i < bds_disc->size; i++)
        disc_dim[i] = (int)gsl_vector_int_get(bds_disc, i);

    model->crosstab = alloc_tab(bds_disc->size, disc_dim);

    model->dp_weight = dp_weight;
    model->mu_a = gsl_vector_alloc(n_cont);
    model->mu_b = gsl_vector_alloc(n_cont);
    model->xo_a = gsl_vector_alloc(bds_orde->size);
    model->xo_b = gsl_vector_alloc(bds_orde->size);
    model->sigma_a = gsl_vector_alloc(n_cont);
    model->sigma_b = gsl_vector_alloc(n_cont);
    model->kappa_a = gsl_vector_alloc(bds_orde->size);
    model->kappa_b = gsl_vector_alloc(bds_orde->size);
    model->lambda_a = lambda_a;
    model->lambda_b = lambda_b;

    model->bds_disc = bds_disc;
    model->n_disc = bds_disc->size;
    model->bds_orde = bds_orde;
    model->n_orde = bds_orde->size;
    model->n_cont = n_cont;

    // this needs to be in sync with missingness_pattern
    int mask = 1;
    model->mask_missing_disc_data = 0;
    model->mask_missing_orde_data = 0;
    model->mask_missing_cont_data = 0;
    for (i = 0; i < model->n_disc; i++) {
        model->mask_missing_disc_data |= mask;
        mask = mask << 1;
    }
    for (i = 0; i < model->n_orde; i++) {
        model->mask_missing_orde_data |= mask;
        mask = mask << 1;
    }
    for (i = 0; i < model->n_cont; i++) {
        model->mask_missing_cont_data |= mask;
        mask = mask << 1;
    }

    model->mis_pat = malloc(sizeof(int) * max_rows);

    return model;
}

// compute the missingness pattern for the given row of data
int 
missingness_pattern(const banmi_model_t *model, int row) {
    int pattern = 0, mask = 1, i, ix[2];
    ix[0] = row;

    for (i = 0; i < model->n_disc; i++) {
        ix[1] = i;
        if (tab_get(model->disc, ix) < 0)
            pattern |= mask;

        mask = mask << 1;
    }

    for (i = 0; i < model->n_orde; i++) {
        if (gsl_matrix_int_get(model->orde, row, i) < 0)
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

    order_blocks(model->disc->dat, model->mis_pat, model->n_rows, model->n_disc);
    order_blocks(model->disc_imp->dat, model->mis_pat, model->n_rows, model->n_disc);
    order_rows_int(model->orde, model->mis_pat, model->n_rows);
    order_rows_int(model->orde_imp, model->mis_pat, model->n_rows);
    order_rows(model->cont, model->mis_pat, model->n_rows);
    order_rows(model->cont_imp, model->mis_pat, model->n_rows);
    order(model->mis_pat, model->mis_pat, model->n_rows);

    model->n_complete = n_complete;
    model->n_missing = model->n_rows - n_complete;

    return n_complete;
}

// Compute crosstab, mu_a, and mu_b from the data. Recall that lambda_a,
// lambda_b, kappa_a, kappa_b, sigma_a, sigma_b, and dp_weight are set
// explicitly by the user in new_banmi_model.
//
// This function should be called after sort_data_by_missingness_pattern, since
// that function sets model->n_complete.
//
void
init_hyperparameters(banmi_model_t *model) {
    int i, j, disc_ix[2], crosstab_ix[model->n_disc];

    // compute the crosstab of the complete rows

    for (i = 0; i < model->crosstab->size; i++)
        model->crosstab->dat[i] = 1;

    for (i = 0; i < model->n_complete; i++) {
        disc_ix[0] = i;

        for (j = 0; j < model->n_disc; j++) {
            disc_ix[1] = j;
            crosstab_ix[j] = tab_get(model->disc, disc_ix);
        }

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

    int bd, wo;
    for (j = 0; j < model->n_orde; j++) {
        insert_index = 0;
        bd = gsl_vector_int_get(model->bds_orde, j);
        for (i = 0; i < model->n_rows; i++) {
            if ((wo = gsl_matrix_int_get(model->orde, i, j)) >= 0)
                temp[insert_index++] = (wo+0.0) / bd;
        }

        mean = gsl_stats_mean(temp, 1, insert_index);
        var = gsl_stats_variance(temp, 1, insert_index);

        gsl_vector_set(model->xo_a, j,
                       mean * (mean * (1 - mean) / var - 1));
        gsl_vector_set(model->xo_b, j,
                       (1 - mean) * (mean * (1 - mean) / var - 1));

        // set note regarding sigma_a, sigma_b above
        mean_shape = 1.06 * pow(var, 0.5) * pow(insert_index, -0.2);
        gsl_vector_set(model->kappa_a, j, (insert_index+0.0) / 2.0);
        gsl_vector_set(model->kappa_b, j, 2.0 * mean_shape / insert_index);
    }

}

// Set initial values for missing values in the imputed data set. This function
// should be called after init_hyperparams.
//
void
init_missing_values(gsl_rng *rng, banmi_model_t *model) {
    int i, j, last_pattern, this_pattern, n_missing, flat_ix, disc_ix[2], 
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
            flat_ix = sample(rng, margtab->size, margtab->dat);
            tab_array_index(marg_ix, margtab, flat_ix);
            disc_ix[0] = i;
            for (j = 0; j < margtab->n; j++) {
                disc_ix[1] = missing_ix[j];
                tab_set(model->disc_imp, disc_ix, marg_ix[j]);
            }
        }

        double this_orde;
        if (this_pattern & model->mask_missing_orde_data) {
            for (j = 0; j < model->n_orde; j++) {
                if (gsl_matrix_int_get(model->orde, i, j) < 0) {
                    this_orde = gsl_ran_beta(rng, gsl_vector_get(model->xo_a, j), 
                                                  gsl_vector_get(model->xo_b, j));
                    this_orde *= gsl_vector_int_get(model->bds_orde, j);
                    gsl_matrix_int_set(model->orde_imp, i, j, (int)this_orde);
                }
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
kernel(const gsl_vector_int *bds_disc, const gsl_vector_int *bds_orde, 
       const gsl_vector *u, const gsl_vector_int *wo, const int *w, 
       const gsl_vector *mu, const gsl_vector_int *xo, const int *x, 
       const gsl_vector *sigma, const gsl_vector *kappa, const gsl_vector *lambda) 
{
    int i;
    double result = 1.0;

    for (i = 0; i < bds_disc->size; i++) {
        if (w[i] == x[i])
            result *= 1 - gsl_vector_get(lambda, i);
        else
            result *= gsl_vector_get(lambda, i) / (gsl_vector_int_get(bds_disc, i) - 1);
    }

    for (i = 0; i < wo->size; i++) {
        result *= gsl_ran_gaussian_pdf(gsl_vector_int_get(wo, i) - gsl_vector_int_get(xo, i),
                                       gsl_vector_get(kappa, i));
    }

    for (i = 0; i < u->size; i++) {
        result *= gsl_ran_gaussian_pdf(gsl_vector_get(u, i) - gsl_vector_get(mu, i),
                                       gsl_vector_get(sigma, i));
    }

    return result;
}

// Initialize the means/modes of the mixture model. This function should be 
// called after init_hyperparams.
//
void
init_latent_variables(gsl_rng *rng, banmi_model_t *model) {
    int i, j, flat_ix, x_ix[2], crosstab_ix[model->n_disc];

    for (i = 0; i < model->n_rows; i++) {
        for (j = 0; j < model->n_cont; j++) {
            gsl_matrix_set(model->mu, i, j,
                           gsl_ran_beta(rng, gsl_vector_get(model->mu_a, j), 
                                             gsl_vector_get(model->mu_b, j)));
        }

        for (j = 0; j < model->n_orde; j++) {
            double beta = gsl_ran_beta(rng, gsl_vector_get(model->xo_a, j),
                                            gsl_vector_get(model->xo_b, j));
            gsl_matrix_int_set(model->xo, i, j,
                               (int)(beta * gsl_vector_int_get(model->bds_orde, j)));
        }

        x_ix[0] = i;
        flat_ix = sample(rng, model->crosstab->size, model->crosstab->dat);
        tab_array_index(crosstab_ix, model->crosstab, flat_ix);
        for (j = 0; j < model->n_disc; j++) {
            x_ix[1] = j;
            tab_set(model->x, x_ix, crosstab_ix[j]);
        }
    }

    double s;
    for (i = 0; i < model->n_cont; i++) {
        s = sqrt(1.0 / gsl_ran_gamma(rng, gsl_vector_get(model->sigma_a, i), 
                                          gsl_vector_get(model->sigma_b, i)));
        gsl_vector_set(model->sigma, i, s);
    }

    for (i = 0; i < model->n_orde; i++) {
        s = sqrt(1.0 / gsl_ran_gamma(rng, gsl_vector_get(model->kappa_a, i), 
                                          gsl_vector_get(model->kappa_b, i)));
        gsl_vector_set(model->kappa, i, s);
    }

    for (i = 0; i < model->n_disc; i++) {
        s = gsl_ran_beta(rng, model->lambda_a, model->lambda_b);
        gsl_vector_set(model->lambda, i, s);
    }
}

void
draw_new_latent_variables(gsl_rng *rng, banmi_model_t *model) {

    int i, j, k, flat_ix, choice, x_ix[model->n_disc], y_ix[model->n_disc], 
        crosstab_ix[model->n_disc];
    double weight[model->n_rows];

    gsl_vector *u_row = gsl_vector_alloc(model->n_cont);
    gsl_vector *mu_row = gsl_vector_alloc(model->n_cont);
    gsl_vector_int *wo_row = gsl_vector_int_alloc(model->n_orde);
    gsl_vector_int *xo_row = gsl_vector_int_alloc(model->n_orde);
    gsl_vector *a_row = gsl_vector_alloc(model->n_cont);
    gsl_vector_int *i_row = gsl_vector_int_alloc(model->n_orde);

    // draw new means/modes

    for (i = 0; i < model->n_rows; i++) {
        gsl_matrix_get_row(u_row, model->cont_imp, i);
        gsl_matrix_get_row(mu_row, model->mu, i);
        gsl_matrix_int_get_row(wo_row, model->orde_imp, i);
        gsl_matrix_int_get_row(xo_row, model->xo, i);

        weight[0] = model->dp_weight * kernel(model->bds_disc,
                                              model->bds_orde,
                                              u_row,
                                              wo_row,
                                              model->disc_imp->dat + i*model->n_disc,
                                              mu_row,
                                              xo_row,
                                              model->x->dat + i*model->n_disc,
                                              model->sigma,
                                              model->kappa,
                                              model->lambda);

        k = 1; // insertion index
        for (j = 0; j < model->n_rows; j++) {
            if (j == i)
                continue;

            gsl_matrix_get_row(u_row, model->cont_imp, i);
            gsl_matrix_get_row(mu_row, model->mu, j);
            gsl_matrix_int_get_row(wo_row, model->orde_imp, i);
            gsl_matrix_int_get_row(xo_row, model->xo, j);

            weight[k] = kernel(model->bds_disc,
                               model->bds_orde,
                               u_row,
                               wo_row,
                               model->disc_imp->dat + i*model->n_disc,
                               mu_row,
                               xo_row,
                               model->x->dat + j*model->n_disc,
                               model->sigma,
                               model->kappa,
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

            for (j = 0; j < model->n_orde; j++) {
                double this_xo = gsl_ran_beta(rng, gsl_vector_get(model->xo_a, j), 
                                                   gsl_vector_get(model->xo_b, j)); 
                this_xo *= gsl_vector_int_get(model->bds_orde, j);
                gsl_matrix_int_set(model->xo, i, j, (int)this_xo);
            }

            x_ix[0] = i;
            flat_ix = sample(rng, model->crosstab->size, model->crosstab->dat);
            tab_array_index(crosstab_ix, model->crosstab, flat_ix);
            for (j = 0; j < model->n_disc; j++) {
                x_ix[1] = j;
                tab_set(model->x, x_ix, crosstab_ix[j]);
            }
        } else {
            // set (x[i], mu[i]) to (x[choice], mu[choice])
            if (choice <= i)
                choice--;

            x_ix[0] = i;
            y_ix[0] = choice;
            for (j = 0; j < model->n_disc; j++) {
                x_ix[1] = j; y_ix[1] = j;
                tab_set(model->x, x_ix, tab_get(model->x, y_ix));
            }

            gsl_matrix_int_get_row(i_row, model->xo, choice);
            gsl_matrix_int_set_row(model->xo, i, i_row);

            gsl_matrix_get_row(a_row, model->mu, choice);
            gsl_matrix_set_row(model->mu, i, a_row);
        }

    }

    gsl_vector_free(u_row);
    gsl_vector_free(mu_row);
    gsl_vector_int_free(wo_row);
    gsl_vector_int_free(xo_row);
    gsl_vector_free(a_row);
    gsl_vector_int_free(i_row);

    // draw new sigma, kappa, and lambda

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

    double kappa_a_post, kappa_b_post;
    for (j = 0; j < model->n_orde; j++) {
        kappa_a_post = gsl_vector_get(model->kappa_a, j) + model->n_rows / 2.0;
        kappa_b_post = gsl_vector_get(model->kappa_b, j);

        for (i = 0; i < model->n_rows; i++) {
            diff = gsl_matrix_int_get(model->orde_imp, i, j) -
                   gsl_matrix_int_get(model->xo, i, j);
            kappa_b_post += pow(diff, 2.0) / 2.0;
        }
        kappa_b_post = 1.0 / kappa_b_post;

        gsl_vector_set(model->kappa, j,
                       sqrt(1.0 / gsl_ran_gamma(rng, kappa_a_post, kappa_b_post)));
    }

    double lambda_a_post, lambda_b_post;
    int disc_ix[2];
    for (j = 1; j < model->n_disc; j++) {
        lambda_a_post = model->lambda_a;
        lambda_b_post = model->lambda_b;
        disc_ix[1] = j;

        for (i = 0; i < model->n_rows; i++) {
            disc_ix[0] = i;
            if (tab_get(model->disc_imp, disc_ix) != tab_get(model->x, disc_ix))
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
    int i, j, xo, choice, disc_ix[2];
    double u, mu, sigma, kappa, base, frac;

    for (i = model->n_complete; i < model->n_rows; i++) {
        disc_ix[0] = i;

        if (model->mis_pat[i] & model->mask_missing_cont_data) {
            // there are missing discrete values
            for (j = 0; j < model->n_disc; j++) {
                disc_ix[1] = j;
                if (tab_get(model->disc, disc_ix) < 0) {
                    u = gsl_rng_uniform(rng);
                    if (u > 1 - gsl_vector_get(model->lambda, j)) {
                        tab_set(model->disc_imp, disc_ix, tab_get(model->x, disc_ix));
                    } else {
                        choice = (int) (gsl_rng_uniform(rng) * (gsl_vector_int_get(model->bds_disc, j) - 1));

                        if (choice >= tab_get(model->x, disc_ix))
                            choice++;

                        tab_set(model->disc_imp, disc_ix, choice);
                    }
                }
            }
        }

        if (model->mis_pat[i] & model->mask_missing_orde_data) {
            // there are missing ordered values
            for (j = 0; j < model->n_orde; j++) {
                kappa = gsl_vector_get(model->kappa, j);
                if (gsl_matrix_int_get(model->orde, i, j) < 0) {
                    xo = gsl_matrix_int_get(model->xo, i, j);
                    frac = modf(xo + gsl_ran_gaussian(rng, kappa), &base);

                    if (frac > 0.5)
                        base += 1;

                    gsl_matrix_int_set(model->orde_imp, i, j, (int)base);
                                   
                }
            }
        }

        if (model->mis_pat[i] & model->mask_missing_cont_data) {
            // there are missing continuous values
            for (j = 0; j < model->n_cont; j++) {
                sigma = gsl_vector_get(model->sigma, j);
                if (gsl_matrix_get(model->cont, i, j) < 0) {
                    mu = gsl_matrix_get(model->mu, i, j);
                    gsl_matrix_set(model->cont_imp, i, j,
                                   mu + gsl_ran_gaussian(rng, sigma));
                }
            }
        }
    }
}

void 
banmi_data_augmentation(gsl_rng *rng, banmi_model_t *model, int n_iter) {

    // initialize the model
    sort_data_by_missingness_pattern(model);
    init_hyperparameters(model);
    init_missing_values(rng, model);
    init_latent_variables(rng, model);

    // perform the data augmentation

    int t;
    for (t = 0; t < n_iter; t++) {
        draw_new_latent_variables(rng, model);
        draw_new_missing_values(rng, model);
    }
}
