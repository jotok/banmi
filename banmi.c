#include "banmi.h"

banmi_model_t*
new_banmi_model(int max_rows, gsl_vector *bds_disc,
                double dp_weight, double sigma_a, double sigma_b, 
                double lambda_a, double lambda_b) 
{
    banmi_model_t *model = malloc(sizeof(banmi_model_t));

    int data_dim[2];
    data_dim[0] = max_rows; data_dim[1] = bds_disc->size;

    model->disc = alloc_tab(2, data_dim);
    model->disc_imp = alloc_tab(2, data_dim);
    model->cont = malloc(sizeof(double) * max_rows);
    model->cont_imp = malloc(sizeof(double) * max_rows);

    model->x = alloc_tab(2, data_dim);
    model->mu = malloc(sizeof(double) * max_rows);

    // TODO update banmi_util to take vectors?
    int disc_dim[bds_disc->size];
    int i;
    for (i = 0; i < bds_disc->size; i++)
        disc_dim[i] = (int)gsl_vector_get(bds_disc, i);

    model->crosstab = alloc_tab(bds_disc->size, disc_dim);

    model->dp_weight = dp_weight;
    model->sigma_a = sigma_a;
    model->sigma_b = sigma_b;
    model->lambda_a = lambda_a;
    model->lambda_b = lambda_b;

    model->bds_disc = bds_disc;
    model->n_disc = bds_disc->size;
    model->mask_missing_cont_data = (int)pow(2, (double)bds_disc->size);

    model->mis_pat = malloc(sizeof(int) * max_rows);

    return model;
}

// compute the missingness pattern for the given row of data
int 
missingness_pattern(const banmi_model_t *model, int row) {
    int pattern = 0, i, ix[2];
    ix[0] = row;

    for (i = 0; i < model->n_disc; i++) {
        ix[1] = i;
        if (tab_get(model->disc, ix) < 0)
            pattern |= 1 << i;
    }

    if (model->cont[row] < 0)
        pattern |= 1 << model->n_disc;

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

    order_d(model->cont, model->mis_pat, model->n_rows);
    order_d(model->cont_imp, model->mis_pat, model->n_rows);
    order_blocks(model->disc->dat, model->mis_pat, model->n_rows, model->n_disc);
    order_blocks(model->disc_imp->dat, model->mis_pat, model->n_rows, model->n_disc);
    order(model->mis_pat, model->mis_pat, model->n_rows);

    model->n_complete = n_complete;
    model->n_missing = model->n_rows - n_complete;

    return n_complete;
}

// Compute crosstab, mu_a, and mu_b from the data. Recall that lambda_a, lambda_b,
// and dp_weight are set explicitly by the user in new_banmi_model.
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

    double cont_mean = 0.0;
    j = 0;
    for (i = 0; i < model->n_rows; i++) {
        if (model->cont[i] >= 0) {
            cont_mean += model->cont[i];
            j++;
        }
    }
    cont_mean /= j;

    double cont_var = 0;
    for (i = 0; i < model->n_rows; i++) {
        if (model->cont[i] >= 0)
            cont_var += pow(model->cont[i] - cont_mean, 2.0);
    }
    cont_var /= j - 1;

    model->mu_a = cont_mean * (cont_mean * (1 - cont_mean) / cont_var - 1);
    model->mu_b = (1 - cont_mean) * (cont_mean * (1 - cont_mean) / cont_var - 1);
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

        // skip this step if missing continuous data only
        if (this_pattern == model->mask_missing_cont_data)
            continue;

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

    // fill in missing continuous values

    for (i = model->n_complete; i < model->n_rows; i++) {
        if (model->cont[i] < 0)
            model->cont_imp[i] = gsl_ran_beta(rng, model->mu_a, model->mu_b);
    }
}

double 
kernel(gsl_vector *bds_disc, double u, const int *w, double mu, const int *x, double sigma, double lambda) 
{
    int i, n_agr = 0;
    for (i = 0; i < bds_disc->size; i++) {
        if (w[i] == x[i])
            n_agr++;
    }
    int n_dis = bds_disc->size - n_agr;

    double result = pow(1 - lambda, n_agr) * pow(lambda, n_dis) * gsl_ran_gaussian_pdf(u - mu, sigma);
    for (i = 0; i < bds_disc->size; i++) {
        if (w[i] != x[i])
            result *= 1.0 / (gsl_vector_get(bds_disc, i) - 1);
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
        model->mu[i] = gsl_ran_beta(rng, model->mu_a, model->mu_b);

        x_ix[0] = i;
        flat_ix = sample(rng, model->crosstab->size, model->crosstab->dat);
        tab_array_index(crosstab_ix, model->crosstab, flat_ix);
        for (j = 0; j < model->n_disc; j++) {
            x_ix[1] = j;
            tab_set(model->x, x_ix, crosstab_ix[j]);
        }
    }

    model->sigma = sqrt(1.0 / gsl_ran_gamma(rng, model->sigma_a, model->sigma_b));
    model->lambda = gsl_ran_beta(rng, model->lambda_a, model->lambda_b);
}

void
draw_new_latent_variables(gsl_rng *rng, banmi_model_t *model) {

    int i, j, k, flat_ix, choice, x_ix[model->n_disc], y_ix[model->n_disc], 
        crosstab_ix[model->n_disc];
    double weight[model->n_rows];

    // draw new means/modes

    for (i = 0; i < model->n_rows; i++) {
        weight[0] = model->dp_weight * kernel(model->bds_disc,
                                                  model->cont_imp[i],
                                                  model->disc_imp->dat + i*model->n_disc,
                                                  model->mu[i],
                                                  model->x->dat + i*model->n_disc,
                                                  model->sigma,
                                                  model->lambda);
        k = 1; // insertion index
        for (j = 0; j < model->n_rows; j++) {
            if (j == i)
                continue;

            weight[k] = kernel(model->bds_disc,
                               model->cont_imp[i],
                               model->disc_imp->dat + i*model->n_disc,
                               model->mu[j],
                               model->x->dat + j*model->n_disc,
                               model->sigma,
                               model->lambda);

            k++;
        }

        choice = sample_d(rng, model->n_rows, weight);

        if (choice == 0) {
            // sample from G_0
            model->mu[i] = gsl_ran_beta(rng, model->mu_a, model->mu_b);

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

            model->mu[i] = model->mu[choice];
        }
    }

    // draw new sigma and lambda

    double sigma_a_post = model->sigma_a + model->n_rows / 2.0;
    double sigma_b_post = model->sigma_b;
    for (i = 0; i < model->n_rows; i++) {
        sigma_b_post += pow(model->cont_imp[i] - model->mu[i], 2.0) / 2.0;
    }
    sigma_b_post = 1.0 / sigma_b_post;

    model->sigma = sqrt(1.0 / gsl_ran_gamma(rng, sigma_a_post, sigma_b_post));

    double lambda_a_post = model->lambda_a;
    double lambda_b_post = model->lambda_b;
    int disc_ix[model->n_disc];

    for (i = 0; i < model->n_rows; i++) {
        disc_ix[0] = i;
        for (j = 0; j < model->n_disc; j++) {
            disc_ix[1] = j;
            if (tab_get(model->disc_imp, disc_ix) != tab_get(model->x, disc_ix))
                lambda_a_post++;
            else
                lambda_b_post++;
        }
    }

    model->lambda = gsl_ran_beta(rng, lambda_a_post, lambda_b_post);
}

void
draw_new_missing_values(gsl_rng *rng, banmi_model_t *model) {
    int i, j, choice, disc_ix[2];
    double u;

    for (i = model->n_complete; i < model->n_rows; i++) {

        disc_ix[0] = i;

        if (model->mis_pat[i] != model->mask_missing_cont_data) {
            // there are missing discrete values
            for (j = 0; j < model->n_disc; j++) {
                disc_ix[1] = j;
                if (tab_get(model->disc, disc_ix) < 0) {
                    u = gsl_rng_uniform(rng);
                    if (u > 1 - model->lambda) {
                        tab_set(model->disc_imp, disc_ix, tab_get(model->x, disc_ix));
                    } else {
                        choice = (int) (gsl_rng_uniform(rng) * (gsl_vector_get(model->bds_disc, j) - 1));

                        if (choice >= tab_get(model->x, disc_ix))
                            choice++;

                        tab_set(model->disc_imp, disc_ix, choice);
                    }
                }
            }
        }

        if (model->mis_pat[i] >= model->mask_missing_cont_data) {
            // the continuous value is missing
            model->cont_imp[i] = model->mu[i] + gsl_ran_gaussian(rng, model->sigma);
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
