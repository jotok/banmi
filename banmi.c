#include "banmi.h"

banmi_model_t*
alloc_banmi_model() {
    banmi_model_t *model = malloc(sizeof(banmi_model_t));

    int *data_dim = (int[]) {MaxRows, NDiscrete};

    model->disc = alloc_tab(2, data_dim);
    model->disc_imp = alloc_tab(2, data_dim);
    model->cont = malloc(sizeof(double) * MaxRows);
    model->cont_imp = malloc(sizeof(double) * MaxRows);

    model->x = alloc_tab(2, data_dim);
    model->mu = malloc(sizeof(double) * MaxRows);

    model->crosstab = alloc_tab(NDiscrete, DiscDim);

    model->mis_pat = malloc(sizeof(int) * MaxRows);

    return model;
}

// Open flas1.txt and read it into the data fields in the model. flas1.txt is
// a space-delimited file with 280 lines, including a header. Each line 
// contains 14 values of which this program will use 6. Additionally, set the
// n_complete field in the model to the number of rows read.
//
// From the FLAS description file:
//   'LAN2' 1=Spanish, 0=other.
//   'LAN3' 1=German, 0=other.
//   'LAN4' 1=Russian, 0=other.
//   'AGE' age group (1=less than 20, 2=20+).
//   'PRI' number of prior foreign language courses (1=none, 2=1-2, 3=3+).
//   'SATM' Scholastic Aptitude Test, math score
//
// Subtract 1 from the AGE and PRI so that their values can be used as array 
// indices
//
int
read_data_from_file(banmi_model_t *model) {

    FILE *f = fopen("data/flas1.txt", "r");
    while (fgetc(f) != '\n') ; // advance past the header line

    int i = 0, disc_ix[2], lan2, lan3, lan4, age, pri, sex, mlat, flas, satv, 
        satm, eng, grd;
    double hgpa, cgpa;
    while (i < MaxRows) {
        fscanf(f, "%d %d %d %d %d %d %d %d %d %d %d %lf %lf %d",
               &lan2, &lan3, &lan4, &age, &pri, &sex, &mlat, &flas, &satv, 
               &satm, &eng, &hgpa, &cgpa, &grd);

        if (feof(f)) break;

        disc_ix[0] = i;
        disc_ix[1] = 0; 
        tab_set(model->disc, disc_ix, lan2);
        tab_set(model->disc_imp, disc_ix, lan2);
        disc_ix[1] = 1; 
        tab_set(model->disc, disc_ix, lan3);
        tab_set(model->disc_imp, disc_ix, lan3);
        disc_ix[1] = 2; 
        tab_set(model->disc, disc_ix, lan4);
        tab_set(model->disc_imp, disc_ix, lan4);
        disc_ix[1] = 3; 
        tab_set(model->disc, disc_ix, (age > 0)? age - 1 : age);
        tab_set(model->disc_imp, disc_ix, (age > 0)? age - 1 : age);
        disc_ix[1] = 4; 
        tab_set(model->disc, disc_ix, (pri > 0)? pri - 1 : pri);
        tab_set(model->disc_imp, disc_ix, (pri > 0)? pri - 1 : pri);

        model->cont[i] = (double)satm;
        model->cont_imp[i] = (double)satm;

        i++;
    } 

    fclose(f);

    model->n_rows = i;
    return i;
}

// compute the missingness pattern for the given row of data
int 
missingness_pattern(const banmi_model_t *model, int row) {
    int pattern = 0, i, ix[2];
    ix[0] = row;

    for (i = 0; i < NDiscrete; i++) {
        ix[1] = i;
        if (tab_get(model->disc, ix) < 0)
            pattern |= 1 << i;
    }

    if (model->cont[row] < 0)
        pattern |= 1 << NDiscrete;

    return pattern;
}

// convert a missingness pattern to an array of margin indices.
// return the number of margins with missing data.
// Example:
//   [0 -1  0 -1  1]  => [1  3]
int 
missingness_pattern_to_index(int *ix, int pattern) {
    int i, j = 0, mask = 1;
    for (i = 0; i < NDiscrete; i++) {
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
    order_blocks(model->disc->dat, model->mis_pat, model->n_rows, NDiscrete);
    order_blocks(model->disc_imp->dat, model->mis_pat, model->n_rows, NDiscrete);
    order(model->mis_pat, model->mis_pat, model->n_rows);

    model->n_complete = n_complete;
    model->n_missing = model->n_rows - n_complete;

    return n_complete;
}

// Compute crosstab, mu_a, and mu_b from the data. Set lambda_a, lambda_b,
// and dp_weight according to constants in banmi.h.
//
// This function should be called after sort_data_by_missingness_pattern, since
// that function sets model->n_complete.
//
void
init_hyperparameters(banmi_model_t *model) {
    int i, j, disc_ix[2], crosstab_ix[NDiscrete];

    // compute the crosstab of the complete rows

    for (i = 0; i < model->crosstab->size; i++)
        model->crosstab->dat[i] = 1;

    for (i = 0; i < model->n_complete; i++) {
        disc_ix[0] = i;

        for (j = 0; j < NDiscrete; j++) {
            disc_ix[1] = j;
            crosstab_ix[j] = tab_get(model->disc, disc_ix);
        }

        tab_set(model->crosstab, crosstab_ix, 
                tab_get(model->crosstab, crosstab_ix) + 1);
    }

    // take the distribution of means for SATM to be (scaled) beta
    // recall, SATM takes values between 0 and 800. Note that I'm using all
    // given SATM values, not just the ones in complete data rows.

    double satm_mean = 0.0;
    j = 0;
    for (i = 0; i < model->n_rows; i++) {
        if (model->cont[i] >= 0) {
            satm_mean += model->cont[i]/800;
            j++;
        }
    }
    satm_mean /= j;

    double satm_var = 0;
    for (i = 0; i < model->n_rows; i++) {
        if (model->cont[i] >= 0)
            satm_var += pow(model->cont[i]/800 - satm_mean, 2.0);
    }
    satm_var /= j - 1;

    model->mu_a = satm_mean * (satm_mean * (1 - satm_mean) / satm_var - 1);
    model->mu_b = (1 - satm_mean) * (satm_mean * (1 - satm_mean) / satm_var - 1);

    model->dp_weight = DPWeight;
    model->sigma_a = SigmaA;        
    model->sigma_b = SigmaB;  
    model->lambda_a = LambdaA;        
    model->lambda_b = LambdaB;        
}

// Set initial values for missing values in the imputed data set. This function
// should be called after init_hyperparams.
//
void
init_missing_values(gsl_rng *rng, banmi_model_t *model) {
    int i, j, last_pattern, this_pattern, n_missing, flat_ix, disc_ix[2], 
        missing_ix[NDiscrete], marg_ix[NDiscrete];
    tab_t *margtab = NULL;

    // fill in missing discrete values

    last_pattern = 0;
    for (i = model->n_complete; i < model->n_rows; i++) {

        this_pattern = model->mis_pat[i];

        // skip this step if missing continuous data only
        if (this_pattern == MissingContinuousData)
            continue;

        if (this_pattern != last_pattern) {
            // rebuild the marginal count table
            last_pattern = this_pattern;
            n_missing = missingness_pattern_to_index(missing_ix, this_pattern);

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
            model->cont_imp[i] = gsl_ran_beta(rng, model->mu_a, model->mu_b) * 800;
    }
}

double 
kernel(double u, const int *w, double mu, const int *x, double sigma, double lambda) 
{
    int i, n_agr = 0;
    for (i = 0; i < NDiscrete; i++) {
        if (w[i] == x[i])
            n_agr++;
    }
    int n_dis = NDiscrete - n_agr;

    double result = pow(1 - lambda, n_agr) * pow(lambda, n_dis) * gsl_ran_gaussian_pdf(u - mu, sigma);
    for (i = 0; i < NDiscrete; i++) {
        if (w[i] != x[i])
            result *= 1.0 / (DiscDim[i] - 1);
    }

    return result;
}

// Initialize the means/modes of the mixture model. This function should be 
// called after init_hyperparams.
//
void
init_latent_variables(gsl_rng *rng, banmi_model_t *model) {
    int i, j, flat_ix, x_ix[2], crosstab_ix[NDiscrete];

    for (i = 0; i < model->n_rows; i++) {
        model->mu[i] = gsl_ran_beta(rng, model->mu_a, model->mu_b) * 800.0;

        x_ix[0] = i;
        flat_ix = sample(rng, model->crosstab->size, model->crosstab->dat);
        tab_array_index(crosstab_ix, model->crosstab, flat_ix);
        for (j = 0; j < NDiscrete; j++) {
            x_ix[1] = j;
            tab_set(model->x, x_ix, crosstab_ix[j]);
        }
    }

    model->sigma = sqrt(1.0 / gsl_ran_gamma(rng, model->sigma_a, model->sigma_b));
    model->lambda = gsl_ran_beta(rng, model->lambda_a, model->lambda_b);
}

void
draw_new_latent_variables(gsl_rng *rng, banmi_model_t *model) {

    int i, j, k, flat_ix, choice, x_ix[NDiscrete], y_ix[NDiscrete], 
        crosstab_ix[NDiscrete];
    double weight[model->n_rows];

    // draw new means/modes

    for (i = 0; i < model->n_rows; i++) {
        weight[0] = model->dp_weight * kernel(model->cont_imp[i],
                                              model->disc_imp->dat + i*NDiscrete,
                                              model->mu[i],
                                              model->x->dat + i*NDiscrete,
                                              model->sigma,
                                              model->lambda);
        k = 1; // insertion index
        for (j = 0; j < model->n_rows; j++) {
            if (j == i)
                continue;

            weight[k] = kernel(model->cont_imp[i],
                               model->disc_imp->dat + i*NDiscrete,
                               model->mu[j],
                               model->x->dat + j*NDiscrete,
                               model->sigma,
                               model->lambda);

            k++;
        }

        choice = sample_d(rng, model->n_rows, weight);

        if (choice == 0) {
            // sample from G_0
            model->mu[i] = gsl_ran_beta(rng, model->mu_a, model->mu_b) * 800.0;

            x_ix[0] = i;
            flat_ix = sample(rng, model->crosstab->size, model->crosstab->dat);
            tab_array_index(crosstab_ix, model->crosstab, flat_ix);
            for (j = 0; j < NDiscrete; j++) {
                x_ix[1] = j;
                tab_set(model->x, x_ix, crosstab_ix[j]);
            }
        } else {
            // set (x[i], mu[i]) to (x[choice], mu[choice])
            if (choice <= i)
                choice--;

            x_ix[0] = i;
            y_ix[0] = choice;
            for (j = 0; j < NDiscrete; j++) {
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
    int disc_ix[NDiscrete];

    for (i = 0; i < model->n_rows; i++) {
        disc_ix[0] = i;
        for (j = 0; j < NDiscrete; j++) {
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

        if (model->mis_pat[i] != MissingContinuousData) {
            // there are missing discrete values
            for (j = 0; j < NDiscrete; j++) {
                disc_ix[1] = j;
                if (tab_get(model->disc, disc_ix) < 0) {
                    u = gsl_rng_uniform(rng);
                    if (u > 1 - model->lambda) {
                        tab_set(model->disc_imp, disc_ix, tab_get(model->x, disc_ix));
                    } else {
                        choice = (int) (gsl_rng_uniform(rng) * (DiscDim[j] - 1));

                        if (choice >= tab_get(model->x, disc_ix))
                            choice++;

                        tab_set(model->disc_imp, disc_ix, choice);
                    }
                }
            }
        }

        if (model->mis_pat[i] >= MissingContinuousData) {
            // the continuous value is missing
            model->cont_imp[i] = model->mu[i] + gsl_ran_gaussian(rng, model->sigma);
        }
    }
}

int 
main() {

    // initialize the RNG
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);

    // initialize the model
    banmi_model_t *model = alloc_banmi_model();
    read_data_from_file(model);
    sort_data_by_missingness_pattern(model);
    init_hyperparameters(model);
    init_missing_values(rng, model);
    init_latent_variables(rng, model);

    // perform the data augmentation

    int t;
    for (t = 0; t < NIter; t++) {
        draw_new_latent_variables(rng, model);
        draw_new_missing_values(rng, model);
    }

    // show the result
    int i, disc_ix[2], lan2, lan3, lan4, age, pri;
    for (i = 0; i < model->n_rows; i++) {
        disc_ix[0] = i;
        disc_ix[1] = 0; lan2 = tab_get(model->disc_imp, disc_ix);
        disc_ix[1] = 1; lan3 = tab_get(model->disc_imp, disc_ix);
        disc_ix[1] = 2; lan4 = tab_get(model->disc_imp, disc_ix);
        disc_ix[1] = 3; age = tab_get(model->disc_imp, disc_ix);
        disc_ix[1] = 4; pri = tab_get(model->disc_imp, disc_ix);

        printf("%2d %2d %2d %2d %2d %6.2f\n", 
               lan2, lan3, lan4, age, pri, model->cont_imp[i]);
    }

    return 0;
}
