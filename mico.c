#include "mico.h"

// Allocate a new mico model
//
mico_model_t*
new_mico_model(int max_rows, int n_cols, int n_comp, double sigma_a, double sigma_b) {
    mico_model_t *model = malloc(sizeof(mico_model_t));

    model->y = gsl_matrix_alloc(max_rows, n_cols);
    model->yi = gsl_matrix_alloc(max_rows, n_cols);

    model->xi = gsl_matrix_alloc(n_comp, n_cols);

    model->mu = gsl_vector_alloc(n_cols);
    model->theta = gsl_vector_alloc(n_cols);

    model->n_comp = n_comp;
    model->sigma_a = sigma_a;
    model->sigma_b = sigma_b;
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

// Sort the rows of data by missingness pattern and set the n_complete and
// n_missing fields. Returns the number of complete rows.
//
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
init_parameters(mico_model_t *model) {
    // place a prior on the marginal variance so that the expected value is
    // the variance of the observed values, and the weight is the number of
    // complete values.
    int i, j, insert_index = 0;
    double this_y, mean, var, alpha, beta, temp[model->n_rows];

    for (j = 0; j < model->n_cols; j++) {
        insert_index = 0;
        for (i = 0; i < model->n_rows; i++) {
            if (!isnan(this_y = gsl_matrix_get(model->y, i, j)))
                temp[insert_index++] = this_y;
        }

        mean = gsl_stats_mean(temp, 1, insert_index);
        var = gsl_stats_variance(temp, 1, insert_index);

        alpha = mean;
        beta = insert_index+0.0;
        gsl_vector_set(model->mu_a, j, alpha);
        gsl_vector_set(model->mu_b, j, beta);

        alpha = (insert_index+0.0) / 2.0;
        beta = 1.0 / (var * (alpha - 2));
        gsl_vector_set(model->theta_a, j, alpha);
        gsl_vector_set(model->theta_b, j, beta);
    }

    // draw values
}

// Set initial values in the imputed data set.
//
static void
init_missing_values(gsl_rng *rng, mico_model_t *model) {
}

void
mico_data_augmentation(gsl_rng *rng, mico_model_t *model, int n_iter) {

    if (!model->is_initialized) {
        sort_data_by_missingness_pattern(model);
        init_parameters(model);
        init_missing_values(rng, model);
    }
}
