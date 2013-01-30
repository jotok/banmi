#include "mico.h"

// Allocate a new mico model
//
mico_model_t*
new_mico_model(int max_rows, int n_cols, int n_comp) {
    mico_model_t *model = malloc(sizeof(mico_model_t));

    model->y = gsl_matrix_alloc(max_rows, n_cols);
    model->yi = gsl_matrix_alloc(max_rows, n_cols);

    model->mu = gsl_matrix_alloc(n_comp, n_cols);
    model->theta = gsl_vector_alloc(n_cols);

    model->n_comp = n_comp;
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

    order_blocks_d(model->y->data, model->mis_pat, model->n_rows, model->y->tda);
    order_blocks_d(model->yi->data, model->mis_pat, model->n_rows, model->yi->tda);
    order(model->mis_pat, model->mis_pat, model->n_rows);

    model->n_complete = n_complete;
    model->n_missing = model->n_rows - n_complete;

    return n_complete;
}
