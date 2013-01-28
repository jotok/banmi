#include <libguile.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include "banmi.h"
#include "banmi_util.h"

static scm_t_bits thit_model_tag;
static SCM thit_error;
static gsl_rng *rng;

static SCM
thit_new_model(SCM s_max_rows, SCM s_bds_disc, SCM s_n_cont, SCM s_dp_weight,
               SCM s_lambda_a, SCM s_lambda_b) 
{
    int max_rows = scm_to_int(s_max_rows);
    int n_cont = scm_to_int(s_n_cont);
    double dp_weight = scm_to_double(s_dp_weight);
    double lambda_a = scm_to_double(s_lambda_a);
    double lambda_b = scm_to_double(s_lambda_b);

    int n_disc = scm_to_int(scm_length(s_bds_disc));
    gsl_vector_int *bds_disc =  gsl_vector_int_alloc(n_disc);

    int i, b;
    for (i = 0; i < n_disc; i++) {
        b = scm_to_int(scm_list_ref(s_bds_disc, scm_from_int(i)));
        gsl_vector_int_set(bds_disc, i, b);
    }

    banmi_model_t *model = new_banmi_model(max_rows, bds_disc, n_cont, dp_weight,
                                           lambda_a, lambda_b);

    SCM smob;
    SCM_NEWSMOB(smob, thit_model_tag, model);

    return smob;
}

static size_t
thit_free_model(SCM s_model) {
    banmi_model_t *model = (banmi_model_t*) SCM_SMOB_DATA(s_model);

    gsl_matrix_int_free(model->disc);
    gsl_matrix_int_free(model->disc_imp);
    gsl_matrix_free(model->cont);
    gsl_matrix_free(model->cont_imp);
    gsl_matrix_int_free(model->x);
    gsl_matrix_free(model->mu);
    gsl_vector_free(model->lambda);
    gsl_vector_free(model->sigma);
    free_tab(model->crosstab);
    gsl_vector_free(model->mu_a);
    gsl_vector_free(model->mu_b);
    gsl_vector_free(model->sigma_a);
    gsl_vector_free(model->sigma_b);
    gsl_vector_int_free(model->bds_disc);
    free(model->mis_pat);

    free(model);
    return 0;
}

SCM
thit_scm_from_vector(gsl_vector *v) {
    SCM s_vector = scm_c_make_vector(v->size, SCM_BOOL_F);

    int i;
    for (i = 0; i < v->size; i++)
        scm_vector_set_x(s_vector, scm_from_int(i),
                         scm_from_double(gsl_vector_get(v, i)));

    return s_vector;
}

SCM
thit_scm_from_vector_int(gsl_vector_int *v) {
    SCM s_vector = scm_c_make_vector(v->size, SCM_BOOL_F);

    int i;
    for (i = 0; i < v->size; i++)
        scm_vector_set_x(s_vector, scm_from_int(i),
                         scm_from_int(gsl_vector_int_get(v, i)));

    return s_vector;
}

SCM
thit_scm_from_banmi_data(gsl_matrix_int *disc, gsl_matrix *cont, int n_rows) {
    SCM s_rows = scm_c_make_vector(n_rows, SCM_BOOL_F);

    int i, j;
    for (i = 0; i < n_rows; i++) {
        SCM this_row = scm_c_make_vector(disc->size2 + cont->size2, SCM_BOOL_F);

        for (j = 0; j < disc->size2; j++) {
            scm_vector_set_x(this_row, scm_from_int(j),
                             scm_from_int(gsl_matrix_int_get(disc, i, j)));
        }

        for (j = 0; j < cont->size2; j++) {
            scm_vector_set_x(this_row, scm_from_int(j + disc->size2),
                             scm_from_double(gsl_matrix_get(cont, i, j)));
        }

        scm_vector_set_x(s_rows, scm_from_int(i), this_row);
    }

    return s_rows;
}


SCM
thit_get_lambda(SCM s_model) {
    scm_assert_smob_type(thit_model_tag, s_model);
    banmi_model_t *model = (banmi_model_t*)SCM_SMOB_DATA(s_model);
    return thit_scm_from_vector(model->lambda);
}

SCM
thit_get_sigma(SCM s_model) {
    scm_assert_smob_type(thit_model_tag, s_model);
    banmi_model_t *model = (banmi_model_t*)SCM_SMOB_DATA(s_model);
    return thit_scm_from_vector(model->sigma);
}

SCM
thit_get_data(SCM s_model) {
    scm_assert_smob_type(thit_model_tag, s_model);
    banmi_model_t *model = (banmi_model_t*)SCM_SMOB_DATA(s_model);
    return thit_scm_from_banmi_data(model->disc, model->cont, model->n_rows);
}

SCM
thit_get_imputed_data(SCM s_model) {
    scm_assert_smob_type(thit_model_tag, s_model);
    banmi_model_t *model = (banmi_model_t*)SCM_SMOB_DATA(s_model);
    return thit_scm_from_banmi_data(model->disc_imp, model->cont_imp, model->n_rows);
}

// Load one of data into the model. The vararg should be a list of values with
// length equal to the number of discrete columns plus the number of continuous
// columns. The first values in the vararg are taken to be discrete, followed
// by the continuous values.
//
SCM
thit_load_row_x(SCM s_model, SCM s_varargs) {
    scm_assert_smob_type(thit_model_tag, s_model);
    banmi_model_t *model = (banmi_model_t*)SCM_SMOB_DATA(s_model);

    int n_col = model->n_disc + model->n_cont;
    int n_args = scm_to_int(scm_length(s_varargs));

    if (model->n_rows == model->disc->size1)
        scm_error_scm(thit_error, 
                      scm_from_locale_string("load-row!"),
                      scm_from_locale_string("The model is full, can't add more rows."),
                      scm_list_2(scm_from_int(n_col), scm_from_int(n_args)),
                      SCM_BOOL_F);

    if (n_args != n_col)
        scm_error_scm(thit_error, 
                      scm_from_locale_string("load-row!"),
                      scm_from_locale_string("Expected ~A values, got ~A."),
                      scm_list_2(scm_from_int(n_col), scm_from_int(n_args)),
                      SCM_BOOL_F);

    int j, ival;
    for (j = 0; j < model->n_disc; j++) {
        ival = scm_to_int(scm_list_ref(s_varargs, scm_from_int(j)));
        gsl_matrix_int_set(model->disc, model->n_rows, j, ival);
    }

    double dval;
    for (j = 0; j < model->n_cont; j++) {
        dval = scm_to_double(scm_list_ref(s_varargs, scm_from_int(j + model->n_disc)));
        gsl_matrix_set(model->cont, model->n_rows, j, dval);
    }

    model->n_rows++;

    return SCM_BOOL_T;
}

SCM
thit_data_augmentation_x(SCM s_model, SCM s_n_iter) {
    scm_assert_smob_type(thit_model_tag, s_model);
    banmi_model_t *model = (banmi_model_t*)SCM_SMOB_DATA(s_model);
    int n_iter = scm_to_int(s_n_iter);
    banmi_data_augmentation(rng, model, n_iter);
    return SCM_BOOL_T;
}

SCM
thit_count_unique_modes(SCM s_model) {
    scm_assert_smob_type(thit_model_tag, s_model);
    banmi_model_t *model = (banmi_model_t*)SCM_SMOB_DATA(s_model);
    return scm_from_int(banmi_count_unique_modes(model));
}

void
banmi_thit(void) {
    thit_model_tag = scm_make_smob_type("banmi_model", sizeof(banmi_model_t*));
    scm_set_smob_free(thit_model_tag, thit_free_model);

    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, time(NULL));

    thit_error = scm_from_locale_symbol("thit-error");

    scm_c_define_gsubr("new-banmi-model", 6, 0, 0, thit_new_model);
    scm_c_define_gsubr("banmi-get-lambda", 1, 0, 0, thit_get_lambda);
    scm_c_define_gsubr("banmi-get-sigma", 1, 0, 0, thit_get_sigma);
    scm_c_define_gsubr("banmi-get-data", 1, 0, 0, thit_get_data);
    scm_c_define_gsubr("banmi-get-imputed-data", 1, 0, 0, thit_get_imputed_data);
    scm_c_define_gsubr("banmi-load-row!", 1, 0, 1, thit_load_row_x);
    scm_c_define_gsubr("banmi-data-augmentation!", 2, 0, 0, thit_data_augmentation_x);
    scm_c_define_gsubr("banmi-count-unique-modes", 1, 0, 0, thit_count_unique_modes);
}
