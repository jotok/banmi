#include <libguile.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "banmi.h"
#include "banmi_util.h"

static scm_t_bits model_tag;

static SCM
new_model(SCM s_max_rows, SCM s_bds_disc, SCM s_n_cont, SCM s_dp_weight,
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
    SCM_NEWSMOB(smob, model_tag, model);

    return smob;
}

static size_t
free_model(SCM s_model) {
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
get_discrete_value(SCM s_model, SCM s_row, SCM s_col) {
    scm_assert_smob_type(model_tag, s_model);
    banmi_model_t *model = (banmi_model_t*)SCM_SMOB_DATA(s_model);
       
    return scm_from_int(gsl_matrix_int_get(model->disc,
                                           scm_to_int(s_row),
                                           scm_to_int(s_col)));
}

SCM
set_discrete_value_x(SCM s_model, SCM s_row, SCM s_col, SCM s_val) {
    scm_assert_smob_type(model_tag, s_model);
    banmi_model_t *model = (banmi_model_t*)SCM_SMOB_DATA(s_model);

    gsl_matrix_int_set(model->disc, scm_to_int(s_row),
                                    scm_to_int(s_col),
                                    scm_to_int(s_val));

    return SCM_BOOL_T;
}

SCM
get_continuous_value(SCM s_model, SCM s_row, SCM s_col) {
    scm_assert_smob_type(model_tag, s_model);
    banmi_model_t *model = (banmi_model_t*)SCM_SMOB_DATA(s_model);

    return scm_from_double(gsl_matrix_get(model->cont,
                                          scm_to_int(s_row),
                                          scm_to_int(s_col)));
}

SCM
set_continuous_value_x(SCM s_model, SCM s_row, SCM s_col, SCM s_val) {
    scm_assert_smob_type(model_tag, s_model);
    banmi_model_t *model = (banmi_model_t*)SCM_SMOB_DATA(s_model);

    gsl_matrix_set(model->cont, scm_to_int(s_row),
                                scm_to_int(s_col),
                                scm_to_double(s_val));

    return SCM_BOOL_T;
}

void
init_banmi(void) {
    model_tag = scm_make_smob_type("banmi_model", sizeof(banmi_model_t*));
    scm_set_smob_free(model_tag, free_model);

    scm_c_define_gsubr("new-model", 6, 0, 0, new_model);
    scm_c_define_gsubr("get-discrete-value", 3, 0, 0, get_discrete_value);
    scm_c_define_gsubr("set-discrete-value!", 4, 0, 0, set_discrete_value_x);
    scm_c_define_gsubr("get-continuous-value", 3, 0, 0, get_continuous_value);
    scm_c_define_gsubr("set-continuous-value!", 4, 0, 0, set_continuous_value_x);
}
