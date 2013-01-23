#include <stdio.h>
#include <time.h>
#include <gsl/gsl_rng.h>

#include "banmi.h"

#define MaxRows 300
#define NDiscrete 5
#define NContinuous 4
#define NIter 50

#define DPWeight 50        // weight paramter to the DP
#define LambdaA 2.0        // alpha parameter to lambda prior
#define LambdaB 8.0        // beta parameter to lambda prior

static double BdsDisc[5] = {2, 2, 2, 2, 3};

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

    int i = 0, lan2, lan3, lan4, age, pri, sex, mlat, flas, satv, 
        satm, eng, grd;
    double hgpa, cgpa;
    while (i < MaxRows) {
        fscanf(f, "%d %d %d %d %d %d %d %d %d %d %d %lf %lf %d",
               &lan2, &lan3, &lan4, &age, &pri, &sex, &mlat, &flas, &satv, 
               &satm, &eng, &hgpa, &cgpa, &grd);

        if (feof(f)) break;

        gsl_matrix_int_set(model->disc, i, 0, (lan2 >= 0) ? lan2 : NAN);
        gsl_matrix_int_set(model->disc, i, 1, (lan3 >= 0) ? lan3 : NAN);
        gsl_matrix_int_set(model->disc, i, 2, (lan4 >= 0) ? lan4 : NAN);
        gsl_matrix_int_set(model->disc, i, 3, (age >= 0)? age - 1 : NAN);
        gsl_matrix_int_set(model->disc, i, 4, (pri >= 0)? pri - 1 : NAN);

        gsl_matrix_set(model->cont, i, 0, (satv >= 0) ? banmi_from_ordered_value(satv, 800) : NAN);
        gsl_matrix_set(model->cont, i, 1, (satm > 0) ? banmi_from_ordered_value(satm, 800) : NAN);
        gsl_matrix_set(model->cont, i, 2, (hgpa >= 0) ? (hgpa+0.0)/4.0 : NAN);
        gsl_matrix_set(model->cont, i, 3, (cgpa >= 0) ? (cgpa+0.0)/4.0 : NAN);

        i++;
    } 

    fclose(f);

    model->n_rows = i;
    return i;
}

int
count_unique_modes(banmi_model_t *model) {
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

int
main() {
    gsl_vector_int *bds_disc = gsl_vector_int_alloc(NDiscrete);
    int i;
    for (i = 0; i < NDiscrete; i++)
        gsl_vector_int_set(bds_disc, i, BdsDisc[i]);

    banmi_model_t *model = new_banmi_model(MaxRows, bds_disc, NContinuous,
                                           DPWeight, LambdaA, LambdaB);
    read_data_from_file(model);

    gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, time(NULL));
    banmi_data_augmentation(rng, model, NIter);

    // show the result
    for (i = 0; i < model->n_rows; i++) {

        printf("%2d %2d %2d %2d %2d %4d %4d %4.2f %4.2f\n", 
               gsl_matrix_int_get(model->disc_imp, i, 0),
               gsl_matrix_int_get(model->disc_imp, i, 1),
               gsl_matrix_int_get(model->disc_imp, i, 2),
               gsl_matrix_int_get(model->disc_imp, i, 3) + 1,
               gsl_matrix_int_get(model->disc_imp, i, 4) + 1,
               banmi_to_ordered_value(gsl_matrix_get(model->cont_imp, i, 0), 800),
               banmi_to_ordered_value(gsl_matrix_get(model->cont_imp, i, 1), 800),
               gsl_matrix_get(model->cont_imp, i, 2) * 4.0,
               gsl_matrix_get(model->cont_imp, i, 3) * 4.0);
    }

    printf("\n");
    printf("lambda %4.2f %4.2f %4.2f %4.2f %4.2f\n",
           gsl_vector_get(model->lambda, 0), 
           gsl_vector_get(model->lambda, 1), 
           gsl_vector_get(model->lambda, 2), 
           gsl_vector_get(model->lambda, 3), 
           gsl_vector_get(model->lambda, 4));
    printf("sigma %4.2f %4.2f %4.2f %4.2f\n",
           gsl_vector_get(model->sigma, 0), 
           gsl_vector_get(model->sigma, 1), 
           gsl_vector_get(model->sigma, 2), 
           gsl_vector_get(model->sigma, 3));
    printf("unique modes: %d\n", count_unique_modes(model));

    return 0;
}
