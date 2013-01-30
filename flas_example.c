#include <stdio.h>
#include <time.h>

#include "banmi.h"

#define MaxRows 300
#define NDiscrete 7
#define NContinuous 7
#define NIter 250
#define NImpute 10

#define DPWeight 50        // weight paramter to the DP
#define InitCrosstab 0.1   // minimum weight in crosstab cell
#define LambdaA 2.0        // alpha parameter to lambda prior
#define LambdaB 8.0        // beta parameter to lambda prior

static double BdsDisc[7] = {2, 2, 2, 2, 3, 2, 2};

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

    gsl_vector_int *disc = gsl_vector_int_alloc(model->n_disc);
    gsl_vector *cont = gsl_vector_alloc(model->n_cont);

    while (i < MaxRows) {
        fscanf(f, "%d %d %d %d %d %d %d %d %d %d %d %lf %lf %d",
               &lan2, &lan3, &lan4, &age, &pri, &sex, &mlat, &flas, &satv, 
               &satm, &eng, &hgpa, &cgpa, &grd);

        if (feof(f)) break;

        gsl_vector_int_set(disc, 0, lan2);
        gsl_vector_int_set(disc, 1, lan3);
        gsl_vector_int_set(disc, 2, lan4);
        gsl_vector_int_set(disc, 3, (age >= 0)? age - 1 : -1);
        gsl_vector_int_set(disc, 4, (pri >= 0)? pri - 1 : -1);
        gsl_vector_int_set(disc, 5, sex);
        gsl_vector_int_set(disc, 6, (grd >= 0)? grd - 1 : -1);

        gsl_vector_set(cont, 0, (mlat >= 0) ? banmi_from_ordered_value(mlat, 110) : -1);
        gsl_vector_set(cont, 1, (flas >= 0) ? banmi_from_ordered_value(flas, 40) : -1);
        gsl_vector_set(cont, 2, (satv >= 0) ? banmi_from_ordered_value(satv, 800) : -1);
        gsl_vector_set(cont, 3, (satm > 0) ? banmi_from_ordered_value(satm, 800) : -1);
        gsl_vector_set(cont, 4, (eng > 0) ? banmi_from_ordered_value(eng, 120) : -1);
        gsl_vector_set(cont, 5, (hgpa >= 0) ? (hgpa+0.0)/4.0 : -1);
        gsl_vector_set(cont, 6, (cgpa >= 0) ? (cgpa+0.0)/4.0 : -1);

        banmi_load_row(model, disc, cont);

        i++;
    } 

    gsl_vector_int_free(disc);
    gsl_vector_free(cont);
    fclose(f);

    return i;
}

int
main() {
    gsl_vector_int *bds_disc = gsl_vector_int_alloc(NDiscrete);
    int i_impute, i;
    for (i = 0; i < NDiscrete; i++)
        gsl_vector_int_set(bds_disc, i, BdsDisc[i]);

    banmi_model_t *model = new_banmi_model(MaxRows, bds_disc, NContinuous,
                                           DPWeight, InitCrosstab, LambdaA, LambdaB);
    read_data_from_file(model);

    gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, time(NULL));
    printf("imputation lan2 lan3 lan4 age pri sex grd mlat flas satv satm eng hgpa cgpa\n");

    for (i_impute = 0; i_impute < NImpute; i_impute++) {

        banmi_data_augmentation(rng, model, NIter);


        // show the result
        for (i = 0; i < model->n_rows; i++) {

            printf("%d %d %d %d %d %d %d %d %d %d %d %d %d %.2f %.2f\n", 
                   i_impute + 1,
                   gsl_matrix_int_get(model->disc_imp, i, 0),
                   gsl_matrix_int_get(model->disc_imp, i, 1),
                   gsl_matrix_int_get(model->disc_imp, i, 2),
                   gsl_matrix_int_get(model->disc_imp, i, 3) + 1,
                   gsl_matrix_int_get(model->disc_imp, i, 4) + 1,
                   gsl_matrix_int_get(model->disc_imp, i, 5),
                   gsl_matrix_int_get(model->disc_imp, i, 6) + 1,
                   banmi_to_ordered_value(gsl_matrix_get(model->cont_imp, i, 0), 110),
                   banmi_to_ordered_value(gsl_matrix_get(model->cont_imp, i, 1), 40),
                   banmi_to_ordered_value(gsl_matrix_get(model->cont_imp, i, 2), 800),
                   banmi_to_ordered_value(gsl_matrix_get(model->cont_imp, i, 3), 800),
                   banmi_to_ordered_value(gsl_matrix_get(model->cont_imp, i, 4), 120),
                   gsl_matrix_get(model->cont_imp, i, 5) * 4.0,
                   gsl_matrix_get(model->cont_imp, i, 6) * 4.0);
        }
    }

    return 0;
}
