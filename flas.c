#include <stdio.h>
#include <gsl/gsl_rng.h>

#include "banmi.h"

#define MaxRows 300
#define NDiscrete 5
#define NOrdered 2
#define NContinuous 2
#define NIter 50

#define DPWeight 50        // weight paramter to the DP
#define LambdaA 2.0        // alpha parameter to lambda prior
#define LambdaB 3.0        // beta parameter to lambda prior

static double BdsDisc[5] = {2, 2, 2, 2, 3};
static int BdsOrde[2] = {800, 800};

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

        gsl_matrix_int_set(model->orde, i, 0, satv);
        gsl_matrix_int_set(model->orde_imp, i, 0, satv);
        gsl_matrix_int_set(model->orde, i, 1, satm);
        gsl_matrix_int_set(model->orde_imp, i, 1, satm);

        gsl_matrix_set(model->cont, i, 0, (hgpa+0.0)/4.0);
        gsl_matrix_set(model->cont_imp, i, 0, (hgpa+0.0)/4.0);
        gsl_matrix_set(model->cont, i, 1, (cgpa+0.0)/4.0);
        gsl_matrix_set(model->cont_imp, i, 1, (cgpa+0.0)/4.0);

        i++;
    } 

    fclose(f);

    model->n_rows = i;
    return i;
}

int
main() {
    gsl_vector_int *bds_disc = gsl_vector_int_alloc(NDiscrete);
    int i;
    for (i = 0; i < NDiscrete; i++)
        gsl_vector_int_set(bds_disc, i, BdsDisc[i]);

    gsl_vector_int *bds_orde = gsl_vector_int_alloc(NOrdered);
    for (i = 0; i < NOrdered; i++)
        gsl_vector_int_set(bds_orde, i, BdsOrde[i]);

    banmi_model_t *model = new_banmi_model(MaxRows, bds_disc, bds_orde, NContinuous,
                                           DPWeight,
                                           LambdaA, LambdaB);
    read_data_from_file(model);

    gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
    banmi_data_augmentation(rng, model, NIter);

    // show the result
    int disc_ix[2], lan2, lan3, lan4, age, pri;
    for (i = 0; i < model->n_rows; i++) {
        disc_ix[0] = i;
        disc_ix[1] = 0; lan2 = tab_get(model->disc_imp, disc_ix);
        disc_ix[1] = 1; lan3 = tab_get(model->disc_imp, disc_ix);
        disc_ix[1] = 2; lan4 = tab_get(model->disc_imp, disc_ix);
        disc_ix[1] = 3; age = tab_get(model->disc_imp, disc_ix) + 1;
        disc_ix[1] = 4; pri = tab_get(model->disc_imp, disc_ix) + 1;

        printf("%2d %2d %2d %2d %2d %4d %4d %4.2f %4.2f\n", 
               lan2, lan3, lan4, age, pri, 
               gsl_matrix_int_get(model->orde_imp, i, 0),
               gsl_matrix_int_get(model->orde_imp, i, 1),
               gsl_matrix_get(model->cont_imp, i, 0) * 4.0,
               gsl_matrix_get(model->cont_imp, i, 1) * 4.0);
    }

    printf("\n");
    printf("lambda %4.2f %4.2f %4.2f %4.2f %4.2f\n",
           gsl_vector_get(model->lambda, 0), 
           gsl_vector_get(model->lambda, 1), 
           gsl_vector_get(model->lambda, 2), 
           gsl_vector_get(model->lambda, 3), 
           gsl_vector_get(model->lambda, 4));
    printf("kappa %4.2f %4.2f\n",
           gsl_vector_get(model->kappa, 0), 
           gsl_vector_get(model->kappa, 1));
    printf("sigma %4.2f %4.2f\n",
           gsl_vector_get(model->sigma, 0), 
           gsl_vector_get(model->sigma, 1));

    return 0;
}
