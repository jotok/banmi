#include "banmi_util.h"

tab_t* alloc_tab(int n, const int *dim) {
    tab_t *t = malloc(sizeof(tab_t));
    t->n = n;
    t->size = 1;
    t->dim = malloc(n * sizeof(int));
    t->step = malloc(n * sizeof(int));

    int this_step = 1;

    int i;
    for (i = 0; i < n; i++) {
        t->dim[i] = dim[i];
        t->size *= dim[i];
        t->step[n - 1 - i] = this_step;

        this_step *= dim[n - 1 - i];
    }

    t->dat = malloc(t->size * sizeof(int));

    return t;
}

void free_tab(tab_t *t) {
    free(t->dat);
    free(t->step);
    free(t->dim);
    free(t);
}

void tab_array_index(int *arr_ix, const tab_t *t, int flat_ix) {
    // assumes arr_ix is properly alloc'd
    int i, s;
    for (i = 0; i < t->n; i++) {
        s = t->step[i];
        arr_ix[i] = flat_ix / s;
        flat_ix = flat_ix % s;
    }
}

int tab_flat_index(int *flat_index, const tab_t *t, const int *arr_ix) {
    *flat_index = 0;

    int i;
    for (i = 0; i < t->n; i++) {
        *flat_index += arr_ix[i] * t->step[i];
    }

    return *flat_index;
}

int tab_get(const tab_t *t, const int *arr_ix) {
    int flat_ix;
    tab_flat_index(&flat_ix, t, arr_ix);
    return t->dat[flat_ix];
}

void tab_set(tab_t *t, const int *arr_ix, int val) {
    int flat_ix;
    tab_flat_index(&flat_ix, t, arr_ix);
    t->dat[flat_ix] = val;
}

tab_t* 
new_marginal_tab(const tab_t *t, const int *margins, const int n_margins) {
    int i, j, k, m_dim[n_margins], m_ix[n_margins], t_ix[t->n];

    for (i = 0; i < n_margins; i++) {
        m_dim[i] = t->dim[margins[i]];
    }

    tab_t *m = alloc_tab(n_margins, m_dim);

    for (i = 0; i < m->size; i++) 
        m->dat[i] = 0.0;

    for (i = 0; i < t->size; i++) {
        tab_array_index(t_ix, t, i);
        for (j = 0; j < n_margins; j++) {
            m_ix[j] = t_ix[margins[j]];
        }
        tab_flat_index(&k, m, m_ix);
        m->dat[k] += tab_get(t, t_ix);
    }

    return m;
}

// Compute the conditional distribution over the given margins conditioned on the given
// values.
//
// Params:
//   t               The base table for this computation
//   free_margins    The indices of the free margins
//   conditioned_on  The values assigned to the non-free margins
//   n_free          Number of free margins
//
// Note that the length of conditioned_on should be t->n - n_margins.
// 
// TODO: it might be more sensible to pass the indices of the conditioned margins as an
// argument, but the current version started as a copy/paste of new_marginal_tab.
//
tab_t*
new_conditional_tab(const tab_t *t, const int *free_margins, const int *conditioned_on, 
                    const int n_free) {

    int i, j, k, flat_ix, m_dim[n_free], m_ix[n_free], t_ix[t->n];

    int n_fixed = t->n - n_free;
    int fixed_margins[n_fixed];

    for (i = 0; i < n_free; i++) {
        m_dim[i] = t->dim[free_margins[i]];
    }

    tab_t *m = alloc_tab(n_free, m_dim);

    for (i = 0; i < m->size; i++) 
        m->dat[i] = 0.0;


    // get the indices of values, i.e. non-margin indices
    k = 0;
    bool is_free;
    for (i = 0; i < t->n; i++) {
        is_free = false;
        for (j = 0; j < m->n; j++) {
            if (i == free_margins[j]) {
                is_free = true;
                break;
            }
        }

        if (is_free)
            continue;

        fixed_margins[k++] = i;
    }

    bool skip_row;
    for (i = 0; i < t->size; i++) {
        skip_row = false;
        tab_array_index(t_ix, t, i);

        for (j = 0; j < n_fixed; j++) {
            if (t_ix[fixed_margins[j]] != conditioned_on[j]) {
                skip_row = true;
                break;
            }
        }

        if (skip_row)
            continue;

        for (j = 0; j < n_free; j++) {
            m_ix[j] = t_ix[free_margins[j]];
        }
        tab_flat_index(&flat_ix, m, m_ix);
        m->dat[flat_ix] += tab_get(t, t_ix);
    }

    return m;
}

#define Swap(a, b, temp) { temp = a; a = b; b = temp; }

int _partition(int *values, int *by, int left, int right, int pivot) {
    int itemp, pivot_value = by[pivot];
    double dtemp;

    Swap(by[pivot], by[right], itemp)
    Swap(values[pivot], values[right], dtemp)

    int i, store_index = left;
    for (i = left; i < right; i++) {
        if (by[i] < pivot_value) {
            Swap(by[i], by[store_index], itemp)
            Swap(values[i], values[store_index], dtemp)
            store_index++;
        }
    }

    Swap(by[store_index], by[right], itemp);
    Swap(values[store_index], values[right], dtemp)
    return store_index;
}

void _order(int *values, int *by, int left, int right) {
    if (left < right) {
        int pivot = (right + left) / 2;
        pivot = _partition(values, by, left, right, pivot);
        _order(values, by, left, pivot - 1);
        _order(values, by, pivot + 1, right);
    }
}

void order(int *values, const int *by, int n) {
    int i, *by_copy = malloc(n * sizeof(int));

    for (i = 0; i < n; i++) 
        by_copy[i] = by[i];

    _order(values, by_copy, 0, n-1);

    free(by_copy);
}

int _partition_d(double *values, int *by, int left, int right, int pivot) {
    int itemp, pivot_value = by[pivot];
    double dtemp;

    Swap(by[pivot], by[right], itemp)
    Swap(values[pivot], values[right], dtemp)

    int i, store_index = left;
    for (i = left; i < right; i++) {
        if (by[i] < pivot_value) {
            Swap(by[i], by[store_index], itemp)
            Swap(values[i], values[store_index], dtemp)
            store_index++;
        }
    }

    Swap(by[store_index], by[right], itemp);
    Swap(values[store_index], values[right], dtemp)
    return store_index;
}

void _order_d(double *values, int *by, int left, int right) {
    if (left < right) {
        int pivot = (right + left) / 2;
        pivot = _partition_d(values, by, left, right, pivot);
        _order_d(values, by, left, pivot - 1);
        _order_d(values, by, pivot + 1, right);
    }
}

void order_d(double *values, const int *by, int n) {
    int i, *by_copy = malloc(n * sizeof(int));

    for (i = 0; i < n; i++) 
        by_copy[i] = by[i];

    _order_d(values, by_copy, 0, n-1);

    free(by_copy);
}

int _partition_blocks(int *values, int *by, int left, int right, int pivot, int block_size) {
    int i, temp, pivot_value = by[pivot];
    
    Swap(by[pivot], by[right], temp)
    for (i = 0; i < block_size; i++) 
        Swap(values[pivot*block_size + i], values[right*block_size + i], temp)

    int j, store_index = left;
    for (i = left; i < right; i++) {
        if (by[i] < pivot_value) {
            Swap(by[i], by[store_index], temp)
            for (j = 0; j < block_size; j++)
                Swap(values[i*block_size + j], values[store_index*block_size + j], temp)
            store_index++;
        }
    }

    Swap(by[store_index], by[right], temp)
    for (i = 0; i < block_size; i++) 
        Swap(values[store_index*block_size + i], values[right*block_size + i], temp)

    return store_index;
}

void _order_blocks(int *values, int *by, int left, int right, int block_size) {
    if (left < right) {
        int pivot = (right + left) / 2;
        pivot = _partition_blocks(values, by, left, right, pivot, block_size);
        _order_blocks(values, by, left, pivot - 1, block_size);
        _order_blocks(values, by, pivot + 1, right, block_size);
    }
}

void order_blocks(int *values, const int *by, int n_blocks, int block_size) {
    int i, *by_copy = malloc(n_blocks * sizeof(int));

    for (i = 0; i < n_blocks; i++) 
        by_copy[i] = by[i];

    _order_blocks(values, by_copy, 0, n_blocks-1, block_size);

    free(by_copy);
}

int _partition_blocks_d(double *values, int *by, int left, int right, int pivot, int block_size) {
    int i, pivot_value = by[pivot];
    double temp;
    
    Swap(by[pivot], by[right], temp)
    for (i = 0; i < block_size; i++) 
        Swap(values[pivot*block_size + i], values[right*block_size + i], temp)

    int j, store_index = left;
    for (i = left; i < right; i++) {
        if (by[i] < pivot_value) {
            Swap(by[i], by[store_index], temp)
            for (j = 0; j < block_size; j++)
                Swap(values[i*block_size + j], values[store_index*block_size + j], temp)
            store_index++;
        }
    }

    Swap(by[store_index], by[right], temp)
    for (i = 0; i < block_size; i++) 
        Swap(values[store_index*block_size + i], values[right*block_size + i], temp)

    return store_index;
}

void _order_blocks_d(double *values, int *by, int left, int right, int block_size) {
    if (left < right) {
        int pivot = (right + left) / 2;
        pivot = _partition_blocks_d(values, by, left, right, pivot, block_size);
        _order_blocks_d(values, by, left, pivot - 1, block_size);
        _order_blocks_d(values, by, pivot + 1, right, block_size);
    }
}

void order_blocks_d(double *values, const int *by, int n_blocks, int block_size) {
    int i, *by_copy = malloc(n_blocks * sizeof(int));

    for (i = 0; i < n_blocks; i++) 
        by_copy[i] = by[i];

    _order_blocks_d(values, by_copy, 0, n_blocks-1, block_size);

    free(by_copy);
}

int sample(gsl_rng *rng, int n, const int *weight) {
    double cumprob = 0.0;
    int i, total = 0;

    for (i = 0; i < n; i++)
        total += weight[i];

    double u = gsl_rng_uniform(rng);

    for (i = 0; i < n; i++) {
        cumprob += (weight[i] + 0.0) / total;
        if (cumprob > u)
            break;
    }

    return i;
}

// is there a better solution than copy/paste?
int sample_d(gsl_rng *rng, int n, const double *weight) {
    double cumprob = 0.0, total = 0.0;
    int i;

    for (i = 0; i < n; i++)
        total += weight[i];

    double u = gsl_rng_uniform(rng);

    for (i = 0; i < n; i++) {
        cumprob += weight[i] / total;
        if (cumprob > u)
            break;
    }

    return i;
}

int
binary_search_d(double u, const double *cumweights, int left, int right) {
    if (left >= right - 1 || u < cumweights[left + 1])
        return left;

    int mid = (left + right) / 2;
    if (u < cumweights[mid])
        return binary_search_d(u, cumweights, left, mid);
    else
        return binary_search_d(u, cumweights, mid, right);
}

int
cumsample_d(gsl_rng *rng, int n, const double *cumweights) {
    double u = gsl_rng_uniform(rng) * cumweights[n-1];
    return binary_search_d(u, cumweights, 0, n);
}
