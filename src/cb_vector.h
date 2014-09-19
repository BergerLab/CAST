#ifndef __CABLAST_VECTOR_H__
#define __CABLAST_VECTOR_H__

//Apparently this is required to make pthread_rwlock* stuff available.
#define __USE_UNIX98

#include <pthread.h>

struct cb_vector_array {
    int size;
    int capacity;
    int th;
    void **data;
};

struct cb_vector {
    int capacity;
    struct cb_vector_array *data;

    pthread_rwlock_t lock_add;
};

struct cb_vector *cb_vector_init(int capacity);
void cb_vector_free(struct cb_vector *vector);
void cb_vector_free_no_data(struct cb_vector *vector);
void *cb_vector_get(struct cb_vector *vector, int index);
void cb_vector_append(struct cb_vector *vector, void *data);

#endif
