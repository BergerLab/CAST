#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "cb_vector.h"

/*Takes in a capacity and initializesa cb_vector_array struct with that
  capacity*/
struct cb_vector_array *cb_vector_array_init(int capacity){
    struct cb_vector_array *array = malloc(sizeof(*array));
    assert(array);

    array->data = malloc(capacity*sizeof(*(array->data)));
    assert(array->data);

    array->capacity = capacity;
    array->size     = 0;
    array->th       = 0;

    return array;
}

//Free a cb_vector_array struct and its data
void cb_vector_array_free(struct cb_vector_array *array){
    //Wait until no thread is accessing the array
    while (__sync_fetch_and_or(&(array->th), 0) > 0) {}

    for (int i = 0; i < array->size; i++)
        free(array->data[i]);
    cb_vector_array_free_no_data(array);
}

//Free a cb_vector_array struct without freeing its data
void cb_vector_array_free_no_data(struct cb_vector_array *array){
    free(array->data);
    free(array);
}


/*Takes in a cb_vector_array struct and an index and returns the element at
  that index or NULL if no element exists at that index.*/
void *cb_vector_array_get(struct cb_vector_array *array, int index){
    void *el = NULL;

    /*When getting the data, increment the number of threads accessing the
      array*/
    __sync_fetch_and_add(&(array->th), 1);
    el = (index >= array->size || index < 0) ? NULL : array->data[index];
    __sync_fetch_and_sub(&(array->th), 1);

    return el;
}

//Append an element to the array
void cb_vector_array_append(struct cb_vector_array *array, void *data){
    array->data[array->size] = data;
    (array->size)++;
}


//Initialize a cb_vector struct
struct cb_vector *cb_vector_init(int capacity){
    int errno;

    struct cb_vector *vector = malloc(sizeof(*vector));
    assert(vector);

    vector->size     = 0;
    vector->capacity = capacity;
    vector->data     = cb_vector_array_init(capacity);

    if (0 != (errno = pthread_rwlock_init(&vector->lock_add, NULL))) {
        fprintf(stderr, "Could not create rwlock. Errno: %d\n", errno);
        exit(1);
    }

    return vector;
}

//Free a cb_vector struct and its data
void cb_vector_free(struct cb_vector *vector){
    int errno;

    //Destroy the lock
    if (0 != (errno = pthread_rwlock_destroy(&vector->lock_add))) {
        fprintf(stderr, "Could not destroy rwlock. Errno: %d\n", errno);
        exit(1);
    }

    cb_vector_array_free(vector->data);
}

//Get a cb_vector_struct's array of data, used for manually freeing the data
void **cb_vector_get_data(struct cb_vector *vector){
    return vector->data->data;
}

//Free a cb_vector struct without freeing its data
void cb_vector_free_no_data(struct cb_vector *vector){
    int errno;

    //Destroy the lock
    if (0 != (errno = pthread_rwlock_destroy(&vector->lock_add))) {
        fprintf(stderr, "Could not destroy rwlock. Errno: %d\n", errno);
        exit(1);
    }

    cb_vector_array_free_no_data(vector->data);
}

/*Takes in a vector and an index and gets the element in the vector at that
  index, returning NULL if there is no element at that index.*/
void *cb_vector_get(struct cb_vector *vector, int index){
    return cb_vector_array_get(vector->data, index);
}

/*Takes in a cb_vector struct and creates a new cb_vector_array struct with
 *twice the capacity of the vector's old array.  The vector's data pointer
 *is then atomically set to the new array and the old array is freed once no
 *threads are accessing it.
 */
void cb_vector_expand(struct cb_vector *vector){
    struct cb_vector_array *old_array=vector->data, 
                           *new_array=cb_vector_array_init(vector->capacity*2);
    void **data = vector->data->data, **new_data = new_array->data;
    int size = vector->size;

    new_array->size = size;

    for (int i = 0; i < size; i++)
        new_data[i] = data[i];

    __sync_val_compare_and_swap(&vector->data, vector->data, new_array);
    cb_vector_array_free(old_array);
}

/*Appends an element to the vector in a thread-safe manner, expanding the
  vector if it has reached its capacity.*/
void cb_vector_append(struct cb_vector *vector, void *data){
    pthread_rwlock_wrlock(&vector->lock_add);
    if (vector->size > vector->capacity)
        cb_vector_expand(vector);
    cb_vector_array_append(vector->data, data);
    pthread_rwlock_unlock(&vector->lock_add);
}
