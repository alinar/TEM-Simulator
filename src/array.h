/***********************************************************************
 * Copyright 2008-2010, Hans Rullgard, Stockholm University and 
 * Lars-Goran Ofverstedt, Karolinska Institute
 *
 * This file is part of TEM Simulator.
 *
 * TEM Simulator is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * TEM Simulator is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with TEM Simulator.  If not, see <http://www.gnu.org/licenses/>.
 *
 ***********************************************************************/

#ifndef ARRAY_HEADER
#define ARRAY_HEADER

/***********************************************************************
 * The array struct represents a 3-dimensional array of doubles
 * which is stored in dynamically allocated memory. 
 * When an array has been declared, one of the functions init_array,
 * init_array_shared or init_array_prealloc must be called before the
 * array is used. 
 * When the array has been used, free_array should be called to release 
 * allocated memory.
 ***********************************************************************/

typedef double   array_data_type;
typedef long int array_index_type;

typedef struct {
  array_index_type size[3];     /* Size of the array */
  array_index_type nalloc;      /* Size of allocated memory */
  int              prealloc;    /* Memory is preallocated, array should not free */
  int              *ref_count;  /* Reference counter for shared memory */
  array_data_type  *data;       /* Pointer to memory storing contents of array */
} array;


/***********************************************************************
 * Function:  init_array
 * Purpose:   Sets the size of the array and allocates memory to store 
 *            contents. The memory is released when free_array is called 
 *            on the last array using that memory.
 * Arguments: a - pointer to the array to initialize.
 *            n1, n2, n3 - Size of the array.
 * Return:    0 on success, nonzero on failure.
 * Note:      Calling init_array on an array which is already initialized,
 *            can cause memory leaks. 
 */

int init_array(array *a, 
               array_index_type n1, 
               array_index_type n2, 
               array_index_type n3);

/***********************************************************************
 * Function:  init_array_alloc
 * Purpose:   Initialize array with extra allocated memory.The memory is 
 *            released when free_array is called on the last array using 
 *            that memory.
 * Arguments: a - pointer to the array to initialize.
 *            n1, n2, n3 - Size of the array.
 *            nalloc - Amount of memory to allocate. Must be at least n1*n2*n3.
 * Return:    0 on success, nonzero on failure.
 * Note:      Calling init_array on an array which is already initialized,
 *            can cause memory leaks. 
 */

int init_array_alloc(array *a, 
                     array_index_type n1, 
                     array_index_type n2, 
                     array_index_type n3,
                     array_index_type nalloc);

/***********************************************************************
 * Function:  init_array_shared
 * Purpose:   Initialize array a using the same memory as another array b.
 * Arguments: a - pointer to the array to initialize.
 *            n1, n2, n3 - Size of the array. b must have enough memory 
 *            allocated to accomodate the data of a.
 *            b - pointer to an already initialized array.
 * Return:    0 on success, nonzero on failure.
 * Note:      Calling init_array_shared on an array which is already 
 *            initialized, can cause memory leaks. 
 */

int init_array_shared(array *a, 
                      array_index_type n1, 
                      array_index_type n2, 
                      array_index_type n3, 
                      array *b);

/***********************************************************************
 * Function:  init_array_prealloc
 * Purpose:   Sets the size of the array and sets data storage to 
 *            previously allocated memory. The memory is never freed by 
 *            array functions.
 * Arguments: a - pointer to the array to initialize.
 *            n1, n2, n3 - size of array.
 *            x - pointer to allocated memory. Must point to a block of 
 *            memory containing at least n1*n2*n3*sizeof(array_data_type) 
 *            bytes.
 * Return:    0 on success, nonzero on failure.
 * Note:      Calling init_array_prealloc on an array which is already 
 *            initialized, can cause memory leaks. 
 */

int init_array_prealloc(array *a, 
                        array_index_type n1, 
                        array_index_type n2, 
                        array_index_type n3, 
                        void *x);


/***********************************************************************
 * Function:  free_array
 * Purpose:   Free memory allocated by the initialization functions. 
 *            After free_array has been called, the array can be 
 *            reinitialized by calling init_array, init_array_shared or 
 *            init_array_prealloc. Making several subseqent calls to
 *            free_array has no further effect.
 * Arguments: a - pointer to the array to free. 
 */

void free_array(array *a);


/***********************************************************************
 * Function:  array_initialized
 * Purpose:   Check if array is ready to use. Mainly to be used by other
 *            array functions for debugging. 
 * Arguments: a - pointer to the array to check.
 * Return:    1 if array is OK to use, 0 otherwise. If the array is not 
 *            OK to use a warning is printed.
 */

int array_initialized(const array *a);


/***********************************************************************
 * Function:  nele_array
 * Purpose:   Return the total number of elements in array.
 * Arguments: a - pointer to array.
 * Return:    Total number of elements in a. 
 */

array_index_type nele_array(const array *a);


/***********************************************************************
 * Function:  array_index_in_range
 * Purpose:   Check if a triple is a valid index into array.
 * Arguments: a - Pointer to array.
 *            i, j, k - Zero-based indices to check
 * Return:    Nonzero if (i,j,k) is a valid index into a.
 */

int array_index_in_range(const array *a, 
                         array_index_type i, 
                         array_index_type j, 
                         array_index_type k);


/***********************************************************************
 * Function:  set_array_entry
 * Purpose:   Set array entry to given value.
 * Arguments: a - Pointer to array
 *            i, j, k - Zero-based indices of entry to set.
 *            x - New value of entry.
 * Note:      If i, j, k are not valid indices into a (as indicated by
 *            array_index_in_range) set_array_entry does nothing.
 */

void set_array_entry(array *a, 
                     array_index_type i, 
                     array_index_type j, 
                     array_index_type k, 
                     array_data_type x);


/***********************************************************************
 * Function:  add_to_array_entry
 * Purpose:   Add number to array entry.
 * Arguments: a - Pointer to array
 *            i, j, k - Zero-based indices of entry to add to.
 *            x - Value to add to entry.
 * Note:      If i, j, k are not valid indices into a (as indicated by
 *            array_index_in_range) add_to_array_entry does nothing.
 */

void add_to_array_entry(array *a, 
                        array_index_type i, 
                        array_index_type j, 
                        array_index_type k, 
                        array_data_type x);


/***********************************************************************
 * Function:  get_array_entry
 * Purpose:   Get value of array entry.
 * Arguments: a - Pointer to array
 *            i, j, k - Zero-based indices of entry to get.
 * Return:    Value of array entry. If i, j, k are not valid indices into
 *            a the function returns zero.
 */

array_data_type get_array_entry(const array *a, 
                                array_index_type i, 
                                array_index_type j, 
                                array_index_type k);


/***********************************************************************
 * Function:  get_array_entry_ptr
 * Purpose:   Get pointer to array entry.
 * Arguments: a - Pointer to array
 *            i, j, k - Zero-based indices of entry to get.
 * Return:    Pointer to array entry. If i, j, k are not valid indices 
 *            into a the function returns a null pointer.
 */

array_data_type *get_array_entry_ptr(array *a, 
                                     array_index_type i, 
                                     array_index_type j, 
                                     array_index_type k);


/***********************************************************************
 * Function:  fill_array
 * Purpose:   Set all entries in array to the same value.
 * Arguments: a - Pointer to array.
 *            x - Value which is set in all entries.
 */

int fill_array(array *a, 
               array_data_type x);


/***********************************************************************
 * Function:  same_size_array
 * Purpose:   Check if arrays have the same size.
 * Arguments: a, b - Pointers to arrays.
 * Return:    Nonzero if arrays have the same size.
 */

int same_size_array(const array *a, 
                    const array *b);


/***********************************************************************
 * Function:  copy_array
 * Purpose:   Copy contents of array a to array b.
 * Arguments: a - Pointer to array to copy from.
 *            b - Pointer to array to copy to.
 * Return:    Zero on success, nonzero on failure. Fails if the arrays do
 *            not have the same size.
 */

int copy_array(const array *a, 
               array *b);


/***********************************************************************
 * Function:  add_array_const
 * Purpose:   Add same value to all entries of an array.
 * Arguments: a - Pointer to array.
 *            x - Value to add to each entry.
 * Return:    Zero on success, nonzero on failure.
 */

int add_array_const(array *a, 
                    array_data_type x);


/***********************************************************************
 * Function:  add_array
 * Purpose:   Add multiple of array a to array b.
 * Arguments: a, b - Pointers to arrays.
 *            x - Scalar which a is multiplied by before adding to b.
 * Return:    Zero on success, nonzero on failure. Fails if the arrays do
 *            not have the same size.
 */

int add_array(const array *a, 
              array *b, 
              array_data_type x);


/***********************************************************************
 * Function:  add_array_offset
 * Purpose:   Add array a to array b with offset. Arrays need not have 
 *            the same size. a(i1, j1, k1) is added to b(i1+i, j1+j, k1+k)
 *            for all indices (i1, j1, k1) such that indices to both 
 *            arrays are in range.
 * Arguments: a, b - Pointers to arrays.
 *            i, j, k - Offset.
 * Return:    Zero on success, nonzero on failure.
 */

int add_array_offset(array *a, 
                     array *b, 
                     array_index_type i, 
                     array_index_type j, 
                     array_index_type k);


/***********************************************************************
 * Function:  mult_array_const
 * Purpose:   Multiply all array entries by a constant.
 * Arguments: a - Pointer to array.
 *            x - Value each entry in a is multiplied by.
 * Return:    Zero on success, nonzero on failure.
 */

int mult_array_const(array *a, 
                     array_data_type x);


/***********************************************************************
 * Function:  mult_array
 * Purpose:   Multiply each element of b by corresponding element of a.
 * Arguments: a, b - Pointers to arrays.
 * Return:    Zero on success, nonzero on failure. Fails if the arrays do
 *            not have the same size. 
 */

int mult_array(const array *a, 
               array *b);


/***********************************************************************
 * Function:  norm_sq_array
 * Purpose:   Compute the sum of the squares of the entries of array.
 * Arguments: a - Pointer to array.
 * Return:    Sum of squares of entries of a. 
 */

array_data_type norm_sq_array(const array *a);


/***********************************************************************
 * Function:  max_array
 * Purpose:   Compute the maximum entry of array.
 * Arguments: a - Pointer to array.
 * Return:    Maximum entry of a. 
 */

array_data_type max_array(const array *a);


/***********************************************************************
 * Function:  min_array
 * Purpose:   Compute the minimum entry of array.
 * Arguments: a - Pointer to array.
 * Return:    Minimum entry of a. 
 */

array_data_type min_array(const array *a);


/***********************************************************************
 * Function:  mean_array
 * Purpose:   Compute the mean value of entries of array.
 * Arguments: a - Pointer to array.
 * Return:    Mean value of entries of a. 
 */

array_data_type mean_array(const array *a);


/***********************************************************************
 * Function:  boundary_mean_array
 * Purpose:   Compute the mean value of entries along the boundary of 
 *            array. Boundary elements are those for which at least one
 *            of the indices is either 0 or its maximal value.
 * Arguments: a - Pointer to array.
 * Return:    Mean value of entries along boundary of a. 
 */

array_data_type boundary_mean_array(const array *a);


/***********************************************************************
 * Function:  laplace_array
 * Purpose:   Compute the discrete laplacian of array a and put the 
 *            result in array b. For the computation of the laplacian
 *            at the boundaries, a is treated as padded with zeros.
 * Arguments: a, b - Pointers to arrays.
 * Return:    Zero on success, nonzero on failure. Fails if the arrays do
 *            not have the same size. 
 */

int laplace_array(const array *a, 
                  array *b);

#endif
