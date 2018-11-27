#cython: language_level=3
#cython: boundscheck=False
#cython: cdivision=False
#cython: wraparound=False

"""
cLev - a Cython extension to calculate the sequence levenshtein distance.

I found the following refrences particularly helpful while hacking this up:

    https://github.com/gfairchild/pyxDamerauLevenshtein
    https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance
    http://hackmap.blogspot.com/2008/04/levenshtein-in-cython.html
"""

from libc.stdlib cimport malloc, free

#------------------------------------------------------------------------------

cdef inline Py_ssize_t MIN3(Py_ssize_t a, Py_ssize_t b, Py_ssize_t c):
    """
    A faster c function to calc the min of 3 values. Note that we use Py_ssize_t
    for the types. Essentially, Py_ssize_t is a signed int that plays nicely
    with Python's internals and accomidates various architectures.
    See: http://stackoverflow.com/q/20987390
    """
    cdef Py_ssize_t m = a
    if m > b: m = b
    if m > c: m = c
    return m

#------------------------------------------------------------------------------

cpdef Py_ssize_t dist(str s1, str s2, bint debug=False):
    """
    Dangerously optimized version of the sequence levenshtein distance.

    Input:
        s1, s2 - strings of interest
        debug - print useful information
    Output:
        The Sequence Levenshtein Distance as defined in:
        doi:10.1186/1471-2105-14-272
    Tests:
        assert dist("CAGG", "CGTC") == 2
        assert dist("CGTC", "CAGG") == 2
        assert dist("TTCC", "TCCATGCATA") == 1
        assert dist("TCCATGCATA", "TTCC") == 1
        assert dist("ACAC", "TCCATGCATA") == 2
        assert dist("CGAA", "TCCATGCATA") == 3
        assert dist("TAGG", "TCCATGCATA") == 3
    """
    # reduce the search space by trimming any common prefixes
    cdef Py_ssize_t len_s1 = len(s1)
    cdef Py_ssize_t len_s2 = len(s2)
    cdef Py_ssize_t idx = 0
    while idx < len_s1 and idx < len_s2 and s1[idx] == s2[idx]:
        idx += 1
    s1 = s1[idx:]
    s2 = s2[idx:]

    if debug == True:
        print('s1: {}\ts2: {}'.format(s1, s2))

    # check for perfect substrings
    if not s1:
        return len_s2
    if not s2:
        return len_s1

    # make sure the longer string goes down the column and recalculate lengths
    # note we are overwriting old variables here!
    len_s1 = len(s1)
    len_s2 = len(s2)
    cdef Py_ssize_t tmp_len
    if len_s1 < len_s2:
        tmp_len = len_s1
        tmp_str = s1        # retains python object
        # assign 2 to 1
        len_s1 = len_s2
        s1 = s2
        # restore 1 into 2
        len_s2 = tmp_len
        s2 = tmp_str

        if debug == True:
            print('s1: {}\ts2: {}'.format(s1, s2))
            print('len_s1: {}\tlen_s2: {}'.format(len_s1, len_s2))

    # initiallize your matrix with malloc. DONT FORGET TO FREE WHEN FINISHED!
    cdef Py_ssize_t row = len_s1 + 1
    cdef Py_ssize_t col = len_s2 + 1
    cdef Py_ssize_t *mat = <Py_ssize_t *>malloc(row * col * sizeof(Py_ssize_t))
    if not mat:
        raise MemoryError()

    # wrap the following in a try-finally for saftey
    # cython will automatically optimize any `for i in range` if i has a ctype
    # we have to define all vars outside the try-finally block
    cdef Py_ssize_t r, c
    cdef Py_ssize_t dist, min_dist, test_r, test_c
    try:
        # initiallize the first row and col
        for c in range(0, col):
            mat[c] = c                      # M[0][0:c] = 0*col + c
        for r in range(0, row):
            mat[r*col] = r                  # M[0:r][0] = r*col + 0

        # begin the actual calculation
        # remember we start from index 1 since 0 has already been initiallized
        # row-wise indexing can be calc'd with M[x][y] = x*c + y, where c = #cols
        for r in range(1, row):
            for c in range(1, col):
                if s1[r-1] == s2[c-1]:
                    mat[r*col + c] = mat[(r-1)*col + c-1]
                else:
                    mat[r*col + c] = MIN3(
                           mat[(r-1)*col + c] + 1,    # deletion
                           mat[r*col + c-1] + 1,      # intertion
                           mat[(r-1)*col + c-1] + 1)  # substitution

        dist = mat[len_s1*col + len_s2]

        # tilo's additional operations, likely optimizable
        min_dist = dist
        for r in range(0, row):
            test_r = mat[r*col + len_s2]              # m[i][last] = i*c + last
            if test_r < min_dist: min_dist = test_r
        for c in range(0, col):
            test_c = mat[len_s1*col + c]              # m[last][i] = last*c + i
            if test_c < min_dist: min_dist = test_c

        # useful debugging information
        if debug == True:
            print('Dynamic Programming Matrix', flush=True)
            for r in range(0, row):
                for c in range(0, col):
                    print(mat[r*col + c], end=" ", flush=True)
                print("", flush=True)

            print('Row Slice')
            for r in range(0, row):
                print(mat[r*col + len_s2], end=" ", flush=True)
            print("", flush=True)

            print('Col Slice')
            for c in range(0, col):
                print(mat[len_s1*col + c], end=" ", flush=True)
            print("", flush=True)

    # free malloc'd memory and return
    finally:
        free(mat)
    return min_dist
