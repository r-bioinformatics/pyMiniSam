#cython:
cimport cython

from libc.stdio cimport FILE
from libc.stdio cimport *
from libc.string cimport memset

cdef unsigned int BYTES_4 = 4

cdef char cigar_convert(unsigned int a):
    if a == 0:
        return 'M'
    elif a == 1:
        return 'I'
    elif a == 2:
        return 'D'
    elif a == 3:
        return 'N'
    elif a == 4:
        return 'S'
    elif a == 5:
        return 'H'
    elif a == 6:
        return 'P'
    elif a == 7:
        return '='
    elif a == 8:
        return 'X'

cdef char seq_convert(unsigned int a):
    if a == 0:
        return '='
    elif a == 1:
        return 'A'
    elif a == 2:
        return 'C'
    elif a == 3:
        return 'M'
    elif a == 4:
        return 'G'
    elif a == 5:
        return 'R'
    elif a == 6:
        return 'S'
    elif a == 7:
        return 'V'
    elif a == 8:
        return 'T'
    elif a == 9:
        return 'W'
    elif a == 10:
        return 'Y'
    elif a == 11:
        return 'H'
    elif a == 12:
        return 'K'
    elif a == 13:
        return 'D'
    elif a == 14:
        return 'B'
    elif a == 15:
        return 'N'



cdef unsigned int get_bytes_from_list_core(bytes o, unsigned int start, unsigned int n_bytes):
    cdef unsigned int myInt1 = 0
    cdef unsigned int ONE = 1
    cdef unsigned int ZERO = 0

    cdef unsigned int i = 0
    while i < n_bytes:
        myInt1 += o[start+i] << (i*8)
        i += 1
    return myInt1

cpdef tuple get_bytes_from_list(bytes o, unsigned int start, unsigned int n_bytes):
    myInt1 = get_bytes_from_list_core(o, start, n_bytes)
    return myInt1, start + n_bytes


cpdef unsigned int get_bits_as_int_from_bam(object bam, unsigned int n_bytes):
    return int.from_bytes(bam.read(n_bytes), byteorder='little')

cpdef tuple get_string_from_list(bytes o, unsigned int start, unsigned int n_bytes):
    return o[start:start + n_bytes].decode('utf-8')[:-1], start + n_bytes


cpdef tuple get_cigar_from_list(bytes o, unsigned int start, unsigned int n_cigar_ops, unsigned int ul_seq):
    # if n_cigar_ops == 2 and int.from_bytes(o[start], byteorder='little') == l_seq:
    #      return None, start + BYTES_4


    cdef unsigned int op_start, op_end, op, bases
    cdef unsigned int seq

    cdef unsigned int n = 0
    cigar_ops = []
    while n < n_cigar_ops:
        op_start = start + (BYTES_4 * n)
        op_end = op_start + BYTES_4
        seq = int.from_bytes(o[op_start:op_end], byteorder='little')
        bases = seq >> 4
        op = seq & 15
        cigar_ops.append(str(bases))
        cigar_ops.append(chr(cigar_convert(op)))
        n += 1

    return "".join(cigar_ops), start + (n_cigar_ops * BYTES_4)


cpdef tuple get_quality_from_list(bytes o, unsigned int start, unsigned int n_bytes):

    cdef char quality[1000]
    cdef unsigned int ONE = 1
    cdef unsigned int ZERO = 0
    memset(quality, 0, n_bytes)
    cdef unsigned int x = 0
    while x < n_bytes:
        quality[x] = o[start + x] + 33
        x += 1
    return quality[:n_bytes], start + n_bytes


cpdef tuple get_seq_from_list(bytes o, unsigned int start, unsigned int length):

    cdef unsigned int b, b1, b2, p=0
    cdef char c1, c2
    cdef char sequence[1000]
    # memset(sequence, 0, length)
    cdef unsigned int ONE = 1
    cdef unsigned int ZERO = 0

    cdef unsigned int n = 0
    while n < length:

        b = o[start + n]
        b1 = b >> 4
        b2 = b & 15
        c1 = seq_convert(b1)
        c2 = seq_convert(b2)
        if c1 != '=':
            sequence[p] = c1
        else:
            break
        p += 1
        if c2 != '=':
            sequence[p] = c2
        else:
            break
        p += 1
        n+=1

    return sequence[:p], start + length


cpdef tuple get_extra_flags_from_bam(bam, unsigned int start):
    return None, start