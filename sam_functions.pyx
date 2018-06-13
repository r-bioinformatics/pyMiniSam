from libc.stdio cimport FILE
from libc.stdio cimport *
from libc.string cimport memset
cimport cython

BYTES_4 = 4

CIGAR_CONVERT = {
    0: 'M',
    1: 'I',
    2: 'D',
    3: 'N',
    4: 'S',
    5: 'H',
    6: 'P',
    7: '=',
    8: 'X',
}

SEQ_CONVERT = {
    0: '=',
    1: 'A',
    2: 'C',
    3: 'M',
    4: 'G',
    5: 'R',
    6: 'S',
    7: 'V',
    8: 'T',
    9: 'W',
    10: 'Y',
    11: 'H',
    12: 'K',
    13: 'D',
    14: 'B',
    15: 'N',
}


cdef char seq_convert(int a):
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



cpdef get_bytes_from_list(object o, int start, int n_bytes):
    cdef unsigned int myInt1 = 0
    for i in range(n_bytes):
        myInt1 += o[start+i] << (i*8)
    return myInt1, start + n_bytes


cpdef unsigned int get_bits_as_int_from_bam(object bam, int n_bytes):
    return int.from_bytes(bam.read(n_bytes), byteorder='little')

cpdef get_string_from_list(object o, int start, int n_bytes):
    return o[start:start + n_bytes].decode('utf-8')[:-1], start + n_bytes


cpdef get_cigar_from_list(object o, int start, int n_cigar_ops, int l_seq):
    # if n_cigar_ops == 2 and int.from_bytes(o[start], byteorder='little') == l_seq:
    #      return None, start + BYTES_4

    cigar_ops = []

    for n in range(n_cigar_ops):
        op_start = (BYTES_4 * n) + start
        op_end = op_start + BYTES_4
        seq = int.from_bytes(o[op_start:op_end], byteorder='little')
        bases = seq >> 4
        op = seq & 15
        cigar_ops.append(str(bases) + CIGAR_CONVERT[op])

    return cigar_ops, start + (n_cigar_ops * BYTES_4)


cpdef get_quality_from_list(object o, start, n_bytes):

    cdef char quality[1000]
    memset(quality, 0, n_bytes)
    for x in range(n_bytes):
        quality[x] = o[start + x] + 33
    return quality[:n_bytes], start + n_bytes


cpdef get_seq_from_list(object o, int start, int length):

    cdef int b, b1, b2, p=0
    cdef char c1, c2
    cdef char sequence[1000]
    memset(sequence, 0, length)
    for n in range(length):

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

    return sequence[:p], start + length


cpdef get_extra_flags_from_bam(bam, start):
    return None, start