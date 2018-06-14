#cython:
cimport cython

from libc.stdio cimport FILE
from libc.stdio cimport *
from libc.string cimport memset
from cpython cimport array

from AlignedRead import AlignedRead

DEF MAX_STR_SIZE = 1000
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



cdef unsigned int get_bytes_from_list_core(bytes o, unsigned int start, int n_bytes):
    cdef unsigned int myInt1 = 0
    cdef int i = 0
    while i < n_bytes:
        myInt1 += o[start+i] << (i*8)
        i += 1
    return myInt1

cdef unsigned int get_bytes_from_list(bytes o, unsigned int start, int n_bytes):
    cdef unsigned int myInt1
    myInt1 = get_bytes_from_list_core(o, start, n_bytes)
    return myInt1


cpdef unsigned int get_bits_as_int_from_bam(object bam,  int n_bytes):
    cdef bytes temp
    cdef unsigned int myInt1
    cdef int i = 0
    temp = bam.read(n_bytes)
    if len(temp) < n_bytes:
        return 0

    while i < n_bytes:
        myInt1 += temp[i] << (i*8)
        i += 1
    return myInt1

cdef str get_string_from_list(bytes o, unsigned int start, int n_bytes):
    if not o:
        return None, None

    cdef char[MAX_STR_SIZE] stringiness
    memset(stringiness, 0, MAX_STR_SIZE)
    cdef int n = 0
    while n < n_bytes:
        stringiness[n] = o[start + n]
        n += 1

    return stringiness.decode('utf-8')


cdef str get_cigar_from_list(bytes o, unsigned int start, unsigned int n_cigar_ops, int ul_seq):
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

    return "".join(cigar_ops)


cdef str get_quality_from_list(bytes o, unsigned int start, int n_bytes):

    cdef char quality[MAX_STR_SIZE]
    memset(quality, 0, MAX_STR_SIZE)
    cdef int x = 0
    while x < n_bytes:
        quality[x] = o[start + x] + 33
        x += 1
    return quality[:n_bytes].decode('utf-8')


cdef str get_seq_from_list(bytes o, unsigned int start, unsigned int length):

    cdef unsigned int b, b1, b2, p=0
    cdef char c1, c2
    cdef char sequence[MAX_STR_SIZE]
    # memset(sequence, 0, length)

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

    return sequence[:p].decode('utf-8')


cpdef object get_extra_flags_from_bam(bam, unsigned int start):
    return None, start

cpdef object get_read(object bam, dict references):

    cdef bytes get_record
    cdef unsigned int pos = 0
    cdef unsigned int size_to_read, ref_id, coord, l_read_name, mapq, bin, n_cigar_op, flag, l_seq, next_refID
    cdef unsigned int next_pos, tlen
    cdef str read_name, cigar, seq, qual



    size_to_read = get_bits_as_int_from_bam(bam, 4)
    if size_to_read == 0:
        return None
    get_record = bam.read(size_to_read)

    ref_id = get_bytes_from_list(get_record, pos, 4)
    pos += 4
    coord = get_bytes_from_list(get_record, pos, 4)
    pos += 4
    l_read_name = get_bytes_from_list(get_record, pos, 1)
    pos += 1
    mapq = get_bytes_from_list(get_record, pos, 1)
    pos +=1
    bin = get_bytes_from_list(get_record, pos, 2)
    pos += 2
    n_cigar_op = get_bytes_from_list(get_record, pos, 2)
    pos += 2
    flag = get_bytes_from_list(get_record, pos, 2)
    pos += 2
    l_seq = get_bytes_from_list(get_record, pos, 4)
    pos += 4
    next_refID = get_bytes_from_list(get_record, pos, 4)
    pos += 4
    next_pos = get_bytes_from_list(get_record, pos, 4)
    pos += 4
    tlen = get_bytes_from_list(get_record, pos, 4)
    pos += 4
    read_name = get_string_from_list(get_record, pos, l_read_name)
    pos += l_read_name
    cigar = get_cigar_from_list(get_record, pos, n_cigar_op, l_seq)
    pos += (n_cigar_op * BYTES_4)
    seq = get_seq_from_list(get_record, pos, (l_seq + 1) // 2)
    pos += (l_seq +1)//2
    qual = get_quality_from_list(get_record, pos, l_seq)
    pos += l_seq
    #
    # extra_flag, pos = get_extra_flags_from_bam(get_record, pos)
    #
    #
    aligned_read = AlignedRead(references[ref_id]['name'], coord + 1,
                               mapq, bin,
                               flag, next_refID,
                               next_pos, tlen,
                               read_name, cigar,
                               seq, qual)

    return aligned_read
    # assert(pos == size_to_read)