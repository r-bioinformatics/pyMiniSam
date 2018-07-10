#cython: languagelevel=3, boundscheck=False, wraparound=False, nonecheck=False

cimport cython

from libc.stdio cimport FILE
from libc.stdio cimport *
from libc.stdlib cimport malloc, free
from libc.string cimport strcpy
from libc.string cimport memset

DEF MAX_STR_SIZE = 100
cdef unsigned int BYTES_4 = 4
cdef int BITS_PER_BYTE = 8
cdef char QUAL_OFFSET = 33
cdef char INT_TO_CHAR = 48
cdef int ONE_INT = 1
cdef int ZERO = 0

@cython.optimize.use_switch(True)
cdef char cigar_convert(unsigned int a) nogil:
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

@cython.optimize.use_switch(True)
cdef char seq_convert(unsigned int a) nogil:
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


cdef struct cigar_struct:
    int end
    char* cigar_string
    int* positions
    int pos_len

cdef unsigned int get_bytes_from_list_core(bytes o, int start, int n_bytes):
    cdef long myInt1 = 0
    cdef int i = 0
    while i < n_bytes:
        myInt1 += <long>o[start+i] << (i*BITS_PER_BYTE)
        i += 1
    return <unsigned int>myInt1

cdef int get_bytes_from_list(bytes o, int start, int n_bytes):
    cdef int myInt1
    myInt1 = get_bytes_from_list_core(o, start, n_bytes)
    return myInt1


cpdef int get_bits_as_int_from_bam(object bam, int n_bytes):
    cdef bytes temp
    cdef int myInt1
    cdef int i = 0
    temp = bam.read(n_bytes)
    if temp == b'':
        return 0
    while i < n_bytes:
        myInt1 += (<int>temp[i]) << (i*BITS_PER_BYTE)
        i += 1
    return myInt1

cdef char* get_string_from_list(bytes o, int start, int n_bytes):

    cdef int n = 0
    cdef char* stringiness = <char *>malloc(n_bytes)
    while n < n_bytes:
        stringiness[n] = <char>o[start + n]
        n += 1

    #return makestring(stringiness, n_bytes)
    return stringiness

@cython.cdivision(True)
cdef int itoa_dest(int i, char* dest, int dest_pos) nogil:
    cdef char rev[MAX_STR_SIZE]
    cdef int digit_count = 0
    cdef int digit_reverse = 0
    while i >0:
        digit = i % 10
        rev[digit_count] = <char>digit + INT_TO_CHAR
        digit_count += 1
        i /= 10

    digit_count -= 1
    while digit_count >= 0:

        dest[dest_pos + digit_reverse] = rev[digit_count]
        digit_count -= 1
        digit_reverse += 1

    return dest_pos + digit_reverse


cdef cigar_struct get_cigar_from_list(bytes o, int start, int n_cigar_ops, int l_seq, int coord):
    # if n_cigar_ops == 2 and int.from_bytes(o[start], byteorder='little') == l_seq:
    #      return None, start + BYTES_4

    cdef cigar_struct result
    cdef int op_start, op_end, bases
    cdef int seq
    cdef char* str_pointer
    cdef int cigar_string_pos = 0
    cdef int sum_pos = 0
    cdef int n = 0  # op_count
    cdef int p = 0  # positions index
    cdef int i = 0  # iterator
    cdef char op, op_type
    result.cigar_string = <char *>malloc((l_seq+1)* sizeof(char))
    result.positions = <int *>malloc((l_seq*+1) * sizeof(int))
    memset(result.positions, 0, (l_seq+1)*sizeof(int))


    while n < n_cigar_ops:
        op_start = start + (BYTES_4 * n)
        op_end = op_start + BYTES_4
        seq = get_bytes_from_list_core(o, op_start, BYTES_4)
        bases = seq >> 4
        op = seq & 15

        cigar_string_pos = itoa_dest(bases, result.cigar_string, cigar_string_pos)
        op_type = cigar_convert(op)
        result.cigar_string[cigar_string_pos] = op_type

        if op_type == 'M': # or op_type == 'D' or op_type == '=' :
            for i in range(bases):
                result.positions[p] = coord + i + sum_pos
                p += 1

        if op_type != 'H' and op_type != 'S':
            sum_pos = sum_pos + bases

        cigar_string_pos += 1
        n += 1


    result.pos_len = p
    result.end = sum_pos + coord
    result.cigar_string[cigar_string_pos] = '\0'

    return result


cdef char* get_quality_from_list(bytes o, int start, int n_bytes):

    cdef char * quality = <char *> malloc(n_bytes+1)
    cdef int x = 0
    while x < n_bytes:
        quality[x] = <char>(o[start + x]) + QUAL_OFFSET
        x += 1
    quality[x] ='\0'
    #return makestring(quality[:n_bytes], n_bytes)
    return quality


cdef char* get_seq_from_list(bytes o, int start, int length):

    cdef int b, b1, b2, p=0
    cdef char c1, c2
    cdef char* sequence = <char *>malloc(((length+1)*2)+ 1)
    memset(sequence, 0, ((length+1)*2)+ 1)

    cdef int n = 0
    while n < length:

        b = <int>o[start + n]
        b1 = b >> 4
        b2 = b & 15
        c1 = seq_convert(b1)
        c2 = seq_convert(b2)
        if c1 != '=':
            sequence[p] = c1
        else:
            sequence[p] = '\0'
        p += 1
        if c2 != '=':
            sequence[p] = c2
        else:
            sequence[p] = '\0'
        p += 1
        n+=1
    # return makestring(sequence[:p], p)
    return sequence


cpdef object get_extra_flags_from_bam(bam, int start):
    return None, start


cpdef dict get_read(object bam, dict references):

    cdef bytes get_record
    cdef unsigned int pos = 0
    cdef unsigned int size_to_read
    cdef unsigned int ref_id
    cdef unsigned int coord
    cdef unsigned int l_read_name
    cdef unsigned int mapq
    cdef unsigned int bin_num
    cdef unsigned int n_cigar_op
    cdef unsigned int flag
    cdef unsigned int l_seq
    cdef unsigned int next_refID
    cdef unsigned int next_pos, tlen
    cdef char* read_name
    cdef char* cigar
    cdef char* seq
    cdef char* qual
    cdef int p = 0
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
    bin_num = get_bytes_from_list(get_record, pos, 2)
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
    cigar_struct = get_cigar_from_list(get_record, pos, n_cigar_op, l_seq, coord)
    pos += (n_cigar_op * BYTES_4)
    seq = get_seq_from_list(get_record, pos, (l_seq + 1) // 2)
    pos += (l_seq +1)//2
    qual = get_quality_from_list(get_record, pos, l_seq)
    pos += l_seq
    #
    # extra_flag, pos = get_extra_flags_from_bam(get_record, pos)
    #
    #

    aligned_read = {'reference_name': references[ref_id]['name'],
                    'start': coord,
                    'end': cigar_struct.end,
                    'positions': [cigar_struct.positions[p] for p in range(cigar_struct.pos_len)],
                    'mapq': mapq,
                    'bin': bin_num,
                    'flag': flag,
                    'next_refID': next_refID,
                    'next_pos': next_pos,
                    'tlen': tlen,
                    'read_name': read_name.decode('utf-8'),
                    'cigar': cigar_struct.cigar_string.decode('utf-8'),
                    'seq': seq.decode('utf-8'),
                    'qual': qual.decode('utf-8')}



    # free(read_name)
    # free(cigar_struct.positions)
    # free(cigar_struct.cigar_string)
    #
    # free(seq)
    # free(qual)

    return aligned_read

    # assert(pos == size_to_read)
