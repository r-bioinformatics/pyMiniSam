#cython: languagelevel=3, boundscheck=False, wraparound=False, nonecheck=False

cimport cython

from libc.stdio cimport FILE
from libc.stdio cimport *
from libc.stdlib cimport malloc, free
from libc.string cimport strcpy

from AlignedRead import AlignedRead

DEF MAX_STR_SIZE = 1000
cdef unsigned int BYTES_4 = 4
cdef int BITS_PER_BYTE = 8
cdef char QUAL_OFFSET = 33
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

cdef char* makestring(char* src, int length):
    cdef char * c_string = <char *> malloc((length + 1) * sizeof(char))
    if not c_string:
        raise MemoryError()
    strcpy(c_string, src)
    return c_string

cpdef int get_bits_as_int_from_bam(object bam, int n_bytes):
    cdef bytes temp
    cdef int myInt1
    cdef int i = 0
    temp = bam.read(n_bytes)
    if temp == b'':
        return 0
    # if not temp or len(temp) < n_bytes:
    #     return 0

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


cdef char* get_cigar_from_list(bytes o, int start, int n_cigar_ops, int ul_seq):
    # if n_cigar_ops == 2 and int.from_bytes(o[start], byteorder='little') == l_seq:
    #      return None, start + BYTES_4


    cdef int op_start, op_end, op, bases
    cdef int seq
    cdef char *result

    cdef int n = 0
    cigar_ops = []
    while n < n_cigar_ops:
        op_start = start + (BYTES_4 * n)
        op_end = op_start + BYTES_4
        seq = get_bytes_from_list_core(o, op_start, BYTES_4)
        bases = seq >> 4
        op = seq & 15
        cigar_ops.append(str(bases))
        cigar_ops.append(chr(cigar_convert(op)))
        n += 1

    return makestring("".join(cigar_ops).encode('utf-8'), len(cigar_ops))
    #return "".join(cigar_ops).encode('utf-8')


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
    cdef char* sequence = <char *>malloc((length+1)*2)
    # memset(sequence, 0, length)

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
            break
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

    aligned_read = {'ref': references[ref_id]['name'],
                    'coord': coord + 1,
                    'mapq': mapq,
                    'bin': bin_num,
                    'flag': flag,
                    'next_refID': next_refID,
                    'next_pos': next_pos,
                    'tlen': tlen,
                    'read_name': read_name.decode('utf-8'),
                    'cigar': cigar,
                    'seq': seq.decode('utf-8'),
                    'qual': qual.decode('utf-8')}


    # aligned_read = AlignedRead(references[ref_id]['name'],
    #                            coord + 1,
    #                            mapq,
    #                            bin_num,
    #                            flag,
    #                            next_refID,
    #                            next_pos,
    #                            tlen,
    #                            read_name,
    #                            cigar,
    #                            seq,
    #                            qual)

    free(read_name)
    # free(cigar)
    free(seq)
    free(qual)

    return aligned_read
    # assert(pos == size_to_read)
