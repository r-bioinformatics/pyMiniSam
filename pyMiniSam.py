import gzip
import argparse
from time import time

# import pyximport; pyximport.install()
from sam_functions import get_bits_as_int_from_bam
from sam_functions import get_string_from_list
from sam_functions import get_quality_from_list
from sam_functions import get_seq_from_list
from sam_functions import get_cigar_from_list
from sam_functions import get_extra_flags_from_bam
from sam_functions import get_bytes_from_list


def parse_arguments(parser=None):
    if not parser:
        parser = argparse.ArgumentParser()

    parser.add_argument('--filename', help="bam to read", required=True)

    args = parser.parse_args()
    return args


class AlignedRead(object):

    def __init__(self, ref, coord,
                 mapq=None,
                 bin=None,
                 flag=None,
                 next_refID=None,
                 next_pos=None,
                 tlen=None,
                 read_name=None,
                 cigar=None,
                 seq=None,
                 qual=None):
        self.ref = ref
        self.coord = coord
        self.mapq = mapq 
        self.bin = bin
        self.flag = flag
        self.next_refID = next_refID
        self.next_pos = next_pos
        self.tlen = tlen
        self.read_name = read_name
        self.cigar = cigar
        self.seq = seq
        self.qual = qual

    def __repr__(self):
        return f"{self.read_name}\t{self.flag}\t{self.ref}\t{self.coord}\t{self.mapq}\t{self.seq}\t{self.qual}\t{self.cigar}"


class pyMiniSam(object):

    def __init__(self, args):
        self.filename = args.filename
        self.references = {}

    def read_compressed_file(self):
        count = 0
        time0 = time()
        with gzip.open(self.filename, "rb") as bam:

            # read headers
            data = bam.read(4)
            header_len = get_bits_as_int_from_bam(bam, 4)
            header_text = bam.read(header_len).decode("utf-8")[:-1]
            reference_sequences_count = get_bits_as_int_from_bam(bam, 4)

            for x in range(reference_sequences_count):
                lname = get_bits_as_int_from_bam(bam, 4)
                ref_name = bam.read(lname).decode("utf-8")[:-1]
                ref_lengh = get_bits_as_int_from_bam(bam, 4)
                self.references[x] = {'name': ref_name, 'length': ref_lengh}

            while True:
                size_to_read = get_bits_as_int_from_bam(bam, 4)
                if size_to_read == 0:
                    break
                count += 1
                pos = 0
                get_record = bam.read(size_to_read)

                ref_id, pos = get_bytes_from_list(get_record, pos, 4)  # at position 0
                coord, pos = get_bytes_from_list(get_record, pos, 4)   # at position 4
                l_read_name, pos = get_bytes_from_list(get_record, pos, 1)
                mapq, pos = get_bytes_from_list(get_record, pos, 1)
                bin, pos = get_bytes_from_list(get_record, pos, 2)
                n_cigar_op, pos = get_bytes_from_list(get_record, pos, 2)
                flag, pos = get_bytes_from_list(get_record, pos, 2)
                l_seq, pos = get_bytes_from_list(get_record, pos, 4)
                next_refID, pos = get_bytes_from_list(get_record, pos, 4)
                next_pos, pos = get_bytes_from_list(get_record, pos, 4)
                tlen, pos = get_bytes_from_list(get_record, pos, 4)

                seq = None
                qual = None

                read_name, pos = get_string_from_list(get_record, pos, l_read_name)
                cigar, pos = get_cigar_from_list(get_record, pos, n_cigar_op, l_seq)
                seq, pos = get_seq_from_list(get_record, pos, (l_seq+1)//2)
                # qual, pos = get_quality_from_list(get_record, pos, l_seq)
                #
                # extra_flag, pos = get_extra_flags_from_bam(get_record, pos)
                #
                #
                aligned_read = AlignedRead(self.references[ref_id]['name'], coord+1,
                                           mapq, bin,
                                           flag, next_refID,
                                           next_pos, tlen,
                                           read_name, cigar,
                                           seq, qual)

                # print(aligned_read)
                # return
                #assert(pos == size_to_read)

        print(f"{count} records read in {time()-time0} seconds")


def main():
    args = parse_arguments()
    default_class = pyMiniSam(args)
    default_class.read_compressed_file()


if __name__ == '__main__':
    main()