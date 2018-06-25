import gzip
import argparse
import io
from time import time
from subprocess import Popen, PIPE

# import pyximport; pyximport.install()
from sam_functions import get_bits_as_int_from_bam
from sam_functions import get_read


def parse_arguments(parser=None):
    if not parser:
        parser = argparse.ArgumentParser()

    parser.add_argument('--filename', help="bam to read", required=True)

    args = parser.parse_args()
    return args


class pyMiniSam(object):

    def __init__(self, args):
        self.filename = args.filename
        self.references = {}

    def read_compressed_file(self):
        count = 0
        time0 = time()

        p = Popen(["gunzip", "-c", self.filename], stdout=PIPE)
        bam = io.BytesIO(p.communicate()[0])
        if p.returncode != 0:

            print("Something went wrong")
            return None
            # read headers
        data = bam.read(4)
        header_len = get_bits_as_int_from_bam(bam, 4)
        header_text = bam.read(header_len).decode("utf-8")[:-1]
        reference_sequences_count = get_bits_as_int_from_bam(bam, 4)

        for x in range(reference_sequences_count):
            lname = get_bits_as_int_from_bam(bam, 4)
            ref_name = bam.read(lname)[:-1].decode('utf-8')
            ref_lengh = get_bits_as_int_from_bam(bam, 4)
            self.references[x] = {'name': ref_name, 'length': ref_lengh}

        while True:
            read = get_read(bam, self.references)
            if read is None:
                break
            print(read)
            count += 1

        print(f"{count} records read in {time()-time0} seconds")


def main():
    args = parse_arguments()
    default_class = pyMiniSam(args)
    default_class.read_compressed_file()


if __name__ == '__main__':
    main()