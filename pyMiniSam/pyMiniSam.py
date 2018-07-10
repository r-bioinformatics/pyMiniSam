import argparse
from time import time
from subprocess import Popen, PIPE

import pyximport
pyximport.install()
from pyMiniSam.sam_functions import get_bits_as_int_from_bam
from pyMiniSam.sam_functions import get_read


def parse_arguments(parser=None):
    if not parser:
        parser = argparse.ArgumentParser()

    parser.add_argument('--filename', help="bam to read", required=True)

    args = parser.parse_args()
    return args


class pyMiniSam(object):

    def __init__(self, filename):
        self.filename = filename
        self.references = {}

        self.p = Popen(["gunzip", "-c", self.filename], stdout=PIPE)
        self.pipe = self.p.stdout
        data = self.pipe.read(4)
        header_len = get_bits_as_int_from_bam(self.pipe, 4)
        header_text = self.pipe.read(header_len).decode("utf-8")[:-1]
        self.reference_sequences_count = get_bits_as_int_from_bam(self.pipe, 4)

        for x in range(self.reference_sequences_count):
            lname = get_bits_as_int_from_bam(self.pipe, 4)
            ref_name = self.pipe.read(lname)[:-1]
            ref_lengh = get_bits_as_int_from_bam(self.pipe, 4)
            self.references[x] = {'name': ref_name.decode('utf-8'), 'length': ref_lengh}

    def get_reads(self):
        count = 0
        time0 = time()

        while True:
            read = get_read(self.pipe, self.references)
            if read is None:
                break
            yield read
            count += 1

        print(f"{count} records read in {time()-time0} seconds")
        self.close()

    def close(self):
        self.p.terminate()


def main():
    args = parse_arguments()
    pyminisam = pyMiniSam(args.filename)
    reads = pyminisam.get_reads()
    for read in reads:
        print(read)
        # print(read['ref'])



if __name__ == '__main__':
    main()