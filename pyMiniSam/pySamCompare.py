
import argparse
import pysam


def parse_arguments(parser=None):
    if not parser:
        parser = argparse.ArgumentParser()

    parser.add_argument('--filename', help="bam to read", required=True)

    args = parser.parse_args()
    return args


class BamParser(object):

    def __init__(self, args):
        self.filename = args.filename

    def read_compressed_file(self):
        file_handle = pysam.AlignmentFile(self.filename, 'rb')
        for read in file_handle.fetch(until_eof=True):
            print(read)

        file_handle.close()


def main():
    args = parse_arguments()
    default_class = BamParser(args)
    default_class.read_compressed_file()


if __name__ == '__main__':
    main()