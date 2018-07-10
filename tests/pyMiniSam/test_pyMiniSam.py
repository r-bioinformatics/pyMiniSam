import unittest
from pyMiniSam.pyMiniSam import pyMiniSam
import pysam


FILENAME = "/Users/anthony/temp/subset.bam"


def pysam_to_dict(read):
    return {'reference_start': read.reference_start,
            'reference_end': read.reference_end,
            'positions': read.positions,
            'reference_name': read.reference_name,
            'cigar': read.cigarstring}

def pyminisam_to_dict(read):
    return {'reference_start': read['start'],
            'reference_end': read['end'],
            'positions': read['positions'],
            'reference_name': read['reference_name'],
            'cigar': read['cigar']}


class TestMiniSam(unittest.TestCase):

    @staticmethod
    def read_bam_filefirstline_pySam(fh):

        for read in fh.fetch(until_eof=True):
            yield pysam_to_dict(read)

    @staticmethod
    def read_bam_filefirstline_pyMiniSam(fh):

        for read in fh.get_reads():
            yield pyminisam_to_dict(read)

    def test_get_reads(self):
        pysam_file_handle = pysam.AlignmentFile(FILENAME, 'rb')
        pyminisam_file_handle = pyMiniSam(FILENAME)

        for read1, read2 in zip(self.read_bam_filefirstline_pyMiniSam(pyminisam_file_handle),
                                self.read_bam_filefirstline_pySam(pysam_file_handle)):

            if read1['reference_end'] != read2['reference_end']:
                print(read1['cigar'], read2['cigar'])

            self.assertEqual(read1['reference_start'], read2['reference_start'])
            self.assertEqual(read1['reference_end'], read2['reference_end'])
            self.assertEqual(read1['positions'], read2['positions'])
            self.assertEqual(read1['reference_name'], read2['reference_name'])

        pysam_file_handle.close()
        pyminisam_file_handle.close()