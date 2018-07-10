import unittest
from pyMiniSam.pyMiniSam import pyMiniSam
import pysam


FILENAME = "/Users/anthony/temp/subset.bam"


def read_to_dict(read):
    return {'reference_start': read.reference_start,
            'reference_end': read.reference_end,
            'positions': read.positions,
            'reference_name': read.reference_name}

def short_dict(read):
    return {'reference_start': read['start'],
            'reference_end': read['end'],
            'positions': read['positions'],
            'reference_name': read['reference_name']}


class TestMiniSam(unittest.TestCase):

    @staticmethod
    def read_bam_filefirstline_pySam():
        file_handle = pysam.AlignmentFile(FILENAME, 'rb')
        for read in file_handle.fetch(until_eof=True):
            file_handle.close()
            return read_to_dict(read)

    @staticmethod
    def read_bam_filefirstline_pyMiniSam():
        reads = pyMiniSam(FILENAME)
        for read in reads.get_reads():
            reads.close()
            return short_dict(read)

    def test_get_reads(self):
        read1 = self.read_bam_filefirstline_pyMiniSam()
        read2 = self.read_bam_filefirstline_pySam()
        print(f"pyMiniSam = {read1}")
        print(f"pySam = {read2}")
        self.assertEqual(read1['reference_start'], read2['reference_start'])
        self.assertEqual(read1['reference_end'], read2['reference_end'])
        self.assertEqual(read1['positions'], read2['positions'])
        self.assertEqual(read1['reference_name'], read2['reference_name'])
