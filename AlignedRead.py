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
        return f"{self.read_name}\t{self.flag}\t{self.ref}\t" \
               f"{self.coord}\t{self.mapq}\t{self.seq}\t{self.qual}\t{self.cigar}"