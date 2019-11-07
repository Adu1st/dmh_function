#!/usr/local/python3.6/bin/python
from intervaltree import IntervalTree

class GRange:
    def __init__(self, chromosome, strand, start, end):
        self.chromosome = chromosome
        self.strand = strand
        self.start = start
        self.end = end
    def get_start_end(self):
        return self.start, self.end

class Gene(GRange):
    def __init__(self, chromosome, strand, start, end):
        super().__init__(chromosome, strand, start, end)
        self.transcript = {}
    def __setitem__(self, transcript_id, transcript):
        if isinstance(transcript, Transcript):
            self.transcript.setdefault(transcript_id, transcript)
            self._update_boundary()
        else:
            raise ValueError('Must be Trancsript!')
    def __getitem__(self, transcript_id):
        return self.transcript[transcript_id]
    def __delitem__(self, transcript_id):
        del self.transcript[transcript_id]
    def __str__(self):
        return 'Chr: %s\nStrand: %s\nStart: %d\nEnd: %d\n%d transcripts: %s' % (self.chromosome, self.strand, self.start, self.end, len(self.transcript), ','.join(self.transcript))
    def __repr__(self):
        return 'Chr: %s\nStrand: %s\nStart: %d\nEnd: %d\n%d transcripts: %s' % (self.chromosome, self.strand, self.start, self.end, len(self.transcript), ','.join(self.transcript))

    def add_transcript(self, transcript_id, start=0, end=0):
        self.transcript.setdefault(transcript_id, Transcript(self.chromosome, self.strand, start, end))
        self._update_boundary()
    def _update_boundary(self):
        boundary = [t.get_start_end() for t in self.transcript.values()]
        self.start = min(x[0] for x in boundary)
        self.end = max(x[1] for x in boundary)
    def get_length(self):
        gene_tree = IntervalTree()
        for t in self.transcript.values():
            for e in t.exon:
                gene_tree.addi(e[0], e[1])
        gene_tree.merge_overlaps()
        return sum(x.end-x.begin+1 for x in gene_tree)

class Transcript(GRange):
    def __init__(self, chromosome, strand, start, end):
        super().__init__(chromosome, strand, start, end)
        self.exon = []
        self.CDS = []
    def __str__(self):
        return 'Chr: %s\nStrand: %s\nStrat: %d\nEnd: %d\nLength: %d\nCoding length: %d\nExons: %d\n%s' \
            % (self.chromosome, self.strand, self.start, self.end, self.get_length(), self.get_length('coding'), len(self.exon), ' '.join(str(x) for x in self.exon))
    def __repr__(self):
        return 'Chr: %s\nStrand: %s\nStrat: %d\nEnd: %d\nLength: %d\nCoding length: %d\nExons: %d\n%s' \
            % (self.chromosome, self.strand, self.start, self.end, self.get_length(), self.get_length('coding'), len(self.exon), ' '.join(str(x) for x in self.exon))
    def add_element(self, start, end, element_type='exon'):
        if element_type == 'exon':
            self.exon.append((start, end))
            self.exon.sort(key=lambda x: x[0])
            self._update_boundary()
        elif element_type == 'CDS':
            self.CDS.append((start, end))
            self.CDS.sort(key=lambda x: x[0])
    def _update_boundary(self):
        self.start = self.exon[0][0]
        self.end = self.exon[-1][1]
    def get_length(self, length_type='transcript'):
        if length_type == 'transcript':
            return sum(x[1]-x[0]+1 for x in self.exon)
        elif length_type == 'coding':
            if self.CDS:
                return sum(x[1]-x[0]+1 for x in self.CDS)
            else:
                return 0
        else:
            raise ValueError('Parameter length_type should be either transcript or coding!')

