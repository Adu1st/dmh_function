#!/usr/local/python3.6/bin/python3
from itertools import groupby, accumulate
from re import findall
from collections import Counter
from io import IOBase

def read_fasta(fasta_file):
    def _read_fasta(fin):
        def _parse_fa_header(header):
            return header.strip('>\n').split()[0]
        fa_reader = groupby(fin, lambda line: line.startswith('>'))
        for fa_rec in fa_reader:
            if fa_rec[0]:
                header_list = list(fa_rec[1])
                if len(header_list) == 1:
                    header = _parse_fa_header(header_list[0])
                    seq_iter = next(fa_reader)
                    seq = ''.join(line.strip() for line in seq_iter[1])
                    yield header, seq
                else:
                    for header in header_list[:-1]:
                        yield _parse_fa_header(header), ''
                    header = _parse_fa_header(header_list[-1])
                    seq_iter = next(fa_reader)
                    seq = ''.join(line.strip() for line in seq_iter[1])
                    yield header, seq
    if isinstance(fasta_file, IOBase):
        _read_fasta(fasta_file)
    elif isinstance(fasta_file, str):
        with open(fasta_file) as fin:
            _read_fasta(fin)
    else:
        raise TypeError('Input fasta file should be a string or file handle!')

def reverse_seq(seq):
    rev_nuc = str.maketrans('ATGCNatgc', 'TACGNtacg')
    return seq[::-1].translate(rev_nuc)

def RNA2DNA(seq):
    rna2dna = str.maketrans('Uu', 'Tt')
    return seq.translate(rna2dna)
def DNA2RNA(seq):
    dna2rna = str.maketrans('Tt', 'Uu')
    return seq.translate(dna2rna)

def split_seq(seq, wrap=60):
    return findall('.{0,%d}' % wrap, seq)

def translate(cds):
    codon_table = {'TTT': 'F', 'CTT': 'L', 'ATT': 'I', 'GTT': 'V',
        'TTC': 'F', 'CTC': 'L', 'ATC': 'I', 'GTC': 'V',
        'TTA': 'L', 'CTA': 'L', 'ATA': 'I', 'GTA': 'V',
        'TTG': 'L', 'CTG': 'L', 'ATG': 'M', 'GTG': 'V',
        'TCT': 'S', 'CCT': 'P', 'ACT': 'T', 'GCT': 'A',
        'TCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A',
        'TCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A',
        'TCG': 'S', 'CCG': 'P', 'ACG': 'T', 'GCG': 'A',
        'TAT': 'Y', 'CAT': 'H', 'AAT': 'N', 'GAT': 'D',
        'TAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D',
        'TAA': '', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E',
        'TAG': '', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E',
        'TGT': 'C', 'CGT': 'R', 'AGT': 'S', 'GGT': 'G',
        'TGC': 'C', 'CGC': 'R', 'AGC': 'S', 'GGC': 'G',
        'TGA': '', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G',
        'TGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'}
    if not cds//3:
        raise ValueError('Length of input sequence is not divisible by 3!')
    prot = ''
    cds = RNA2DNA(cds)
    cds = cds.upper()
    for i, codon in enumerate(split_seq(cds, 3)):
        prot += codon_table.get(codon, 'X')
    return prot

def _calc_assembly_NL_cutoff(genome_size, quantile):
    return genome_size / (quantile/100)
def _calc_NLx(seq_len_list, quantile):
    len_list = sorted(seq_len_list, reverse=True)
    cutoff = _calc_assembly_NL_cutoff(sum(len_list), quantile)
    for i, len_acc in enumerate(accumulate(len_list)):
        if len_acc >= cutoff:
            return i, len_list[i]
def calc_assembly_Nx(seq_len_list, quantile=50):
    return _calc_NLx(seq_len_list, quantile)[1]
def calc_assembly_Lx(seq_len_list, quantile=50):
    return _calc_NLx(seq_len_list, quantile)[0] + 1

def calc_GC(seq):
    seq_stat = Counter(seq.upper())
    if isinstance(tmpSeq, str):
        GC = seq_stat['G'] + seq_stat['C']
        return GC / (GC + seq_stat['A'] + seq_stat['T'])
    else:
        raise ValueError('Input sequence must be string!')

