from re import search, findall
from itertools import groupby, accumulate
from collections import Counter
import numpy as np
import pandas as pd

def valid_sequence(seq, rna=False):
    """
    Test the legaltiy of input seq. If not, raise a error, otherwise return None.
        seq: input sequence.
        rna: Boolen. If True, the input sequence should be RNA consisting of AUGC. Otherwise should be DNA consisting of ATGC. Default is False.
    """
    if not isinstance(seq, str):
        raise TypeError('Input sequence should be string!')
    nuc = set('ATGC')
    if not (set(seq.upper()) <= nuc):
        raise ValueError('Input sequence contains letters other than ATGC!')

def split_seq(seq, wrap=60):
    """Split input seq into pieces with no longer than wrap."""
    valid_sequence(seq)
    return findall('.{1,%d}' % wrap, seq)

def reverse_seq(seq, rna=False):
    """
    Return the reverse complement sequence of input seq.
        rna: Boolen. If True, the input sequence should be RNA. Otherwise should be DNA. Default is False.
    """
    valid_sequence(seq, rna=rna)
    if rna:
        rev_nuc = str.maketrans('AUGCaugc', 'UACGuacg')
    else:
        rev_nuc = str.maketrans('ATGCatgc', 'TACGtacg')
    return seq[::-1].translate(rev_nuc)

def calc_GC(seq):
    """
    Calculate GC content of input seq.
    """
    valid_sequence(seq)
    seq_stat = Counter(seq.upper())
    GC = (seq_stat['G'] + seq_stat['C']) / len(seq)
    return GC

def ePCR(target_seq, primer5, primer3):
    """
    Test the PCR product by input target_seq, primer5 and primer3.
    Return a tuple, first element is boolen whether there exists products by giving primer pair, second element is the sequence of product (empty while not exist).
    Raise ValueError when multiple location found in target_seq.
    """
    valid_sequence(target_seq)
    ix = False
    pcr_seq = ''
    for seq in (target_seq, reverse_seq(target_seq)):
        cnt5, cnt3 = seq.count(primer5), seq.count(primer3)
        if cnt5 > 1 or cnt3 > 1:
            raise ValueError('Multiple location of primer pair in target sequence. Primer5: %d, Primer3: %d.' % (cnt5, cnt3))
        pos5 = seq.find(primer5)
        pos3 = reverse_seq(seq).find(primer3)
        pos3 = len(seq) - pos3 - 1 if pos3 != -1 else -1
        if pos3 >= pos5 != -1:
            pcr_seq = seq[pos5:(pos3+1)]
            ix = True
    return ix, pcr_seq

def translate(cds):
    """
    Translate input cds into protein sequence. Input cds should be divided by 3 and without IUPAC.
    """
    RNA2DNA = lambda seq: seq.upper().replace('U', 'T')
    cds = RNA2DNA(cds)
    valid_sequence(cds)
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
    prot = ''
    if len(cds) % 3:
        raise ValueError('The CDS sequence seems to have wrong length which cannot be divided by 3!')
    prot = ''.join((codon_table[codon] for codon in split_seq(cds, 3)))
    return prot

def read_fasta(fasta_file):
    """
    Read fasta file.
    Return a generator of many tuples. First element is sequence header (between > and the first space). Second element is sequence.
    """
    parser_header = lambda header: header.strip('\n>').split()[0]
    with open(fasta_file) as fin:
        fa_reader = groupby(fin, lambda line: line.startswith('>'))
        for fa_rec in fa_reader:
            if fa_rec[0]:
                header_list = list(fa_rec[1])
                if len(header_list) == 1:
                    header = parser_header(header_list[0])
                    seq_iter = next(fa_reader)
                    seq = ''.join(line.strip() for line in seq_iter[1])
                    yield header, seq
                else:
                    for header in header_list[:-1]:
                        yield parser_header(header), ''
                    header = parser_header(header_list[-1])
                    seq_iter = next(fa_reader)
                    seq = ''.join(line.strip() for line in seq_iter[1])
                    yield header, seq

def read_gtf_file(gtf_file, selected_type=('CDS', 'exon')):
    """
    Read GTF file. Only return records whose types are in selected_type (default is CDS and exon).
    Return a pandas.DataFrame with 7 columns:
        1: seq_name
        2: seq_strand
        3: seq_start
        4: seq_end
        5: seq_type
        6: gene_id
        7: transcript_id
    """
    def generate_gtf_selected_stream(gtf_file, selected_type=None):
        with open(gtf_file) as fin:
            for line in fin:
                if not line.startswith('#'):
                    cnts = line.rstrip('\n').split('\t')
                    if (not selected_type) or (cnts[2] in selected_type):
                        yield cnts
    def parse_gtf_cnts(cnts):
        if not isinstance(cnts, list):
            raise ValueError('Input data should be a list but %s' % type(cnts))
        if len(cnts) != 9:
            raise ValueError('The length of input list should be 9 but %d' % len(cnts))
        def parse_gtf_features(col9):
            tid_p = r'transcript_id "([^;]*)"'
            gid_p = r'gene_id "([^;]*)"'
            try:
                tid = search(tid_p, col9).group(1)
            except:
                raise ValueError('%s lost transcript_id!' % col9)
            try:
                gid = search(gid_p, col9).group(1)
            except:
                gid = ''
            return tid, gid
        tid, gid = parse_gtf_features(cnts[8])
        out_rec = {
                'seq_name': cnts[0],
                'seq_strand': cnts[6],
                'seq_start': int(cnts[3]),
                'seq_end': int(cnts[4]),
                'seq_type': cnts[2],
                'gene_id': gid,
                'transcript_id': tid
        }
        return out_rec
    gtf_iter = generate_gtf_selected_stream(gtf_file, selected_type=selected_type)
    anno_db = [parse_gtf_cnts(x) for x in gtf_iter]
    return pd.DataFrame(anno_db, columns=['seq_name', 'seq_strand', 'seq_start', 'seq_end', 'seq_type', 'gene_id', 'transcript_id'])

def calc_assembly_NL(seq_len_list, quantile=50):
    """
    Calculate N50 or L90 to evaluate quanlity of genome assembly.
    Input seq_len_list is a list or generator of length of all sequences.
    Input quantile 50 for L50 and N50, 90 for L90 and N90.
    Return a tuple with first element L-quantile and second element N-quantile.
    """
    len_list = sorted(seq_len_list, reverse=True)
    genome_size = sum(len_list)
    cutoff = genome_size*(quantile/100)
    for i, accu_len in enumerate(accumulate(len_list)):
        if len_acc >= cutoff:
            return i, len_list[i]

def get_sequence(seq_id, anno_db, genome, seq_type='CDS', exon_split='', flank_len=0):
    """
    Get sequence of gene by seq_id. The anno_db and genome is required. Default seq_type is CDS. Output sequence can be splited at junction loci by exon_split.
    """
    def get_sequence_from_genome_by_anno_db(df, genome):
        tmp_seq = genome[df['seq_name']]
        tmp_seq_start = df['seq_start']-1
        tmp_seq_end = df['seq_end']
        df['seq'] = tmp_seq[tmp_seq_start:tmp_seq_end]
        return df
    gene_db = anno_db.query(f'seq_type == "{seq_type}" and transcript_id == "{seq_id}"')
    gene_db = gene_db.sort_values(by='seq_start')
    gene_db = gene_db.apply(get_sequence_from_genome_by_anno_db, axis=1, genome=genome)
    seq_name = gene_db['seq_name'].unique()
    seq_strand = gene_db['seq_strand'].unique()
    if len(seq_name) > 1 or len(seq_strand) > 1:
        raise TypeError('Different seq_name or seq_strand for the same gene!')
    else:
        seq_name = seq_name[0]
        seq_strand = seq_strand[0]
    if isinstance(exon_split, str):
        cds_seq = exon_split.join(gene_db['seq'])
    else:
        raise TypeError('Input exon_split should be a string but %s' % type(exon_split))
    if flank_len:
        gene_start = gene_db.iloc[0]['seq_start']
        gene_end = gene_db.iloc[-1]['seq_end']
        tmp_seq = genome[seq_name]
        flank5 = tmp_seq[(gene_start-flank_len-1):(gene_start-1)].lower()
        flank3 = tmp_seq[gene_end:(gene_end+flank_len)].lower()
    else:
        flank5 = ''
        flank3 = ''
    cds_seq = flank5 + cds_seq + flank3
    if np.all(gene_db['seq_strand'] == '-'):
        cds_seq = reverse_seq(cds_seq)
    elif np.all(gene_db['seq_strand'] == '+'):
        pass
    else:
        raise ValueError('Different strand in the elements of %s' % seq_id)
    return cds_seq
