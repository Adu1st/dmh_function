#!/usr/local/python3.6/bin/python
from intervaltree import IntervalTree
from re import search
import numpy as np
import pandas as pd

def parse_features(features, filetype='gtf'):
    if filetype == 'gtf':
        fpattern = r'([^\s]*) +([^;]*)'
    elif filetype == 'gff3':
        fpattern = r'([^\s]*) *=* ([^;]*)'
    else:
        raise ValueError('Parameter filetype should be gtf or gff3')
    for cnts in features.strip().split(';'):
        cnts = cnts.strip()
        p = search(fpattern, cnts)
        value_id = p.group(1)
        value_data = p.group(2).strip().strip('"')
        if value_data.isdigit():
            value_data = int(value_data)
        else:
            try:
                value_data = float(value_data)
            except:
                pass
        yield value_id, value_data
def parse_annotation_line(line, filetype='gtf'):
    cnts = line.strip().split('\t')
    seq_chr = cnts[0]
    seq_type = cnts[2]
    seq_start = int(cnts[3])
    seq_end = int(cnts[4])
    seq_strand = cnts[6]
    feature_dict = dict(parse_features(cnts[8], filetype=filetype))
    if filetype == 'gtf':
        seq_id = feature_dict['transcript_id']
        seq_parent = feature_dict['gene_id']
    elif filetype == 'gff3':
        seq_id = feature_dict['ID']
        seq_parent = feature_dict['Parent']
    else:
        raise ValueError('Parameter filetype should be gtf or gff3')
    return seq_chr, seq_strand, seq_start, seq_end, seq_type, seq_id, seq_parent
def read_gtf_file(gtf_file, selected_type=('exon', 'CDS')):
    anno_type = np.dtype([
        ('seq_chr', 'U'),
        ('seq_strand', 'S2'),
        ('seq_start', np.uint32),
        ('seq_end', np.uint32),
        ('seq_type', 'U'),
        ('gene_id', 'U'),
        ('transcript_id', 'U')
        ])
    def generate_gtf_stream(gtf_file, selected_type):
        with open(gtf_file) as fin:
            for line in fin:
                cnts = parse_annotation_line(line)
                if cnts[4] in selected_type:
                    yield cnts
    anno_db = np.array(
            [x for x in generate_gtf_stream(gtf_file, selected_type)],
            dtype=anno_type
            )
    return pd.DataFrame(anno_db)
