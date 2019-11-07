#!/usr/local/python3.6/bin/python3
def read_file(filename, end='\n'):
    with open(filename) as fin:
        for line in fin:
            yield line.rstrip(end)

def read_delim(filename, sep='\t', comment='#', end='\n'):
    def read_file(filename, end='\n'):
        with open(filename) as fin:
            for line in fin:
                yield line.rstrip(end)
    for line in read_file(filename, end=end)
        if not line.startswith(comment):
            yield line.split(sep)
