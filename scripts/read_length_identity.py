#!/usr/bin/env python

import sys
import collections
import statistics


def main():
    read_lengths = {}
    read_alignments = collections.defaultdict(list)

    paf_filename = sys.argv[1]
    with open(paf_filename, 'rt') as paf:
        for line in paf:
            paf_parts = line.strip().split('\t')
            if len(paf_parts) < 11:
                continue
            read_name = paf_parts[0]
            read_length = int(paf_parts[1])
            read_lengths[read_name] = read_length
            start = int(paf_parts[2])
            end = int(paf_parts[3])
            identity = 100.0 * int(paf_parts[9]) / int(paf_parts[10])
            read_alignments[read_name].append((start, end, identity))

    print('\t'.join(['Name', 'Length', 'Identity']))
    for read_name, read_length in read_lengths.items():
        alignments = read_alignments[read_name]
        identity_by_base = [0.0] * read_length
        for start, end, identity in alignments:
            for i in range(start, end):
                if identity > identity_by_base[i]:
                    identity_by_base[i] = identity
        whole_read_identity = statistics.mean(identity_by_base)
        print('\t'.join([read_name, str(read_length), str(whole_read_identity)]))


if __name__ == '__main__':
    main()