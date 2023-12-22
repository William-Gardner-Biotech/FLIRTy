import argparse
import re
"""
Python program that will tranform a primer bed file into
an amplicon bed file. It relies on finding a common name between two primers
and directionality.
"""

parser = argparse.ArgumentParser('I/O files and optional slicing conditions')

parser.add_argument('--input', '-i', required=True, type=str, default = 'data/SHIV_env_primers.bed',
                    help='Path to the input primer folder')

parser.add_argument('--output', '-o', required=False, default='data/generated_amplicons.bed',
                    help='File name for your amplicon bed file')

parser.add_argument('--slice', '-s', required=False, choices=['in', 'out'], default='out',
                    help='How much of the amplicon to include based on primers, \n\
                        in = inner coordinates which will exclude the primer sequences,\n\
                        out = outer coordinates used to allow primers included in amplicon')

args = parser.parse_args()

pos = []
neg = []

with open(args.input, 'r') as primers:
    for k in primers:
        k = k.strip('\n')
        print(k)
        split_k = k.split('\t')
        if split_k[5] == '+':
            # k[1:] removes the pesky seq name for regex ease
            pos.append('\t'.join(split_k[1:]))
        elif split_k[5] == '-':
            neg.append('\t'.join(split_k[1:]))
        else:
            exit('Unknown character found in primer bed file')
        ref_seq_name = split_k[0]
    for x in pos:
        print(x)
    
    if len(pos) != len(neg):
        exit('Uneven number of primers given')
    if args.slice == 'out':
        # Recursive function to find and create the amplicons
        OG_len = len(pos)
        def build_out(positives, negatives, original_length, ref_seq, ampli_bed):
            # Make the default bigger than the largest organism's genome (Japanese canopy plant)
            minimum = 149000000000
            for count, feat in enumerate(positives):
                num1 = re.search(r'\b\d+\b', feat)
                num1 = int(num1.group(0))
                if num1 < minimum: 
                    minimum = num1
                    index = count
            pos1 = minimum
            positives.pop(index)
            name1 = original_length-len(positives)

            minimum = 149000000000
            for count, feat in enumerate(negatives):
                matches = re.findall(r'\b\d+\b', feat)
                num2 = matches[1]
                if int(num2) < minimum: 
                    minimum = int(num2)
                    index = count
            pos2 = minimum
            negatives.pop(index)
            name2 = original_length-len(negatives)

            # Quick check that both lists reduced evenly
            if name1 != name2:
                exit('Error with matching the primers')

            bed_line = f"{ref_seq}\t{pos1}\t{pos2}\tAMP_{name1}\t+"
            ampli_bed.append(bed_line)

            if positives:
                return build_out(positives, negatives, OG_len, ref_seq, ampli_bed)
            else:
                return ampli_bed
            
        list_of_amplicons = build_out(pos, neg, OG_len, ref_seq_name, ampli_bed=[])

    for line in list_of_amplicons:
        print(line)
        
            