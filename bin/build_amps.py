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
        # Recursive function to find and create the amplicons WITH primers
        OG_len = len(pos)
        # TODO add a check for one unmatched FWD and one unmatched REV that would technically be even but not correspond to one another
        def build_out(positives, negatives, original_length, ref_seq, ampli_bed):
            """
            Recursive function that takes the list of FWD and REV primers and then looks for the minimums of each list.
            Once found they are considered a pair because tiled amplicons will follow the logic of first primer set are furthest left on the sequence and therefore minimums.
            We pop these out of the lists and add them to a growing list of the amplicon file.
            """
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
            # Simple way to make the names unique and ascend as it goes along
            name1 = original_length-len(positives)

            minimum = 149000000000
            for count, feat in enumerate(negatives):
                # RE will find the second match as this is where the matching minimum is located on the negatives
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

            # bed formatting the string
            bed_line = f"{ref_seq}\t{pos1}\t{pos2}\tAMP_{name1}\t+"
            ampli_bed.append(bed_line)

            # Check when the list has been popped out so we can return a value
            if positives:
                return build_out(positives, negatives, OG_len, ref_seq, ampli_bed)
            else:
                return ampli_bed
            
        list_of_amplicons = build_out(pos, neg, OG_len, ref_seq_name, ampli_bed=[])
        bulk_string = '\n'.join(list_of_amplicons)

    elif args.slice == 'in':
        # Recursive function to find and create the amplicons WITHOUT primers
        OG_len = len(pos)
        def build_in(positives, negatives, original_length, ref_seq, ampli_bed):
            """
            Recursive function that takes the list of FWD and REV primers and then looks for the minimums of each list.
            Using one simple swap from the first to the second coordinate we get the smaller of the two numbers respectively finding the smaller amplicon.
            Once found they are considered a pair because tiled amplicons will follow the logic of first primer set are furthest left on the sequence and therefore minimums.
            We pop these out of the lists and add them to a growing list of the amplicon file.
            """
            # Make the default bigger than the largest organism's genome (Japanese canopy plant)
            minimum = 149000000000
            for count, feat in enumerate(positives):
                # RE will find the second match as this is where the matching minimum is located on the negatives
                matches = re.findall(r'\b\d+\b', feat)
                num1 = matches[1]
                if int(num1) < minimum: 
                    minimum = int(num1)
                    index = count
            pos1 = minimum
            positives.pop(index)
            # Simple way to make the names unique and ascend as it goes along
            name1 = original_length-len(positives)

            minimum = 149000000000
            for count, feat in enumerate(negatives):
                matches = re.findall(r'\b\d+\b', feat)
                num2 = matches[0]
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
                return build_in(positives, negatives, OG_len, ref_seq, ampli_bed)
            else:
                return ampli_bed
            
        list_of_amplicons = build_in(pos, neg, OG_len, ref_seq_name, ampli_bed=[])
        # Combine the list by newlines to make a writeable bed format
        bulk_string = '\n'.join(list_of_amplicons)

    with open(args.output, 'w') as bedfile:
        bedfile.write(bulk_string)
