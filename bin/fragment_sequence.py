from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse

# This function will cut the sequence up into a bunch of smaller fasta sequences and combine into a multifasta file
def chop(seqID, sequence, chop_length:int):

    if len(sequence) < chop_length:
        exit(f'sequence too short to be chopped by {chop_length} segments')

    # Find the number of chops to perform
    clean_partition = False
    while clean_partition == False:
        print(f"Chopping with {int(chop_length)}")
        # Conditional to check that our final seq is at least 90% the length of other sequences
        if len(sequence) % chop_length >= chop_length*0.7 and len(sequence) % chop_length <= chop_length*0.9:
            clean_partition = True
            print(f'Found ideal chop length: {int(chop_length)}')
        elif chop_length < 20:
            exit(f'Failed to find ideal chop length')
        else:
            # Reduce the size of the segments until we get a good length to maintain the final seq is long enough
            chop_length *= 0.9

    number_of_segments = int(len(sequence)//(chop_length))
    # NOTE: Lots of int() rounding of floats, could be causing numerical shifts

    # Length, number of segments, length of segments
    print(f'{len(sequence)}, {number_of_segments}, {int(chop_length)}')

    new_seq_list = []
    # Remove float normalizing error
    chop_length = int(chop_length)
    # For loop that will chop the seqs down
    for x in range(0, number_of_segments+1): 
        print(f'{x}, {sequence[:(chop_length)]}')
        new_seq = Seq(sequence[:(chop_length)])

        # Append a description to each new sequence that tracks the starting position relative to the reference
        descript = f"Subsegment: {x+1}, Relative_Pos: {chop_length*(x)}"
        new_seq_record = SeqRecord(new_seq, id=f'{seqID}_subsegment_{x+1}', description=descript)
        sequence = sequence[int(chop_length)::]
        new_seq_list.append(new_seq_record)

    return new_seq_list

def main():
    """
    A program that takes a single fasta file and breaks it apart into many sequence fragments of a chosen length.\n\
    The program seeks to split the sequence evenly to ensure the final segment is still useful.\n\
    It uses basic arithmetic to ensure that the final segment will be within 70-90% of the other segments length.\n\
    To accomplish this ideal segment length the input fragment size may be modified by reducing it's size by 90% each iteration.\n\
    This program is designed for large sequences and whose fragments should not be smaller than 100 bp.
    """
    parser = argparse.ArgumentParser(
        description="Program designed to take a fasta file containing one sequence and fragment it into segments based upon user's requested size\n\
            NOTE: input fastas should be single sequence and fragmenting below 100bp is not advised",
            formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('--input', '-i', required=False, type=str, default='data/example.fasta',
                    help='Path to the input fasta file \
                        RUN > python3 bin/fragment_sequence.py -i data/example.fasta')
    
    parser.add_argument('--output', '-o', required=False, type=str, default='data/example_fragmented.fasta',
                    help='Name/path of the output fasta file \
                        RUN > python3 bin/fragment_sequence.py -i data/example.fasta -o data/example_fragmented.fasta')
    
    parser.add_argument('--frag_size', '-f', required=False, type=int, default=500,
                        help='Desired size of fragments\n\
                            RUN > python3 bin/fragment_sequence.py -i data/example.fasta -o data/example_fragmented.fasta -f 500')
    
    args = parser.parse_args()

    fasta_sequences = SeqIO.parse(open(args.input), 'fasta')

    # Get the first element from the iterator
    first_record = next(fasta_sequences, None)

    if first_record:
        # Check if there is only one sequence
        second_record = next(fasta_sequences, None)
        if second_record is not None:
            # More than one sequence, print a message
            exit("Error: Multiple sequences found. Processing only the first sequence.")
        else:
            # Only one sequence, proceed
            name, sequence = first_record.id, str(first_record.seq)
            #print(f'{name} {sequence}')
            chopped_segments = chop(name, sequence, args.frag_size)
            total_frag_length = sum(len(str(record.seq)) for record in chopped_segments)
            if total_frag_length == len(sequence):
                SeqIO.write(chopped_segments, args.output, "fasta")
            else:
                exit('PROCESS TERMINATED: Lost nucleotides to fragmenting error')

    else:
        exit("Error: No sequences found in the input file.")

if __name__ == "__main__":
    main()