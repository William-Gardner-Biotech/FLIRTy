import sys

# 1: primer string, 2: amplicon bed file, 3: output amplicon indv bed file path

with open(sys.argv[2], 'r') as amplicons:
    primer = sys.argv[1]
    primer = primer.split('\t')
    print(primer[1], primer[2])
    for raw_amp in amplicons:
        amp = raw_amp.split('\t')
        if primer[1] == amp[1] or primer[1] == amp[2] or primer[2] == amp[1] or primer[2] == amp[2]:
            print(f"Amp: {amp[1], amp[2]}")
            with open(sys.argv[3], 'w') as bed_file:
                bed_file.write(raw_amp)