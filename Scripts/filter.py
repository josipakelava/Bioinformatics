with open('SRR1031159_full.fasta', 'r') as f:
    line = f.readline()
    while line:
        print(line.strip())
        print(f.readline().strip())
        f.readline()
        f.readline()
        line = f.readline()
