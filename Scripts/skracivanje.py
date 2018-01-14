with open('SRR1031159_1_full.fasta', 'r') as f:
    head = f.readline().strip()
    while head:
        seq = f.readline().strip()
        if seq.find('N') == -1:
            print(head)
            print(seq)
        head = f.readline().strip()
