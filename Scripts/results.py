ls = [[] for i in range(16)]

nxt = lambda x : float(x.readline().strip().split()[-1])
for i in range(1, 11):
    with open('SRR1031159_1_out/out'+str(i)+('.txt'), 'r') as f:
        f.readline()
        for i in range(16):
            ls[i].append(nxt(f))

txt = [
    "Basic init time: ",
    "Basic no. of kmers: ",
    "Basic query time: ",
    "Basic FPR: ",
    "One sided init time: ",
    "1 Sided query time: ",
    "One sided FPR: ",
    "Two sided init time: ",
    "Two sided potential edges: ",
    "Two sided edges: ",
    "2 Sided query time: ",
    "Two sided FPR: ",
    "Best match init time: ",
    "Best match no. of kmers: ",
    "Best fit query time: ",
    "Best fit FPR: "
]

for l, t in zip(ls, txt):
    print(t, sum(l) / len(l))
