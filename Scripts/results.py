ls = [[] for i in range(20)]

nxt = lambda x : float(x.readline().strip().split()[-1])
for i in range(1, 11):
    with open('../TestData/1e4_out/out'+str(i)+('.txt'), 'r') as f:
        f.readline()
        for i in range(20):
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
    "Hitting set init time: ",
    "Hitting set no. of kmers: ",
    "Hitting set query time: ",
    "Hitting set FPR: "
]

for l, t in zip(ls, txt):
    print(t, sum(l) / len(l))
