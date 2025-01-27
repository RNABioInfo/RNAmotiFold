# This script prints out the diff between the json and api sequences from the bgsu as of 3.88, hardcoded cause im lazy.
# There are some significant differences and Im considering if I should always keep the bulged and unbulged versions cause for UAA/GAN it
# seems like the unbuled is only part ? This is in reference to the realization that the json only contains unbulged versions while the api
# returns the full sequences. In the Automated Classification paper they mention that the unbulged verisons are the "core" of the motif
# and any nucleotide that does not form any interactions is not part of this core, which is whack cause then I'd have to decide if I want the
# unbulged nucleotides on a motif to motif basis? Shouldn't be an issue since I can put that into the json and just put an if unbulged -> get json_sequences
# condition into the update motif collection. But do I really want that ? Does it really make a difference ?
# Analysis for our paper was done including the unbulged sequences.
Tj = [
    "UUCGAGA",
    "UGAGAAU",
    "UUUGA",
    "UUUGU",
    "UGAAA",
    "UUCGA",
    "UUCAA",
    "UUAGA",
    "UUCGAUU",
    "UGAGAGU",
    "AGUAU",
    "UUCGACU",
    "UUCAAAU",
    "UGAAAAG",
    "UUCGAAU",
    "UUCAAUU",
    "UUCGAGC",
    "UUCAG",
    "UGGAA",
    "UGGGA",
    "UGCAACU",
    "UGAGA",
    "UUCGAUC",
    "UGAAACU",
    "UUCGAGU",
    "UGCAA",
    "CUUGA",
    "UUCAAGU",
    "AUUGA",
    "UGUAA",
]
Ta = [
    "UUCAGAU",
    "CUUGAAAC",
    "UGGAAC",
    "UUUGAAAA",
    "UGCAACU",
    "UGUAAUU",
    "UUCAAGU",
    "UUAGAUU",
    "UGCAACCC",
    "UUCGAGA",
    "UGGAAU",
    "UGAGAU",
    "UGCAAGU",
    "AUUGAAGC",
    "UUUGUAGA",
    "UGAGAAU",
    "UUCGAGC",
    "UGCAACUC",
    "UGGGAA",
    "UGAAAUU",
    "UGGAAA",
    "UGAGAGU",
    "UGAGAAC",
    "UUCGAUC",
    "UUCAAUU",
    "UUCGAAU",
    "UUCGAGU",
    "UGCAAGC",
    "AGUAAAU",
    "UGAAAAG",
    "UGGGAU",
    "UUCAAAU",
    "UGAAACU",
    "UUUGUAUA",
    "UUCGAUU",
    "UUCGACU",
]
Gj = [
    "UUAA",
    "GAAA",
    "GGAA",
    "GCAU",
    "GCAA",
    "UCAC",
    "GAGA",
    "GGGA",
    "CAAA",
    "GCGA",
    "GUAA",
    "UAGC",
    "UAAC",
    "GCAG",
    "GUGA",
    "AAAA",
    "UGAA",
    "GUAG",
]
Ga = [
    "UAAC",
    "UAGC",
    "GCAA",
    "GAUAAC",
    "GCAU",
    "GUAAG",
    "CAAA",
    "GGAAC",
    "GUAG",
    "GCAG",
    "UCAC",
    "GAAAG",
    "AAAA",
    "GAAA",
    "GGGA",
    "GUAA",
    "UGAA",
    "GGCGA",
    "GGAA",
    "UUCAA",
    "GAGA",
    "GUGA",
    "GAUAAG",
    "GAAAU",
    "GCGA",
    "GCUAAC",
]
Uj = []
Ua = ["CAAG", "GUGA", "CAGAU", "UUUG", "UACG", "UCCG", "GCUA", "UUCG"]
Lj = [
    "CUCGAAA",
    "UUCCCAA",
    "UUGGCAA",
    "UUGAAAA",
    "CUCAUAA",
    "CUGAUAA",
    "UUGGUAG",
    "CUGGAAA",
    "UUGGUAA",
    "UUGAGGU",
    "CUGUAAA",
    "UUCACGU",
    "CUCCAGA",
    "CUGAAAA",
    "GUCGUCG",
    "CUGUAGA",
    "CUUACAA",
    "UUCGGGA",
    "UUGCCAA",
    "CUGUCAC",
    "UUGAAGA",
]
La = [
    "CUGAUAA",
    "CUGUCAC",
    "CUUACAA",
    "UUAGGUAG",
    "CUGAAAA",
    "UUGCCAA",
    "UUCCCAA",
    "UUCACGU",
    "UUGGUAA",
    "CUCAUAA",
    "UUGAAAA",
    "CUGGAAA",
    "UUCGGGA",
    "CUGUAAA",
    "CUCGAAA",
    "CUCCAGA",
    "UUGAAGA",
    "UUGGCAA",
    "CUGUAGA",
    "GUCGUCG",
    "UUGAGGU",
]
Aj = ["UAA$GA", "UAA$UA", "UUA$GA", "CAA$GA"]
Aa = [
    "GAU$UAA",
    "UAA$GAA",
    "UUA$GAA",
    "CAA$GAAA",
    "UAA$UAU",
    "UAA$GAG",
    "UAA$GAU",
    "UAA$GAAA",
    "GAA$UAA",
    "UAA$GAAAU",
]
Mj = [
    "GU$GA",
    "AG$GG",
    "CU$CU",
    "UG$AC",
    "UG$GA",
    "UC$UC",
    "CA$GA",
    "CU$UU",
    "UA$GA",
    "GC$UU",
    "CA$CG",
    "AA$AG",
    "AU$AU",
    "GA$UA",
    "AG$AA",
    "CC$GC",
    "CU$UA",
    "GUUC",
    "AA$CA",
    "UC$UU",
    "CC$AC",
    "CG$GA",
    "GU$UU",
    "CA$CA",
    "GG$GA",
    "AA$GG",
    "UU$UU",
    "CC$CC",
    "CG$CA",
    "CC$CU",
    "GA$GA",
    "UU$UC",
    "CU$UC",
    "AC$AC",
]
Ma = [
    "UU$UC",
    "UC$UU",
    "CA$CCUUG",
    "AAGGA$AG",
    "CU$CU",
    "CG$GA",
    "UC$UC",
    "CC$AC",
    "GGCA$UG",
    "CU$CC",
    "AA$GG",
    "AC$AC",
    "CC$CC",
    "UU$UU",
    "UA$CU",
    "GCAU$AU",
    "CU$UU",
    "CC$GC",
    "UU$GU",
    "GU$GA",
    "UA$GA",
    "GA$GA",
    "AC$CC",
    "AGGGU$GG",
    "UC$CU",
    "CA$CA",
    "UU$GC",
    "GG$GUA",
    "GA$CA",
    "GGG$AA",
    "GG$AG",
    "GAA$GA",
    "AAG$AA",
    "GUUC",
    "CA$CG",
    "UG$AC",
    "CA$AA",
    "AGA$GG",
    "CU$UA",
]
Cj = [
    "CAC$AA",
    "GAC",
    "CCC",
    "AAU",
    "CAA$AA",
    "CGC",
    "GAA",
    "UCAA$A",
    "CAA",
    "CAC",
    "CAAU$A",
    "CUC",
    "GGC",
]
Ca = [
    "UCAA$GA",
    "GAC",
    "CAC$C",
    "AAU$CC",
    "CGC$U",
    "CAA",
    "CGAA$AA",
    "GAA$G",
    "CUC$C",
    "CUC$U",
    "CAC$G",
    "CC$CCC",
    "GGC",
    "CAC$AA",
    "U$CAC",
    "CAA$G",
    "CAAU$GA",
    "CAC$A",
]
Sj = [
    "GGGUA$AAGA",
    "UCAGUA$GAACC",
    "UAGUA$GAAC",
    "AGGUA$GAAUG",
    "GUA$GAG",
    "UUAGUA$GAACC",
    "GAGUA$GAAA",
    "GGGUA$GAGA",
    "AGGUA$AAAUA",
    "UCAGUA$GAACU",
    "AUGAGUA$GAAAGG",
    "UAAGUA$GAACU",
    "AGAACU$UUAUAC",
    "GAGUAC$AAAA",
    "GGUA$GAC",
    "UAAGUA$GAAAU",
    "CGAGUA$GAAAC",
    "CAGUA$GAAC",
    "UAGUA$GAAG",
    "AAGAGUA$GAAGCA",
    "UUAGUA$GAAAC",
    "GUUA$GAG",
    "GUA$AAA",
    "GAGUAC$GAAA",
    "UAGUA$GAAA",
    "GAAGCAACG$AAGUA",
    "GAGUA$AAAA",
    "CAGUA$GACC",
    "UGAGUA$GAAAU",
    "CAGUAGAA$AAGAAC",
    "CAGUAC$CGACC",
    "GAGUAU$GAAA",
    "AGUA$GAA",
    "UGAACU$UCAUAA",
    "GGGUA$GAAUA",
    "GUA$GAA",
    "UAGUA$GAAU",
]
Sa = [
    "UAGUA$GAAC",
    "UGAGUA$GAAAU",
    "CCGAGUA$GAAAC",
    "GUA$GAA",
    "AGGUA$AAAUA",
    "AUGUA$GAA",
    "AAAA$GAGUA",
    "CGACC$CAGUAC",
    "GAGUA$GAAA",
    "AAGAAC$CAGUAGAA",
    "GAGUAC$AAAA",
    "GAGUAC$GAAA",
    "GGGUA$GAAUA",
    "AAGAGUA$GAAGGCA",
    "GAC$GGUA",
    "AUGAGUA$GAAAGG",
    "GAACU$UCAGUA",
    "GGAGUA$GAGA",
    "UUAGUA$GAACC",
    "UAAGUA$GAAAU",
    "GAGUAU$GAAA",
    "AGGUA$GAAUG",
    "CAGUA$GAAC",
    "GGAGUA$AAGA",
    "AGAACU$UUAGUAC",
    "GACC$CAGUA",
    "GAACU$UAAGUA",
    "GUA$AAA",
    "UGAACU$UCAAUAA",
    "UAGUA$GAAU",
    "AUAGUA$GAAG",
    "AAGUA$GAAGCAACG",
    "UUAGUA$GAAAC",
    "UCAGUA$GAACC",
    "GUA$GAG",
    "GAG$GUUA",
    "UAGUA$GAAA",
]
Kj = [
    "GUGAU$GGA",
    "CUGACA$CGA",
    "GAGACA$CGA",
    "AUGAU$GGA",
    "AAGA$GA",
    "GAGAC$CAA",
    "AUGAA$GGA",
    "UUGAU$UCA",
    "GAGAU$GGA",
    "CCGACA$CGA",
    "C$ACCGGGG",
    "GAGAA$GGA",
    "GAGAC$ACA",
    "AAGAU$UGA",
    "UGGAA$GGA",
    "CGGAU$UGA",
    "GUGAA$GGA",
    "GGAA$GGAU",
    "GUGAG$GUA",
    "UUGAU$GGA",
    "GAGAU$AUA",
    "UAGAA$GGA",
    "GCAA$GA",
    "GAGAU$UGA",
    "UUAAC$ACA",
    "GAGAA$AGA",
    "GUGAC$CGA",
    "UAGACA$CGA",
    "CGAGAA$UGAA",
    "GAGAA$CGA",
    "GCGAUA$AGUA",
    "UGGAU$UGA",
    "UUGAA$GGA",
    "UAGAUU$UUGAUU",
    "GGAA$AGAA",
    "AGGAU$GGA",
    "GUUAAA$UAA",
    "AAGAC$ACA",
    "GAAA$GA",
    "GUUGAA$GGA",
    "GAGA$GA",
    "UUGAC$CGA",
    "GGGAA$GA",
    "GGGAA$GGA",
]
Ka = [
    "GGA$GAGGAU",
    "GGA$GAUGAA",
    "GGAAU$GAGAU",
    "GGA$GUUGAA",
    "UCA$UUGGAU",
    "GAAGAU$UGA",
    "UAUGAA$GGA",
    "GAUGA$GA",
    "UUUGAA$GGA",
    "GAAGAA$GGA",
    "GGA$UAUGAA",
    "CGA$GUUUGAC",
    "UGUGAU$UGA",
    "UGCGAA$GGA",
    "UGUGAA$GGA",
    "AUCGAU$GGA",
    "AAUGA$GA",
    "UAA$GUUAUAA",
    "GGAG$CGACCUUGAAAUAC",
    "AUUA$GACGAU",
    "UGA$GAAGAU",
    "UGA$CGUGAU",
    "GGA$GUGGAU",
    "UAAGACA$CGA",
    "GA$GAGAAA",
    "GGA$GAAGAA",
    "GGA$AGAGAU",
    "GUUGAA$GGA",
    "CAA$GAUGAC",
    "AUCA$GACGAC",
    "GGAA$AGGAA",
    "GGUGAA$GA",
    "UGA$UGUGAU",
    "AAUGAU$UGA",
    "GA$AAUGA",
    "AUCA$UUCAAC",
    "GUA$GUAGAG",
    "UAAGGAUU$UUGAUU",
    "GA$GCGAAA",
    "GGAGAA$GGA",
    "AGA$GACGAA",
    "CCAGACA$CGA",
    "AUGGAA$GGA",
    "UGAA$CGAUGAA",
    "GAAGACA$CGA",
    "CUAGACA$CGA",
    "GGA$UUUGAU",
    "GACGAU$GGA",
    "CGAUGAA$UGAA",
    "GACGAA$CGA",
    "GAUGAA$GGA",
    "AGUA$GCAGAUA",
    "CGA$UUUGAC",
    "AAGACA$AAUGAC",
]

print("T-Loops:")
print(set(Tj) - set(Ta))
print(set(Ta) - set(Tj))
print("GNRA:")
print(set(Gj) - set(Ga))
print(set(Ga) - set(Gj))
print("UNCG:")
print(set(Uj) - set(Ua))
print(set(Ua) - set(Uj))
print("Anticodon:")
print(set(Lj) - set(La))
print(set(La) - set(Lj))
print("UAA/GAN:")
print(set(Aj) - set(Aa))
print(set(Aa) - set(Aj))
print("Tandems:")
print(set(Mj) - set(Ma))
print(set(Ma) - set(Mj))
print("C-Loops:")
print(set(Cj) - set(Ca))
print(set(Ca) - set(Cj))
print("Sarcin-Ricin:")
print(set(Sj) - set(Sa))
print(set(Sa) - set(Sj))
print("Kink-Turn:")
print(set(Kj) - set(Ka))
print(set(Ka) - set(Kj))
