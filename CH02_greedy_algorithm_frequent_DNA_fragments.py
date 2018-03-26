#finds list of t k-length substrings of t strings
#minimizes mismatches - Laplace version
#outputs t k-length strings (list made into \n strings)
def greedy_motif_search_Laplace(k, t, Dna):
    best_motifs = []
    for i in range(len(Dna)):
        best_motifs.append(Dna[i][0:k])
    motifs = []
    for j in range(len(Dna[0])-k):
        motifs = (Dna[0][j:j+k]).split()
        profile = []
        for l in range(t-1):
            profile = create_profile_Laplace(motifs)
            new_motif = profile_most_probable(Dna[l+1],k,profile)
            motifs.extend(new_motif)
        if score(motifs) < score(best_motifs):
            best_motifs = motifs
    return '\n'.join(best_motifs)
    
#finds list of t k-length substrings of t strings
#minimizes mismatches
#outputs t k-length strings (list made into \n strings)
def greedy_motif_search(k, t, Dna):
    best_motifs = []
    for i in range(len(Dna)):
        best_motifs.append(Dna[i][0:k])
    motifs = []
    for j in range(len(Dna[0])-k):
        motifs = (Dna[0][j:j+k]).split()
        profile = []
        for l in range(t-1):
            profile = create_profile(motifs)
            new_motif = profile_most_probable(Dna[l+1],k,profile)
            motifs.extend(new_motif)
        if score(motifs) < score(best_motifs):
            best_motifs = motifs
    return '\n'.join(best_motifs)

#scores the number of mismatches in an alignment
#output is an integer
def score(motifs):
    dna_lib = {'A': 0,'C':1,'G':2,'T':3}
    score = 0
    for n in range(len(motifs[0])):
        nucleotide_set = [0,0,0,0]
        for l in range(len(motifs)):
            nucleotide_set[dna_lib[motifs[l][n]]] += 1
        this_score = len(motifs)-max(nucleotide_set)
        score = score + this_score
    return score
        
#based on a list of string of equal length
#calculates probability profiles using Laplace's rule of succession
#returns 4xlen(strings) matrix of probabilities         
def create_profile_Laplace(motifs):
    profile = []
    nucleotide_set = [1,1,1,1]
    dna_lib = {'A': 0,'C':1,'G':2,'T':3}
    for j in range(len(motifs[0])):
        nucleotide_set = [1,1,1,1]
        profile_part = []
        for i in range(len(motifs)):
            nucleotide_set[dna_lib[motifs[i][j]]] += 1
        for k in range(4):
            profile_part.append(nucleotide_set[k]/(len(motifs)+4))
        profile.append(profile_part)
    return profile
    
#based on a list of string of equal length
#calculates probability profiles
#returns 4xlen(strings) matrix of probabilities         
def create_profile(motifs):
    profile = []
    nucleotide_set = [0,0,0,0]
    dna_lib = {'A': 0,'C':1,'G':2,'T':3}
    for j in range(len(motifs[0])):
        nucleotide_set = [0,0,0,0]
        profile_part = []
        for i in range(len(motifs)):
            nucleotide_set[dna_lib[motifs[i][j]]] += 1
        for k in range(4):
            profile_part.append(nucleotide_set[k]/len(motifs))
        profile.append(profile_part)
    return profile

#find the best k-length match along a string 'text'
#given probability profile
#output is a list with the first best match
def profile_most_probable(text, k, profile):
    dna_lib = {'A': 0,'C':1,'G':2,'T':3}
    probabilities = []
    for i in range(len(text)-k+1):
        probability = 1
        counter = 0
        for base in text[i:i+k]:
            base_number = dna_lib[base]
            probability = probability * profile[counter][base_number]
            counter += 1
        probabilities.append(probability)
    if len(probabilities) > 0: max_prob = max(probabilities)
    else: return text
    answer = []
    for j in range(len(probabilities)):
        if probabilities[j] == max_prob: #answer.append(text[j:j+k])
            return text[j:j+k].split()
            
with open('CH02_data.txt') as f:
    firstline = f.readline().split(' ')
    a=int(firstline[0])
    b=int(firstline[1])
    d=[]
    restlines = f.readlines()
    for i in range(len(restlines)):
        d.append(restlines[i].strip('\n'))

print(greedy_motif_search_Laplace(a,b,d))
#print(greedy_motif_search(a,b,d))
#print(score(['AAAT','AAAA']))
#print(create_profile(['AAAAA','AAAAT']))
#print(profile_most_probable('ATTCAACAGACCCATCTGATCCCCTGGACAGCTGGCGCACTGGACCGATCCTAGTTTCTCCTCATAAGTGCCACTATTCTCAAACGGAAATAAGTCTTGCTGGGTAAGACCCGCGAGCCCTGTCTAGTTCGTTCCGCGGTTGGCCCACCGTCGCAGTAGTCCCCCCGGTGCTAAAACATGCGAGCTACATCTCTCAAAGCTTGGCATCATGAACACGTGCCTCTCCCTTTTTGACCGAAGTCTCTGCCGTGGTCTAACATGAATTGAGAGCAATATGAACCCACTTTTAGCGGCGTCTAGCTGAGGCGATGCACCATTAGCTCCACGGTGACCCCTCACCTGGGCCACGACTGCTATCCGTCCTTATCTATTCGGTCGCCGTGGAACTATGGACGGTGGGAGGGACCACCCGAATAAAAACGATCAAGCCATTCTCCCTGAGACAACTACTGGTTCCCAGCCCATCCTCAGTTTCCGGTTGGTCTAGCCTAGATGAATGTATGCCGCGGGATCTTTAACGATGCCGTGTGACTAATCTTAATGACCCCGGTATGGCAGTCGTTTGGTAACAACCTTATTGATAGGAGAGATTACGCGGTCTGCCTTGACTTTGATGAGGCTCAGCGATGGTCTGGGCCCGATGGTCGGCGGGTATTGCGGCCTGCAGCCAGGGGTTGGGCCCGAAGGCACTGTCGAGTCGCTAAAAAATGTGCCAAGCTATCTACTTTTCATAATGTGGGTACCCCTGACGGTCGATATCTCACTCAAACCAGCGTTGGTGATTGGCGCTGCGAACGCGGTATAATAGAATACCGTCCGACTTTCAGTTTTTGATGTGAATGTCGGGTCTACTTTAACACCAATAATGCTCGCGGCTGAACGCTAGCAGCTGTCACGAATCCGCTAGGGTATACTGCCCGCGAGCGGCTCTTAAACTTTTTGCCGGAGGCTCGCGCTACTACAGCCGTGTGTGGACTAGACGGCCTTCCTAGTTGCTCAT',14,
#[[0.338, 0.211, 0.324, 0.225, 0.296, 0.225, 0.254, 0.225, 0.254, 0.38, 0.211, 0.211, 0.225, 0.155],
#[0.169, 0.254, 0.268, 0.211, 0.239, 0.38, 0.197, 0.225, 0.225, 0.183, 0.211, 0.324, 0.31, 0.254],
#[0.282, 0.31, 0.211, 0.296, 0.197, 0.155, 0.31, 0.225, 0.239, 0.268, 0.268, 0.239, 0.211, 0.324],
#[0.211, 0.225, 0.197, 0.268, 0.268, 0.239, 0.239, 0.324, 0.282, 0.169, 0.31, 0.225, 0.254, 0.268]]))