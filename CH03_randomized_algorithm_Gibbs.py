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
            
#find the most probable set of k-mer strings 
#from longer set of strings Dna, given profile probabilities
#output is a list of k-mer strings
def motifs_from_profile(profile, Dna):
    k = len(profile)
    best_motifs = []
    for i in range(len(Dna)):
        best_string = profile_most_probable(Dna[i], k, profile)
        best_motifs.extend(best_string)
    return best_motifs

#search for a better list of k-mers in a list of Dna strings
#by randomized search - need to run this many times for a good fit!
#output is a list of k-mers and the score of this motif
def randomized_motif_search(Dna, k, t):
    motifs = []
    for i in range(len(Dna)):
        position = random.randint(0,(len(Dna[0])-k))
        motifs.append(Dna[i][position:(position+k)])
    best_motifs = motifs
    while True:
        profile = create_profile_Laplace(motifs)
        motifs = motifs_from_profile(profile, Dna)
        if score(motifs) < score(best_motifs):
            best_motifs = motifs
        else: return (best_motifs,score(best_motifs))

#search for a better list of k-mers in a list of Dna strings
#by randomized Gibbs search - need to run this many times for a good fit!
#output is a list of k-mers and the score of this motif
def Gibbs_sampler(Dna, k, t, N):
    motifs = []
    for i in range(len(Dna)):
        position = random.randint(0,(len(Dna[0])-k))
        motifs.append(Dna[i][position:(position+k)])
    best_motifs = motifs
    for j in range(N):
        a = random.randint(0,t-1)
        del motifs[a]
        profile = create_profile_Laplace(motifs)
        motif_a = profile_random_Gibbs(Dna[a], k, profile)
        motifs.insert(a,motif_a)
        if score(motifs) < score(best_motifs):
            best_motifs = motifs
    return (best_motifs,score(best_motifs))

#choose a random k-mer from a string
#output is the good k-mer
def profile_random_Gibbs(text, k, profile):
    pees = []
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
    choice = random_Gibbs(probabilities)
    return text[choice:choice+k]
    
#input is a list of p2s that do not have to sum to one
#inside is a list of p's that do sum to one
#output is a random integer between 0 and N-1
def random_Gibbs(pees):
    pees2 = []
    for i in range(len(pees)):
        pees2.append(pees[i]/sum(pees))
    die = random.uniform(0,1)
    threshold = 0
    for j in range(len(pees2)):
        threshold = threshold + pees2[j]
        if die < threshold: return j
        
#run the randomized_motif_search algorithm multiple times
#and save the output to two lists
#then choose the best (one or several) motifs based on the scores
import random
import time
t0 = time.time()
with open('CH03_data_1.txt') as f:
    firstline = f.readline().split(' ')
    k=int(firstline[0])
    t=int(firstline[1])
    d=[]
    restlines = f.readlines()
    for i in range(len(restlines)):
        d.append(restlines[i].strip('\n'))
tests_motifs = []
tests_scores = []
for i in range(2):
    test = randomized_motif_search(d,k,t)
    tests_motifs.append(test[0])
    tests_scores.append(test[1])
for j in range(len(tests_scores)):
    if tests_scores[j] == min(tests_scores):
        print(min(tests_scores))
        print('\n'.join(tests_motifs[j]),'\n')
t1 = time.time()
print('time: ',t1-t0)

#run the profile_random_Gibbs algorithm multiple times
#and save the output to two lists
#then choose the best (one or several) motifs based on the scores
t0 = time.time()
with open('CH03_data_2.txt') as f:
    firstline = f.readline().split(' ')
    k=int(firstline[0])
    t=int(firstline[1])
    N=int(firstline[2])
    d=[]
    restlines = f.readlines()
    for i in range(len(restlines)):
        d.append(restlines[i].strip('\n'))
tests_motifs = []
tests_scores = []
for i in range(2000):
    test = Gibbs_sampler(d,k,t,N)
    tests_motifs.append(test[0])
    tests_scores.append(test[1])
for j in range(len(tests_scores)):
    if tests_scores[j] == min(tests_scores):
        print(tests_scores[j])
        print('\n'.join(tests_motifs[j]),'\n')
t1 = time.time()
print('time: ',t1-t0)
