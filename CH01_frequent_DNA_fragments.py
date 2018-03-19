#Code to find frequent short DNA fragments in genomes
#in search of DnaA-boxes in bacteria
#applied to the Salmonella enterica genome
#as part of Bioinformatic algorithms course

import matplotlib.pyplot as plt

#find the G:C skew along a DNA string
#output is a list
def skew_along(dna):
    count = 0
    out = [0]
    for t in dna:
        if t == 'C': 
            count = count - 1
        if t == 'G':
            count = count + 1
        out.append(count)
    return out
    
#find the position(s) where the G:C skew is at its minimum
#output is a string with spaces
def min_skew(dna):
    skew = skew_along(dna)
    m = min(skew)
    out = []
    for i in range(len(skew)):
        if skew[i] == m: out.append(i)
    return(' '.join(str(x) for x in out))
    
#calculates the number of nucleotides between two DNA strings
#output is an integer
def Hamming_distance(p,q):
    hd = 0
    for i in range(len(p)):
        if p[i] != q[i]: hd += 1
    return hd
    
#finds the number of times a DNA word and its neighbors occur in a DNA string
#a neighbor deviates max d nucleotides from the DNA word
#output is an integer
def pattern_matching(text, pattern, d):
    #d is the max number of mismatches
    positions = []
    for i in range(len(text) - len(pattern)+1):
        if Hamming_distance(text[i:i+len(pattern)], pattern) <= d:
            positions.append(str(i))
    return(len(positions))
    #return (' '.join(positions))
    
#produces a list of DNA strings that deviate max d nucleotides from input DNA string
#RECURSIVE
#output is a list
def neighbors(pattern,d):
    if d == 0: return pattern
    if len(pattern) == 1: return ['A','C','G','T']
    neighborhood = []
    suffix_neighbors = neighbors(pattern[1:],d)
    for t in suffix_neighbors:
        if Hamming_distance(pattern[1:],t) < d:
            for n in ['A','C','G','T']:
                neighborhood.append(n + t)
        else: neighborhood.append(pattern[0] + t)
    return neighborhood
    
#algorithm to turn a DNA string into a number RECURSIVE
#output is a number
def pattern_to_number(pattern):
    if pattern == '':
        return(0)
    symbol = pattern[-1]
    prefix = pattern[0:-1]
    return(4*pattern_to_number(prefix)+symbol_to_number(symbol))
    
#algorithm to turn a nucleotide into a number <4
#output is an integer
def symbol_to_number(letter):
    letters = {'A':0,'C':1,'G':2,'T':3}
    return letters[letter]

#turns a number <4 into a nucleotide
#returns a single letter string
def number_to_symbol(r):
    letters = {0:'A',1:'C',2:'G',3:'T'}
    return letters[r]

#algorithm to turn number into a DNA string RECURSIVE
#output is a string    
def number_to_pattern(index,k):
    if k == 1:
        return number_to_symbol(index)
    remainder = index%4
    prefix_index = (index - remainder)/4
    symbol = number_to_symbol(remainder)
    prefix_pattern = number_to_pattern(prefix_index,k-1)
    return (prefix_pattern + symbol)
    
#returns reverse complement of a DNA string
#output is a string
def rc(string):
    base_complement = {'A':'T','T':'A','G':'C','C':'G'}
    string_rc = ''
    for i in range(len(string)):
        string_rc = string_rc + base_complement[string[i]]
    return(string_rc[::-1])
    
#algorithm to find frequent k-mer words in DNA with max d mismatches
#output is a string with spaces
def frequent_words_with_mismatches(text, k, d):
    frequent_patterns = []
    close = []
    frequency_array = []
    for i in range(pow(4,k)):
        close.append(0)
        frequency_array.append(0)
    for i in range(len(text)-k+1):
        neighborhood = neighbors(text[i:i+k],d)
        for t in neighborhood:
            index = pattern_to_number(t)
            close[index] = 1
    for i in range(pow(4,k)):
        if close[i] == 1:
            pattern = number_to_pattern(i,k)
            frequency_array[i] = pattern_matching(text,pattern,d)
    maxCount = max(frequency_array)
    for i in range(pow(4,k)):
        if frequency_array[i] == maxCount:
            pattern = number_to_pattern(i,k)
            frequent_patterns.append(pattern)
    return ' '.join(frequent_patterns)
    
#algorithm to find frequent k-mer words in DNA with max d mismatches
#also reverse complement of k-mer is searches for
#output is a string with spaces
def frequent_words_with_mismatches_and_rc(text, k, d):
    frequent_patterns = []
    close = []
    frequency_array = []
    for i in range(pow(4,k)):
        close.append(0)
        frequency_array.append(0)
    for i in range(len(text)-k+1):
        neighborhood = neighbors(text[i:i+k],d)
        #neighborhood.extend(neighbors(rc(text[i:i+k]),d))
        for t in neighborhood:
            index = pattern_to_number(t)
            close[index] = 1
    for i in range(pow(4,k)):
        if close[i] == 1:
            pattern = number_to_pattern(i,k)
            frequency_array[i] = pattern_matching(text,pattern,d) + pattern_matching(text,rc(pattern),d)
    maxCount = max(frequency_array)
    for i in range(pow(4,k)):
        if frequency_array[i] == maxCount:
            pattern = number_to_pattern(i,k)
            frequent_patterns.append(pattern)
    return ' '.join(frequent_patterns)

#faster algorithm to find frequent k-mer words in DNA with max d mismatches
#output is a list of words
def frequent_words_with_mismatches_by_sorting(text, k, d):
    frequent_patterns = []
    neighborhoods = []
    index = []
    count = []
    for i in range(len(text)-k+1):
        neighborhoods.extend(neighbors(text[i:i+k],d))
    #neighborhoods_set = list(set(neighborhoods))
    for i in range(len(neighborhoods)):
        pattern = neighborhoods[i]
        index.append(pattern_to_number(pattern))
        count.append(1)
    index.sort()
    for i in range(len(neighborhoods)-1):
        if index[i] == index[i+1]:
            count[i+1] = count[i] + 1
    max_count = max(count)
    for i in range(len(neighborhoods)):
        if count[i] == max_count:
            pattern = number_to_pattern(index[i],k)
            frequent_patterns.append(pattern)
    return frequent_patterns

#read the genome from text file into a single long string
Senterica_genome = ''
with open('Salmonella_enterica.txt') as f:
    data = f.readlines()[1:]
f.close()
for l in data:
    Senterica_genome = Senterica_genome + l.strip('\n')
print(len(Senterica_genome)) #total genome length

#plot the G:C skew along the genome
y = skew_along(Senterica_genome)
x = []
for i in range(len(y)):
    x.append(i)
plt.scatter(x,y)
plt.show()
print(min_skew(Senterica_genome)) #find the minimum
#answer: 3764856 3764858

#find frequent 9-mers in the 1000 bp region around the G:C minimum
suspect_region = Senterica_genome[(3764856-500):(3764858+500)]
print(frequent_words_with_mismatches_by_sorting(suspect_region,9,1))
#answer: ['TCCAAATAA', 'TGAAAATCA', 'TTATCCACA']

#find the positions of TTATCCACA along the genome
answer = []
for i in range(len(Senterica_genome)-9+1):
    if Senterica_genome[i:i+9] == 'TTATCCACA':
        answer.append(i)
print(answer)
print(len(answer)) #how many times is it in the genome?
#answer: 24

#plot the occurrences of the reverse complement of TTATCCACA along the genome
answer2 = []
for i in range(len(Senterica_genome)-9+1):
    if Senterica_genome[i:i+9] == rc('TTATCCACA'):
        answer2.append(i)
print(answer2)
print(len(answer2))
plt.scatter(answer2,range(len(answer2)))
plt.show()