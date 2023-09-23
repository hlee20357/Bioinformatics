#countwithpseudocunts defintiion 
# Input:  A set of kmers Motifs
# Output: CountWithPseudocounts(Motifs). Returns the count matrix of Motifs with pseudocounts as a dictionary of lists. 
def CountWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    Countwithpseudo = {} 
    #establishes empty dictionary 
    for character in 'ACGT': 
        Countwithpseudo[character] = []
         #makes the dictionary with A, C, G, T as the key values 
        for i in range(k): 
            Countwithpseudo[character].append(1)
            #prints out 1's equal to the length of the Motif. This helps establish the scoring matrix because it starts with the pseudocount already.
            #Adds 1 to each column for A,G, C and T. 
            #Essentially, by starting out with 1, it does the Laplace pseudocount already. 
    for i in range(t): 
    #iterating through all motifs eventually, one by one. 
        for j in range(k): 
        #iterating through the length of one motif at a time 
            nucleotide = Motifs[i][j]
            #store nucleotide based off of a Motif's certain position within the matrix. 
            # For example, Motif[1][2] would access the second string's 3rd position. 
            Countwithpseudo[nucleotide][j] +=1
            #Adds 1 to the matrix for a given character and position. 
            # The [symbol] aspect tells you which key to access, and the [j] aspect tells you which index to change within the key. 
            # Another way to look at it is that the [symbol] tells you the row to access (aka access the row meant for the particualr nucleotide) and the [j] tells you what column to access along the row. 
            
    return Countwithpseudo 

#Profile with Pseudocounts definition 
#Outputs the profile of a list of strings (called Motifs). It uses Laplace's pseudocounts to make the profile. 
def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = CountWithPseudocounts(Motifs) 
    # output variable

    for symbol in 'ACGT': 
    #access the symbol within the dictionary 
        for j in range(k): 
        #access the specific position within the nucleotide 
            profile[symbol][j] = (profile[symbol][j])/(t + 4)
            #The adjustment here is adding 4 to t because +1 is added to all nucleotide counts, and since there are 4 nucleotides, 4 is added. 

    return profile


#Motifs definition 
#Input: Takes a profile matrix PROFILE corresponding to a list of strings Dna as input and returns a list of Profile-most probable k-mers in each string from Dna.
#Note: Need Pr(text, profile) to calculate probabilities and Profilemostprobablekmer (text, k, profile) definitions 
def Motifs(Profile, Dna):
    t = len(Dna) 
    k = 4
    #test code had a k-mer of 4
    Motifs = []
    for i in range(t): 
        text = Dna[i]
        most_probable_kmer_in_string = ProfileMostProbableKmer(text,k,Profile)
        Motifs.append(most_probable_kmer_in_string)
    return Motifs
    # insert your code here

# Insert your ProfileMostProbablePattern(Text, k, Profile) and Pr(Pattern, Profile) functions here.

def ProfileMostProbableKmer(text, k, profile):
    prob1 = Pr(text[0:k],profile)
    #initialize prob1 so that you can use as a comparison 
    x = len(text) 
    most_probable_k_mer = text[0:k]
    for z in range(1,x-k+1): 
    #iterate through the whole length of text 
    #Need to add the +1 to account for the range function being exclusive. Also, you can only search up to the x-k position.  
        prob2 = Pr(text[z:z+k], profile) 
        #calculate prob2. Basically slides a window down text, trying out every k-mer other than the first one to compare it to the first k-mers probability. 
        if prob1 < prob2:
        #compare prob1 to prob 2. If prob1 is less, replace prob 1 with prob 2. If prob 1 is more than prob2, prob 1 stays the same. 
            prob1 = prob2
            most_probable_k_mer = text[z:z+k]
    
    return most_probable_k_mer

def Pr(Text, Profile):
    prob = 1
    x = len(Text) 
    for i in range(x): 
    #goes through every index of the Text 
        nucleotide = Text[i]
        #obtains the nucleotide for a specific index of Text 
        prob = prob * Profile[nucleotide][i]
        #Profile[nucleotide][i] accesses the float number within the profile. It first accesses which row based off of the [nucleotide] and which column based off of the [i]
        #keeps multiplying prob unti the length of Text is satisfied 
    return prob



# RandomMotifs definition
#Uses random.randint to choose a random k-mer from each of t different strings Dna, and returns a list of t strings. 

import random 
# Input:  A list of strings Dna, and integers k and t
# Output: RandomMotifs(Dna, k, t). A list of randomly selected k-mers from each string. 
# HINT:   You might not actually need to use t since t = len(Dna), but you may find it convenient
def RandomMotifs(Dna, k, t):
    k = 3
    #Example used to test code had a k-mer of length 3. 
    t = len(Dna) 
    list = []
    for i in range(t): 
        string_access = Dna[i]
        x = random.randint(0,7)
        random_kmer = string_access[x:x+3]
        list.append(random_kmer)
    return list

def Score(Motifs):
    consensus_string = Consensus(Motifs) 
    #stores output of definition Consensus(Motifs) as variable consensus_string
    count_matrix = CountWithPseudocounts(Motifs)
    #stores output of definition Count(Motifs) as variable count_matrix 
    k = len(consensus_string)
    #getting the nucleotide length of motifs
    t= len(Motifs)
    #getting the number of motifs that are present
    score = 0
    #initialize 
    for x in range(k): 
    #for loop to go through the length of the consensus string
        nucleotide = consensus_string[x] 
        #getting the nucleotide from the x index of consensus string 
        value = count_matrix[nucleotide][x]
         #getting the value associated with the nucleotide of the consensus string 
        no_match_total = t-(value)
        #subtracting the (value associated with the nucleotide of the consensus string) from the max value that can be in a column. In this case, the max value is the number of motifs that are present. 
        score = score + no_match_total
        #add to the score. 
    return score


def Consensus(Motifs):
    k = len(Motifs[0]) 
    #sets k equal to the length of all motifs 
    count = CountWithPseudocounts(Motifs) 
    #pulls up the output of Count(motifs) and stores it as a variable called count
    consensus = "" 
    #open string. Will concatenate frequent Symbol with time 
    for j in range(k): 
    #outer loop so this will change the column that the analysis will be done on. Goes index 0, 1, 2, etc. 
        m = 0 
        #initalize for the start of going through a column. Will always go back to this when switching from column to column 
        frequentSymbol = "" 
        #Initalize for the column. Will always go back to this when switching from column to column 
        for symbol in "ACGT": 
        #For a specific column, it looks at all the nucleotides 
            if count[symbol][j] > m:
                m = count[symbol][j] 
                #For a given column, the first round will be A's count. Then, the second round will compare the count of A (stored as m) to the count of C. 
                frequentSymbol = symbol 
                #Makes variable frequentSymbol as the symbol A, C, G, or T if the specific character satisified the if statement. 
        consensus += frequentSymbol 
        #once the for loop is done for ACGT for one column, the frequent symbol is added to the consensus string 
    return (consensus)


# Input:  Positive integers k and t, followed by a list of strings Dna
# Output: RandomizedMotifSearch(Dna, k, t)
def RandomizedMotifSearch(Dna, k, t):
    M = RandomMotifs(Dna, k, t)
    #stores the output of defintiion random motifs into variable M 
    BestMotifs = M
    #Initially makes the Best Motifs equal to the random motifs. Once the whlie loop is accessed, it never goes back to this. 
    while True:
        Profile = ProfileWithPseudocounts(M)
        #Makes the profile from the input M. In the first round, it is the RandomMotifs(Dna,k,t). However, when you go through the while loop the 2nd time, you will make the profile from M, which is making the Profile from (Motifs(Profile)) 
        # Once this line is completed in line to, Profile(Motifs(Profile)) has been completed. 
        M = Motifs(Profile, Dna)
        #Stores the motifs made from the profile and the list of strings in Dna in variable M. 
        # Reupdates what variable M is for in the first round of the while loop, which would be now M = Motifs(Profile). 
        # In the second while loop iteration, it will be Motifs(Profile(Motifs((profile))))
        if Score(M) < Score(BestMotifs):
        #if the score of the motifs made from the profile does not the score of Best Motifs, then M is updated to be the best motif. We want the lowest score possible. 
            BestMotifs = M
        else:
            return BestMotifs 
            # returns the BestMotifs that exist. 
            # Will stop the infinite while loop. 
            # This only occcurs if the profile made does not exceed the current BestMotifs 
