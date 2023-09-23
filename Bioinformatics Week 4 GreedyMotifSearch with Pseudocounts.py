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


#Updated GreedyMotifSearch with Pseudocounts definition
# Input:  A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
# Output: GreedyMotifSearch(Dna, k, t)
#The only difference between this GreedyMotifSearch and GreedyMotifSearch with Pseudocounts code-wise is the usage of ProfilewithPseudocounts definition.

def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    n = len(Dna[0])
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
        #sets BestMotifs equal to the first k-mer from each string in dna. These strings will serve as the best-scoring motifs found thus far. 
        # This helps initialize what is considered BestMotifs. It won't be this at the very end, but it will help start the process of comparison. 
    for i in range(n-k+1):
        Motifs =[]
        #initalize motifs. Resetting every single time because we are just doing the same iterative process over and over again but with a different k-mer every single time. 
        #This i  leads to every k-mer possible being tested along the length of the first string in Dna.
        #This is because the k-mer chosen in Dna[0] sets up the profile that leads to the selection of other k-mers in Dna[1] and beyond. 
        Motifs.append(Dna[0][i:i+k])
        # print (Motifs)
        # print ('Empty')
        #Test code to understand the flow
        for j in range(1,t): 
        #this for loop must be satisified before doing +1 on the "for i in range..." loop 
        #the loop will occur t times since range function is exclusive. 
        #print ('inside for x loop')
        #print (x)
        #print (Motifs[0:x])
        #Motifs.append('new addition')
        #print ('')
        #Test code to see the flow of code 
            P = ProfileWithPseudocounts(Motifs[0:j])
            #determines the profile for a set of Motifs. This profile changes as j grows. 
            # For a given k-mer, it'll make a profile from j motifs. 
            #Makes a profile with pseudocounts
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
            #adds the profile most probable pattern to the Motif list. 
            #So for example, for a given k-mer from the first string: 
            #Round 1: Profile will be done on Motifs [0:1], which is Motifs[0]. The profile will be done on Motifs[0]. Then, consideringt the profile of Motifs[0], Motifs[1] will be added to Motifs, and it will be equal to the Profile-most porbable k-mer in Dna[1], the second string in Dna. j @ 1 is done. Now j is at 2. 
            #Round 2: Now j is at 2. Profile will be done on Motifs[0:2], aka the Profile will be made on Motifs [0] and Motifs[1], where Motifs[1] came from the second index in Dna. Using the profile made from Motifs[0] and Motifs[1], Motifs[2] will be made, and it will be eqeual to the profile-most probable k-mer in Dna[2]. Now Motifs has 3 indices from 0 to 2. j@2 is done. Now move onto j@3
            #The cycle conintues until j is up. 
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
            #This if loop is only accessed after the for loop "for j in range (1,t)" is completed. 
            #After selecting a k-mer from each string in Dna to obtain a collection of strings Motifs, GreedyMotifSearch checks whether Motifs outscores the current best scoring collection of motifs, BestMotifs. If the new set of Motifs beats the current set of BestMotifs in scoring, BestMotifs becomes the new set of discovered motifs. 
            #Otherwise,  you move to the outer for loop. The k-mer changes. And what we say in rounds 1 and 2 and beyond will be repeated. 
#after everyghint is run, the BestMotifs is returned. 
    return (BestMotifs)

def Consensus(Motifs):
    # insert your code here
    k = len(Motifs[0]) 
    #sets k equal to the length of all motifs 
    count = CountWithPseudocounts(Motifs) 
    #pulls up the output of CountWithPseudocounts(motifs) and stores it as a variable called count
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

# Then copy your ProfileMostProbableKmer(Text, k, Profile) and Pr(Text, Profile) functions here.
# Write your ProfileMostProbableKmer() function here.
# The profile matrix assumes that the first row corresponds to A, the second corresponds to C,
# the third corresponds to G, and the fourth corresponds to T.
# You should represent the profile matrix as a dictionary whose keys are 'A', 'C', 'G', and 'T' and whose values are lists of floats

def ProfileMostProbableKmer(text, k, profile):
    prob1 = Pr(text[0:k],profile)
    #initialize prob1 so that you can use as a comparison 
    x = len(text) 
    most_probable_k_mer = text[0:k]
    for z in range(1,x-k+1): 
    #iterate through the whole length of text 
    #Need to add the +1 to account for the range function being exclusive. Also, you can only search up to the x-k position.  
        prob2 = Pr(text[z:z+k], profile) 
    #calculate prob2
        if prob1 < prob2:
#       compare prob1 to prob 2. If prob1 is less, replace prob 1 with prob 2. If prob 1 is more than prob2, prob 1 stays the same. 
            prob1 = prob2
            most_probable_k_mer = text[z:z+k]
    
    return most_probable_k_mer

def Pr(Text, Profile):
    # insert your code here
    prob = 1
    x = len(Text) 
    for i in range(x): 
    #goes through every index of the Text 
        nucleotide = Text[i]
        #obtains the nucleotide for a specific index of Text 
        prob = prob * Profile[nucleotide][i]
        #ProfileWithPseudocounts [nucleotide][i] accesses the float number within the profile. It first accesses which row based off of the [nucleotide] and which column based off of the [i]
        #keeps multiplying prob unti the length of Text is satisfied 
    return prob

def Score(Motifs):
    # Insert code here
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

