#Probability Definition
# Input:  String Text and matrix Profile (A,C, G, T make up the rows, and the columns are the positions along the string Text)
# Output: Probability of the string to occur based off of the probabilities present in the matrix profile
def Pr(Text, Profile):
    prob = 1
    x = len(Text) 
    for i in range(x): 
    #goes through every index of the Text 
        nucleotide = Text[i]
        #obtains the nucleotide for a specific index of Text 
        prob = prob * Profile[nucleotide][i]
        #Profile[nucleotide][i] accesses the float number within the profile. It first accesses which row based off of the [nucleotide] and which column based off of the [i]
        #keeps multiplying prob until the length of Text is satisfied 
    return prob

#ProfileMostProbableKmer definition
# The profile matrix assumes that the first row corresponds to A, the second corresponds to C, the third corresponds to G, and the fourth corresponds to T.
# You should represent the profile matrix as a dictionary whose keys are 'A', 'C', 'G', and 'T' and whose values are lists of floats
def ProfileMostProbableKmer(text, k, profile):
    prob1 = Pr(text[0:k],profile)
    #initialize prob1 so that you can use as a comparison 
    x = len(text) 
    most_probable_k_mer = text[0:k]
    #initialize the most probable k_mer to act as a way to start the comparison between k-mers. 
    for z in range(1,x-k+1): 
    #iterate through the whole length of text. You start at 1 because you've already completed dealing with index 0 (index 0 was dealt with before the for loop)
    #Need to add the +1 to account for the range function being exclusive. Also, you can only search up to the x-k position.  
        prob2 = Pr(text[z:z+k], profile) 
        #calculate prob2
        if prob1 < prob2:
        #compare prob1 to prob 2. If prob1 is less, replace prob 1 with prob 2. If prob 1 is more than prob2, prob 1 stays the same. 
            prob1 = prob2
            most_probable_k_mer = text[z:z+k]
            #if prob2 has a higher probability than the original prob1, then the most_probable_kmer gets updated to be the substring that won. 
    return most_probable_k_mer

#Count(Motifs) definition
# Takes a list of strings Motifs as input and returns the count matrix of Motifs (as a dictionary of lists). 
def Count(Motifs):
    count = {} 
    # initializing the count dictionary
    #Motifs = ['AACCGGTT']
    #test code 
    k = len(Motifs[0]) 
    #k = length of the 1st string in Motifs. Okay to use since all the motifs are of the same length
    for symbol in "ACGT":
        count[symbol] = [] 
        #establishes key value pairing within the dictionary count. Basically establishes 4 keys of A, C, G, and T in the dictionary, and all of them have open lists as values 
        for i in range(k): 
            count[symbol].append(0) 
            #prints out 0's equal to the length of the Motif. This helps establish the scoring matrix because you can add to the 0's. 
    t = len(Motifs) 
    for i in range(t): 
    #iterating through all the strings in Motifs. Overarching for loop so you'll go through all the strings in Motifs. 
        for j in range(k): 
        #iterating through the length of one motif at a time
        #print (j)
        #print ('ABove is j')
        #test code to see the flow of the code 
            symbol = Motifs[i][j] 
            #accesses motif's certain position as if it was a matrix. Accesses the matrix's row and column position (in that order). For example, Motif[1][2] would access the second string's 3rd position
            #The access obtains the nucleotide at that certain position
        #print (symbol)
        #print ('This is the symbol')
        #Test code to see the flow of the code 
            count[symbol][j] += 1 
            #Adds 1 to the matrix for a given character and position. The [symbol] aspect tells you which key to access, and the [j] aspect tells you which index to change within the key. 
            #Another way to look at it is that the [symbol] tells you the row to access (aka access the row meant for the particualr nucleotide) and the [j] tells you what column to access along the row. 
            #  print (count)
            #print ('Above is count')
            #Test code to see the flow of the code 
    return count

#Profile of Motifs definition 
#Uses the Count definition and makes a profile matrix from the count matrix. 
#In a profile matrix, each intger within the matrix is dividied by the number of strings that make up Motifs. 
def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = Count(Motifs) 
    #Makes the profile the same as the output as Count(Motifs). 
    #That means the profile is a dictionary with keys of A, C, G, and T with values of how many of that nucleotide was at that position
    for symbol in 'ACGT': 
    #To access the symbol within the profile AKA accessing the row 
            for j in range(k): 
            #to access the index along the particular nucleotide AKA accessing the column 
                profile[symbol][j] = profile[symbol][j]/t 
                #This results in a profile matrix Profile(Motifs) for which element (i,j) is the frequency of the i-th nucleotide in the j-th column of the motif matrix (i.e., the number of occurrences of the i-th nucleotide divided by t, the number of nucleotides in the column)
                #You are just updating profile[symbol][j] to represent a profile matrix rather than a count matrix. 

    return profile

#Entropy definition outputs the entropy of the motifs in a list. 
Import math 
#Motifs = ['TCGGGGGTTTTT', 'CCGGTGACTTAC', 'ACGGGGATTTTC', 'TTGGGGACTTTT', 'AAGGGGACTTCC','TTGGGGACTTCC', 'TCGGGGATTCAT','TCGGGGATTCCT', 'TAGGGGAACTAC', 'TCGGGTATAACC']
#Above two lines are used to test definition
def Entropy(Motifs):
    n = len(Motifs[0])
    #sets up your ability to work through each index position 
    probability_matrix = Profile(Motifs)
    total = 0 
    for x in range(n): 
        for symbol in 'ACGT': 
        #this double for loop means you go through each column's ACGT associated values before moving onto the next column 
            prob_value = probability_matrix[symbol][x]      
            #obtain probability distribution value 
            if prob_value == 0: 
                continue
            #need this because log2(0) leads to errors. Thus, if the code ever reaches this situation, you need to move onto the next letter 
            ind_entropy = (prob_value)*(math.log2(prob_value))
            #formula to calculate individual entropy
            total = total + ind_entropy
            #total is the sum of all the individual entropies 
    total = total *(-1)
    #at the end, entropy is defined as the negative of the total sum of individual entropies 
    return total

#Consensus definition
#For a given set of Motifs, the output is the consensus string
#It can do that by determining the max value within a column, and the associated nucleotide within that column becomes part of the consensus string. 
def Consensus(Motifs):
    k = len(Motifs[0]) 
    #sets k equal to the length of all motifs 
    count = Count(Motifs) 
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
                #Makes variable frequentSymbol as the symbol A, C, G, or T if the specific character satisified the if statement (of having the largest value)
        consensus += frequentSymbol 
        #once the for loop is done for ACGT for one column, the frequent symbol is added to the consensus string 
    return (consensus)

#Score definition 
#Since score is defined as the total number of lowercase letters present in Motifs, you can approach the task by calculating the no match total in a column. 
#No match total in a column is calculated by (number of strings in Motifs)-(number of the most present nucleotide in that column of Motifs)
def Score(Motifs):
    consensus_string = Consensus(Motifs) 
    #stores output of definition Consensus(Motifs) as variable consensus_string
    count_matrix = Count(Motifs)
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

#GreedyMotifSearch definiition
# Input:  A list of strings DNA to analzyze, and integers k and t (where t is the number of kmers in Dna)
# Output: GreedyMotifSearch(Dna, k, t)
def GreedyMotifSearch(Dna, k, t):
    n = len(Dna[0])
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
        #sets BestMotifs equal ot the first k-mer from each string in dna. 
        #These strings will serve as the best-scoring motifs found thus far. This helps initialize what is considered BestMotifs. 
        #It won't be this first selected k-mer at the very end, but it will help start the process of comparison. 
    for i in range(n-k+1):
        Motifs =[]
        #initalize motifs. Resetting every single time because we are just doing the same iterative process over and over again but with a different k-mer every single time. 
        #This i  leads to every k-mer possible being tested along the length of the string 
        Motifs.append(Dna[0][i:i+k])
        #    print (Motifs)
        #    print ('Empty')
        #Code above sets up one motif to represent string 0 (the first string present in DNA)
        #This for loop is only accessed only after the [for j in range(1,t: ] is complete. 
        for j in range(1,t): 
        #this for loop must be satisified before doing +1 on the "for i in range..." loop 
        #the loop will occur t times since range function is exclusive. 
        #print ('inside for x loop')
        #print (x)
        #print (Motifs[0:x])
        #Motifs.append('new addition')
        #print ('')
            P = Profile(Motifs[0:j])
            #determines the profile for a set of Motifs. This profile changes as j grows. For a given k-mer, it'll make a profile from j motifs. 
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
            #adds the profile most probable pattern to the Motif list. 
            #Round 1: Profile will be done on Motifs [0:1], which is Motifs[0]. The profile will be done on Motifs[0]. Then, consideringt the profile of Motifs[0], Motifs[1] will be added to Motifs, and it will be equal to the Profile-most porbable k-mer in Dna[1], the second string in Dna. j @ 1 is done. Now j is at 2. 
            #Round 2: Now j is at 2. Profile will be done on Motifs[0:2], aka the Profile will be made on Motifs [0] and Motifs[1], where Motifs[1] came from the second index in Dna. Using the profile made from Motifs[0] and Motifs[1], Motifs[2] will be made, and it will be eqeual to the profile-most probable k-mer in Dna[2]. Now Motifs has 3 indices from 0 to 2. j@2 is done. Now move onto j@3
            #The cycle conintues until j is up. 
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
            #This if loop is only accessed after the for loop "for j in range (1,t)" is completed. 
            #After selecting a k-mer from each string in Dna to obtain a collection of strings Motifs, GreedyMotifSearch checks whether Motifs outscores the current best scoring collection of motifs, BestMotifs. If the new set of Motifs beats the current set of BestMotifs in scoring, BestMotifs becomes the new set of discovered motifs. 
            #Otherwise,  you move to the outer for loop. The k-mer changes. And what we say in rounds 1 and 2 and beyond will be repeated. 

    #after everything is run, the BestMotifs is returned. 
    return (BestMotifs)

#Because the GreedyMotifDefinition has a lot of flow going on, I made an example code to see how it flows 
#for z in range (6): 
#    print ('hello')
#for i in range (5): 
#    print ('hi')
#    for x in range (3): 
#        print ('sigh')
#    if 2 ==2: 
#       print ('yes)
#In this example, 6 'hello's are printed. Then a 'hi' then 3 'sigh' then a 'yes' are printed. then another 'hi' then 3 'sigh' then a 'yes' are printed. Etc. 
#Inner loops must be satisfied before changing the iterating number on the outer for loop. If loop is only accessed after the for loop preceding it is finished 