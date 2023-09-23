
#Definition Symbol Array 
#Input is the string Genome and the symbol (aka a nucleotide) 
#Output is an array that calculates how many nucleotides within a window that is half of the length of the string Genome 
def SymbolArray(Genome, symbol):
    array = {} 
    #establishes empty array 
    n = len(Genome) 
    #n = length of Genome 
    ExtendedGenome = Genome + Genome[0:n//2] 
    #this extendeds the genome to account for a potential window that includes both of the begninning and end of the DNA (if the DNA was linear)
    for i in range(n): 
    #range (n) goes up to (length of Genome) -1 
        array[i] = PatternCount(ExtendedGenome[i:i+(n//2)], symbol) 
        #used to set the array index starting at 0 to a certain value. The ExtendedGenome variable is put into the "Text" argument for the PatternCount, and the symbol variable is put into the "Pattern" argument for the PatternCount
        #print (i)
        #print (array[i])
        #print ('Empty space')
        #print statements above used to determine flow of the code 
        
    return array
    
def PatternCount(Text, Pattern): 
    count = 0 
    #print (range(len(Text)-len(Pattern)+1))
    for i in range(len(Text)-len(Pattern)+1):
        #print ('Coming from PatternCount')
        #print (i)
        #print ('Separation coming from PatternCount2')
        #print statements above used to determine flow of the code 
        if Text[i:i+len(Pattern)] == Pattern: 
            count = count+1
            #print ('Is this accessed?')
            #print ('count: ', count)
            #print statements above used to determine flow of code. 
    return count

SymbolArray('AAAAGGGG', 'A')
#Line above to be used to test the code. 

#For SymbolArray('AAAAGGGG', 'A'), n = 8. Extended Genome states to make a new string that consists of parts 'AAAAGGGG' + 'AAAA' (since brackets are stating from 0 to 4, but the 4 is exclusive so it's actually index positions 0 to 3) 
#i in range(n) would mean i goes from 0 to 7 since range is exclusive of the last number. That means that there will be seven array values with indices from 0 to 7. 
#for the first cycle, the starting array @ 0th index = Pattern Count function. 
#For the patternCount function, it's taking ExtendedGenome[i: i+ n(//2)] and symbol as inputs. 
#in the for loop, the range will be from length of ExtendedGenome[i: i+ n(//2)] - len of symmbol + 1. Since len of symbol =1, then the -1 and +1 will cancel out. 
#len of ExtendedGenome[i:i+(n//2)] is 4. Therefore, the range in the PatternCount definition is 4. It counts from 0 to 3. 
#Therefore, in the first iteration, the for loop in PatternCount will change from 0 to 3 every single time. Now to the if loop. 
#in the if loop, it's checking every single character in the window established by ExtendedGenome[i:i+(n//2)] is equal to the symbol. If it does, add 1 to the count.  
#For the first iteration of the SymbolArray function, the if loop in PatternCount will check indices 0 to 3. That's why in the first iteration of Symbol Array, array[i] is 4. 
#In the second iteration of SybmolArray, you are checking indices 1 to 4. The for loop in PatternCount ensures that you iterate through indices 1 to 4 (the for loop iterates 4 times, but the range remains the same as to 0 to 4), and the if loop is the code that checks if index 1 is equal to pattern, index 2 is equal to pattern, and so forth in the established window. 

# Symbolarray definition is slow because it iterates through len(Genome) times! 
# Then, to compute PatternCount(ExtendedGenome[i+i(n//2)], symbol) we individually compare symbol against n//2 symbols 
# In combination, that means there’s (n^2)//2 comparisons to execute 

def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]

    # look at the first half of Genome to compute first array value. This will be the first thing that's run. It will run before the for loop. 
    array[0] = PatternCount(Genome[0:n//2], symbol)

    for i in range(1, n):
        # start by setting the current array value equal to the previous array value
        #Range function is exclusive after the commma so it goes to position n-1, which is what we want. If it went to the n, then that's a repeat as if you looked at the window from position 0. Because the range is inclusive for the first number presented, you start by looking at position 1 of ExtendedGenome. 
        array[i] = array[i-1] 
        #establishes that the current array index's value is equal to the curretn array index's value of -1. Therefore array[5] = array [4], etc. This is done first before the if statements because the if statements are their to adjust the index value for the paritcular i. In the case of 'AAAAGGGG' with ExtendedGenome 'AAAAGGGGAAAA', if i = 4, array[4] = array[3] = 1 BEFORE the if loops. 
        #print ('This is', i, 'aka i before the if loops')
        #print ('This is ', array[i], 'for array[i]')
        #print ('This is ', array[i-1], 'for array[i-1]')
        #Print statements above used to see flow of the code. 

        # the current array value can differ from the previous array value by at most 1
        #print (ExtendedGenome[i-1])
        #print (ExtendedGenome[i+(n//2)-1])
        #Print statements above to see flow of code. 
        if ExtendedGenome[i-1] == symbol: 
        #because you shifted one position to the right, you lose your leftmost position. If that leftmost position was your symbol, you have to subtract one because you lost it in the current window. You are down one occurrence with the shifted window.  
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol: 
        #This calculates the end of the window on the right. If this base pair is equal to symbol, that means you gained one in the window. Aka the +1. 
            array[i] = array[i]+1
        #print ('This is ', array[i], 'after the if loops')
        #print ('Loop has finished. Onto the next i')
        #print ()
        #Print statements above to see flow of code 
    return array
    
def PatternCount(Text, Pattern): 
    count = 0 
    print (range(len(Text)-len(Pattern)+1))
    for i in range(len(Text)-len(Pattern)+1):
        #print ('Coming from PatternCount')
        #print (i)
        #print ('Separation coming from PatternCount2')
        #print statements above used to see flow of code. 
        if Text[i:i+len(Pattern)] == Pattern: 
            count = count+1
            #print ('Is this accessed?')
            #print ('count: ', count)
            #print ('Space')
            #print statements above used to see flow of code. 
    return count

print ((FasterSymbolArray(('AAAAGGGG'), 'A')))
#Print statement above is used to test code. 

#SkewArray definition 
#skew(K): #G-#C for the first k nucleotides of a genome 
# Input:  A String Genome
# Output: The skew array of Genome as a list.
def SkewArray(Genome):
    count = 0
    list=[0]
    for char in Genome:
        if char == 'G':
            count = count+1
        if char == 'C':
            count = count-1
        list.append(count)
    return list

#Determines minimum_skew position(s) present in the skew array
def MinimumSkew(Genome): 
    positions = []
    lowest_value = min(SkewArray(Genome)) 
    #You need SkewArray(Genome) as the recipient of the min function because the SkewArray(Genome) function outputs the list that you want to find the min from. 
    n = len(SkewArray(Genome))
    for i in range(n): 
        if (SkewArray(Genome))[i] == lowest_value:
            positions.append(i)
    return positions


#Hamming Distance finds the number of nucleotide differences between two strings p and q. 
def HammingDistance(p, q):
    count = 0 
    n = len(p)
    for i in range(n):
        if p[i] != q[i]:
            count = count + 1
    return count 

#Approximate Pattern Matching Problem Code: 
#Find all approximate occurrences of a pattern in a string. 
#Input: Strings Pattern and Text along with integer d(allowable amount of mismatches between Pattern and Text)
# output is all starting positions where Pattern appears as a substring of Text with at most d mismatches 
def ApproximatePatternMatching(Text, Pattern, d):
    positions = [] 
    # initializing list of positions
    n = len(Text)
    m = len(Pattern)
    for i in range (n-m+1): 
        if HammingDistance(Text[i:i+m],Pattern) <=d: 
            positions.append(i)
    return positions

#ApproximatePatternCount definition
# Input:  Strings Pattern and Text, and an integer d
# Output: The number of times Pattern appears in Text with at most d mismatches
def ApproximatePatternCount(Pattern, Text, d):
    count = 0 
    # initialize count variable
    n = len(Text)
    m = len(Pattern)
    for i in range (n-m+1): 
    #Goes through all the positions in which pattern can fit in the text
        if HammingDistance(Text[i:i+m],Pattern) <=d: 
            count = count+1
    return count

