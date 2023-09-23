#Pattern Count 
#Input are two strings text and pattern 
#Output is the number of times pattern appears as a substring in text 
def PatternCount(Text, Pattern): 
    count = 0 
    for i in range(len(Text)-len(Pattern)+1): 
        if Text[i:i+len(Pattern)] == Pattern: 
            count = count+1
    return count
#in this case, I starts at 0 and goes until whatever the end is specified. 
#range is exclusive, meaning, it doesn’t count reach the last number. For example, if the Text is length 10 and Pattern is length 3, then range in this case goes to 8, so I will go to 0 to 7 
#When you do index searching, the left side of the colon is inclusive while the right side of the colon is exclusive. Thus, Text[i:i+ len(Pattern)], when I = 0, is looking at indices 0 to 2 although the index searching goes from 0 to 3. 

#FrequencyMap definition
#Input is string text and integer k 
#Output is a dictionary with keys being all k-mers that exist in the text and the value being how many times that k-mer appears in the text 
def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        if Pattern not in freq: 
            freq[Pattern] = 0
        if Pattern in freq: 
            freq[Pattern] = freq[Pattern]+1
    return freq

#FreuqentWords
#Using FreuqencyMap, it identifies the max value from the FrequencyMap and outputs a list of words that appear the most often in String Text. 
def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text,k)
    m = max(freq.values())
    #gets the max values in the FrequencyMap dictionary
    for key in freq: 
        if freq[key]==m:
            words.append(key)
            #checks if the key has the max value, and if so, appends the key to the words list. 
    return words

#Reverse definition 
#Takes a pattern and reverses the order of the pattern 
def Reverse(Pattern):
    x = len(Pattern)
    #print (x)
    reverse = '' 
    #have to set it here because this is a constant. It establishes the empty string in which the loop will act on 
    for i in range(len(Pattern)): 
    #goes from 0 to 7 for 'ABCDEFGH' because len[Patern] = 8 
        #print (i)
        #print (Pattern[len(Pattern)-i-1]) 
        #have to do -1 because the first [i] variable is 0. 
        #print (type((Pattern[len(Pattern)-i-1]))) 
        #verifies that the output is a string
        reverse = reverse + str((Pattern[len(Pattern)-i-1])) 
        #makes sure that you concatentate a string
        #print (reverse)
    return reverse
#print (Reverse('ABCDEFGH'))
#print statement is test code 


#Complement definition
#Takes a pattern and produces its complement
def Complement(Pattern):
    Comp = ''
    for char in (Pattern): 
    #goes through every character in Pattern argument
        if char == 'A':
            Comp = Comp + 'T'
        if char == 'T':
            Comp = Comp + 'A'
        if char == 'G':
            Comp = Comp + 'C'
        if char == 'C':
            Comp = Comp + 'G'
    return Comp
    
#print (Complement('ATGCATGC')) 
#test sequence

#Reverse Complement definition uses Reverse and Complement definitions
def ReverseComplement(Pattern):
    Pattern = Reverse(Pattern) 
    # reverse all letters in a string
    Pattern = Complement(Pattern) 
    # complement each letter in a string
    return Pattern

#PatternMatching definition
#Takes in strings pattern and genome. Basically, you slide a window down the whole genome in increments of k. Then, if the window matches the pattern, then you state the first position where the pattern was found in the genome. 
def PatternMatching(Pattern,Genome):
    Positions = []
    n = len(Genome)
    k = len(Pattern)
    for i in range (n-k+1): 
        Window = Genome[i:i+k]
        if Window == Pattern: 
            Positions.append(i)
    return Positions
    
#print (PatternMatching('ATAT', 'GATATATGCATATACTT'))
#test code 

