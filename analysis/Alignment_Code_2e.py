#Graphing functions
import matplotlib.pyplot as plt
import numpy as np

#Gloabal lists used throughout multiple functions
means_yap1 = []
means_skn7 = []

#Function to initialize the geneNames and sequences lists taken from 
# pathoyeastract.org and set up the motifs
def buildLists():
    #Select Data file
    f = open('../Data/genes_and_promoter.txt', 'r')
    
    #Lists of gene names and their promoter sequences
    geneNames = []
    sequences = []     

    #Declare a variable to store the sequences
    promoter=''
    for line in f.readlines():
        #At every new gene, begin a new sequence for the promoter and store the names
        if line[0] == '>':
            line = line[1:]
            sentence = line.split()
            geneNames.append(sentence[0])
            sequences.append(promoter)
            promoter = ''
        #If the line doesn't start with '>' then add it to the currently building promoter
        else: 
            promoter += line.strip()
    #For last gene
    sequences.append(promoter)
    
    #Remove first entry which is blank to synchronize the indicies of the lists
    del sequences[0]
    
    master = [geneNames, sequences]
    return master
    
#Check if skn7 motif is present in a given promoter
def skn7_Test(promoter):
    motif = 'G_C__GGCC'
    motif2 = 'G_C__GCCC'
    for i in range(len(promoter)-len(motif)+1):
        sample = promoter[i:i+len(motif)].upper()
        if ((motif[0]==sample[0]) and (motif[2] == sample[2]) and 
            ((motif[5:8]==sample[5:8]) or (motif2[5:8] == sample[5:8]))):
            return True
    return False

#Check if yap1 motif is present in a given promoter
def yap1_Test(promoter):
    motif = 'TTAC'
    for i in range(len(promoter)-len(motif)):
        sample = promoter[i:i+len(motif)].upper()
        if motif==sample:
            return True
    return False

#Build global means lists with percentages for given motifs
def run():
    master = buildLists()
    geneNames = master[0]
    sequences = master[1]
    
    count1 = 0
    #Find probabilities of group 1
    for i in range(12):
        if (yap1_Test(sequences[i]) == True):
            count1 = count1+1
    means_yap1.append((count1/12)*100)
    count1 = 0
    for i in range(12):
        if (skn7_Test(sequences[i]) == True):
            count1 = count1+1
    means_skn7.append((count1/12)*100)    
    
    count2=0
    #Find probabilities of group 2
    for i in range(18):
        if (yap1_Test(sequences[i+12]) == True):
            count2 = count2+1
    means_yap1.append((count2/18)*100)
    count2 = 0
    for i in range(18):
        if (skn7_Test(sequences[i+12]) == True):
            count2 = count2+1
    means_skn7.append((count2/18)*100)       
    
    count3=0
    #Find probabilities of group 3
    for i in range(8):
        if (yap1_Test(sequences[i+30]) == True):
            count3 = count3+1
    means_yap1.append((count3/8)*100)
    count3 = 0
    for i in range(8):
        if (skn7_Test(sequences[i+30]) == True):
            count3 = count3+1
    means_skn7.append((count3/8)*100)

#Plot final means lists
def plot():
    run()
    n_groups = 3
    fig, ax = plt.subplots()
    index = np.arange(n_groups)
    bar_width = 0.35
    opacity = 0.8
     
    rects1 = plt.bar(index, means_yap1, bar_width,
                     alpha=opacity,
                     color='black',
                     label='Yap1 Consensus within 1000bp upstream')
     
    rects2 = plt.bar(index + bar_width, means_skn7, bar_width,
                     alpha=opacity,
                     color='white',
                     label='Skn7 Consensus within 1000bp upstream',
                     ec = 'black')
     
    plt.ylabel('Percent of the group')
    plt.title('Identified Consensus Sites')
    plt.xticks(index + bar_width, ('1', '2', '3'))
    plt.legend(loc = 'lower left')
     
    plt.tight_layout()
    plt.show()    

plot()

