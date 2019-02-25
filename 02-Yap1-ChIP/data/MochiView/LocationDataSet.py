# Reorganizes the Kuo et al values found in the SeriesMatrix File to be imported as a Location/Data set in Mochiview
# Series Matrix data can be found at:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15818
# Replicates 1, 2, and 3 can be found under samples GSM397447, GSM397448, and GSM397449, respectively

def removeSpaces(line):
    ret = ''
    words = line.split('\t')
    for word in words:
        if len(word) >0:
            ret += word + '\t'
    return ret.rstrip('\t')

def removeE(line, min, max):
    ret = ''
    newWords = []
    newWord2 = []
    words = line.split('\t')
    for word in words:
        newWords.append(word.replace('\n',''))
        if word.startswith('-') and float(word) < min:
            min = float(word)
        elif '.' in word and float(word) > max:
            max = float(word)

    for index,word in enumerate(newWords):
        if 'E-' in word:
            if word[0] == '-':
                word = '-0.000' + word[1]
            else:
                word = '0.000' + word[0]
            newWord2.append(word)
        else:
            newWord2.append(word)


    for word in newWord2:
        if len(word) >0:
            ret += word + '\t'
    return ret + '\n', min, max

keep = []
with open('SeriesMatrix') as file:
    min = 0
    max = 0
    for line in file.readlines():
        line = line.replace(':', '\t')
        line = line.replace('-', '\t', 1)
        line = line.replace('"', '')

        line = removeSpaces(line)
        line, min, max = removeE(line, min, max)
        keep.append(line)

    print('Minimum value: '+ str(min) + ' ' + 'Maximum Value: '+str(max))
    file.close()

total = ''
for key in keep:
    total += key

total = 'SEQ_NAME\tSTART\tEND\tREP1\tREP2\tREP3\n' + total

with open('LocationDataSet.txt', 'w') as file:
    file.write(total)
    file.close()
