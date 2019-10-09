# DEFINE KEY PHRASES HERE
key_phrases = ['cell wall', 'hyph']




import _pickle as pickle
import pandas as pd

table = pd.read_csv('ClassifiedTable.txt', sep='\t', header=(0), index_col=0)


clust5 = table[table.loc[:, 'Prediction'] == 5]
clust5 = list(clust5.index)

clust10 = table[table.loc[:, 'Prediction'] == 10]
clust10 = list(clust10.index)

clust11 = table[table.loc[:, 'Prediction'] == 11]
clust11 = list(clust11.index)

pickleIn = open('GOanno.pickle', 'rb')
geneAnno = pickle.load(pickleIn)


counter_list = []
anno_list = []

targets = clust5 + clust10 + clust11

for gene in targets:
    for anno in geneAnno[gene]:
        for phrase in key_phrases:
            if phrase in anno:
                counter_list.append(gene)
                anno_list.append(anno)            

print(str(len(set(counter_list))) + ' genes found in CgYap1 target list related to keywords: ', str(key_phrases).replace("'", "")[1:-1])
print('Genes: ', set(counter_list))
print('Annotations: ',set(anno_list))
