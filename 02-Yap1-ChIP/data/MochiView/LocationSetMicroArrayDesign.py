# This file modifies the MicroArray design values found at the link below to a format which can be imported as a Location
# set in MochiView v1.46

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc=GPL8477&id=18401&db=GeoDb_blob33

keep = ''

with open('MicroArrayDesignRAW.txt') as file:
    for line in file.readlines():
        if not line.startswith('chr0'):
            continue
        line = line.split('\t')
        line[0] = line[0].replace(':','\t')
        line[0] = line[0].replace('-','\t')
        start = line[0].split('\t')[1]
        end = line[0].split('\t')[2]
        line[0] = line[0] + '\t+\t'
        keep = keep + '\n' + line[0] + line[1]+'\t'+start+'\t'+end+'\t'+'1'+'\t'+start+'\t'+end
    keep = keep.replace('\n','',1)
    keep = 'SEQ_NAME	START	END	STRAND	FEATURE_NAME	TXN_START	TXN_END	EXON_COUNT	EXON_STARTS	EXON_ENDS'+'\n' + keep
    print(keep)
    file.close()
    
with open('locationSetMicroArrayDesign.txt','w') as file:
    file.write(keep)
    file.close()
