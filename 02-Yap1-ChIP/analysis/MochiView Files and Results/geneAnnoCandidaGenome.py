# File modified from:
# http://candidagenome.org/download/gff/C_glabrata_CBS138/archive/C_glabrata_CBS138_version_s01-m01-r01_features.gff

# This modification allows the .gff file returned to be compatible as a location set in MochiView

import re


keep= ''

with open ('geneAnnoCandidaGenome.txt') as file:
    for line in file.readlines():
        if re.search('Chr._C_glabrata_CBS138',line):
            replace = 'c' + line[1:3] + '0' + line[3]
            line = re.sub('Chr._C_glabrata_CBS138',replace, line)
            keep = keep + line
        elif re.search('mito_C_glabrata_CBS138',line):
            replace = 'chrCaglfM'
            line = re.sub('mito_C_glabrata_CBS138',replace, line)
            keep = keep + line
        else:
            keep = keep + line
    file.close()

with open('annoCandidaGenome.gff', 'w') as f:
    f.write(keep)
    f.close()