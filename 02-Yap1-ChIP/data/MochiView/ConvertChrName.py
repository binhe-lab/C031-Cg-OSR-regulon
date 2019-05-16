# This file converts the C_glabrata_CBS138_features_mochiview.txt file into a different format for the
# chromosome names to use with Mochiview in the output file C_glabrata_CBS138_features_mochiview.txt

alphabet = dict()
for i in range(1,27):
    letter = chr(ord('A') + i -1)
    alphabet[i] = letter

keep = []
with open('C_glabrata_CBS138_current_feature_mochi.txt') as file:
    for line in file.readlines():
        if line.startswith('Chr'):
            parts = line.split('\t')
            if len(parts[0]) == 4:
                num = int(line[3])
                chr = 'chr0' + alphabet[num]
                parts[0] = chr

            if len(parts[0]) == 5:
                num = int(line[3:5])
                chr = 'chr0' + alphabet[num]
                parts[0] = chr

            keep.append('\t'.join(parts))
        elif line.startswith('mito'):
            chr = 'chrCaglfM'
            parts = line.split('\t')
            parts[0] = chr
            keep.append('\t'.join(parts))
        else:
            keep.append(line)

    file.close()

with open('C_glabrata_CBS138_features_mochiview.txt', 'w') as file:
    new = ''.join(keep)
    file.write(new)
