motif_yap1 = 'TTAC'
motif_skn7_1 = 'G_C__GGCC'
motif_skn7_2 = 'G_C__GCCC'

promoter = 'GCCTTGCCCTTAGATTACCCCGCCCEA'

def skn7_Test(promoter, motif, motif2):
    for i in range(len(promoter)-len(motif)):
        sample = promoter[i:i+len(motif)]
        if ((motif[0]==sample[0]) and (motif[2] == sample[2]) and ((motif[5:]==sample[5:]) or (motif2[5:] == sample[5:]))):
            return True
    return False

def yap1_Test(promoter, motif):
    for i in range(len(promoter)-len(motif)):
        sample = promoter[i:i+len(motif)]
        if motif==sample:
            return True
    return False

print('Yap1 found: ',yap1_Test(promoter,motif_yap1))

print('Skn7 found: ',skn7_Test(promoter,motif_skn7_1,motif_skn7_2))
