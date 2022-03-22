#!/usr/bin/env python
# coding: utf-8
import os

def load_Alignments(file):
    with open(file, 'r') as fr:
        lines = fr.readlines()
        synteny = []
        for line in lines:
            if line.startswith('## Alignment'):
                try:
                    synteny.append([Alignment, block, block_evalue])
                except:
                    print(f"{file} starts line {line.strip()} !")
                ele = line.split('=')
                align_id = ele[0].split(' ')[2]
                score = float(ele[1].split(' ')[0])
                evalue = float(ele[2].split(' ')[0])
                N = int(ele[3].split(' ')[0])
                position = ele[3].split(' ')[1]  + ' ' + ele[3].split(' ')[2]
                Alignment = [score, evalue, N, position, align_id]
                block = []
                block_evalue = []
            elif not line.startswith('#') and N > 0:
                ele = line.split('\t')
                block.append(ele[1] + '@' + ele[2])
                block_evalue.append(ele[-1].strip())
                N -= 1
        try:
            synteny.append([Alignment, block, block_evalue])
        except:
            print(f'Warning! {file} maybe a empty file!!!!!!!!!!!!!!!1')
            
    return synteny


def compare_list(list1, list2):
    overlap = len(set(list1).intersection(set(list2)))
    return overlap*2


def choose_better_one(Alignment1, Alignment2):
    evalue1, evalue2 = Alignment1[1], Alignment2[1]
    if evalue1 < evalue2:
        return 1    
    else :
        return 2


def generate_blocks(synteny, n):
    n = int(n)
    a = synteny[0]
    e = synteny[2]
    Alignment = f'## Alignment {n}: score={str(a[0])} e_value={str(a[1])} N={str(a[2])} {a[3]}'
    block, i = '', 0
    for b in synteny[1]:
        pair = b.replace('@', '\t')
        block += f'{n}- {int(i)}:\t{pair}\t{e[i]}\n'
        i += 1
    w = Alignment + block
    return w


def generate_Alignments(file1, file2):
    synteny1, synteny2 = load_Alignments(file1), load_Alignments(file2)
    synteny1_bak, synteny2_bak = synteny1[:], synteny2[:]
    basename = os.path.basename(file1)
    with open('../merged_aligns/' + basename.replace('collinearity', 'aligns'), 'w') as fw:
        with open('../merged_aligns/' + basename.replace('collinearity', 'report'), 'w') as fw1:
            n = 0
            fw1.write(file1 + ' Alignments ' + str(len(synteny1)) + ' ' + file2 + ' Alignments ' + str(len(synteny2)) + '\n')
            for s1 in synteny1:
                Alignment1 = s1[0]
                N1 = Alignment1[2]
                block1 = s1[1]
                for s2 in synteny2:
                    N2 = s2[0][2]
                    block2 = s2[1]
                    overlap = compare_list(block1, block2)
                    sum_pairs = N1 + N2
                    if overlap / sum_pairs > 0.6:
                        Alignment2 = s2[0]
                        A = choose_better_one(Alignment1, Alignment2)
                        if A == 1:
                            synteny = s1
                        elif A == 2:
                            synteny = s2
                        w = generate_blocks(synteny, n)
                        fw.write(w)
                        n += 1
                        try:
                            synteny1_bak.remove(s1)
                        except:
                            print('file1 one-to-mutil file2 match ' + str(n))
                            fw1.write('file1 one-to-mutil match ' + str(n))
                        try:
                            synteny2_bak.remove(s2)
                        except:
                            print('file1 mutil-to-one file2 match ' + str(n))
                            fw1.write('file1 mutil-to-one file2 match ' + str(n))
                        fw1.write(Alignment1[4] + '\t-----matches-------\t' + Alignment2[4] + '\n')
            fw1.write(file1 + ' Unique Alignments ' + str(len(synteny1_bak)) + file2 + ' Unique Alignments ' + str(len(synteny2_bak)) + '\n')
            fw1.write('sum overlap: ' + str(n) + '\n')
            m = max(len(synteny1_bak), len(synteny2_bak))
            for i in range(m):
                if i < len(synteny1_bak) and i < len(synteny2_bak):
                    w = generate_blocks(synteny1_bak[i], n)
                    fw.write(w)
                    n += 1
                    w = generate_blocks(synteny2_bak[i], n)
                    fw.write(w)
                    n += 1
                elif i >= len(synteny1_bak):
                    w = generate_blocks(synteny2_bak[i], n)
                    fw.write(w)
                    n += 1
                elif i >= len(synteny2_bak):
                    w = generate_blocks(synteny1_bak[i], n)
                    fw.write(w)
                    n += 1
                i += 1
            fw1.write('sum pairs: ' + str(n) + '\n')


generate_Alignments(file1, file2)
