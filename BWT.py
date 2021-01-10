'''Script d'alignement de read sur un génome
'''
import argparse
import copy
from datetime import datetime


def str_return_seq(seq2):
    '''Reurn of the inverse sequence
    '''
    dico_base = { 'A' : 'T',
                  'T' : 'A',
                  'C' : 'G',
                  'G' : 'C'
                  }
    seq_inv = []
    seq=[]
    seq[:0]= seq2
    for base in range(0,len(seq)):
        seq_inv.append(dico_base[seq[base]])
    seq_inv = "".join(seq_inv)
    seq_inv = seq_inv[::-1]
    return seq_inv

def str_compute_bwt(seq):
    '''Generate burrows-wheeler sequence
    '''
    bwt = ""
    mat_bwt = []
    for pos in range(0,len(seq)+1):
        list_tmp = seq[(len(seq)-pos):len(seq)] + "$" + seq[0:(len(seq)-pos)]
        mat_bwt.append(list_tmp)
    mat_bwt.sort()
    for p in range(len(seq)+1):
        bwt += mat_bwt[p][len(seq)]
    return bwt

'''Fonction de génération alternative de la bwt, sans générer la matrice
La fonction repose sur un algorithme de tri.
On part du principe que la bwt se génère à partir de 2 colonnes,
la première colonne de la matrice triée par ordre alphabétique,
et la dernière, qui sont les lettres précédentes de la première.
Donc si on connaît la 1ère colonne, on peut reconstruire la bwt.
Pour cela, l'idée est ici d'associer les positions de la séquence
avec les lettres, et de trier ces positions.
On décide entre 2 positions laquelle est la meilleure de la manière suivante :
- on prend 2 positions, on regarde les lettres associées
  si l'une des deux est plus petite (ex : A < T) on trie, sinon
  on passe à la lettre suivante jusqu'à avoir déterminé quelle portion 
  de séquence était la meilleure. (nb : on fait attention si on arrive à 
  la fin de la séquence)

L'algorithme de tri utilisé est un algorithme de tri par insertion,
on a donc en mémoire n x 2 (séquence + positions associées), et en temps
de calcul n(n-1)/2, soit O(n^2). On sait qu'on peut diminuer le temps
de calcul en utilisant un tri fusion (nlog(n)), mais on choisit ici d'utiliser
une autre façon de diminuer le temps de calcul.
Étant donné que le temps de calcul dépend de la taille de la séquence,
on choisit de diviser au préalable la séquence en 4 puisqu'on sait que les
positions associées à A < positions associées à C etc...
On définit en fait un dictionnaire de k-mer déjà trié de taille 1.
Donc si on génére un dictionnaire de k-mer trié de taille supérieure
avec une fonction, on peut d'autant plus diminuer le temps de 
génération de bwt.
NB : lorsqu'on teste la fonction ci-dessus sur la séquence complète, on 
a environ 10mn de calcul, lorsqu'on passe au dictionnaire de k-mer de
taille 1, on passe à 3mn pour toute la bwt.
'''

def str_compute_bwt_alt(seq):
    '''Generate burrows-wheeler sequence
    '''
    seq = "$" + seq
    vec_pos = range(len(seq))
    dico_base = {"A" : [],
                 "C" : [],
                 "G" : [],
                 "T" : []}
    for i in range(1, len(seq)):
        dico_base[seq[i]].append(i) #génération des positions

    for key in dico_base.keys():
        dico_base[key] = insertion_sort(seq, dico_base[key]) #tri des listes

    bwt = seq[len(seq)-1]
    for key in dico_base.keys():
        for i in dico_base[key]:
            bwt =  bwt + seq[i-1] #génération de la bwt

    return bwt

def insertion_sort(seq, list_base):
    '''Fonction de tri par insertion
    '''
    list_copy = list(list_base)
    N = len(list_copy)
    for n in range(1,N):
        unsort = list_copy[n]
        sort = n-1
        while sort>=0 and best_of_two(seq, list_copy[sort], unsort):
            list_copy[sort+1] = list_copy[sort]
            sort = sort-1
        list_copy[sort+1] = unsort
        print("\rAvancee tri : {}/{}".format(n, N), end = "")

    return list_copy

def best_of_two(seq, pos1, pos2):
    '''Fonction de comparaison de position
    '''
    for i in range(0,len(seq)):
        if pos2+i >= len(seq) or seq[pos1+i] > seq[pos2+i]:
            return True

        elif seq[pos1+i] < seq[pos2+i]:
            return False

    return

def dict_col_first(seq):
    '''Générer dictionnaire correspondant à la colonne 1 de la matrice bwt
    '''
    dict_letters = {}
    seq = list(seq)
    seq.sort()
    seq = "".join(seq)

    for letter in range(len(seq)):
        if seq[letter] in dict_letters:
            dict_letters[seq[letter]] += 1
        else:
            dict_letters[seq[letter]] = 1

    tmp = copy.deepcopy(dict_letters)
    save_nb = 0
    for cles in dict_letters.keys():
        dict_letters[cles] = save_nb
        save_nb += tmp[cles]

    return dict_letters

def vec_generate_num(seq):
    '''Classer les bases par numéro de l'apparition de la base
    (T0, T1, T2...)
    '''
    dict_letters = {}
    vec = []

    for letter in range(len(seq)):
        if seq[letter] in dict_letters:
            dict_letters[seq[letter]] += 1
        else:
            dict_letters[seq[letter]] = 1

        vec.append(dict_letters[seq[letter]]-1)
    return vec

def str_reconstruct(dico, seq, val_seq):
    '''Reconstruire séquence originale
    '''
    pos = 0
    new_seq = ""
    for i in range(len(seq)):
        new_seq += seq[pos]
        pos = dico[seq[pos]] + val_seq[pos]
    new_seq = new_seq[:-1]
    new_seq = new_seq[::-1]
    return new_seq

def vec_generate_pos(dico, seq, val_seq):
    '''Retrouver position des bases de bt dans seq originale
    '''
    pos = 0
    pos_btw = [0]*(len(seq))
    for i in range(len(seq)):
        pos_btw[pos] = len(seq)-2-i
        pos = dico[seq[pos]] + val_seq[pos]
    return pos_btw

def dict_generate_OCC(dico, seq):
    '''Générer matrice OCC
    Nb bases ATGC à chaque ligne de bwt
    '''
    mat_mot = []
    tmp = {"$":0,"A":0,"C":0,"G":0,"T":0}
    for pos in range(len(seq)):
        tmp[seq[pos]] += 1
        copie = copy.deepcopy(tmp)
        mat_mot.append(copie)

    return mat_mot

def tupple_mapping(occ, dico_cpt, read):
    '''Retrouver read dans la séquence et
       renvoyer position dans la colonne first
    '''
    i = len(read)-1
    begin = dico_cpt[read[i]]
    end = dico_cpt[read[i]] + occ[-1][read[i]] -1
    while i>=1 and begin <= end:
        i = i-1
        begin = dico_cpt[read[i]] + occ[begin-1][read[i]]
        end = dico_cpt[read[i]] + occ[end][read[i]]-1
    return (begin,end,i)

def vec_search_pos_read(tpl_pos, pos):
    '''Retrouver positions réelles dans la séquence
    '''
    begin = tpl_pos[0]
    end = tpl_pos[1]
    if begin == end:
        return pos[begin] + 1
    elif end - begin > 1:
        return -2

    return -1

'''Fonction pour gérer substitution insertion délétion
Pour gérer les erreurs de read, on choisit d'utiliser une fonction
un peu lourde en temps de calcul mais fonctionnelle.
L'idée est, lorsqu'un read n'est pas aligné, de retenir la position
à laquelle le read n'a pas pu être aligné. On effectue ensuite les
modifications envisageables à cette position (substitution, insertion,
délétion) et on réaligne. Si le read ne s'aligne toujours pas, on retient
la modification qui a permis d'aligner le plus de lettres du read, puis on
effectue de nouveau un test des modifications à la nouvelle position 
(ce qui explique pourquoi on a choisi d'effectuer une fonction récursive).
On décide de définir un score de 5% de modifications au maximum pour
chaque read.

L'un des problèmes de la fonction est notamment qu'elle ne discrimine
pas les différentes modifications (une substitution vaut autant qu'une
insertion ou qu'une délétion, le score augmente de 1).
'''

def sub_indel(read, pos, score, trsf_occ, trsf_dico_pos, trsf_pos_btw):
    '''Gestion sub et indel
    '''
    if pos<1 or pos>len(read)-1 or score/len(read) > 0.05:
        return -1

    bases = ["A", "T", "G", "C"]
    list_test = []

    #Substitution
    for base in bases:
        readsub = read[0:pos] + base + read[pos+1:len(read)]

        pos_fin_bis = tupple_mapping(trsf_occ, trsf_dico_pos, readsub)
        pos_vrm_fin_bis = vec_search_pos_read(pos_fin_bis, trsf_pos_btw)

        if pos_vrm_fin_bis != -1:
            return pos_vrm_fin_bis
        list_test.append([readsub, pos_fin_bis[2]])

    #Deletion
    readdel = read[0:pos] + read[pos+1:len(read)]
    pos_fin_bis = tupple_mapping(trsf_occ, trsf_dico_pos, readdel)
    pos_vrm_fin_bis = vec_search_pos_read(pos_fin_bis, trsf_pos_btw)
    if pos_vrm_fin_bis != -1:
        return pos_vrm_fin_bis
    list_test.append([readdel, pos_fin_bis[2]])

    #Insertion
    for base in bases:
        readins = read[0:pos+1] + base + read[pos+1:len(read)]

        pos_fin_bis = tupple_mapping(trsf_occ, trsf_dico_pos, readins)
        pos_vrm_fin_bis = vec_search_pos_read(pos_fin_bis, trsf_pos_btw)

        if pos_vrm_fin_bis != -1:
            return pos_vrm_fin_bis
        list_test.append([readins, pos_fin_bis[2]])

    #Meilleure modification effectuée
    rang_min = 0
    for i in range(1,len(list_test)):
        if list_test[i][1] < list_test[rang_min][1]:
            rang_min = i

    score += 1
    sub_indel(list_test[rang_min][0], list_test[rang_min][1], score, trsf_occ, trsf_dico_pos, trsf_pos_btw)

def seq_read_fasta(filename):
    '''Lecture de fichier fasta
    '''
    vec_seq = []
    seq_tmp = []

    with open(filename, "r") as filin:
        for line in filin.readlines():
            if line.startswith(">"):
                if seq_tmp:
                    seq_tmp = "".join(seq_tmp)
                    vec_seq.append(seq_tmp)
                    seq_tmp = []
            else:
                seq_tmp.append(line.rstrip('\n'))

        seq_tmp = "".join(seq_tmp)
        vec_seq.append(seq_tmp)

    return vec_seq


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Alignment of reads using BWT.')
    parser.add_argument('seq_ref', type=str,
                        help='an input file with sequence of reference (fasta format)')
    parser.add_argument('--read', type=str,
                        help='reads for alignment')
    parser.add_argument('--filout', type=str,
                        help='file for results', default = "results.txt")
    parser.add_argument('--bwt', type=str,
                        help='bwt file')
    parser.add_argument('--bwt_inv', type=str,
                        help='bwt inv file')
    args = parser.parse_args()

    ref = "".join(seq_read_fasta(args.seq_ref))
    ref_inv = str_return_seq(ref)
    motif = seq_read_fasta(args.read)

    '''
    t_i = datetime.now()
    bwt = str_compute_bwt(ref)
    bwt_inv = str_compute_bwt(ref_inv)
    print("Temps bwt= {}".format(datetime.now()-t_i))
    '''

    t_i = datetime.now()

    if args.bwt is not None:
        with open(args.bwt, "r") as file_bwt:
            bwt = file_bwt.read()
    else:
        bwt = str_compute_bwt_alt(ref)
        with open("bwt_tmp.txt", "w") as file_bwt:
            file_bwt.write(bwt)

    if args.bwt_inv is not None:
        with open(args.bwt_inv, "r") as file_bwt:
            bwt_inv = file_bwt.read()
    else:
        bwt_inv = str_compute_bwt_alt(ref_inv)
        with open("bwt_tmp_inv.txt", "w") as file_bwt:
            file_bwt.write(bwt_inv)
    #print(bwt)
    print("\nTemps bwt= {}".format(datetime.now()-t_i))

    t_i = datetime.now()
    dico_pos = dict_col_first(bwt)
    dico_pos_inv = dict_col_first(bwt_inv)
    #print(dico_pos)
    print("Temps dico pos= {}".format(datetime.now()-t_i))

    t_i = datetime.now()
    vec_num = vec_generate_num(bwt)
    vec_num_inv = vec_generate_num(bwt_inv)
    #print(vec_num)
    print("Temps vec_num= {}".format(datetime.now()-t_i))

    t_i = datetime.now()
    seq_fin = str_reconstruct(dico_pos, bwt, vec_num)
    seq_fin_inv = str_reconstruct(dico_pos_inv, bwt_inv, vec_num_inv)
    #print(seq_fin)
    print("Temps seq_fin= {}".format(datetime.now()-t_i))

    t_i = datetime.now()
    pos_btw = vec_generate_pos(dico_pos, bwt, vec_num)
    pos_btw_inv = vec_generate_pos(dico_pos_inv, bwt_inv, vec_num_inv)
    #print(pos_btw)
    print("Temps pos_btw= {}".format(datetime.now()-t_i))

    t_i = datetime.now()
    occ = dict_generate_OCC(dico_pos, bwt)
    occ_inv = dict_generate_OCC(dico_pos_inv, bwt_inv)
    #print(occ)
    print("Temps occ= {}".format(datetime.now()-t_i))

    t_i = datetime.now()

    vec_ecriture = ''

    for num_motif in range(0,len(motif)):

        pos_fin1 = tupple_mapping(occ, dico_pos, "".join(motif[num_motif]))
        #print(pos_fin)
        pos_vrm_fin = vec_search_pos_read(pos_fin1, pos_btw)
        #print(pos_vrm_fin)
        if pos_vrm_fin == -1:
            pos_fin2 = tupple_mapping(occ_inv, dico_pos_inv, "".join(motif[num_motif]))
            #print(pos_fin)
            pos_vrm_fin = vec_search_pos_read(pos_fin2, pos_btw_inv)
            #print(pos_vrm_fin)

        if pos_vrm_fin == -1:
            score = 0
            pos_vrm_fin = sub_indel("".join(motif[num_motif]), pos_fin1[2], score, occ, dico_pos, pos_btw)

        if pos_vrm_fin == -1:
            score = 0
            pos_vrm_fin = sub_indel("".join(motif[num_motif]), pos_fin2[2], score, occ_inv, dico_pos_inv, pos_btw_inv)

        if num_motif%10000 == 0:
            print("\rTemps mapping= {}".format(datetime.now()-t_i), end = "")

        vec_ecriture+="readname{}    {}\n".format(num_motif, pos_vrm_fin)

    print("\nTemps mapping total= {}".format(datetime.now()-t_i))

    with open(args.filout, "w") as filout:
        filout.write(vec_ecriture)
    print("\n")
