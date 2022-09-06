import matplotlib.pyplot as plt
import random as random
from statistics import mode
import re
#from Bio import pairwise2
NUM_BP_FINESTRA = 1000
NUM_SEQUENZE = 5
NUM_BASI_MER = 12
PERCENTUALE= 0.1

def individuazione_sovrapposizioni(pos_in_seq1, pos_fin_seq1, pos_in_seq2, pos_fin_seq2, pos_in_sov, pos_fin_sov):
    i = 0
    j = 0
    while i < len(pos_in_seq1) or j < len(pos_in_seq1):
        if pos_in_seq1[i] == pos_in_seq2[j] and pos_fin_seq1[i] == pos_fin_seq2[j]:
            pos_in_sov.append(pos_in_seq1[i])
            pos_fin_sov.append(pos_fin_seq1[i])
        else:
            if pos_in_seq1[i] < pos_in_seq2[j] and pos_fin_seq1[i] < pos_fin_seq2[j]:
                if pos_in_seq2[j] < pos_fin_seq1[i]:
                    pos_in_sov.append(pos_in_seq2[j])
                    pos_fin_sov.append(pos_fin_seq1[i])
                else:
                    break
            elif pos_in_seq1[i] > pos_in_seq2[j] and pos_fin_seq1[i] > pos_fin_seq2[j]:
                if pos_in_seq1[i] < pos_fin_seq2[j]:
                    pos_in_sov.append(pos_in_seq1[i])
                    pos_fin_sov.append(pos_fin_seq2[j])
                else:
                    break
        if i == len(pos_in_seq1):
            i = i - 1
            j = j + 1
        elif j == len(pos_in_seq2):
            j = j - 1
            i = i + 1
        else:
            i = i + 1
            j = j + 1

def match_pattern(app1, app2):
    if app1 == app2:
        return True
    else:
        caratteri_diversi = 0
        j = 0
        while j < len(app1) and j < len(app2):
            if app1[j] != app2[j]:
                caratteri_diversi = caratteri_diversi + 1
            j = j + 1
        if len(app1) != len(app2):
            if len(app1) > len(app2):
                caratteri_diversi = caratteri_diversi + (len(app1) - len(app2))
            else:
                caratteri_diversi = caratteri_diversi + (len(app2) - len(app1))
        if caratteri_diversi <= (distanza_frequente/3):
            return True
        else:
            return False

f = open('prova.fna')
app = f.readline()
cromosoma = f.read() #cromosoma senza la prima riga
cromosoma = re.sub("\n","",cromosoma)
i = 0
j = NUM_BP_FINESTRA
finestra = []
finestre_mers = []

#suddivisione del cromosoma in finestre da 1kbp = 1000 basi azotate
while j <= len(cromosoma):
    finestra.append(cromosoma[i:j])
    if j == len(cromosoma):
        break
    i = i + NUM_BP_FINESTRA
    if i < (len(cromosoma) - NUM_BP_FINESTRA):
        j = j + NUM_BP_FINESTRA
    else:
        j = len(cromosoma)
cont = 0

# suddivisione delle finestre in 12-mers
while cont < len(finestra):
    i = 0
    j = NUM_BASI_MER
    mers = []
    while i <= (len(finestra[cont]) - NUM_BASI_MER):
        mers.append(finestra[cont][i:j])
        i = i + 1
        if i < (len(finestra[cont]) - NUM_BASI_MER):
            j = j + 1
        else:
            j = len(finestra[cont])
            mers.append(finestra[cont][i:j])
            break

    finestre_mers.append(mers)
    cont = cont + 1
i=0
j=1
cont = 0
matches = []
matches_finestra = []

#calcolo matches dei mers per ogni finestra
while cont < len(finestra):
    cont_macthes = 0
    match = True
    i=0
    j=1
    while i < (len(finestre_mers[cont]) - 1):
        while j < len(finestre_mers[cont]):
            if finestre_mers[cont][i] != finestre_mers[cont][j]:
                match = False
            else:
                cont_macthes = cont_macthes + 1
                matches.append(cont_macthes)
            j = j + 1
        i = i + 1
        j = i + 1
    matches_finestra.append(cont_macthes)
    cont = cont + 1

#ricalcolo finestre in base al numero di mathces
cont = 0
nuove_finestre = []
pos = []
while cont < len(finestra):
    dieci_percento = len(finestre_mers[cont]) * PERCENTUALE
    if matches_finestra[cont] > dieci_percento:
        pos.append(cont)
    cont = cont + 1

cont = 0
i = 0
pos_precedente = 0
n_finestra = ""
entrato = False
while cont < len(finestra):
    if i < len(pos):
        if cont == pos[i]:
            if pos[i] == (pos_precedente + 1):
                n_finestra = n_finestra + finestra[cont]
                entrato = True
            else:
                if entrato:
                    nuove_finestre.append(n_finestra)
                    entrato = False
                    n_finestra = finestra[cont]
                    nuove_finestre.append(n_finestra)
                else:
                    n_finestra = finestra[cont]
                    nuove_finestre.append(n_finestra)
            pos_precedente = pos[i]
            i = i + 1
        cont = cont + 1
    else:
        break

#calcolo distanze tra i matches dei vari mers
cont = 0
distanze = []
matches = []
cont_macthes = 0
while cont < len(nuove_finestre):
    cont_macthes = 0
    match = True
    i=0
    j=1
    while i < (len(nuove_finestre[cont]) - NUM_BASI_MER):
        while j < len(nuove_finestre[cont]):
            pos_finale1 = i + NUM_BASI_MER
            pos_finale2 = j + NUM_BASI_MER
            if nuove_finestre[cont][i:pos_finale1] != nuove_finestre[cont][j:pos_finale2]:
                match = False
            else:
                matches.append(cont_macthes)
                cont_macthes = cont_macthes +1
                distanze.append(j - i)
            j = j + 1
        i = i + 1
        j = i + 1
    cont = cont + 1

#creazione istogramma similarita' media di ciscuna finestra
figura, grafico= plt.subplots()
grafico.bar(x = matches, height = distanze, fc = "orange")
grafico.set_xlabel("Matches")
grafico.set_ylabel("Distanza")
grafico.set_title("Istogramma distanze")
plt.show()

#trovare distanza più frequente
distanza_frequente = mode(distanze)

print("Distanza con frequenza maggiore: ", distanza_frequente)

#generazione 5 sequenze random di lunghezza N pari alla distanza più frequente
i = 0
sequenze = []
sequenze_regioni = []
cont = 0
while cont < len(nuove_finestre):
    sequenze = []
    i = 0
    while i < NUM_SEQUENZE:
        j = 0
        limite = len(nuove_finestre[cont]) - distanza_frequente
        numero_casuale = random.randrange(0, limite)
        seq = []
        while j < distanza_frequente:
                seq.append(nuove_finestre[cont][numero_casuale])
                j = j + 1
                numero_casuale = numero_casuale + 1
        i = i + 1
        sequenze.append("".join(seq))
    sequenze_regioni.append(sequenze)
    cont = cont + 1
print("Sequenze generate", sequenze_regioni)

#match delle sequenze generate con la regione di interesse
cont = 0
#matches_sequenze_regioni = []
#match_sequenze = []
pos_iniziali = []
pos_finali = []
pos_in_se1 = []
pos_fin_se1 = []
pos_in_se2 = []
pos_fin_se2 = []
pos_in_se3 = []
pos_fin_se3 = []
pos_in_se4 = []
pos_fin_se4 = []
pos_in_se5 = []
pos_fin_se5 = []

while cont < len(nuove_finestre):
    #match_sequenze = []
    k = 0
    posizioni_iniziali_sequenze = []
    posizioni_finali_sequenze= []
    while k < len(sequenze_regioni[cont]):
        i = 0
        pos_iniziali = []
        pos_finali = []
        while i < len(nuove_finestre[cont]):
            controllo = match_pattern(sequenze_regioni[cont][k],  nuove_finestre[cont][i:(i+distanza_frequente)])
            if controllo:
                pos_iniziali.append(i)
                pos_finali.append(i+distanza_frequente)
            i = i + 1
        if k == 0:
            pos_in_se1.append(pos_iniziali)
            pos_fin_se1.append(pos_finali)
        elif k == 1:
            pos_in_se2.append(pos_iniziali)
            pos_fin_se2.append(pos_finali)
        elif k == 2:
            pos_in_se3.append(pos_iniziali)
            pos_fin_se3.append(pos_finali)
        elif k == 3:
            pos_in_se4.append(pos_iniziali)
            pos_fin_se4.append(pos_finali)
        elif k == 4:
            pos_in_se5.append(pos_iniziali)
            pos_fin_se5.append(pos_finali)
        k = k + 1
    cont = cont + 1

#individuazione overlap
sovrapposizioni = []
pos_in_sov = []
pos_fin_sov = []
posizioni_iniziali_sovrapposizioni = []
posizioni_finali_sovrapposizioni = []
cont = 0
while cont < len(nuove_finestre):
    pos_in_sov = []
    pos_fin_sov = []
    individuazione_sovrapposizioni(pos_in_se1[cont], pos_fin_se1[cont], pos_in_se2[cont], pos_fin_se2[cont], pos_in_sov, pos_fin_sov)
    individuazione_sovrapposizioni(pos_in_se1[cont], pos_fin_se1[cont], pos_in_se3[cont], pos_fin_se3[cont], pos_in_sov, pos_fin_sov)
    individuazione_sovrapposizioni(pos_in_se1[cont], pos_fin_se1[cont], pos_in_se4[cont], pos_fin_se4[cont], pos_in_sov, pos_fin_sov)
    individuazione_sovrapposizioni(pos_in_se1[cont], pos_fin_se1[cont], pos_in_se5[cont], pos_fin_se5[cont], pos_in_sov, pos_fin_sov)

    individuazione_sovrapposizioni(pos_in_se2[cont], pos_fin_se2[cont], pos_in_se3[cont], pos_fin_se3[cont], pos_in_sov, pos_fin_sov)
    individuazione_sovrapposizioni(pos_in_se2[cont], pos_fin_se2[cont], pos_in_se4[cont], pos_fin_se4[cont], pos_in_sov, pos_fin_sov)
    individuazione_sovrapposizioni(pos_in_se2[cont], pos_fin_se2[cont], pos_in_se5[cont], pos_fin_se5[cont], pos_in_sov, pos_fin_sov)

    individuazione_sovrapposizioni(pos_in_se3[cont], pos_fin_se3[cont], pos_in_se4[cont], pos_fin_se4[cont], pos_in_sov, pos_fin_sov)
    individuazione_sovrapposizioni(pos_in_se3[cont], pos_fin_se3[cont], pos_in_se5[cont], pos_fin_se5[cont], pos_in_sov, pos_fin_sov)

    individuazione_sovrapposizioni(pos_in_se4[cont], pos_fin_se4[cont], pos_in_se5[cont], pos_fin_se5[cont], pos_in_sov, pos_fin_sov)
    posizioni_iniziali_sovrapposizioni.append(pos_in_sov)
    posizioni_finali_sovrapposizioni.append(pos_fin_sov)
    cont = cont + 1

cont = 0
sovrapposizioni_regione = []
while cont < len(nuove_finestre):
    i = 0
    sovrapposizioni = []
    while i < len(posizioni_iniziali_sovrapposizioni[cont]):
        sovrapposizioni.append(nuove_finestre[cont][posizioni_iniziali_sovrapposizioni[cont][i]:posizioni_finali_sovrapposizioni[cont][i]])
        i = i + 1
    sovrapposizioni_regione.append(sovrapposizioni)
    cont = cont + 1
print(sovrapposizioni_regione)
#gestione overlap
cont = 0
sovrapposizioni_salvate = []
while cont < len(nuove_finestre):
    i = 0
    sovrapposizioni = []
    while i < len(sovrapposizioni_regione[cont]):
        if len(sovrapposizioni_regione[cont][i]) < 10:
            mezzo = int(len(sovrapposizioni_regione[cont][i])/2)
            sovrapposizioni.append(sovrapposizioni_regione[cont][i][0:mezzo])
            fine = len(sovrapposizioni_regione[cont][i])
            sovrapposizioni.append(sovrapposizioni_regione[cont][i][mezzo:fine])
        else:
            sovrapposizioni.append(sovrapposizioni_regione[cont][i])
        i = i + 1
    sovrapposizioni_salvate.append(sovrapposizioni)
    cont = cont + 1
print(sovrapposizioni_salvate)

#allineamento sequenze
'''
cont = 0
while cont < len(sovrapposizioni_regione):
    i = 0
    while i < (len(sovrapposizioni_regione[cont]) - 1):
        j = i + 1
        allineamenti = pairwise2.align.globalxx(sovrapposizioni_regione[cont][i], sovrapposizioni_regione[cont][j])
        #for al in allineamenti:
            #print(pairwise2.format_alignment(*al))
       #i = i + 1
    cont = cont + 1
'''
cont = 0
while cont < len(sovrapposizioni_regione):
    fi = open("sequenze_sovrapposizioni.fna", "w")
    i = 0
    while i < len(sovrapposizioni_regione[cont]):
        fi.write('>' + str(i+1) + '\n' + str(sovrapposizioni_regione[cont][i]) + '\n')
        i = i + 1
    cont = cont + 1
f1 = open('sequenze_sovrapposizioni_allineate.fna')
sovrapposizioni_allineate = f1.read()
print (sovrapposizioni_allineate)
