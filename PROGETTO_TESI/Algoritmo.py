import matplotlib.pyplot as plt
#import random as random
import re
NUM_BP_FINESTRA = 1000
NUM_SEQUENZE = 5
NUM_BASI_MER = 12
PERCENTUALE= 0.1

f = open('D:/TESI/ncbi-genomes-2022-08-10/prova.fna')
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

pos_iniziale = 0
distanze = []
i=0
j=1
cont = 0
cont_distanze = []
matches = []
matches_finestra = []

#calcolo match per i mer per ogni finestra tenendo conto delle distanze
while cont < len(finestra):
    c_distanze = 0
    cont_macthes = 0
    match = True
    while i < (len(finestre_mers[cont]) - 1):
        while j < len(finestre_mers[cont]):
            pos_iniziale = j
            if finestre_mers[cont][i] != finestre_mers[cont][j]:
                match = False
            else:
                distanze.append(pos_iniziale - i)
                cont_macthes = cont_macthes + 1
                matches.append(cont_macthes)
                c_distanze = c_distanze + 1

            j = j + 1
        i = i + 1
        j = i + 1
    cont_distanze.append(c_distanze)
    cont = cont + 1
    matches_finestra.append(cont_macthes)

#ricalcolo finestre in base alle distanze
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
print(distanze[0:20])
#creazione istogramma similarita' media di ciscuna finestra
figura, grafico= plt.subplots()
grafico.bar(x = matches, height = distanze, fc = "orange")
grafico.set_xlabel("Matches")
grafico.set_ylabel("Distanza")
grafico.set_title("Istogramma distanze")
plt.show()
