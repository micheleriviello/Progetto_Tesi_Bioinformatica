import matplotlib.pyplot as plt
import random as random
import re
import multiprocessing
import concurrent.futures
from ClasseRisultato import RisultatoCalcoloMatch
import time
NUMERO_CORE = multiprocessing.cpu_count()
NUM_BP_FINESTRA = 1000
NUM_SEQUENZE = 5
NUM_BASI_MER = 12
PERCENTUALE = 0.1

def match_pattern(app1, app2, distanza_frequente):
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

def calcola_match(cont, i, nuove_finestre):
    risultato = RisultatoCalcoloMatch()
    risultato.matches = []
    risultato.num_matches = 0
    risultato.frequenza_dist = {}
    cont = int(cont)
    j = int(i + NUM_BASI_MER)
    while j < (len(nuove_finestre[cont]) - NUM_BASI_MER):
        pos_finale1 = i + NUM_BASI_MER
        pos_finale2 = j + NUM_BASI_MER
        if nuove_finestre[cont][i:pos_finale1] == nuove_finestre[cont][j:pos_finale2]:
            risultato.matches.append(risultato.num_matches)
            risultato.num_matches = risultato.num_matches + 1
            distanza = j - (pos_finale1 - 1)
            if distanza in risultato.frequenza_dist:
                risultato.frequenza_dist[distanza] += 1
            else:
                risultato.frequenza_dist[distanza] = 1
        j = j + 1
    return risultato

if __name__ == "__main__":
    f = open('prova.fna')
    app = f.readline()
    cromosoma = f.read() #cromosoma senza la prima riga
    cromosoma = re.sub("\n","",cromosoma)
    i = 0
    j = NUM_BP_FINESTRA
    finestra = []
    finestre_mers = []
    start = time.time()
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
    del cromosoma
    ft = open('filetempi.txt', 'w')
    end = time.time()
    ft.write("Tempo suddivisione finestre: " + str(end-start) + "\n")
    start = time.time()
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
    end = time.time()
    ft.write("Tempo suddivisione finestre in 12-mers: " + str(end-start) + "\n")
    start = time.time()
    #calcolo matches dei mers per ogni finestra
    while cont < len(finestra):
        cont_matches = 0
        match = True
        i=0
        j=1
        while i < (len(finestre_mers[cont]) - 1):
            while j < len(finestre_mers[cont]):
                if finestre_mers[cont][i] != finestre_mers[cont][j]:
                    match = False
                else:
                    cont_matches = cont_matches + 1
                    matches.append(cont_matches)
                j = j + 1
            i = i + 1
            j = i + 1
        matches_finestra.append(cont_matches)
        cont = cont + 1
    end = time.time()
    ft.write("Tempo matches finestre: " + str(end-start) + "\n")
    start = time.time()
    #ricalcolo finestre in base al numero di mathces
    cont = 0
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
    nuove_finestre = []
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
                    else:
                        n_finestra = finestra[cont]
                pos_precedente = pos[i]
                i = i + 1
            cont = cont + 1
        else:
            break
    nuove_finestre.append(n_finestra)
    print(nuove_finestre)
    end = time.time()
    ft.write("Tempo ricalcolo finestre: " + str(end-start) + "\n")
    del mers
    del finestre_mers
    del finestra
    del matches_finestra
    #calcolo distanze tra i matches dei vari mers
    cont = 0
    diz_dis_freq = {}
    matches = []
    cont_matches = 0
    executor = concurrent.futures.ProcessPoolExecutor()
    while cont < len(nuove_finestre):
        start = time.time()
        cont_matches = 0
        irange = NUMERO_CORE
        i = range(0,irange)
        while i[NUMERO_CORE - 1] < (len(nuove_finestre[cont]) - NUM_BASI_MER):
            print("I: ", i, "SU: ", (len(nuove_finestre[cont]) - NUM_BASI_MER), "CONT: ", cont)
            irange = irange + NUMERO_CORE
            results = [executor.submit(calcola_match, cont, k, nuove_finestre) for k in i]
            for f in concurrent.futures.as_completed(results):
                matches = matches + f.result().matches
                cont_matches = cont_matches + f.result().num_matches
                for chiave in f.result().frequenza_dist.keys():
                    if chiave not in diz_dis_freq:
                        diz_dis_freq[chiave] = f.result().frequenza_dist[chiave]
                    else:
                        diz_dis_freq[chiave] = diz_dis_freq[chiave] + f.result().frequenza_dist[chiave]
            i = range(irange, irange + multiprocessing.cpu_count())
        if(i[0] < (len(nuove_finestre[cont]) - NUM_BASI_MER)):
            i = range(i[0], (len(nuove_finestre[cont]) - NUM_BASI_MER))
            print("I: ", i, "SU: ", (len(nuove_finestre[cont]) - NUM_BASI_MER), "CONT: ", cont)
            results = [executor.submit(calcola_match, cont, k, nuove_finestre) for k in i]
            for f in concurrent.futures.as_completed(results):
                matches = matches + f.result().matches
                cont_matches = cont_matches + f.result().num_matches
                for chiave in f.result().frequenza_dist.keys():
                    if chiave not in diz_dis_freq:
                        diz_dis_freq[chiave] = f.result().frequenza_dist[chiave]
                    else:
                        diz_dis_freq[chiave] = diz_dis_freq[chiave] + f.result().frequenza_dist[chiave]
        end = time.time()
        ft.write("Tempo calcolo distanze match finestra " + str(cont+1) + ": " + str(end-start) + "\n")
        cont = cont + 1

    #creazione istogramma similarita' media di ciscuna finestra
    figura, grafico = plt.subplots()
    grafico.bar(x=list(diz_dis_freq.keys()), height=list(diz_dis_freq.values()), fc="orange")
    grafico.set_xlabel("Valori distanza")
    grafico.set_ylabel("Valori frequenza")
    grafico.set_title("Istogramma frequenza distanze")
    plt.show()

    #trovare distanza più frequente
    distanza_frequente = max(diz_dis_freq, key=diz_dis_freq.get)
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
    start = time.time()
    #match delle sequenze generate con la regione di interesse
    cont = 0
    pos_in = []
    pos_fin = []
    posizioni_iniziali_match = []
    posizioni_finali_match = []
    while cont < len(nuove_finestre):
        k = 0
        while k < NUM_SEQUENZE:
            i = 0
            while i < (len(nuove_finestre[cont]) - distanza_frequente):
                p_i = i
                p_f = i + distanza_frequente
                if match_pattern(sequenze_regioni[cont][k], nuove_finestre[cont][p_i:p_f], distanza_frequente):
                    pos_in.append(p_i)
                    pos_fin.append(p_f)
                i = i + 1
            k = k + 1
        posizioni_iniziali_match.append(pos_in)
        posizioni_finali_match.append(pos_fin)
        cont = cont + 1
    end = time.time()
    ft.write("Tempo match tra sequenze e regioni di interesse: " + str(end-start) + "\n")
    #individuazione e gestione delle sovrapposizioni
    cont = 0

    while cont < len(nuove_finestre):
        start = time.time()
        match = []
        f1 = open('sequenze_allineate'+str(cont+1)+'.fna', "w")
        i = 0
        while i < (len(posizioni_iniziali_match[cont])-1):
            j = i + 1
            p_i = 0
            p_f = 0
            trovata_s = False
            while j < len(posizioni_iniziali_match[cont]):
                if posizioni_iniziali_match[cont][i] == posizioni_iniziali_match[cont][j] and posizioni_finali_match[cont][i] == posizioni_finali_match[cont][j]:
                    p_i = posizioni_iniziali_match[cont][i]
                    p_f = posizioni_finali_match[cont][i]
                elif posizioni_finali_match[cont][i] > posizioni_iniziali_match[cont][j] and posizioni_iniziali_match[cont][j] > posizioni_iniziali_match[cont][i]:
                    p_i = posizioni_iniziali_match[cont][j]
                    p_f = posizioni_finali_match[cont][i]
                elif posizioni_finali_match[cont][j] > posizioni_iniziali_match[cont][i] and posizioni_iniziali_match[cont][i] > posizioni_iniziali_match[cont][j]:
                    p_i = posizioni_iniziali_match[cont][i]
                    p_f = posizioni_finali_match[cont][j]
                else:
                    p_i = -1
                    p_f = -1
                if p_i != -1:
                    trovata_s = True
                    if (p_f - p_i) < 10:
                        mezzo = int((p_i + p_f)/2)
                        if p_i < posizioni_iniziali_match[cont][i]:
                            match.append(nuove_finestre[cont][mezzo:posizioni_finali_match[cont][i]])
                            match.append(nuove_finestre[cont][posizioni_iniziali_match[cont][j]:mezzo])
                        else:
                            match.append(nuove_finestre[cont][mezzo:posizioni_finali_match[cont][j]])
                            match.append(nuove_finestre[cont][posizioni_iniziali_match[cont][i]:mezzo])
                    else:
                        if p_i < posizioni_iniziali_match[cont][i]:
                            match.append(nuove_finestre[cont][p_f:posizioni_finali_match[cont][i]])
                            match.append(nuove_finestre[cont][posizioni_iniziali_match[cont][j]:p_f])
                        else:
                            match.append(nuove_finestre[cont][p_i:posizioni_finali_match[cont][j]])
                            match.append(nuove_finestre[cont][posizioni_iniziali_match[cont][i]:p_i])
                j = j + 1
            if trovata_s == False:
                match.append(nuove_finestre[cont][posizioni_iniziali_match[cont][i]:posizioni_finali_match[cont][i]])
            i = i + 1
        k = 0
        while k < len(match):
            if match[k] != '':
                f1.write(">"+str(k+1)+'\n'+match[k]+'\n')
            k = k + 1
        f1.close()
        del match
        end = time.time()
        ft.write("Tempo individuazione e gestione overlap per finestra " + str(cont+1) + ": " + str(end-start) + "\n")
        cont = cont + 1