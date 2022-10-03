import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import random
import re
import numpy as np
import math
import multiprocessing
import concurrent.futures
from ClasseRisultato import RisultatoCalcoloMatch
import time
import json
NUMERO_CORE = multiprocessing.cpu_count()
LIMITE_SINGLE_CORE = 10000
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

def calcola_match(listaRipetizioni):
    risultato = {}

    for i in range(len(listaRipetizioni) - 1):
        for j in range(i, len(listaRipetizioni)):
            if listaRipetizioni[j] > listaRipetizioni[i] + NUM_BASI_MER:
                distanza = listaRipetizioni[j] - listaRipetizioni[i]
                if distanza in risultato:
                    risultato[distanza] += 1
                else:
                    risultato[distanza] = 1

    return risultato

def individua_e_gestisci_sovrapposizioni(cont, posizioni_iniziali_match, posizioni_finali_match, nuove_finestre):
    match = []
    f1 = open('sequenze_allineate'+str(cont+1)+'.fna', "w")
    i = 0
    while i < (len(posizioni_iniziali_match)-1):
        j = i + 1
        p_i = 0
        p_f = 0
        trovata_s = False
        while j < len(posizioni_iniziali_match):
            if posizioni_iniziali_match[i] == posizioni_iniziali_match[j] and posizioni_finali_match[i] == posizioni_finali_match[j]:
                p_i = posizioni_iniziali_match[i]
                p_f = posizioni_finali_match[i]
            elif posizioni_finali_match[i] > posizioni_iniziali_match[j] and posizioni_iniziali_match[j] > posizioni_iniziali_match[i]:
                p_i = posizioni_iniziali_match[j]
                p_f = posizioni_finali_match[i]
            elif posizioni_finali_match[j] > posizioni_iniziali_match[i] and posizioni_iniziali_match[i] > posizioni_iniziali_match[j]:
                p_i = posizioni_iniziali_match[i]
                p_f = posizioni_finali_match[j]
            else:
                p_i = -1
                p_f = -1
            if p_i != -1:
                trovata_s = True
                if (p_f - p_i) < 10:
                    mezzo = int((p_i + p_f)/2)
                    if p_i < posizioni_iniziali_match[i]:
                        match.append(nuove_finestre[mezzo:posizioni_finali_match[i]])
                        match.append(nuove_finestre[posizioni_iniziali_match[j]:mezzo])
                    else:
                        match.append(nuove_finestre[mezzo:posizioni_finali_match[j]])
                        match.append(nuove_finestre[posizioni_iniziali_match[i]:mezzo])
                else:
                    if p_i < posizioni_iniziali_match[i]:
                        match.append(nuove_finestre[p_f:posizioni_finali_match[i]])
                        match.append(nuove_finestre[posizioni_iniziali_match[j]:p_f])
                    else:
                        match.append(nuove_finestre[p_i:posizioni_finali_match[j]])
                        match.append(nuove_finestre[posizioni_iniziali_match[i]:p_i])
            j = j + 1
        if trovata_s == False:
            match.append(nuove_finestre[posizioni_iniziali_match[i]:posizioni_finali_match[i]])
        i = i + 1
    k = 0
    while k < len(match):
        if match[k] != '':
            f1.write(">"+str(k+1)+'\n'+match[k]+'\n')
        k = k + 1
    del match
    f1.close()

if __name__ == "__main__":
    cromosomi = ["cromosoma21"]
    for nomefile in cromosomi:
        cartellaRisultati = "Risultati/"
        ftt = open(cartellaRisultati+'filetempi_overlap_multithread'+nomefile+'.txt', 'w')
        f = open('Genoma_Umano/'+nomefile+'.fna')
        app = f.readline()
        cromosoma = f.read() #cromosoma senza la prima riga
        cromosoma = re.sub("\n","",cromosoma)


        i = 0
        j = NUM_BP_FINESTRA
        finestra = []
        finestre_mers = []
        start = time.time()
        print(len(cromosoma))
        print("Sto qua 1")
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
        print(len(finestra))
        end = time.time()
        ftt.write("Tempo suddivisione finestre: " + str(end-start) + "\n")
        start = time.time()
        # suddivisione delle finestre in 12-mers
        print("Sto qua 2")

        for cont in range(0, len(finestra)):
            print(cont," SU: ",len(finestra), "1")
            lunghezzaFinestraCorrente = len(finestra[cont]) - NUM_BASI_MER
            mers = []
            for i in range(0, lunghezzaFinestraCorrente):
                mers.append(finestra[cont][i : i + NUM_BASI_MER])

            finestre_mers.append(mers)

        i=0
        j=1
        cont = 0
        matches = []
        matches_finestra = []
        end = time.time()
        ftt.write("Tempo suddivisione finestre in 12-mers: " + str(end-start) + "\n")
        start = time.time()
        print("Sto qua 3")
        #calcolo matches dei mers per ogni finestra
        for cont in range(0, len(finestra)):
            print(cont," SU: ",len(finestra))
            diz = {}
            cont_matches = 0
            #match = True
            for i in range(0, len(finestre_mers[cont]) - 1):
                kmer = str(finestre_mers[cont][i: i + NUM_BASI_MER])
                if kmer not in diz.keys():
                    diz[kmer] = 0
                diz[kmer] = diz[kmer] + 1
            for val in diz.values():
                if val > 1:
                    cont_matches += val - 1
            matches_finestra.append(cont_matches)
        end = time.time()



        ftt.write("Tempo matches finestre: " + str(end-start) + "\n")
        start = time.time()
        #ricalcolo finestre in base al numero di mathces
        cont = 0
        pos = []
        print("Sto qua 4")
        fl = open(cartellaRisultati+"Risultati_finestre_pre_raggruppamento_"+nomefile+".txt", "w")
        fl.write('numero_finestra, ' + 'posizione_iniziale, ' + 'ripetizioni' + '\n')
        k = 0
        while cont < len(finestra):
            dieci_percento = len(finestre_mers[cont]) * PERCENTUALE
            if matches_finestra[cont] > dieci_percento:
                fl.write(str(cont+1) + ', ' + str(k) + ', SI\n')
                k = k + NUM_BP_FINESTRA
                pos.append(cont)
            else:
                fl.write(str(cont+1) + ', ' + str(k) + ', NO\n')
                k = k + NUM_BP_FINESTRA
            if k >= len(cromosoma):
                break
            cont = cont + 1
        cont = 0
        fl.close()
        i = 1
        pos_precedente = pos[0]
        n_finestra = finestra[pos[0]]

        entrato = False
        nuove_finestre = []
        ft = open(cartellaRisultati+"Risultati_finestre_raggruppamento_"+nomefile+".txt", "w")
        ft.write('numero_finestra, ' + 'posizione_iniziale, ' + 'numero_basi\n')
        k = 0
        j = 1

        ultima_posizione = 0
        superato = False
        indice = 1
        som = 0
        p_i = 0
        while cont < len(finestra):
            print(cont," SU: ",len(finestra))
            if i < len(pos):
                if cont == pos[i]:
                    if pos[i] == (pos_precedente + 1):
                        n_finestra = n_finestra + finestra[cont]
                        entrato = True
                    else:
                        nuove_finestre.append(n_finestra)
                        ft.write(str(indice) + ', ' + str(p_i) + ', ' + str(len(n_finestra)) + '\n')
                        p_i = pos[i] * NUM_BP_FINESTRA
                        indice = indice + 1
                        if entrato:
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
        ft.write(str(indice) + ', ' + str(p_i) + ', ' + str(len(n_finestra)) + '\n')
        end = time.time()
        ft.close()

        with open(cartellaRisultati+"nuove_finestre_"+nomefile, "w") as fl:
            json.dump(nuove_finestre, fl)

        ftt.write("Tempo ricalcolo finestre: " + str(end-start) + "\n")
        del mers
        del finestre_mers
        del finestra
        del matches_finestra

        #calcolo distanze tra i matches dei vari mers
        executor = concurrent.futures.ProcessPoolExecutor()
        cont = 0
        diz_dis_freq = {}
        matches = []
        cont_matches = 0

        with open(cartellaRisultati+"nuove_finestre_"+nomefile, "r") as fl:
            nuove_finestre = json.load(fl)

        while cont < len(nuove_finestre):
            print("Finestra: ",cont + 1,": ", len(nuove_finestre[cont]))
            cont = cont + 1

        start = time.time()
        for cont in range(0, len(nuove_finestre)):
            diz = {}
            stringaStampa = ""
            for i in range(0, len(nuove_finestre[cont]) - NUM_BASI_MER):
                stringaStampa = "I: ", i, "SU: ", (len(nuove_finestre[cont]) - NUM_BASI_MER), "CONT: ", cont, "SU: ", len(nuove_finestre)
                print(stringaStampa)
                kmer = str(nuove_finestre[cont][i: i + NUM_BASI_MER])
                if kmer not in diz.keys():
                    diz[kmer] = []
                diz[kmer].append(i + NUM_BASI_MER)

            percentuale = 0
            ripetizioniFatte = 0
            bufferPieno = False
            buffer = []
            for listaRipetizioni in diz.values():
                ripetizioniFatte += len(listaRipetizioni)
                percentuale = ripetizioniFatte/((len(nuove_finestre[cont]) - NUM_BASI_MER)/100)
                print(stringaStampa," ", percentuale,"%")
                if(len(listaRipetizioni) > 3000):
                    buffer.append(listaRipetizioni)
                    if(len(buffer) == NUMERO_CORE):
                        print("STO USANDO IL BUFFER")
                        results = [executor.submit(calcola_match, k) for k in buffer]
                        for f in concurrent.futures.as_completed(results):
                            for chiave in f.result().keys():
                                if chiave not in diz_dis_freq:
                                    diz_dis_freq[chiave] = f.result()[chiave]
                                else:
                                    diz_dis_freq[chiave] = diz_dis_freq[chiave] + f.result()[chiave]
                        buffer = []
                else:
                    for i in range(len(listaRipetizioni) - 1):
                        for j in range(i, len(listaRipetizioni)):
                            if listaRipetizioni[j] > listaRipetizioni[i] + NUM_BASI_MER:
                                distanza = listaRipetizioni[j] - listaRipetizioni[i]
                                if distanza in diz_dis_freq:
                                    diz_dis_freq[distanza] += 1
                                else:
                                    diz_dis_freq[distanza] = 1
            if(len(buffer) > 0):
                results = [executor.submit(calcola_match, k) for k in buffer]
                print("STO USANDO IL BUFFER PER L'ULTIMA VOLTA")
                for f in concurrent.futures.as_completed(results):
                    for chiave in f.result().keys():
                        if chiave not in diz_dis_freq:
                            diz_dis_freq[chiave] = f.result()[chiave]
                        else:
                            diz_dis_freq[chiave] = diz_dis_freq[chiave] + f.result()[chiave]
                buffer = []
        end = time.time()
        ftt.write("Tempo calcolo distanze match finestra " + str(cont+1) + ": " + str(end-start) + "\n")

        '''
        for cont in range(0, len(nuove_finestre)):
            start = time.time()
            cont_matches = 0
            chunk = []
            if(len(nuove_finestre[cont]) < LIMITE_SINGLE_CORE):
                for i in range(0, len(nuove_finestre[cont])):
                    print("I: ", i, "SU: ", (len(nuove_finestre[cont]) - NUM_BASI_MER), "CONT: ", cont, "SU: ", len(nuove_finestre), "SINGLECORE")
                    for j in range(i + NUM_BASI_MER, len(nuove_finestre[cont]) - NUM_BASI_MER):
                        pos_finale1 = i + NUM_BASI_MER
                        pos_finale2 = j + NUM_BASI_MER
                        if nuove_finestre[cont][i:pos_finale1] == nuove_finestre[cont][j:pos_finale2]:
                            matches.append(cont_matches)
                            cont_matches = cont_matches + 1
                            distanza = j - (pos_finale1 - 1)
                            if distanza in diz_dis_freq:
                                diz_dis_freq[distanza] += 1
                            else:
                                diz_dis_freq[distanza] = 1
            else:
                irange = NUMERO_CORE
                i = range(0,irange)
                while i[NUMERO_CORE - 1] < (len(nuove_finestre[cont]) - NUM_BASI_MER):
                    irange = irange + NUMERO_CORE
                    chunk.append(i)
                    i = range(irange, irange + NUMERO_CORE)
                if(i[0] < (len(nuove_finestre[cont]) - NUM_BASI_MER)):
                    i = range(i[0], (len(nuove_finestre[cont]) - NUM_BASI_MER))
                    chunk.append(i)
                for i in chunk:
                    print("I: ", i, "SU: ", (len(nuove_finestre[cont]) - NUM_BASI_MER), "CONT: ", cont, "SU: ", len(nuove_finestre))
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
            #ft.write("Tempo calcolo distanze match finestra " + str(cont+1) + ": " + str(end-start) + "\n")
        '''
        #creazione istogramma similarita" media di ciscuna finestra
        '''
        figura, grafico = plt.subplots()
        distanza_frequente = max(diz_dis_freq, key=diz_dis_freq.get)
        plt.rcParams["figure.dpi"] = 1200
        grafico.bar(x=list(diz_dis_freq.keys()), height=list(diz_dis_freq.values()), fc="orange")
        grafico.set_xlabel("Valori distanza")
        grafico.set_ylabel("Valori frequenza")
        grafico.set_title("Istogramma frequenza distanze")
        plt.text(distanza_frequente + 12, diz_dis_freq[distanza_frequente] + 12, "Distanza con frequenza maggiore: "+ str(distanza_frequente))
        cerchio = Circle((distanza_frequente, diz_dis_freq[distanza_frequente]), 3, color="red", lw=3, fill= False, zorder=10)
        plt.gca().add_patch(cerchio)
        plt.savefig(cartellaRisultati+'hist_'+nomefile+'.pdf')
        plt.show()
        '''
        #trovare distanza più frequente
        distanza_frequente = max(diz_dis_freq, key=diz_dis_freq.get)
        print("Distanza con frequenza maggiore: ", distanza_frequente)
        print("Frequenza della distanza: ", diz_dis_freq[distanza_frequente])
        '''
        #generazione 5 sequenze random di lunghezza N pari alla distanza più frequente
        i = 0
        sequenze = []
        sequenze_regioni = []
    
        for cont in range(0, len(nuove_finestre)):
            sequenze = []
            for i in range(0, NUM_SEQUENZE):
                limite = len(nuove_finestre[cont]) - distanza_frequente
                numero_casuale = random.randrange(0, limite)
                seq = []
                for j in range(0, distanza_frequente):
                    seq.append(nuove_finestre[cont][numero_casuale])
                    numero_casuale = numero_casuale + 1
                sequenze.append("".join(seq))
            sequenze_regioni.append(sequenze)
        print("Sequenze generate", sequenze_regioni)
        start = time.time()
        #match delle sequenze generate con la regione di interesse
        pos_in = []
        pos_fin = []
        posizioni_iniziali_match = []
        posizioni_finali_match = []
    
        for cont in range(0, len(nuove_finestre)):
            for k in range(0, NUM_SEQUENZE):
                for i in range(0, (len(nuove_finestre[cont]) - distanza_frequente)):
                    p_i = i
                    p_f = i + distanza_frequente
                    if match_pattern(sequenze_regioni[cont][k], nuove_finestre[cont][p_i:p_f], distanza_frequente):
                        pos_in.append(p_i)
                        pos_fin.append(p_f)
            posizioni_iniziali_match.append(pos_in)
            posizioni_finali_match.append(pos_fin)
        end = time.time()
        #ft.write("Tempo match tra sequenze e regioni di interesse: " + str(end-start) + "\n")
        #individuazione e gestione delle sovrapposizioni
        contrange = NUMERO_CORE
        executor = concurrent.futures.ThreadPoolExecutor()
        start = time.time()
        cont = range(0,contrange)
        while cont[NUMERO_CORE - 1] < len(nuove_finestre):
            [executor.submit(individua_e_gestisci_sovrapposizioni, k, posizioni_iniziali_match[k], posizioni_finali_match[k], nuove_finestre[k]) for k in cont]
            cont = range(contrange, contrange + NUMERO_CORE)
            contrange = contrange + NUMERO_CORE
        if(cont[0] < (len(nuove_finestre))):
            cont = range(cont[0], len(nuove_finestre))
            [executor.submit(individua_e_gestisci_sovrapposizioni, k, posizioni_iniziali_match[k], posizioni_finali_match[k], nuove_finestre[k]) for k in cont]
        end = time.time()
        #ft.write("Tempo individuazione e gestione overlap: " + str(end-start) + "\n")
        '''