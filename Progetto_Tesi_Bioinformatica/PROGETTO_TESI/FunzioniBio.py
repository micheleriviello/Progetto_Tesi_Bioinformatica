import concurrent.futures
import json
import numpy as np
import multiprocessing
import time

'''
---------------------------------------------------------------------
suddivisioneCromosomaFinestre effettua una suddivisione del cromosoma
in finestre composte ciascuna da 1000 basi azotate    
---------------------------------------------------------------------
'''
def suddivisioneCromosomaInFinestre(cromosoma, finestra, NUM_BP_FINESTRA):
    i = 0
    j = NUM_BP_FINESTRA
    while j <= len(cromosoma):
        finestra.append(cromosoma[i:j])
        if j == len(cromosoma):
            break
        i = i + NUM_BP_FINESTRA
        if i < (len(cromosoma) - NUM_BP_FINESTRA):
            j = j + NUM_BP_FINESTRA
        else:
            j = len(cromosoma)
    finestra = np.array(finestra)
'''
---------------------------------------------------------------------
suddivisioneFinestreInMers effettua una suddivisione delle finestre
di ciascun cromosoma in 12-mers composti da 12 basi azotate  
---------------------------------------------------------------------
'''
def suddivisioneFinestreInMers(finestra, finestre_mers, NUM_BASI_MER):
    for cont in range(0, len(finestra)):
        print(cont," SU: ",len(finestra), "1")
        lunghezzaFinestraCorrente = len(finestra[cont]) - NUM_BASI_MER
        mers = []
        for i in range(0, lunghezzaFinestraCorrente + 1):
            if i < lunghezzaFinestraCorrente:
                mers.append(finestra[cont][i : i + NUM_BASI_MER])
            else:
                mers.append(finestra[cont][i : len(finestra[cont])])
        finestre_mers.append(mers)
'''
---------------------------------------------------------------------
calcoloMatch calcola i vari match per ogni kmer di ciascuna finestra 
confrontandoli con il resto della finestra
---------------------------------------------------------------------
'''
def calcoloMatch(finestra, finestre_mers, matches_finestra):
    for cont in range(0, len(finestra)):
        print(cont," SU: ",len(finestra))
        diz = {}
        cont_matches = 0
        for i in range(0, len(finestre_mers[cont])):
            kmer = str(finestre_mers[cont][i])
            if kmer not in diz.keys():
                diz[kmer] = 0
            diz[kmer] = diz[kmer] + 1
        for val in diz.values():
            if val > 1:
                cont_matches += val - 1
        matches_finestra.append(cont_matches)
'''
----------------------------------------------------------------------
ricalcoloFinestre effettua un ricalcolo delle finestre di ciascun 
cromosoma andando a tenere in considerazione quelle che possiedono
un numero di matches superiore al 10% rispetto al numero di mer totali
----------------------------------------------------------------------
'''
def ricalcoloFinestre(cartellaRisultati, nomefile, cromosoma, pos, conta_finestre, visualizzatore_cromosoma, finestra, finestre_mers, matches_finestra, DIVISORE_GRAFICO_RIPETIZIONI, PERCENTUALE, NUM_BP_FINESTRA):
    fl = open(cartellaRisultati+"Test/Risultati_finestre_pre_raggruppamento_"+nomefile+".txt", "w")
    fl.write('numero_finestra, ' + 'posizione_iniziale, ' + 'ripetizioni' + '\n')
    cont = 0
    k = 0
    for _ in range(int(len(finestra)//DIVISORE_GRAFICO_RIPETIZIONI) + 1):
        visualizzatore_cromosoma.append(0)
        conta_finestre.append(0)
    while cont < len(finestra):
        dieci_percento = len(finestre_mers[cont]) * PERCENTUALE
        posizione = int(cont//DIVISORE_GRAFICO_RIPETIZIONI)
        if matches_finestra[cont] > dieci_percento:
            fl.write(str(cont+1) + ', ' + str(k) + ', SI\n')
            visualizzatore_cromosoma[posizione] += 1
            print(posizione,", ",len(finestra))
            conta_finestre[posizione] = posizione
            k = k + NUM_BP_FINESTRA
            pos.append(cont)
        else:
            fl.write(str(cont+1) + ', ' + str(k) + ', NO\n')
            conta_finestre[posizione] = posizione
            k = k + NUM_BP_FINESTRA
        if k >= len(cromosoma):
            break
        cont = cont + 1
    fl.close()
'''
----------------------------------------------------------------------
creazioneNuoveFinestre crea una lista in cui vengono inserite tutte le 
finestre con un numero di ripetizioni maggiore del 10%
----------------------------------------------------------------------
'''
def creazioneNuoveFinestre(cartellaRisultati, nomefile, finestra, pos, nuove_finestre, NUM_BP_FINESTRA):
    entrato = False
    ft = open(cartellaRisultati+"Test/Risultati_finestre_raggruppamento_"+nomefile+".txt", "w")
    ft.write('numero_finestra, ' + 'posizione_iniziale, ' + 'numero_basi\n')
    i = 1
    cont = 0
    pos_precedente = pos[0]
    n_finestra = finestra[pos[0]]
    indice = 1
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
    ft.close()
'''
--------------------------------------------------------------------------
calcoloDistanzeMatch effettua un calcolo delle distanze di tutti i matches
trovati all'interno di ciascuna finestra ricalcolata, salvandone anche il
valore del numero di occorrenze
--------------------------------------------------------------------------
'''
def calcoloDistanzeMatch(cartellaRisultati, nomefile, abilitaFileDebug,ftt, nuove_finestre, NUM_BASI_MER, NUMERO_CORE):
    executor = concurrent.futures.ProcessPoolExecutor()
    cont = 0
    matches = []
    cont_matches = 0

    #with open(cartellaRisultati+"/nuove_finestre_"+nomefile, "r") as fl:
    #    nuove_finestre = json.load(fl)

    #distTemp = open(cartellaRisultati+"/distanze_temp_"+nomefile+".txt", "w")

    start = time.time()
    buffers = {}
    diz_dis_freq = {}
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
        for listaRipetizioni in diz.values():
            ripetizioniFatte += len(listaRipetizioni)
            percentuale = ripetizioniFatte/((len(nuove_finestre[cont]) - NUM_BASI_MER)/100)
            print(stringaStampa," ", percentuale,"%")
            if len(listaRipetizioni) > 2000 and len(listaRipetizioni) < 100000:
                posizione = 2000 * (len(listaRipetizioni)//2000)
                print("STO USANDO IL BUFFER ",posizione,"?")
                if posizione not in buffers.keys():
                    buffers[posizione] = []
                listaRipetizioni = np.array(listaRipetizioni)
                buffers[posizione].append(listaRipetizioni)
                gestisci_buffer(buffers[posizione], executor, diz_dis_freq, NUM_BASI_MER)

            elif len(listaRipetizioni) > 100000:
                buffer = []
                for k in range(NUMERO_CORE):
                    nuovaLista = []
                    buffer.append(nuovaLista)
                for i in range(0, len(listaRipetizioni) - 1, NUMERO_CORE):
                    if i + NUMERO_CORE > len(listaRipetizioni) - 1:
                        for j in range(len(listaRipetizioni) - 1 - i):
                            buffer[j].append(i + j)
                        break
                    for k in range(NUMERO_CORE):
                        buffer[k].append(i + k)
                print("STO ESEGUENDO SU: ", len(listaRipetizioni))
                listaRipetizioni = np.array(listaRipetizioni)
                for i in range(len(buffer)):
                    print(len(buffer[i]))
                    buffer[i] = np.array(buffer[i])
                results = [executor.submit(calcola_match_singolo, listaRipetizioni, k, NUM_BASI_MER) for k in buffer]
                for f in concurrent.futures.as_completed(results):
                    for dist in f.result().keys():
                        if dist not in diz_dis_freq.keys():
                            diz_dis_freq[dist] = 0
                        diz_dis_freq[dist] += f.result()[dist]
                    #distTemp.write(json.dumps(f.result())+"\n")
                buffer.clear()
            else:
                listaRipetizioni = np.array(listaRipetizioni)
                distanze = calcola_distanze_match(listaRipetizioni, NUM_BASI_MER)
                for dist in distanze.keys():
                    if dist not in diz_dis_freq.keys():
                        diz_dis_freq[dist] = 0
                    diz_dis_freq[dist] += distanze[dist]
                '''
                for i in range(len(listaRipetizioni) - 1):
                    for j in range(i, len(listaRipetizioni)):
                        if listaRipetizioni[j] > listaRipetizioni[i] + NUM_BASI_MER:
                            distanza = listaRipetizioni[j] - listaRipetizioni[i]
                            if distanza in temp:
                                temp[distanza] += 1
                            else:
                                temp[distanza] = 1
                '''
                #distTemp.write(json.dumps(temp)+"\n")

    nuoviBuffer = []
    nuovoBuffer = []
    chiaviBuffer = list(buffers.keys())
    i = len(chiaviBuffer) - 1
    while i >= 0:
        for j in buffers[chiaviBuffer[i]]:
            nuovoBuffer.append(j)
            if len(nuovoBuffer) == NUMERO_CORE:
                nuoviBuffer.append(nuovoBuffer)
                nuovoBuffer = []
        i = i - 1
    nuoviBuffer.append(nuovoBuffer)

    del nuovoBuffer
    print(len(nuoviBuffer), "PRIMA: ", len(buffers))
    del buffers

    for buffer in nuoviBuffer:
        print(len(buffer))
        svuota_buffer(buffer, executor, diz_dis_freq, NUM_BASI_MER)

    #distTemp.close()
    #distTemp = open(cartellaRisultati+"/distanze_temp_"+nomefile+".txt", "r")
    '''
    dizionari = distTemp.readlines()

    
    #Prendi dizionario dal disco e metti tutto insieme
    for dizionario in dizionari:
        diz = {}
        diz = json.loads(dizionario)
        for chiave in diz.keys():
            if chiave not in diz_dis_freq:
                diz_dis_freq[chiave] = diz[chiave]
            else:
                diz_dis_freq[chiave] = diz_dis_freq[chiave] + diz[chiave]
    '''

    end = time.time()
    if abilitaFileDebug:
        ftt.write("Tempo calcolo distanze match finestra " + str(cont+1) + ": " + str(end-start) + "\n")
    #ft.write("Tempo calcolo distanze match finestra " + str(cont+1) + ": " + str(end-start) + "\n")


    with open(cartellaRisultati+"SavedFiles/distanza_frequenze"+nomefile, "w") as fl:
        json.dump(diz_dis_freq, fl)

#Da cancellare
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

'''
----------------------------------------------------------------------------------------------
calcola_distanze_match calcola la distanza tra i match di uno stesso k-mer.               
Prende in input una lista contenente le posizioni finali di tutte le occorrenze           
del k-mer presenti all'interno del file contenente il cromosoma.                          
Restituisce in output un dizionario contenente come chiavi le diverse distanze calcolate  
e come valori la frequenza della distanza (il numero di volte in cui è stata trovata)     
----------------------------------------------------------------------------------------------
'''
def calcola_distanze_match(listaRipetizioni, NUM_BASI_MER):
    risultato = {}

    for i in range(1, len(listaRipetizioni)):
        rip = listaRipetizioni[:i]
        rip = rip[listaRipetizioni[i] - NUM_BASI_MER - rip <= 10000]
        distanze = (listaRipetizioni[i] - NUM_BASI_MER) - rip
        distanze = distanze[distanze > 0]
        #distanze = distanze[distanze < 10000]

        for distanza in distanze:
            distanza = int(distanza)
            risultato[distanza] = risultato.get(distanza, 0) + 1

    return risultato


#Da cancellare
def calcola_match_singolo(listaRipetizioni, r, NUM_BASI_MER):
    risultato = {}

    for i in r:
        rip = listaRipetizioni[:i]
        rip = rip[listaRipetizioni[i] - NUM_BASI_MER - rip <= 10000]
        distanze = (listaRipetizioni[i] - NUM_BASI_MER) - rip
        distanze = distanze[distanze > 0]
        #distanze = distanze[distanze < 10000]

        for distanza in distanze:
            distanza = int(distanza)
            if distanza not in risultato.keys():
                risultato[distanza] = 0
            risultato[distanza] += 1
    return risultato


'''
----------------------------------------------------------------------------------------------
gestisci_buffer si occupa della gestione di una lista di liste(buffer) di posizioni finali 
(che sono poi gestite attraverso calcola_match). Prende in input il buffer, l'esecutore per
il multiprocessing (concurrent.futures.ProcessPoolExecutor()), un file in scrittura per scrivere
i risultati su disco (per evitare di riempire la ram) e il numero core della macchina su
cui viene eseguito il codice
----------------------------------------------------------------------------------------------
'''
def gestisci_buffer(buffer, executor, distanze, NUM_BASI_MER):
    numero_core = multiprocessing.cpu_count()
    if(len(buffer) == numero_core):
        print("SI")
        for k in buffer:
            print(len(k))
        results = [executor.submit(calcola_distanze_match, k, NUM_BASI_MER) for k in buffer]
        for f in concurrent.futures.as_completed(results):
            for dist in f.result().keys():
                distanze[dist] = distanze.get(dist, 0) + f.result()[dist]
            '''
            print("RISULTATO")
            for chiave in f.result().keys():
                if chiave not in diz_dis_freq:
                    diz_dis_freq[chiave] = f.result()[chiave]
                else:
                    diz_dis_freq[chiave] = diz_dis_freq[chiave] + f.result()[chiave]
            '''
            #Dump dizionario
        buffer.clear()
'''
----------------------------------------------------------------------------------------------
svuota_buffer si occupa di gestire le ultime liste di posizioni finali di k-mer presenti
all'interno del file del cromosoma che non sono state gestite da gestisci_buffer poichè il 
buffer per poter essere gestito da gestisci_buffer deve essere pieno(il numero di liste
delle posizioni finali deve essere uguale al numero di core della macchina).
Prende in input il buffer, l'esecutore(concurrent.futures.ProcessPoolExecutor()) e un file 
in scrittura per scrivere i risultati su disco (per evitare di riempire la ram)
----------------------------------------------------------------------------------------------
'''
def svuota_buffer(buffer, executor, distanze, NUM_BASI_MER):
    if(len(buffer) > 0):
        results = [executor.submit(calcola_distanze_match, k, NUM_BASI_MER) for k in buffer]
        print("STO USANDO IL BUFFER PER L'ULTIMA VOLTA")
        for f in concurrent.futures.as_completed(results):
            for dist in f.result().keys():
                distanze[dist] = distanze.get(dist, 0) + f.result()[dist]
            '''
            for chiave in f.result().keys():
                if chiave not in diz_dis_freq:
                    diz_dis_freq[chiave] = f.result()[chiave]
                else:
                    diz_dis_freq[chiave] = diz_dis_freq[chiave] + f.result()[chiave]
            '''
            #Dump dizionario
        buffer.clear()

'''
----------------------------------------------------------------------------------------------
individua_e_gestisci_sovrapposizioni: da aggiustare e commentare
----------------------------------------------------------------------------------------------
'''
def individua_e_gestisci_sovrapposizioni(cont, posizioni_iniziali_match, posizioni_finali_match, nuove_finestre):
    match = []
    f1 = open('sequenze_da_allineare'+str(cont+1)+'.fna', "w")
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