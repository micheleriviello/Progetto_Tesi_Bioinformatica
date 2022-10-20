import matplotlib.pyplot as plt
from pathlib import Path
#import random
import re
import numpy as np
#import math
import multiprocessing
#from ClasseRisultato import RisultatoCalcoloMatch
import time
#import _pickle as pickle
import json
import FunzioniBio

NUMERO_CORE = multiprocessing.cpu_count()
LIMITE_SINGLE_CORE = 10000
NUM_BP_FINESTRA = 1000
NUM_SEQUENZE = 5
NUM_BASI_MER = 12
PERCENTUALE = 0.1
DIVISORE_GRAFICO_RIPETIZIONI = 1000
MAX_LIMITE_X = 10000

if __name__ == "__main__":
    limitex = 750
    limitey = 250000
    if limitex > MAX_LIMITE_X:
        limitex = MAX_LIMITE_X
    cromosomi = []
    for i in range(24):
        cromosomi.append("cromosoma"+str(i+1))
    nuoveFinestreCalcolate = False
    distanzeCalcolate = True
    graficoFatto = False
    voglioGraficoIllimitato = True
    voglioGraficoLimitato = True
    abilitaFileDebug = False
    saltaGrafico = False
    for nomefile in cromosomi:
        cartellaRisultati = "Risultati/"+nomefile+"/"
        Path(cartellaRisultati).mkdir(parents=True, exist_ok=True)
        Path(cartellaRisultati+"/Test").mkdir(parents=True, exist_ok=True)
        Path(cartellaRisultati+"/SavedFiles").mkdir(parents=True, exist_ok=True)
        if abilitaFileDebug:
            ftt = open(cartellaRisultati+'Test/filetempi_overlap_multithread'+nomefile+'.txt', 'w')
        if(not nuoveFinestreCalcolate):
            visualizzatore_cromosoma = []
            f = open('Genoma_Umano/'+nomefile+'.fna')
            app = f.readline()
            cromosoma = f.read() #cromosoma senza la prima riga
            cromosoma = re.sub("\n","",cromosoma)

            finestra = []
            start = time.time()
            print(len(cromosoma))
            print("Sto qua 1")
            #suddivisione del cromosoma in finestre da 1kbp = 1000 basi azotate
            FunzioniBio.suddivisioneCromosomaInFinestre(cromosoma, finestra, NUM_BP_FINESTRA)
            print(len(finestra))
            end = time.time()
            if abilitaFileDebug:
                ftt.write("Tempo suddivisione finestre: " + str(end-start) + "\n")

            # suddivisione delle finestre in 12-mers
            start = time.time()
            print("Sto qua 2")
            finestre_mers = []
            FunzioniBio.suddivisioneFinestreInMers(finestra, finestre_mers, NUM_BASI_MER)
            end = time.time()
            if abilitaFileDebug:
                ftt.write("Tempo suddivisione finestre in 12-mers: " + str(end-start) + "\n")

            start = time.time()
            matches = []
            matches_finestra = []
            print("Sto qua 3")
            #calcolo matches dei mers per ogni finestra
            FunzioniBio.calcoloMatch(finestra, finestre_mers, matches_finestra)
            end = time.time()
            if abilitaFileDebug:
                ftt.write("Tempo matches finestre: " + str(end-start) + "\n")

            #ricalcolo finestre in base al numero di mathces
            start = time.time()
            print("Sto qua 4")
            pos = []
            conta_finestre = []
            FunzioniBio.ricalcoloFinestre(cartellaRisultati, nomefile, cromosoma, pos, conta_finestre, visualizzatore_cromosoma, finestra, finestre_mers, matches_finestra, DIVISORE_GRAFICO_RIPETIZIONI, PERCENTUALE, NUM_BP_FINESTRA)

            with open(cartellaRisultati+"SavedFiles/visualizzatore_"+nomefile, "w") as fl:
                json.dump(visualizzatore_cromosoma, fl)


            #creazione delle nuove finestre
            nuove_finestre = []
            FunzioniBio.creazioneNuoveFinestre(cartellaRisultati, nomefile, finestra, pos, nuove_finestre, NUM_BP_FINESTRA)

            with open(cartellaRisultati+"SavedFiles/nuove_finestre_"+nomefile, "w") as fl:
                json.dump(nuove_finestre, fl)
            end = time.time()
            if abilitaFileDebug:
                ftt.write("Tempo ricalcolo finestre: " + str(end-start) + "\n")
            #del mers
            del finestre_mers
            del finestra
            del matches_finestra
            del visualizzatore_cromosoma

        if not saltaGrafico:
            with open(cartellaRisultati+"SavedFiles/visualizzatore_"+nomefile, "r") as fl:
                visualizzatore_cromosoma = json.load(fl)

            for i in range(len(visualizzatore_cromosoma)):
                if visualizzatore_cromosoma[i] > DIVISORE_GRAFICO_RIPETIZIONI/8:
                    visualizzatore_cromosoma[i] = 1
                else:
                    visualizzatore_cromosoma[i] = 0

            figura, grafico = plt.subplots()
            figura.set_figwidth(12)
            figura.set_figheight(3)
            grafico.set_ylim([0, 1])
            figura.gca().axes.get_yaxis().set_visible(False)
            plt.figure(dpi=300)
            grafico.bar(x=np.arange(len(visualizzatore_cromosoma)), height=visualizzatore_cromosoma, fc="orange")
            figura.savefig(cartellaRisultati+'centromero_'+nomefile+'.jpg')
            figura.show()

        if(not distanzeCalcolate):
            #calcolo distanze tra i matches dei vari mers
            if nuoveFinestreCalcolate:
                with open(cartellaRisultati+"SavedFiles/nuove_finestre_"+nomefile, "r") as fl:
                    nuove_finestre = json.load(fl)

            FunzioniBio.calcoloDistanzeMatch(cartellaRisultati, nomefile, abilitaFileDebug,ftt, nuove_finestre, NUM_BASI_MER, NUMERO_CORE)


        #creazione istogramma similarita" media di ciscuna finestra
        if not graficoFatto:
            if not saltaGrafico:
                with open(cartellaRisultati+"SavedFiles/distanza_frequenze"+nomefile, "r") as fl:
                    diz_dis_freq = json.load(fl)
                if voglioGraficoLimitato:
                    distanza_frequente = max(diz_dis_freq, key=diz_dis_freq.get)

                    if limitey > diz_dis_freq[distanza_frequente]:
                        limitey = diz_dis_freq[distanza_frequente]

                    x = np.array([int(dist) for dist in diz_dis_freq.keys()])
                    yaxis = [int(dist) for dist in diz_dis_freq.values()]

                    x = x[x <= limitex]
                    yaxis = yaxis[:limitex]
                    y = []
                    for val in yaxis:
                        if val > limitey:
                            val = limitey
                        y.append(val)

                    figura, grafico = plt.subplots()
                    grafico.set_xlim([0, limitex])
                    grafico.set_ylim([0, limitey])

                    plt.figure(dpi=1200)
                    grafico.bar(x=x, height=y, fc="orange")
                    grafico.set_xlabel("Valori distanza")
                    grafico.set_ylabel("Valori frequenza")
                    grafico.set_title("Istogramma frequenza distanze")
                    #figura.text(int(distanza_frequente) + 12, int(diz_dis_freq[distanza_frequente]) + 12, "Distanza con frequenza maggiore: " + str(distanza_frequente))
                    #cerchio = Circle((int(distanza_frequente), int(diz_dis_freq[distanza_frequente])), 3, color="red", lw=3, fill= False, zorder=10)
                    #figura.gca().add_patch(cerchio)
                    figura.savefig(cartellaRisultati+'istogramma_parziale_'+nomefile+'.jpg')
                    figura.show()

                if voglioGraficoIllimitato:
                    x = [int(dist) for dist in diz_dis_freq.keys()]
                    y = [int(dist) for dist in diz_dis_freq.values()]
                    figura, grafico = plt.subplots()
                    distanza_frequente = max(diz_dis_freq, key=diz_dis_freq.get)
                    plt.figure(dpi=1200)
                    grafico.bar(x=x, height=y, fc="orange")
                    grafico.set_xlabel("Valori distanza")
                    grafico.set_ylabel("Valori frequenza")
                    grafico.set_title("Istogramma frequenza distanze")
                    #figura.text(int(distanza_frequente) + 12, int(diz_dis_freq[distanza_frequente]) + 12, "Distanza con frequenza maggiore: " + str(distanza_frequente))
                    #cerchio = Circle((int(distanza_frequente), int(diz_dis_freq[distanza_frequente])), 3, color="red", lw=3, fill= False, zorder=10)
                    #figura.gca().add_patch(cerchio)
                    figura.savefig(cartellaRisultati+'istogramma_'+nomefile+'.jpg')
                    figura.show()


            #trovare distanza più frequente
            distanza_frequente = max(diz_dis_freq, key=diz_dis_freq.get)
            print("Distanza con frequenza maggiore: ", distanza_frequente)
            print("Frequenza della distanza: ", diz_dis_freq[distanza_frequente])
            distanza_frequente = int(distanza_frequente)

        '''
        #generazione 5 sequenze random di lunghezza N pari alla distanza più frequente
        i = 0
        sequenze = []
        sequenze_regioni = []
    
        for cont in range(0, len(nuove_finestre)):
            sequenze = []
            print("STO QUA VITO: " + str(cont) + " su: " + str(len(nuove_finestre)))
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
            print("STO QUA MICHELE: " + str(cont) + " su: " + str(len(nuove_finestre)))
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
            print("STO FINENDO")
            [executor.submit(individua_e_gestisci_sovrapposizioni, k, posizioni_iniziali_match[k], posizioni_finali_match[k], nuove_finestre[k]) for k in cont]
            cont = range(contrange, contrange + NUMERO_CORE)
            contrange = contrange + NUMERO_CORE
        if(cont[0] < (len(nuove_finestre))):
            cont = range(cont[0], len(nuove_finestre))
            [executor.submit(individua_e_gestisci_sovrapposizioni, k, posizioni_iniziali_match[k], posizioni_finali_match[k], nuove_finestre[k]) for k in cont]
        end = time.time()
        #ft.write("Tempo individuazione e gestione overlap: " + str(end-start) + "\n")
        '''