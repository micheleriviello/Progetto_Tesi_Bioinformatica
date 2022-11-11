import matplotlib.pyplot as plt
import numpy as np
import json
class Istogrammatore:
    '''
    ----------------------------------------------------------------------------------------------
    mostra_istogramma: da aggiustare e commentare
    ----------------------------------------------------------------------------------------------
    '''
    def __init__(self, diz_dis_freq, chiave):
        #distanza_frequente = max(diz_dis_freq, key=diz_dis_freq.get)
        numeroAnnotazioni = 10
        limitex = 10000
        #limitey = 1000000
        x = np.array([int(dist) for dist in diz_dis_freq.keys()])

        numMaxValori = 10
        massimiX = []
        massimiY = []
        #valoreMax = 0
        for i in range(numMaxValori):
            massimiX.append(0)
            massimiY.append(0)
        valoriDizionario = list(diz_dis_freq.values())
        valoriDizionario.sort(reverse=True)
        massimiY = valoriDizionario[0:numMaxValori]
        for i in range(numMaxValori):
            for j in diz_dis_freq.keys():
                if diz_dis_freq[j] == massimiY[i]:
                    massimiX[i] = int(j)
                    break

        limitey = max(massimiY) + (max(massimiY) / 100) * 10

        figura, grafico = plt.subplots(figsize=(9, 7))
        grafico.set_xlim([0, limitex])
        grafico.set_ylim([0, limitey])

        grafico.bar(x=x, height=list(diz_dis_freq.values()), fc="orange", width = 5, snap=False)
        grafico.set_xlabel("Valori distanza")
        grafico.set_ylabel("Valori frequenza")
        grafico.set_title("Istogramma frequenza distanze "+ chiave)

        annotazioni = []
        picchiFissi = {}
        for i in range(numeroAnnotazioni):
            annotation = grafico.annotate(
                text = '',
                xy = (0,0),
                xytext = (-15,15),
                textcoords='offset points',
                bbox = {'boxstyle': 'round', 'fc': 'w'},
                arrowprops= {'arrowstyle': '->'}
            )
            annotation.set_visible(False)
            annotazioni.append(annotation)



        for i in range(numMaxValori):
            if massimiY[i]//1000 > 0:
                valorey = massimiY[i]//1000
                testo = '({:.0f}, {:.0f}k)'
            else:
                valorey = massimiY[i]
                testo = '({:.0f}, {:.0f})'
            if massimiX[i] not in picchiFissi.keys():
                picchiFissi[massimiX[i]] = grafico.annotate(
                    text = testo.format(massimiX[i], valorey),
                    xy = (massimiX[i], massimiY[i]),
                    xytext = (15,15),
                    textcoords='offset points',
                    bbox = {'boxstyle': 'round', 'fc': 'w'},
                    arrowprops= {'arrowstyle': '->'}
                )

        #cursor = Cursor(grafico, horizOn=True, vertOn=True, useblit=True, color = 'r', linewidth= 1)


        def onclick(event):
            raggioAzione = 20
            xcoord = event.xdata
            maxX = 0
            maxY = 0
            if raggioAzione > 0:
                for k in range(int(xcoord) - raggioAzione, int(xcoord) + raggioAzione):
                    if str(k) in diz_dis_freq.keys():
                        if maxY < diz_dis_freq[str(k)]:
                            maxY = diz_dis_freq[str(k)]
                            maxX = k
            else:
                maxY = diz_dis_freq[str(int(xcoord))]
                maxX = int(xcoord)

            if event.button == 1:
                for i in range(numeroAnnotazioni):
                    annotazioni[i].set_visible(False)

                annotazioni[0].xy = (maxX, diz_dis_freq[str(maxX)])
                if diz_dis_freq[str(maxX)]//1000 > 0:
                    valorey = diz_dis_freq[str(maxX)]//1000
                    testo = '({:.0f}, {:.0f}k)'
                else:
                    valorey = diz_dis_freq[str(maxX)]
                    testo = '({:.0f}, {:.0f})'
                text_label = 'P1'+testo.format(maxX, valorey)
                annotazioni[0].set_text(text_label)
                annotazioni[0].set_visible(True)
                for chiave in picchiFissi.keys():
                    picchiFissi[chiave].set_visible(True)
                if maxX in picchiFissi.keys():
                    picchiFissi[maxX].set_visible(False)
            elif event.button == 3:
                picchiTrovati = []
                annotazioneCorrente = 0
                multiploCorrente = maxX
                while multiploCorrente <= 10000 and annotazioneCorrente < numeroAnnotazioni:
                    multiploVicino = multiploCorrente
                    maxY = 0
                    if raggioAzione > 0:
                        for k in range(int(multiploVicino) - raggioAzione, int(multiploVicino) + raggioAzione):
                            if str(k) in diz_dis_freq.keys():
                                if maxY < diz_dis_freq[str(k)]:
                                    maxY = diz_dis_freq[str(k)]
                                    multiploVicino = k

                    if str(multiploVicino) in diz_dis_freq.keys():
                        annotazioni[annotazioneCorrente].xy = (multiploVicino, diz_dis_freq[str(multiploVicino)])
                        if diz_dis_freq[str(multiploVicino)]//1000 > 0:
                            valorey = diz_dis_freq[str(multiploVicino)]//1000
                            testo = '({:.0f}, {:.0f}k)'
                        else:
                            valorey = diz_dis_freq[str(multiploVicino)]
                            testo = '({:.0f}, {:.0f})'
                        text_label = 'P'+str(annotazioneCorrente + 1)+testo.format(multiploVicino, valorey)
                        annotazioni[annotazioneCorrente].set_text(text_label)
                        annotazioni[annotazioneCorrente].set_visible(True)
                        picchiTrovati.append(multiploVicino)
                    else:
                        annotazioni[annotazioneCorrente].set_visible(False)
                    annotazioneCorrente += 1
                    multiploCorrente = multiploCorrente + maxX
                for chiave in picchiFissi.keys():
                    picchiFissi[chiave].set_visible(True)
                for multiplo in picchiTrovati:
                    if multiplo in picchiFissi.keys():
                        picchiFissi[multiplo].set_visible(False)
                while annotazioneCorrente < numeroAnnotazioni:
                    annotazioni[annotazioneCorrente].set_visible(False)
                    annotazioneCorrente += 1
            figura.canvas.draw()

        figura.canvas.mpl_connect('button_press_event', onclick)
        figura.show()