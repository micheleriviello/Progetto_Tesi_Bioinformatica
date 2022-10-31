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

        figura, grafico = plt.subplots(figsize=(15, 10))
        grafico.set_xlim([0, limitex])
        grafico.set_ylim([0, limitey])

        grafico.bar(x=x, height=list(diz_dis_freq.values()), fc="orange", width = 5, snap=False)
        grafico.set_xlabel("Valori distanza")
        grafico.set_ylabel("Valori frequenza")
        grafico.set_title("Istogramma frequenza distanze "+ chiave)

        annotation = grafico.annotate(
            text = '',
            xy = (0,0),
            xytext = (15,15),
            textcoords='offset points',
            bbox = {'boxstyle': 'round', 'fc': 'w'},
            arrowprops= {'arrowstyle': '->'}
        )
        annotation.set_visible(False)


        for i in range(numMaxValori):
            grafico.annotate(
                text = '({:.0f}, {:.0f}k)'.format(massimiX[i], massimiY[i]//1000),
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
            annotation.xy = (maxX,maxY)
            text_label = '({:.0f}, {:.0f}k)'.format(maxX, maxY//1000)
            annotation.set_text(text_label)
            annotation.set_visible(True)
            figura.canvas.draw()

        figura.canvas.mpl_connect('button_press_event', onclick)
        figura.show()