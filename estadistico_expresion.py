#Significancia_estadistica_analisis_de_expresion
#Version 1
#Andrea_Ross

import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy.stats import shapiro


#Prueba de normalidad, según el valor p se corre la prueba T de Student para paramétricos o U de Mann-Whitney para no paramétricos
def test_estadistico(df1, df2, columna):
    p1 = shapiro(df1[columna])[1]
    p2 = shapiro(df2[columna])[1]
    if p1 > 0.05 and p2 > 0.05:
        t_stat, p_value = stats.ttest_ind(df1[columna], df2[columna])
        print(f'Prueba t de Student para {columna}:')
        print('t-statistic:', t_stat)
        print('p-value:', p_value)
        if p_value <= 0.05:
            print("Este valor es estadísticamente significativo")
        print('\n')
    else:
        u_stat, p_value = stats.mannwhitneyu(df1[columna], df2[columna])
        print(f'Prueba U de Mann-Whitney para {columna}:')
        print('U-statistic:', u_stat)
        print('p-value:', p_value)
        if p_value <= 0.05:
            print("Este valor es estadísticamente significativo")
        print('\n')

# Imprimir en la terminal los valores de significancia
df = pd.read_csv('/Users/andreaross/Desktop/EPCM.csv')
df_experiment = df[df['group'] == 'experiment']
df_control = df[df['group'] == 'control']

print('\n')
print('Resultados de las pruebas estadísticas:')
test_estadistico(df_experiment, df_control, 'miR141')
test_estadistico(df_experiment, df_control, 'miR145')
test_estadistico(df_experiment, df_control, 'miR146')
test_estadistico(df_experiment, df_control, 'miR148')
