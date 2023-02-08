#Significancia_estadistica_analisis_de_expresion
#Version 2
#Andrea_Ross

import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy.stats import shapiro


#Prueba de normalidad, según el valor p se corre la prueba T de Student para paramétricos o U de Mann-Whitney para no paramétricos
def test_estadistico(df1, df2, columna):
    p1 = shapiro(df1[columna])[1]
    p2 = shapiro(df2[columna])[1]
    es_normal = False
    if p1 > 0.05 and p2 > 0.05:
        t_stat, p_value = stats.ttest_ind(df1[columna], df2[columna])
        print(f'Prueba t de Student para {columna}:')
        print('t-statistic:', t_stat)
        print('p-value:', p_value)
        if p_value <= 0.05:
            print("Este valor es estadísticamente significativo")
        es_normal = True
        print('\n')
    else:
        u_stat, p_value = stats.mannwhitneyu(df1[columna], df2[columna])
        print(f'Prueba U de Mann-Whitney para {columna}:')
        print('U-statistic:', u_stat)
        print('p-value:', p_value)
        if p_value <= 0.05:
            print("Este valor es estadísticamente significativo")
        print('\n')
    return es_normal

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


#Correlación entre los genes o miRNAS de estudio, Pearson para paramétricos y Spearman para no paramétricos
def analisis_correlacion(df, mirnas, test_normalidad):
    mirnas: ['miR141', 'miR145', 'miR146', 'miR148']
    test_normalidad: [test_estadistico(df_experiment, df_control, 'miR141'),
                test_estadistico(df_experiment, df_control, 'miR145'),
                test_estadistico(df_experiment, df_control, 'miR146'),
                test_estadistico(df_experiment, df_control, 'miR148')]
    
    corr_method = 'pearson' if all(test_normalidad) else 'spearman'
    corr_matrix = df[mirnas].corr(method=corr_method)
    print(f'Análisis de correlación ({corr_method}):')
    print(corr_matrix)
    print('\n')
    
    significant_corrs = corr_matrix[(corr_matrix <= 0.05) & (corr_matrix >= 0)].stack().reset_index()
    reported_pairs = []
    for row in significant_corrs.itertuples(index=False):
        miRNA1, miRNA2 = sorted([row.level_0, row.level_1])
        if (miRNA1, miRNA2) not in reported_pairs:
            print(f"{row[0]} and {row[1]} with a correlation of {row[2]:.3f}")
            reported_pairs.append((miRNA1, miRNA2))
            print('\n')

# Imprimir la matriz de correlación
test_normalidad = [True, True, True, True]
analisis_correlacion(df, ['miR141', 'miR145', 'miR146', 'miR148'], test_normalidad)

