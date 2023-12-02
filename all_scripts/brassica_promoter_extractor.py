import pandas as pd

brapa_tab = pd.read_table('/Volumes/sesame/joerecovery/genomes/brapa/Brassica_rapa.Brapa_1.0.55.chr.gff3', skiprows=17)
print(brapa_tab.head())