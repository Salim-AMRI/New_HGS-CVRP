import os
import matplotlib.pyplot as plt
import numpy as np

# Les données
data = {
    "OX": [
7.279189041	,
7.257250714	,
7.240750714	,
7.209875714	,
7.184806429	,
7.150907857	,
7.124705		
],
    "AOX": [
6.170187599	,
6.148381532	,
6.122922482	,
6.118788814	,
6.099327014	,
6.09610634	,
6.085503129	
        ],
    "GOX_Best": [
6.112340141	,
6.089653086	,
6.065670829	,
6.030568943	,
6.013837921	,
5.994927557	,
5.9867153	
        ]
}



# Création du dossier Fig s'il n'existe pas
output_folder = "Fig"
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Temps en secondes pour chaque point (par exemple, chaque point correspond à 20 secondes)
time = np.arange(0, len(data["OX"]) * 20, 20)  # par pas de 20

# Création du graphique
plt.figure(figsize=(10, 6))
plt.plot(time, data["OX"], label="OX")
plt.plot(time, data["AOX"], label="AOX")
plt.plot(time, data["GOX_Best"], label="GOX_Best")

# Personnalisation du graphique
plt.xlabel("Time (s)")  # Chaque ligne correspond à 20 sec
plt.ylabel("Avg result")  # Chaque colonne donne pour chaque algorithme l'évolution du résultat
plt.legend()
plt.grid(True)

# Enregistrement de la figure
output_path = os.path.join(output_folder, "Fig.png")
plt.savefig(output_path, dpi=300)
plt.show()

print(f"Figure enregistrée dans le dossier : {output_path}")
 
