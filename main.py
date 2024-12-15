# ___   ___ 
#(__ \ (__ )
# / _/  (_ \
#(____((___/

import json
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy.integrate import odeint

donnees = json.load(open(os.path.dirname(__file__) + "/donnees.json", "r")) #Ouvre le fichier en lecteur, et le convertie en json (similaire au dict)

result = input("\n\nAppuyer sur \"Entrée\" entre chaque étape \n1 - Cas réel \n2 - Maquette\n") # Question
if(result == "1"): # Si réponse est égale à 1 (Cas réel)
    # Données pour le cas réel
    donnees["hpenteplot"]= 20
    donnees["S"]= 1.95*1.35
    donnees["mu"]= 0.01
    donnees["m"]=1760
    donnees["hinit"]=20
    donnees["r"]=6
    donnees["hravin"]=1
    donnees["lravin"]=9
    donnees["tmax"]=10
    donnees["pas"]=1001
    donnees["ravintmax"]=0.6
    donnees["ravinpass"]=200
    donnees["receptionLongeur"]= 20

donnees["cos"] = np.cos(np.radians(donnees["alpha"])) # Cos(40)
donnees["sin"] = np.sin(np.radians(donnees["alpha"])) # Sin(40)
donnees["pSC"] = 1/2*(donnees["p"]*donnees["S"]*donnees["C"]) # Pour simplifier les calculs par la suite
donnees["L"] = donnees["hinit"] / donnees["sin"] # Calcul longeur pente avec les valeur initiale

print("Données du livrable : " + str(donnees))
################################################################################################################
# Pente #
################################################################################################################

# Calcul de la vitesse en bas de la pente sans frottement
v_descente_sf = np.sqrt(donnees["g"]*(donnees["sin"])*2*donnees["L"])
print("Vitesse en bas de la pente sans frottement : " + str(round(v_descente_sf, 3)) + " m.s^-1")

# Calcul de la vitesse en bas de la pente avec frottement
v_descente_af = np.sqrt(donnees["g"]*(donnees["sin"]-donnees["mu"]*donnees["cos"])*2*donnees["L"])
print("Vitesse en bas de la pente avec frottement : " + str(round(v_descente_af, 3)) + " m.s^-1")

h = np.linspace(0, donnees["hpenteplot"], donnees["pas"]) # Création d'une liste contenant toutes les hauteurs de la pente entre 0 et 2 avec un pas de 0.01
l = h/donnees["sin"] # Calcul de la longueur de la pente pour chaque hauteur

# Calcul de la vitesse en bas de la pente en fonction de la hauteur sans frottement
evolution_vitesse_sans_frottement = np.sqrt(donnees["g"]*(donnees["sin"])*2*l)
# Calcul de la vitesse en bas de la pente en fonction de la hauteur avec frottement
evolution_vitesse_avec_frottement = np.sqrt(donnees["g"]*(donnees["sin"]-donnees["mu"]*donnees["cos"])*2*l)

# Matplotlib
fig, (ax1, ax2, ax3) = plt.subplots(1, 3) # Création de 3 courbes sur une figure
ax1.plot(h, evolution_vitesse_sans_frottement, "+") # N°1 : Mettre la vitesse en fonction du temps sans frottement, Avec pour symbole "+"
ax1.set_title("Vitesse en fonction de la hauteur \nsans frottement en m.s^-1") # N°1 : Titre du graph

ax2.plot(h, evolution_vitesse_avec_frottement, "tab:orange") # N°2 : Mettre la vitesse en fonction du temps avec frottement.
ax2.set_title("Vitesse en fonction de la hauteur \navec frottement en m.s^-1") # N°2 : Titre du graph

ax3.plot(h, evolution_vitesse_sans_frottement,  "+") # N°3 : Mettre la vitesse en fonction du temps sans frottement, Avec pour symbole "+"
ax3.plot(h, evolution_vitesse_avec_frottement) # N°3 : Mettre la vitesse en fonction du temps avec frottement.
ax3.set_title("Vitesse en fonction de la hauteur \nsans et avec frottement en m.s^-1") # N°3 : Titre du graph

ax1.set(ylabel="Vitesse en m.s^-1") # Pour avoir la notation uniquement sur le premier
for element in [ax1,ax2,ax3]: # Pour chaque axe
    element.set(xlabel="Hauteur (m)") # Défini des noms aux axes 

plt.show() # Affiche les graphs

input("Fin Pente \n\n") # Pause
################################################################################################################
# Plat #
################################################################################################################

#La phase possede les même equation avant et apres le looping on définie donc une fonction
def plat(v0):
    """
    v0 : float,int : Vitesse initiale au début de la phase.
    return : float, int : vitesse a la fin de la phases
    """
    # Calcul de delta vue dans le livrable 2
    a = (1/2)*(-donnees["mu"]*donnees["g"])
    b = v_descente_sf
    c = -donnees["r"]
    delta = b**2 - 4*a*c

    # Temps à la fin de la phases
    t1= (-b-np.sqrt(delta))/ (2*a) # Temps impossible car négatif
    t2= (-b+np.sqrt(delta))/ (2*a) # Temps a la fin de la pente
    print("Temps pour que la voiture atteigne la fin du plat : " + str(round(t2 ,3)) + " s")

    v_plat_af = (-donnees["mu"]*donnees["g"])*t2 + v0 # Calcul de la vitesse en fonction du temps
    return v_plat_af 

v_plat_af1 = plat(v_descente_af) # Appelle de la fonction "plat" avec comme parametre v0 la vitesse avec frottement durant la descente
print("Vitesse a la fin du plan avec frottement : " + str(round(v_plat_af1 , 3)) + " m.s^-1")
    
input("Fin Plat \n\n") # Pause
################################################################################################################
# Looping  #
################################################################################################################

#Equation de mouvement à résoudre avec 'odeint'
def equation(y, t):
    calcul = (donnees["m"]*donnees["g"]*(-np.sin(y[0])-donnees["mu"]*np.cos(y[0]))-y[1]**2*(donnees["mu"]*donnees["m"]*donnees["r"]+donnees["pSC"]*donnees["r"]**2))/(donnees["m"]*donnees["r"])
    return [y[1], calcul]

#vecteur temps
t = np.linspace(donnees["t0"],donnees["tmax"],donnees["pas"])

#Calcul vitesse angulaire (teta prime) en rad/s^-1
omega=v_plat_af1/donnees["r"]

#Résolution
y = odeint(equation, [0,omega], t)
#Savoir quand est-ce qu'on dépasse 2pi, ce qui correspond à la voiture à la fin du looping
i = 0
while y[i][0] < 2*np.pi:
    i = i+1

v_looping_af = y[i-1][1]*donnees["r"]

print("Vitesse initiale : " + str(round(v_plat_af1, 3)) + " m.s^-1 | Vitesse finale : " + str(round(v_looping_af, 3)) + " m.s^-1")

t = np.linspace(donnees["t0"],4,donnees["pas"]) # Création d'une liste contenant le temps entre 0 et 4 avec 1000 valeur

#Tracer des résultats
#Vitesse de la voiture avec frottements en fonction du temps
plt.plot(t[0:i-1],y[0:i-1,1]*donnees["r"])

# Titre plus nom des axes
plt.title("Vitesse de la voiture avec frottements en fonction du temps")
plt.xlabel("Temps en ms") #  Notation de l'axe x
plt.ylabel("Vitesse en m.s^-1") #  Notation de l'axe y
plt.grid()
plt.show() #Affichage

#Vitesse de décrochage sans frottement
v_decrochage_sf =np.sqrt(5*donnees["g"]*donnees["r"])
print("La voiture décroche si la vitesse initiale est inférieure à (sans frottement) : " + str(round(v_decrochage_sf, 3)) + " m.s^-1")

#Vitesse de décrochage avec frottement
test_current_vitesse = 10 #On suppose que la vitesse limite avec frottement est entre 10 et v_decrochage_sf
Rn = 0 # Initialisation de la variable

while (Rn >= 0 and test_current_vitesse >= v_decrochage_sf/2):
    vteta0 = [0,test_current_vitesse/donnees["r"]] # On prend le teta initial lorsque le radian est égale à 0
    Rn = 0 # Réinitialisation de la variable
    vteta = odeint(equation, vteta0, t) # Utilisation de la fonction odeint 
    for index in range(0, donnees["pas"]): # Boucle avec i allant de 0 à pas
        if(Rn < 0): # Si Rn est insuffisant 
            break # Arrêt de la boucle 
        Rn = donnees["g"]*np.cos(vteta[index][0]) + donnees["r"]*vteta[index][1]**2 # Equation de Rn en fonction de teta
    test_current_vitesse-=0.1
        
print("La voiture décroche si la vitesse initiale est inférieure à (avec frottement) : " + str(round(test_current_vitesse+0.1, 3)) + " m.s^-1")

input("Fin Looping \n\n") # Pause
################################################################################################################
# PLAT  #
################################################################################################################
v_plat_af2 = plat(v_looping_af)
print("Vitesse a la fin du plan avec frottement : " + str(round(v_plat_af2 , 3)) + " m.s^-1")
    
input("Fin Plat \n\n")

################################################################################################################
# Ravin  #
################################################################################################################

# equations de mouvements
def x(t,v0): # fonction de notre équation de position sur x
    return v0 * t
def y(t, g, h0): # fonction de notre équation de position sur y
    return -(1/2)*g*t**2

# tableaux
tab_t = np.linspace(donnees["t0"],donnees["ravintmax"],donnees["ravinpass"])
tab_x = x(tab_t,v_plat_af2)
tab_y = y(tab_t,donnees["g"],donnees["hravin"])

fig, (ax1, ax2, ax3) = plt.subplots(1, 3) # Création de 3 courbes sur une figure
ax1.plot(tab_x,tab_y,'r')
ax1.set_title("Trajectoire de la voiture dans le ravin\nsans frottement") #titre du graphe


t = np.linspace(donnees["t0"],donnees["ravintmax"],donnees["ravinpass"])
#définition de l'équation differentielle
def equation(y,t): # fonction qui calcul la dérivée de y à t
    v_initial_x = y[2] # vitesse initial sur x
    v_initial_y = y[3] # vitesse initial sur y
    # Nos equations de mouvements
    X = -donnees["rho"]/(2*donnees["m"])*donnees["SCx"]*v_plat_af2**2 # En x
    Y = -donnees["g"]+donnees["rho"]/(2*donnees["m"])*donnees["SCz"]*v_plat_af2**2 # En y

    return [v_initial_x,v_initial_y,X,Y]

#Résolution
y = odeint(equation, [0,0,v_plat_af2,0], t)

#Tracé des résultats
ax2.set_title("Trajectoire de la voiture dans le ravin\navec frottement") #titre du graphe
ax2.plot(y[:,0], y[:,1]) #toutes les valeurs de la colonne 0 et 1 (de toutes les lignes).
if(result != "1"):
    ax2.set_ylim(-0.46,0.02)

ax3.set_title("Trajectoire de la voiture dans le ravin\nsans et avec frottement") #titre du graphe
ax3.plot(tab_x,tab_y,'r')
ax3.plot(y[:,0], y[:,1]) #toutes les valeurs de la colonne 0 et 1 (de toutes les lignes).

ax1.set(ylabel="Hauteur en m")
for element in [ax1,ax2,ax3]: # Pour chaque axe 
    element.set(xlabel="Longueur en m") # Défini des noms aux axes
    element.grid()
    element.hlines(-donnees["hravin"],donnees["lravin"],donnees["receptionLongeur"],'r') #ligne horizontal qui est la réception de l'autre côté du ravin

plt.show() #Affichage 
