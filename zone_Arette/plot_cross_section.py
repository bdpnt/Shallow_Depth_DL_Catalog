import pygmt
import numpy as np
from pygmt.datasets import load_earth_relief
from datetime import datetime
from math import cos, sin, radians
import os
import xarray as xr
from scipy.ndimage import gaussian_filter

# Parameters
use_err = 'erh' # erv or erh

fichier_seisme = "RESULT/GLOBAL_W_200.txt" # after NLL
FORMAT_fichier = 1 # sortie NLL

# fichier_seisme = "../obs/GLOBAL_W.obs" # before NLL
# FORMAT_fichier = 4 # bulletin OBS

stations_file = "stations/GTSRCE_W_200.txt"

save_file = "RESULT/FIGURES/cross_section/arette_W_200_erh.pdf"

# ---------------------------
# Paramètres de la coupe
# ---------------------------
#lon0, lat0 = -0.8  , 43.08     # perpendiculaire
#lon0,lat0 = -0.6275 + 0.56, 42.95 # pour comp mathieu
lon0,lat0 = -0.6275, 43 # pour comp mathieu

azimut = 0                 # orientation de la coupe en degrés (0 = Nord, 90 = Est)
longueur_coupe = 16          # km
largeur_coupe = 8            # km
prof_coupe = 18             # km
prof_min, prof_max = 0, 15   # km
UNCERT_max_H = 1.5
UNCERT_max_V = 1.5

# ---------------------------
# Fonctions utilitaires
# ---------------------------
R = 111.0  # km par degré approximatif
def dest_point(lon, lat, azimut, dist_km):
    """Calcule un point à partir d'un point (lon,lat), d'un azimut (°N) et d'une distance en km.
       Approximation locale (1° ~ 111 km)."""
    az = radians(azimut)
    dlat = (dist_km * cos(az)) / R
    dlon = (dist_km * sin(az)) / (R * cos(radians(lat)))
    return lon + dlon, lat + dlat

# Calcul des points extrémités de la coupe
lon1, lat1 = lon0, lat0
lon2, lat2 = dest_point(lon0, lat0, azimut, longueur_coupe)

# Calcul des points extrémités de la coupe
lon1a, lat1a = lon0 - largeur_coupe/(R * cos(radians(lat0))), lat0
lon2a, lat2a = dest_point(lon0 - largeur_coupe/(R * cos(radians(lat0))), lat0, azimut, longueur_coupe)

# Calcul des points extrémités de la coupe
lon1b, lat1b = lon0 + largeur_coupe/(R * cos(radians(lat0))), lat0
lon2b, lat2b = dest_point(lon0 + largeur_coupe/(R * cos(radians(lat0))), lat0, azimut, longueur_coupe)

# Région affichée
Region = [lon1a-0.2, lon2b+0.2, lat1a-0.2, lat2b+0.2]
#Region = [-0.75,-0.5,43 ,43.15] # pour comp mathieu
#Region = [-0.9,0.,42.5,43.25]

# ---------------------------
# Chargement du catalogue
# ---------------------------
if FORMAT_fichier == 1:
    data = np.loadtxt(fichier_seisme)
    lon = data[:, 7]
    lat = data[:, 6]
    depth = data[:, 8]
    erv = data[:, 13]
    erh = data[:, 12]
    rms = data[:, 10]
    gap = data[:, 14]
    nbphase = data[:, 11]
    year = [i + 2000 if i < 75 else i + 1900 for i in data[:, 0]]
    month = data[:, 1]
    day = data[:, 2]
    hour = data[:, 3]
    minu = data[:, 4]
    
elif FORMAT_fichier == 2:
    data = np.loadtxt(fichier_seisme)
    lon = data[:, 7]
    lat = data[:, 6]
    depth = data[:, 8]
    erv = np.zeros_like(depth) + UNCERT_max_V - 0.01
    erh = np.zeros_like(depth) + UNCERT_max_H - 0.01
    rms = np.zeros_like(depth)
    gap = data[:, 13]
    nbphase = data[:, 10]
    year = data[:, 0]
    month = data[:, 1]
    day = data[:, 2]
    hour = data[:, 3]
    minu = data[:, 4]

elif FORMAT_fichier == 3:
    data = np.loadtxt(fichier_seisme)
    lon = data[:, 7]
    lat = data[:, 6]
    depth = data[:, 8] 
    erv = np.zeros_like(depth) + UNCERT_max_V - 0.01
    erh = np.zeros_like(depth) + UNCERT_max_H - 0.01
    rms = np.zeros_like(depth)
    mag = data[:, 9]
    year = data[:, 0]
    month = data[:, 1]
    day = data[:, 2]
    hour = data[:, 3]
    minu = data[:, 4]

else:
    with open(fichier_seisme, 'r') as f:
        lines = f.readlines()

    # Count the number of relevant lines
    num_lines = sum(1 for line in lines if line.startswith('# '))

    # Initialize all arrays with zeros
    year = np.zeros(num_lines)
    month = np.zeros(num_lines)
    day = np.zeros(num_lines)
    hour = np.zeros(num_lines)
    minu = np.zeros(num_lines)
    lat = np.zeros(num_lines)
    lon = np.zeros(num_lines)
    depth = np.zeros(num_lines)
    erv = np.full(num_lines, -1)
    erh = np.full(num_lines, -1)
    rms = np.full(num_lines, -1)
    gap = np.full(num_lines, -1)
    nbphase = np.full(num_lines, -1)

    # Fill the arrays
    idx = 0
    for line in lines:
        if line.startswith('# '):
            data = line.rstrip('\n').lstrip('# ').split()
            year[idx] = float(data[0])
            month[idx] = float(data[1])
            day[idx] = float(data[2])
            hour[idx] = float(data[3])
            minu[idx] = float(data[4])
            lat[idx] = float(data[6])
            lon[idx] = float(data[7])
            depth[idx] = float(data[8])
            if data[12] != 'None':
                nbphase[idx] = float(data[12])
            idx += 1

    
# Filtrage
#mask = (erv < UNCERT_max_V) & (erh < UNCERT_max_H) & (rms < 0.5) & (gap < 180) & (nbphase > 6)

#mask = (erv < UNCERT_max_V) & (erh < UNCERT_max_H) & (rms < 0.2) & (gap < 180) & (nbphase > 9)

mask = (erv < UNCERT_max_V) & (erh < UNCERT_max_H) & (rms < 0.5)

lon = lon[mask]
lat = lat[mask]
depth = depth[mask]
erv = erv[mask]
erh = erh[mask]

# ---------------------------
# Figure planimétrique
# ---------------------------
fig = pygmt.Figure()
fig.basemap(region = Region, projection="M6i", frame="a")
fig.coast(shorelines=True, water="lightblue", land="lightgray", resolution="h")

# Relief
grid = load_earth_relief("03s", region=Region)
fig.grdimage(grid=-grid, cmap="gray")

# Failles
fig.plot("../failles/FNP.dat", pen="1.25p", style="f1c/0.25c", fill="black")
fig.plot("../failles/structures_lacan.dat", pen="1.25p", style="f1c/0.25c", fill="black")
# fig.plot("../failles/failles_neotectonic.xy", pen="1.25p", style="f1c/0.25c", fill="red")
fig.plot("../failles/lacan.thrust", pen="1.25p", style="f1c/0.25c", fill="blue")
fig.plot("../failles/lacan.other", pen="1.25p", style="f1c/0.25c", fill="blue")


# Trace de la coupe
fig.plot(x=[lon1, lon2], y=[lat1, lat2], pen="2p,red")
fig.plot(x=[lon1a, lon2a], y=[lat1a, lat2a], pen="0.5p,red")
fig.plot(x=[lon1b, lon2b], y=[lat1b, lat2b], pen="0.5p,red")

#--- Ajout des stations depuis le fichier stations/GTSRCE.txt ---
stations = []
with open(stations_file, "r") as f:
    for line in f:
        if not line.strip():
            continue
        parts = line.split()
        if len(parts) >= 7:
            stas = parts[1]      # code station
            lats = float(parts[3])
            lons = float(parts[4])
            stations.append((stas, lats, lons))

#Tracer les stations sur la carte
for stas, lats, lons in stations:
    fig.plot(
        x=lons,
        y=lats,
        style="t0.5c",     # triangle de 0.4 cm
        fill="red",        # couleur de remplissage
        pen="1p,black",    # contour noir
    )
    fig.text(
        x=lons + 0.01,
        y=lats + 0.005,
        text=stas,
        font="10p,Helvetica-Bold",
        justify="LM",
    )

#--- Ajout des villes comme carrés ---
villes = {
    "Arette": ( -0.717, 43.096 ),
    "Sarrance": ( -0.6008333, 43.0522 ),
    "Oloron-Sainte-Marie": ( -0.6056, 43.1947 )
}

for name, (lon_v, lat_v) in villes.items():
    fig.plot(
        x=lon_v,
        y=lat_v,
        style="s0.4c",     # carré de 0.4 cm
        fill="yellow",     # couleur de remplissage (modifiable)
        pen="1p,black"     # contour noir
    )
    fig.text(
        x=lon_v + 0.01,
        y=lat_v + 0.01,
        text=name,
        font="10p,Helvetica-Bold",
        justify="LM"
    )

#---------------------------
#Séismes colorés par profondeur
#---------------------------
pygmt.makecpt(cmap="viridis", series=[prof_min, prof_max], reverse=True)

# Tracer les points
fig.plot(
    x=lon,
    y=lat,
    style="c0.3c",
    fill=depth,
    cmap=True,
    pen="black"
)

fig.colorbar(frame="af+lProfondeur (km)")

# ---------------------------
# Coupe en profondeur
# ---------------------------
data_cat = np.column_stack((lon, lat, depth, erv, erh))

pygmt.project(
    data=data_cat,
    center=[lon1, lat1],
    endpoint=[lon2, lat2],
    width=[-largeur_coupe, largeur_coupe],
    convention="pz",
    unit=True,
    outfile="cross.dat",
    output_type='file',
)

cross_file = "cross.dat"

# Vérification de la sortie pygmt.project
if not os.path.exists(cross_file) or os.path.getsize(cross_file) == 0:
    print("⚠️ Aucun séisme projeté sur la coupe (cross.dat vide).")
    print("👉 On continue avec seulement la carte planimétrique.")
    plot_coupe = False
else:
    try:
        data = np.loadtxt(cross_file, ndmin=2)
    except Exception as e:
        print(f"⚠️ Erreur de lecture {cross_file} : {e}")
        plot_coupe = False
    else:
        if data.size == 0:
            print("⚠️ Fichier cross.dat vide après lecture.")
            plot_coupe = False
        else:
            plot_coupe = True
            X = data[:, 0]   # distance (km)
            Z = data[:, 1]   # profondeur (km)
            erv = data[:, 2]
            erh = data[:, 3]
            err = np.sqrt(erv ** 2 + erh ** 2)

if plot_coupe == True:
        # Sort data on err
        if FORMAT_fichier != 4 :
            if use_err == 'erv':
                data = np.column_stack((X, Z, erv))

                # Sort by error (descending) so the smallest errors are plotted last
                sorted_data = data[np.argsort(erv)]
                sorted_data = sorted_data[::-1]

                # Unpack the sorted data
                X_sorted, Z_sorted, err_sorted = sorted_data.T
            else:
                data = np.column_stack((X, Z, erh))

                # Sort by error (descending) so the smallest errors are plotted last
                sorted_data = data[np.argsort(erh)]
                sorted_data = sorted_data[::-1]

                # Unpack the sorted data
                X_sorted, Z_sorted, err_sorted = sorted_data.T
        
        fig.shift_origin(yshift="-10c")
        fig.basemap(
                projection="X10/-7",
                region=[0, longueur_coupe, -1, prof_coupe],
                frame=['xafg100+lDistance (km)', 'yafg50+lDepth (km)', "WSen"],
        )

        if FORMAT_fichier != 4 :
            if use_err == 'erv':
                pygmt.makecpt(cmap="magma", series=[0, UNCERT_max_V], reverse=True)
            else:
                pygmt.makecpt(cmap="magma", series=[0, UNCERT_max_H], reverse=True)
            
            fig.plot(
                x=X_sorted, 
                y=Z_sorted, 
                style="c0.15c", 
                fill=err_sorted, 
                cmap=True, 
                pen="0.25p,black", 
                # transparency=30,
            )

            if use_err == 'erv':
                fig.colorbar(
                    frame=['af+lErV'],
                    position="JMR+w5c/0.5c+o0.5c/0c"
                )
            else:
                fig.colorbar(
                    frame=['af+lErH'],
                    position="JMR+w5c/0.5c+o0.5c/0c"
                )

        else:
            fig.plot(
                x=X, 
                y=Z, 
                style="c0.15c", 
                fill="#DF2B2B",
                pen="0.25p,black", 
                # transparency=30,
            )

        
        #dx, dz = 0.05, 0.05  # km

        # Calcul de la densité via un histogramme 2D
        #H, xedges, zedges = np.histogram2d(X, Z,
        #                                   bins=(np.arange(0, longueur_coupe+dx, dx),
        #                                         np.arange(-1, prof_coupe+dz, dz)))

        # Moyenne ou normalisation optionnelle
        #H = H.T  # transpose pour correspondre à la convention (z, x)
        # H = matrice 2D de densité (X en axe 0, Z en axe 1)
        #sigma_x = 5   # lissage horizontal (en nb de cellules)
        #sigma_z = 3   # lissage vertical (en nb de cellules)
        #H_smooth = gaussian_filter(H, sigma=(sigma_x, sigma_z))
        #H_smooth = gaussian_filter(H, sigma=5)  # augmente sigma pour plus de lissage
        #H = np.log1p(H_smooth)  # échelle logarithmique douce
        #H[H < 0.01] = np.nan
        
        #H[H < 0.0002] = np.nan
        # Création d’un DataArray PyGMT (xarray)
        #grid = xr.DataArray(
        #    H,
        #    coords=[("z", (zedges[:-1] + zedges[1:]) / 2),
        #            ("x", (xedges[:-1] + xedges[1:]) / 2)]
        #)

        # --- Affichage de la densité sous la coupe ---
        #fig.shift_origin(yshift="-12c")
        #fig.basemap(
        #    projection="X10/-10",
        #    region=[0, longueur_coupe, -5, prof_coupe],
        #    frame=['xafg100+lDistance (km)', 'yafg50+lDepth (km)', "WSen"],
        #)
 
        #fig.grdimage(grid=grid, cmap="hot", shading=False, nan_transparent=True)
        
        #fig.plot(x=X, y=Z, style="c0.2c", pen="black")
        #fig.colorbar(frame='af+l"Densité (log)"', position="JBC+w8c/0.4c+o0/-0.5c")
        
# ---------------------------
# Export
# ---------------------------
fig.savefig(save_file)
