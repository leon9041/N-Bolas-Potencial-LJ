import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.image as mpimg

# ==========================================
# CONFIGURACIÓN
# ==========================================
RESULTS_DIR = "results"
FILE_POS = os.path.join(RESULTS_DIR, "salida.dat")
FILE_DAT = os.path.join(RESULTS_DIR, "datos_graficas.dat")
IMG_SPRITE = "scripts/pokemini.png" # Asegúrate de tener esta imagen o cambiar el nombre

# Crear carpeta si no existe
os.makedirs(RESULTS_DIR, exist_ok=True)

print(">>> Procesando gráficas y animación con Python <<<")

# ==========================================
# 1. GRÁFICAS DE ENERGÍA Y PRESIÓN
# ==========================================
print(f"-> Leyendo {FILE_DAT}...")
try:
    # Leer datos saltando el header
    # Formato esperado: tiempo, E_cin, E_pot, E_tot, Presion
    data_metrics = np.loadtxt(FILE_DAT, skiprows=1)
    
    if data_metrics.ndim == 1:
        # Si solo hay una linea de datos, la convertimos a matriz 2D
        data_metrics = data_metrics.reshape(1, -1)

    t = data_metrics[:, 0]
    e_cin = data_metrics[:, 1]
    e_pot = data_metrics[:, 2]
    e_tot = data_metrics[:, 3]
    presion = data_metrics[:, 4]

    # Crear figura con 2 subplots compartiendo eje X
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # --- Subplot 1: Energías ---
    ax1.plot(t, e_cin, label='Cinética', color='blue', alpha=0.7, linewidth=1)
    ax1.plot(t, e_pot, label='Potencial (LJ)', color='red', alpha=0.7, linewidth=1)
    ax1.plot(t, e_tot, label='Total (Conservada)', color='black', linestyle='--', linewidth=1.5)
    ax1.set_ylabel('Energía [u.r.]')
    ax1.set_title('Conservación de Energía')
    ax1.legend(loc='upper right')
    ax1.grid(True, alpha=0.5)

    # --- Subplot 2: Presión ---
    ax2.plot(t, presion, label='Presión (Virial)', color='green', linewidth=1)
    ax2.set_ylabel('Presión P*')
    ax2.set_xlabel('Tiempo t*')
    ax2.set_title('Fluctuación de la Presión')
    ax2.legend(loc='upper right')
    ax2.grid(True, alpha=0.5)

    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "energias_presion.png"), dpi=150)
    plt.close()
    print("   [OK] Gráfica guardada: energias_presion.png")

except Exception as e:
    print(f"   [ERROR] No se pudieron graficar energías/presión: {e}")


# ==========================================
# 2. CARGA DE DATOS DE POSICIÓN
# ==========================================
print(f"-> Leyendo {FILE_POS} (esto puede tardar)...")
try:
    raw_data = np.loadtxt(FILE_POS)
    
    # Manejo si solo hay 1 frame
    if raw_data.ndim == 1:
        raw_data = raw_data.reshape(1, -1)

    times = raw_data[:, 0]
    # Resto de columnas: x, y, vx, vy, x, y, vx, vy...
    particles_data = raw_data[:, 1:]
    
    # Calcular N
    N = particles_data.shape[1] // 4
    
    # Reformatear a (Frames, N, 4) -> [x, y, vx, vy]
    data_reshaped = particles_data.reshape(len(times), N, 4)
    print(f"   [INFO] Datos cargados: {len(times)} frames, {N} partículas.")

except Exception as e:
    print(f"   [ERROR] Falló la lectura de posiciones: {e}")
    exit()

# ==========================================
# 3. HISTOGRAMA DE VELOCIDADES (Final)
# ==========================================
print("-> Generando Histograma de Velocidades...")
try:
    last_frame = data_reshaped[-1]
    vx = last_frame[:, 2]
    vy = last_frame[:, 3]
    speeds = np.sqrt(vx**2 + vy**2)

    plt.figure(figsize=(7, 5))
    count, bins, ignored = plt.hist(speeds, bins=20, density=True, 
                                    alpha=0.6, color='purple', edgecolor='black', label='Simulación')

    # Ajuste Maxwell-Boltzmann 2D Teórico
    v_mean = np.mean(speeds)
    if v_mean > 0:
        # En 2D, P(v) = (v / T) * exp(-v^2 / 2T) aprox.
        # sigma_mb relaciona con v_mean: v_mean = sigma * sqrt(pi/2)
        sigma_mb = v_mean / np.sqrt(np.pi / 2)
        v_grid = np.linspace(0, max(speeds)*1.2, 100)
        pdf = (v_grid / sigma_mb**2) * np.exp(-v_grid**2 / (2*sigma_mb**2))
        plt.plot(v_grid, pdf, linewidth=2, color='r', label='Maxwell-Boltzmann (2D)')

    plt.xlabel('Velocidad |v|')
    plt.ylabel('Probabilidad')
    plt.title(f'Distribución de Velocidades (t={times[-1]:.2f})')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(os.path.join(RESULTS_DIR, "histograma_velocidades.png"), dpi=150)
    plt.close()
    print("   [OK] Gráfica guardada: histograma_velocidades.png")

except Exception as e:
    print(f"   [ERROR] Histograma falló: {e}")

# ==========================================
# 4. TRAYECTORIAS (Conexión de puntos)
# ==========================================
print("-> Generando Trayectorias...")
try:
    plt.figure(figsize=(6, 6))
    # Graficar solo las primeras 10 para no saturar
    for i in range(min(20, N)):
        # Extraer x, y de la partícula i a lo largo del tiempo
        x_traj = data_reshaped[:, i, 0]
        y_traj = data_reshaped[:, i, 1]
        plt.plot(x_traj, y_traj, linewidth=0.5, alpha=0.7)

    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.title('Trayectorias (Muestra)')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.grid(True)
    plt.savefig(os.path.join(RESULTS_DIR, "trayectorias.png"), dpi=150)
    plt.close()
    print("   [OK] Gráfica guardada: trayectorias.png")

except Exception as e:
    print(f"   [ERROR] Trayectorias falló: {e}")

# ==========================================
# 5. ANIMACIÓN (GIF)
# ==========================================
print("-> Generando Animación GIF...")

# Configuración de figura
fig, ax = plt.subplots(figsize=(6, 6))
# FONDOS NEGROS PARA ESTILO
fig.patch.set_facecolor('black')
ax.set_facecolor('black')

# IMPORTANTE: LÍMITES FIJOS para ver PBC
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)

# Etiquetas en blanco
ax.set_xlabel("X", color='white')
ax.set_ylabel("Y", color='white')
ax.set_title("Simulación Fluido Lennard-Jones", color='white')
ax.tick_params(axis='x', colors='white')
ax.tick_params(axis='y', colors='white')
# Bordes blancos
for spine in ax.spines.values():
    spine.set_edgecolor('white')

# Cargar Sprite si existe
usar_sprite = False
sprite_img = None
if os.path.exists(IMG_SPRITE):
    try:
        sprite_img = mpimg.imread(IMG_SPRITE)
        usar_sprite = True
    except:
        pass

# Objeto scatter inicial (respaldo si no hay imagen)
scat = ax.scatter([], [], s=60, c='cyan', edgecolors='white', alpha=0.9)
time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes, color='white', fontsize=10)

def init():
    return scat, time_text

def update(frame_idx):
    # Limpiar ejes pero mantener configuración oscura y límites
    ax.clear()
    ax.set_facecolor('black')
    ax.set_xlim(0, 1) # <--- CRÍTICO: Mantener límites fijos
    ax.set_ylim(0, 1) # <--- CRÍTICO: Mantener límites fijos
    
    # Datos del frame actual
    current_data = data_reshaped[frame_idx]
    x_data = current_data[:, 0]
    y_data = current_data[:, 1]
    
    # Título y Tiempo
    t_current = times[frame_idx]
    ax.set_title(f"Simulación LJ - Tiempo: {t_current:.2f}", color='white')

    if usar_sprite:
        # Dibujar Sprites
        for i in range(len(x_data)):
            # Zoom: Ajustar según el tamaño de tu imagen. 0.05 es usualmente bueno.
            imagebox = OffsetImage(sprite_img, zoom=0.2) 
            ab = AnnotationBbox(imagebox, (x_data[i], y_data[i]), frameon=False)
            ax.add_artist(ab)
    else:
        # Dibujar Círculos Neón
        ax.scatter(x_data, y_data, s=50, c='cyan', edgecolors='white', alpha=0.8)

    return []

# Lógica de saltar frames si hay demasiados datos (para que el GIF no pese 500MB)
total_frames = len(times)
skip = 1
if total_frames > 200: skip = 2
if total_frames > 500: skip = 5
if total_frames > 1000: skip = 10
if total_frames > 5000: skip = 50

frames_to_anim = range(0, total_frames, skip)
print(f"   [INFO] Animando {len(frames_to_anim)} cuadros (Salto: {skip})...")

ani = animation.FuncAnimation(fig, update, frames=frames_to_anim, 
                              init_func=init, blit=False, interval=30)

output_gif = os.path.join(RESULTS_DIR, "animacion.gif")
try:
    ani.save(output_gif, writer='pillow', fps=30)
    print(f"   [OK] Animación guardada: {output_gif}")
except Exception as e:
    print(f"   [ERROR] Falló al guardar GIF: {e}")

print(">>> Proceso finalizado exitosamente <<<")