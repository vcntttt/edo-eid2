import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# --- Estilo Tufte: máxima reducción de tinta ---
plt.style.use("default")


def tufte_ax(ax):
    """Aplica estilo Tufte a un eje."""
    # Sin bordes (spines) superiores y derechos
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Ejes más delgados
    ax.spines["left"].set_linewidth(0.5)
    ax.spines["bottom"].set_linewidth(0.5)

    # Quitar ticks innecesarios
    ax.tick_params(axis="both", which="both", length=3, width=0.5)

    # Grid extremadamente suave (o comentado si prefieres sin grid)

    return ax


# --- Parámetros del modelo ---
N = 1000
alpha = 0.003
beta = 0.25

S0, I0, R0 = 999, 1, 0
y0 = [S0, I0, R0]

t0, t_final = 0, 25
t = np.linspace(t0, t_final, 1000)


def sir_nolineal(t, y):
    S, I, R = y
    return [-alpha * S * I, alpha * S * I - beta * I, beta * I]


# --- Resolver con solve_ivp ---
sol = solve_ivp(sir_nolineal, (t0, t_final), y0, t_eval=t, rtol=1e-6, atol=1e-8)

S, I, R = sol.y

# --- Gráfico estilo Tufte ---
fig, ax = plt.subplots(figsize=(8, 4.2))
tufte_ax(ax)

ax.plot(t, S, linewidth=1.2, label="Susceptibles S(t)")
ax.plot(t, I, linewidth=1.2, label="Infectados I(t)")
ax.plot(t, R, linewidth=1.2, label="Recuperados R(t)")

# Etiquetas mínimas, ordenadas
ax.set_xlabel("Tiempo (horas)", fontsize=10)
ax.set_ylabel("Número de equipos", fontsize=10)

# Leyenda sin borde, discreta
legend = ax.legend(frameon=False, fontsize=9)

plt.tight_layout()
plt.savefig("parte2/tufte.png")
plt.show()
