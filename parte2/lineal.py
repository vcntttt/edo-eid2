import numpy as np
import matplotlib.pyplot as plt

# Parámetros del problema
N = 1000
S0 = 999
I0 = 1
R0 = 0

alpha = 0.003
beta = 0.005

# Dominio de tiempo (en horas)
t_max = 1000
t = np.linspace(0, t_max, 1000)

# MODELO A (2 variables)
S_A = S0 * np.exp(-alpha * t)
I_A = I0 + S0 * (1 - np.exp(-alpha * t))


# MODELO B (3 variables)
S_B = S0 * np.exp(-alpha * t)
I_B = I0 * np.exp(-beta * t) + (alpha * S0 / (beta - alpha)) * (
    np.exp(-alpha * t) - np.exp(-beta * t)
)
R_B = N - S_B - I_B

# GRÁFICOS
plt.style.use("default")
fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(t, S_A, label="S_A(t) - susceptibles (Modelo A)")
ax.plot(t, I_A, label="I_A(t) - infectados (Modelo A)")
ax.set_xlabel("Tiempo [h]")
ax.set_ylabel("Número de equipos")
ax.set_title("Modelo A: sistema lineal sin recuperación")
ax.legend(frameon=False)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig("parte2/modelo_a.pdf", format="pdf", bbox_inches="tight")

fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(t, S_B, label="S_B(t) - susceptibles")
ax.plot(t, I_B, label="I_B(t) - infectados")
ax.plot(t, R_B, label="R_B(t) - recuperados")
ax.set_xlabel("Tiempo [h]")
ax.set_ylabel("Número de equipos")
ax.set_title("Modelo B: sistema lineal con recuperación")
ax.legend(frameon=False)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig("parte2/modelo_b.pdf", format="pdf", bbox_inches="tight")
