import numpy as np
import matplotlib.pyplot as plt

# Parámetros del problema
N = 1000
S0 = 999
I0 = 1
R0 = 0

alpha = 0.003
beta = 0.25

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
plt.figure(figsize=(10, 5))
plt.plot(t, S_A, label="S_A(t) - susceptibles (Modelo A)")
plt.plot(t, I_A, label="I_A(t) - infectados (Modelo A)")
plt.xlabel("Tiempo [h]")
plt.ylabel("Número de equipos")
plt.title("Modelo A: sistema lineal sin recuperación")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

plt.figure(figsize=(10, 5))
plt.plot(t, S_B, label="S_B(t) - susceptibles")
plt.plot(t, I_B, label="I_B(t) - infectados")
plt.plot(t, R_B, label="R_B(t) - recuperados")
plt.xlabel("Tiempo [h]")
plt.ylabel("Número de equipos")
plt.title("Modelo B: sistema lineal con recuperación")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
