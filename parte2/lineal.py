import numpy as np
import matplotlib.pyplot as plt
from parameters import alpha, beta, N, S0, I0, t


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
fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(t, S_A, label="S_A(t) - susceptibles (Modelo A)")
ax.plot(t, I_A, label="I_A(t) - infectados (Modelo A)")
ax.set_xlabel("Tiempo [h]")
ax.set_ylabel("Número de equipos")
ax.set_title("Modelo A: sistema lineal sin recuperación")
ax.legend(frameon=False, loc="center right")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig("parte2/modelo_lineal_a.pdf", format="pdf", bbox_inches="tight")

fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(t, S_B, label="S_B(t) - susceptibles", color="C0")
ax.plot(t, R_B, label="R_B(t) - recuperados", color="C2")

# eje secundario para I(t)
ax2 = ax.twinx()
ax2.plot(t, I_B, label="I_B(t) - infectados", color="C1")

ax.set_xlabel("Tiempo [h]")
ax.set_ylabel("S(t), R(t)")
ax2.set_ylabel("I(t)")

# combinar leyendas
lines, labels = ax.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax.legend(lines + lines2, labels + labels2, frameon=False, loc="center right")
ax.set_title("Modelo B: sistema lineal con recuperación")

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax2.spines["top"].set_visible(False)
plt.tight_layout()
plt.savefig("parte2/modelo_lineal_b.pdf", format="pdf", bbox_inches="tight")
plt.show()
