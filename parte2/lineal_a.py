import numpy as np
import matplotlib.pyplot as plt
from parameters import alpha, N, S0, I0, t, tol

# MODELO A (2 variables)
S_A = S0 * np.exp(-alpha * t)
I_A = I0 + S0 * (1 - np.exp(-alpha * t))

# detección de hitos
hitos = {}

# intersección S = I
signo_diff = np.sign(S_A - I_A)
cambios = np.where(np.diff(signo_diff))[0]
if len(cambios) > 0:
    indice = cambios[0]
    hitos["S=I"] = (t[indice], S_A[indice])

# S ≈ 0
indice_s0 = np.where(S_A < tol)[0]
if len(indice_s0) > 0:
    indice = indice_s0[0]
    hitos["S≈0"] = (t[indice], S_A[indice])

# I ≈ N
indice_imax = np.where(I_A > N - tol)[0]
if len(indice_imax) > 0:
    indice = indice_imax[0]
    hitos["I≈N"] = (t[indice], I_A[indice])

# grafico
anotaciones = {
    "S=I": {"txt": "S ≈ I", "col": "purple", "off": (20, 0)},
    "S≈0": {"txt": "S ≈ 0", "col": "black", "off": (10, 20)},
    "I≈N": {"txt": "I ≈ N", "col": "orange", "off": (10, -30)},
}

fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(t, S_A, label="S(t) - Susceptibles")
ax.plot(t, I_A, label="I(t) - Infectados")

for key, (tx, vy) in hitos.items():
    if key in anotaciones:
        c = anotaciones[key]
        ax.scatter(tx, vy, color=c["col"], zorder=5)
        ax.annotate(
            f"{c['txt']}\n({tx:.1f}h, {vy:.0f})",
            xy=(tx, vy),
            xytext=c["off"],
            textcoords="offset points",
            arrowprops=dict(arrowstyle="->", color=c["col"]),
            color=c["col"],
        )

ax.set_xlabel("Tiempo (horas)")
ax.set_ylabel("Número de equipos")
ax.set_title("Modelo A: Sistema lineal sin recuperación")
ax.legend(frameon=False, loc="center right")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig("parte2/lineal_a.pdf", format="pdf", bbox_inches="tight")
plt.show()
