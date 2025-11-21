import numpy as np
import matplotlib.pyplot as plt
from parameters import alpha, beta, N, S0, I0, t, tol


# MODELO B (3 variables)
S_B = S0 * np.exp(-alpha * t)
I_B = I0 * np.exp(-beta * t) + (alpha * S0 / (beta - alpha)) * (
    np.exp(-alpha * t) - np.exp(-beta * t)
)
R_B = N - S_B - I_B

# detección de hitos
hitos = {}

# pico de infectados
indice_pico = np.argmax(I_B)
hitos["pico"] = (t[indice_pico], I_B[indice_pico])

# intersección S = I
signo_diff = np.sign(S_B - I_B)
cambios = np.where(np.diff(signo_diff))[0]
if len(cambios) > 0:
    indice = cambios[0]
    hitos["S=I"] = (t[indice], S_B[indice])

# S ≈ 0
indice_s0 = np.where(S_B < tol)[0]
if len(indice_s0) > 0:
    indice = indice_s0[0]
    hitos["S≈0"] = (t[indice], S_B[indice])

# R ≈ I
mask = R_B > tol
if mask.any():
    indice_ri = np.argmin(np.abs(R_B[mask] - I_B[mask]))
    indice_real = np.where(mask)[0][indice_ri]
    hitos["R=I"] = (t[indice_real], R_B[indice_real])

# R ≈ N
indice_rmax = np.where(R_B > N - tol)[0]
if len(indice_rmax) > 0:
    indice = indice_rmax[0]
    hitos["R≈N"] = (t[indice], R_B[indice])

# I ≈ 0
indice_i0 = np.where(I_B < tol)[0]
if len(indice_i0) > 0:
    indice = indice_i0[0]
    hitos["I≈0"] = (t[indice], I_B[indice])

# grafico
anotaciones = {
    "pico": {"txt": "Pico I(t)", "col": "red", "off": (0, 40)},
    "S=I": {"txt": "S ≈ I", "col": "purple", "off": (10, 30)},
    "S≈0": {"txt": "S ≈ 0", "col": "black", "off": (10, 20)},
    "R=I": {"txt": "R ≈ I", "col": "green", "off": (30, 0)},
    "R≈N": {"txt": "R ≈ N", "col": "brown", "off": (10, -30)},
    "I≈0": {"txt": "I ≈ 0", "col": "orange", "off": (10, 40)},
}

fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(t, S_B, label="S(t) - Susceptibles", color="C0")
ax.plot(t, R_B, label="R(t) - Recuperados", color="C2")

# eje secundario para I(t)
ax2 = ax.twinx()
ax2.plot(t, I_B, label="I(t) - Infectados", color="C1")

for key, (tx, vy) in hitos.items():
    if key in anotaciones:
        c = anotaciones[key]
        if key == "pico" or key == "I≈0":
            ax2.scatter(tx, vy, color=c["col"], zorder=5)
            ax2.annotate(
                f"{c['txt']}\n({tx:.1f}h, {vy:.0f})",
                xy=(tx, vy),
                xytext=c["off"],
                textcoords="offset points",
                arrowprops=dict(arrowstyle="->", color=c["col"]),
                color=c["col"],
            )
        else:
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
ax.set_ylabel("S(t), R(t)")
ax2.set_ylabel("I(t)")

lines, labels = ax.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax.legend(lines + lines2, labels + labels2, frameon=False, loc="center right")
ax.set_title("Modelo B: Sistema lineal con recuperación")

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax2.spines["top"].set_visible(False)
plt.tight_layout()
plt.savefig("parte2/lineal_b.pdf", format="pdf", bbox_inches="tight")
plt.show()
