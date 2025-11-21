# ruff: noqa: E741
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from parameters import alpha, beta, S0, I0, R0, t0, tol

t_final = 25
t = np.linspace(t0, t_final, 1000)
y0 = [S0, I0, R0]

tol_s = tol / 10


def sir_nolineal(t, y):
    S, I, R = y
    dSdt = -alpha * S * I
    dIdt = alpha * S * I - beta * I
    dRdt = beta * I
    return [dSdt, dIdt, dRdt]


# ode45 ~ rk45
sol = solve_ivp(sir_nolineal, (t0, t_final), y0, t_eval=t, rtol=1e-6, atol=1e-8)
t = sol.t
S, I, R = sol.y

# deteccion de hitos
hitos = {}

# pico de infectados
indice = np.argmax(I)
hitos["pico"] = (t[indice], I[indice])

# intersección S = I
# buscamos donde la diferencia cambia de signo
signo_diff = np.sign(S - I)
cambios = np.where(np.diff(signo_diff))[0]

if len(cambios) == 1:
    indice = cambios[0]
    hitos["S=I"] = (t[indice], S[indice])
elif len(cambios) > 1:
    # si hay más de uno, guardamos el primero y el último
    hitos["S=I_1"] = (t[cambios[0]], S[cambios[0]])
    hitos["S=I_2"] = (t[cambios[-1]], S[cambios[-1]])

# S ≈ 0
indice_s0 = np.where(S < tol_s)[0]
if len(indice_s0) > 0:
    indice = indice_s0[0]
    hitos["S≈0"] = (t[indice], S[indice])

# R ≈ I
mask = R > tol  # solo consideramos cuando R ya ha crecido un poco
if mask.any():
    indice_ri = np.argmin(np.abs(R[mask] - I[mask]))
    indice_real = np.where(mask)[0][indice_ri]  # ajuste índice al array original
    hitos["R=I"] = (t[indice_real], R[indice_real])

# R ≈ N
indice_rmax = np.where(R > 1000 - tol)[0]
if len(indice_rmax) > 0:
    indice = indice_rmax[0]
    hitos["R≈1000"] = (t[indice], R[indice])


# graficos
fig, ax = plt.subplots(figsize=(10, 5))

ax.plot(t, S, label="S(t) - Susceptibles", linewidth=2)
ax.plot(t, I, label="I(t) - Infectados", linewidth=2)
ax.plot(t, R, label="R(t) - Recuperados", linewidth=2)

anotaciones = {
    "pico": {"txt": "Pico I(t)", "col": "red", "off": (0, 40)},
    "S=I": {"txt": "S ≈ I", "col": "purple", "off": (10, 30)},
    "S=I_1": {"txt": "S ≈ I (1)", "col": "purple", "off": (10, 30)},
    "S=I_2": {"txt": "S ≈ I (2)", "col": "purple", "off": (10, -30)},
    "S≈0": {"txt": "S ≈ 0", "col": "black", "off": (10, 20)},
    "R=I": {"txt": "R ≈ I", "col": "green", "off": (30, 40)},
    "R≈1000": {"txt": "R ≈ 1000", "col": "brown", "off": (10, -30)},
}

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
ax.set_title("Modelo SIR No Lineal")
ax.legend(frameon=False)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
plt.savefig("parte2/no_lineal.pdf", format="pdf", bbox_inches="tight")
plt.show()
