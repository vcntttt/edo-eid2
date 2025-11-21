# ruff: noqa: E741
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from parameters import alpha, beta, S0, I0, R0, t0, t_final, t

y0 = [S0, I0, R0]


def sir_nolineal(t, y):
    S, I, R = y
    dSdt = -alpha * S * I
    dIdt = alpha * S * I - beta * I
    dRdt = beta * I
    return [dSdt, dIdt, dRdt]


sol = solve_ivp(sir_nolineal, (t0, t_final), y0, t_eval=t, rtol=1e-6, atol=1e-8)

# vectores
t = sol.t
S, I, R = sol.y

# grafiquito
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(t, S, label="S(t) - Susceptibles")
ax.plot(t, I, label="I(t) - Infectados")
ax.plot(t, R, label="R(t) - Recuperados")
ax.set_xlabel("Tiempo (horas)")
ax.set_ylabel("Número de equipos")
ax.legend(frameon=False)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_title("Modelo SIR No Lineal")
plt.tight_layout()
plt.savefig("parte2/no_lineal.pdf", format="pdf", bbox_inches="tight")
plt.show()
