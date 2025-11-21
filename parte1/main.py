import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# =========
# Parámetros base (ajustables)
# =========
I0 = 1.0  # infectados iniciales
N = 1000  # tamaño de red (solo referencia)
r_lab = 3.0  # infecciones/hora en jornada laboral (08:00–18:00 => 0..10 h)
rho = 0.2  # infecciones/hora fuera de jornada
r_peak = 6.0  # pico de campaña (2.4a)
alpha = 0.2  # decaimiento temporal (1/h) (2.4a)
r0_state = 3.0  # r0 para r(I)=r0/(1+cI) (2.4b)
c_state = 0.01  # c para r(I)=r0/(1+cI) (2.4b)
r_log = 0.35  # r logístico (1/h) (2.5)
K = 300.0  # capacidad efectiva logístico (2.5)

# Horizonte y resolución
DT_MIN = 1 / 60  # paso 1 min
T_DAY = 24.0  # horas en un día
T_WORK = 10.0  # 08:00 a 18:00 => 10 h si t=0 equivale a 08:00
T_LINEAR = 72.0  # horizonte para 2.3 y 2.6 (3 días)
T_TEMPORAL = 48.0  # horizonte 2.4(a)
T_STATE = 72.0  # horizonte 2.4(b)
T_LOG = 72.0  # horizonte 2.5

# Salida
OUTDIR = Path("parte1/figs")
OUTDIR.mkdir(parents=True, exist_ok=True)


# =========
# Utilidades
# =========
def r_piecewise(time_h: np.ndarray, r_lab: float, rho: float) -> np.ndarray:
    """
    r(t) por tramos diario:
    - r_lab en [0, T_WORK]
    - rho fuera de ese rango (hasta 24 h), repetido cada día
    """
    tau = np.mod(time_h, T_DAY)  # hora dentro del día actual
    return np.where((tau >= 0) & (tau <= T_WORK), r_lab, rho)


def rk4_scalar(f, y0: float, t: np.ndarray) -> np.ndarray:
    """
    Integrador RK4 para y' = f(y), y(t0)=y0 con malla t (escalar).
    """
    y = np.zeros_like(t)
    y[0] = y0
    for k in range(1, len(t)):
        h = t[k] - t[k - 1]
        yk = y[k - 1]
        k1 = f(yk)
        k2 = f(yk + 0.5 * h * k1)
        k3 = f(yk + 0.5 * h * k2)
        k4 = f(yk + h * k3)
        y[k] = yk + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
    return y


def save_pdf(figpath: Path):
    plt.savefig(figpath, bbox_inches="tight")
    plt.close()


# =========
# 2.2 — r(t) por tramos en 24 h
# =========
day_hours = np.linspace(0, T_DAY, int(T_DAY / DT_MIN) + 1)
r_schedule = np.where((day_hours >= 0) & (day_hours <= T_WORK), r_lab, rho)

plt.figure()
plt.step(day_hours, r_schedule, where="post")
plt.title("2.2 — Tasa de propagación por tramos en 24 h")
plt.xlabel("Tiempo dentro del día (h)  —  0≡08:00, 10≡18:00")
plt.ylabel("r(t) [infecciones/hora]")
plt.xlim(0, 24)
plt.ylim(0, max(r_lab, rho) * 1.2)
plt.grid(True, alpha=0.3)
save_pdf(OUTDIR / "r_schedule_24h.pdf")

# =========
# 2.3 — Modelo lineal por tramos (3 días)
# =========
t_linear = np.arange(0, T_LINEAR + DT_MIN, DT_MIN)
r_t = r_piecewise(t_linear, r_lab, rho)
I_linear = I0 + np.cumsum(r_t) * DT_MIN  # I' = r(t) => integral acumulada

plt.figure()
plt.plot(t_linear, I_linear)
plt.title("2.3 — Modelo lineal por tramos (3 días)")
plt.xlabel("Tiempo (h)")
plt.ylabel("Infectados I(t) [equipos]")
plt.grid(True, alpha=0.3)
save_pdf(OUTDIR / "I_linear_piecewise_3days.pdf")

# =========
# 2.4 — Tasa dependiente del tiempo y del estado
# =========
# (a) r(t) decreciente tras una campaña
t_temporal = np.arange(0, T_TEMPORAL + DT_MIN, DT_MIN)
r_temporal = r_peak * np.exp(-alpha * np.maximum(0.0, t_temporal - 0.0))
I_temporal = I0 + np.cumsum(r_temporal) * DT_MIN

plt.figure()
plt.plot(t_temporal, I_temporal)
plt.title("2.4 — I(t) con tasa temporal decreciente")
plt.xlabel("Tiempo (h)")
plt.ylabel("Infectados I(t) [equipos]")
plt.grid(True, alpha=0.3)
save_pdf(OUTDIR / "I_temporal_decay.pdf")


# (b) r(I) = r0/(1 + c I), integración numérica (RK4) de I' = r(I)
def dI_dt_state(I: float) -> float:  # noqa: E741
    return r0_state / (1.0 + c_state * I)


t_state = np.arange(0, T_STATE + DT_MIN, DT_MIN)
I_state = rk4_scalar(dI_dt_state, I0, t_state)

plt.figure()
plt.plot(t_state, I_state)
plt.title("2.4 — I(t) con tasa dependiente del estado r(I)=r0/(1+cI)")
plt.xlabel("Tiempo (h)")
plt.ylabel("Infectados I(t) [equipos]")
plt.grid(True, alpha=0.3)
save_pdf(OUTDIR / "I_state_dependent_numeric.pdf")

# =========
# 2.5 — Modelo logístico (solución explícita)
# =========
t_log = np.linspace(0, T_LOG, 1000)
A = (K - I0) / I0
I_logistic = K / (1.0 + A * np.exp(-r_log * t_log))

plt.figure()
plt.plot(t_log, I_logistic)
plt.title("2.5 — Modelo logístico: solución explícita")
plt.xlabel("Tiempo (h)")
plt.ylabel("Infectados I(t) [equipos]")
plt.grid(True, alpha=0.3)
save_pdf(OUTDIR / "I_logistic.pdf")

# =========
# 2.6 — Comparación lineal (tramos) vs logístico
# =========
# Interpolamos el logístico al mallado del lineal para comparación limpia
I_logistic_interp = np.interp(t_linear, t_log, I_logistic)

plt.figure()
plt.plot(t_linear, I_linear, label="Lineal (tramos)")
plt.plot(t_linear, I_logistic_interp, label="Logístico")
plt.title("2.6 — Comparación: lineal por tramos vs logístico")
plt.xlabel("Tiempo (h)")
plt.ylabel("Infectados I(t) [equipos]")
plt.legend()
plt.grid(True, alpha=0.3)
save_pdf(OUTDIR / "compare_linear_vs_logistic.pdf")

print("Figuras guardadas en:", OUTDIR.resolve())
