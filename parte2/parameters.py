import numpy as np

# condiciones iniciales
N = 1000  # Número de equipos
S0 = 999  # Número inicial de equipos susceptibles
I0 = 1  # Número inicial de equipos infectados
R0 = 0  # Número inicial de equipos recuperados

# parametros
alpha = 0.003  # tasa de infección
beta = 0.05  # tasa de recuperación

# horas
t0 = 0
t_final = 2500
t = np.linspace(t0, t_final, 1000)

# hitos
tol = 1.0  # tolerancia para "≈ 0" o "≈ N"
