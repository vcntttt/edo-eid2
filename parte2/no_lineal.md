# Funcionamiento de `solve_ivp` y del método Runge–Kutta adaptativo

En esta sección se explica con detalle cómo opera la función `solve_ivp` de SciPy
y cómo funciona internamente el método de Runge–Kutta adaptativo que utiliza por
defecto (RK45, Dormand–Prince).

## 1. ¿Qué es realmente `solve_ivp`?

`solve_ivp` (solve Initial Value Problem) es un integrador general de ecuaciones
diferenciales ordinarias del tipo:

$$
y'(t) = f(t, y(t)), \qquad y(t_0) = y_0.
$$

SciPy no interpreta matemáticas simbólicas.  
No entiende letras como α o β, ni expresiones como “dS/dt”.  
Solo recibe:

- una función Python que define las derivadas,
- un intervalo de tiempo,
- una condición inicial,
- un método numérico,
- tolerancias de error.

Con esta información, avanza paso a paso en el tiempo y aproxima la solución.

## 2. Concepto básico de Runge–Kutta

Los métodos de Runge–Kutta aproximan la solución evaluando varias veces la
función de derivadas $f(t, y)$ en posiciones intermedias del intervalo
$[t, t + h]$.  
El clásico RK4 usa cuatro evaluaciones, pero SciPy no usa RK4 puro.

SciPy usa **RK45 (Dormand–Prince)**, que produce **dos aproximaciones
simultáneas**:

- una de orden 4,
- otra de orden 5.

La diferencia entre ambas permite estimar el error del paso.

## 3. ¿Qué significa "adaptativo"?

Un método adaptativo **no usa un paso fijo h**.  
En cada paso:

1. Calcula dos aproximaciones: $y_4$ y $y_5$.
2. Estima el error local:
   $$
   \text{error} = |y_5 - y_4|.
   $$
3. Compara ese error con las tolerancias solicitadas:
   $$
   |y_5 - y_4| \le \text{atol} + \text{rtol} \cdot |y_5|.
   $$
4. Si el error es pequeño → el paso se **acepta** y se puede **aumentar h**.
5. Si el error es demasiado grande → el paso se **rechaza** y se repite con
   un **h más pequeño**.

Esta estrategia hace que el método sea eficiente y preciso.

## 4. Cómo se ajusta el tamaño del paso

Si el error del paso es razonable, el método aumenta el tamaño del paso usando:

$$
h*{\text{nuevo}} = h*{\text{actual}}
\left( \frac{\text{tolerancia}}{\text{error}} \right)^{1/5}.
$$

Si el error es grande, el valor de la fracción es menor que uno y el paso se
reduce.

Esto permite:

- usar pasos pequeños cerca del pico de infección (donde la solución cambia muy rápido),
- usar pasos grandes donde la solución es suave.

## 5. Tolerancias `rtol` y `atol`

Las tolerancias regulan cuánto error se permite por paso:

### Tolerancia Relativa (`rtol`)

Escala con el tamaño de la solución:

$$
\text{error permitido} \approx \text{rtol} \cdot |y|.
$$

### Tolerancia Absoluta (`atol`)

Controla el error cuando la solución es muy pequeña.

Ambas tolerancias combinadas dan un control flexible del error.

## 6. ¿Usa `solve_ivp` los puntos de `t_eval` para integrar?

No.  
Esto es un punto clave:

**SciPy NO integra usando los puntos de `t_eval`.**

Usa muchos pasos internos adaptativos entre esos puntos y **después interpola**
la solución en los valores de `t_eval`.

Esto significa que:

- puedes pedir 10 puntos o 10 000 en `t_eval`,
- la precisión del método no cambia,
- `t_eval` solo determina dónde quieres ver la solución final.

## 7. ¿Qué pasa en cada iteración interna?

En cada paso, `solve_ivp`:

1. Evalúa tu función `f(t, y)` → calcula derivadas.
2. Calcula varios valores intermedios $k_1, k_2, \dots, k_7$ (Dormand–Prince tiene 7 etapas).
3. Construye:
   - $y_4$: solución de orden 4,
   - $y_5$: solución de orden 5.
4. Estima el error: $ |y_5 - y_4| $.
5. Decide si acepta o rechaza el paso.
6. Ajusta el tamaño del paso.
7. Continúa hasta llegar a $t\_{\text{final}}$.

## 8. ¿Por qué RK45 es ideal para un SIR no lineal?

El modelo SIR no lineal tiene:

- crecimiento rápido inicial,
- un pico abrupto de infectados,
- un desaceleramiento en las colas.

El método adaptativo:

- reduce el paso en la zona del pico,
- aumenta el paso cuando la solución es suave,
- mantiene la precisión sin gastar tiempo de cómputo innecesario.

## 9. Resumen

> El integrador `solve_ivp` utiliza por defecto el método de Runge–Kutta de orden
> 4–5 de Dormand–Prince (RK45), un método adaptativo que estima el error local
> comparando dos aproximaciones de distinto orden. En cada paso, el método ajusta
> automáticamente el tamaño del paso temporal para satisfacer las tolerancias
> relativas (`rtol`) y absolutas (`atol`). Esto permite resolver sistemas no
> lineales como el modelo SIR de manera eficiente, reduciendo el paso en regiones
> donde la solución presenta variaciones rápidas y aumentándolo cuando la solución
> es suave.
