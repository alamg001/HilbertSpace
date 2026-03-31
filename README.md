<div align="center">

<img src="https://www.uni.edu.pe/images/logo_uni.png" width="80"/>

# Doctorado en Ciencias e Ingeniería Estadística

**Antonio Lam García**

Universidad Nacional de Ingeniería

---

<img src="https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEhyXMenkLZGsGeWDp-agKQACtNTPDIhOqLvefpd51lZ-ZH4wvilhNd1amtWj6z_sHGdbd33mQ42q7N7dPfWHHmURTuI9H_90E_06CdaTFlFKBgumTdfzTTnim0SzG-XlVH7a83zWg98DQW_jpTcxgianCjslPrdsMHqSHCNCPKhKEc0kD7HFnR1Qx7I50eK/s943/logo_y_uni.png" width="600" alt="Logo Universidad Nacional de Ingeniería"/>

</div>

---

## CS281B/Stat241B (Primavera 2008) — Teoría del Aprendizaje Estadístico

### Clase 7: Espacios de Hilbert con Núcleo Reproductor

> **Profesor:** Peter Bartlett  
> **Tomador de notas:** Chunhui Gu

---

## 1. Espacios de Hilbert con Núcleo Reproductor

### Espacio de Hilbert y Núcleo

Un producto interno ⟨u, v⟩ puede ser:

1. Un producto punto usual: `⟨u,v⟩ = vᵀw = ∑ᵢ vᵢwᵢ`
2. Un producto de núcleo: `⟨u,v⟩ = k(v,w) = ψ(v)ᵀψ(w)` (donde ψ(u) puede tener dimensiones infinitas)

Sin embargo, un producto interno ⟨·,·⟩ debe satisfacer:

1. **Simetría**

$$\langle u, v \rangle = \langle v, u \rangle \quad \forall\, u, v \in \mathcal{X}$$

2. **Bilinealidad**

$$\langle \alpha u + \beta v,\, w \rangle = \alpha\langle u, w \rangle + \beta\langle v, w \rangle \quad \forall\, u,v,w \in \mathcal{X},\; \forall\, \alpha,\beta \in \mathbb{R}$$

3. **Definición positiva**

$$\langle u, u \rangle \geq 0 \quad \forall\, u \in \mathcal{X}$$

$$\langle u, u \rangle = 0 \iff u = 0$$

---

> **📘 Definición.** Un *Espacio de Hilbert* es un espacio de producto interno completo y separable respecto a la norma definida por el producto interno.

**Ejemplos:**

1. $\mathbb{R}^n$ con $\langle a, b \rangle = a^\top b$
2. Espacio $\ell_2$ de sucesiones sumables cuadráticas: $\langle x, y \rangle = \sum_{i=1}^{\infty} x_i y_i$
3. Espacio $L_2$ de funciones integrables cuadráticas: $\langle f, g \rangle = \int_S f(x)g(x)\,dx$

---

> **📘 Definición.** $k(\cdot,\cdot)$ es un *núcleo reproductor* de un espacio de Hilbert $\mathcal{H}$ si $\forall f \in \mathcal{H}$:
>
> $$f(x) = \langle k(x,\cdot),\, f(\cdot) \rangle$$

---

Un **Espacio de Hilbert con Núcleo Reproducente (RKHS)** es un espacio de Hilbert $\mathcal{H}$ con un núcleo reproductor cuyo recorrido es denso en $\mathcal{H}$. De forma equivalente, es un espacio de Hilbert de funciones con todos los funcionales de evaluación acotados y lineales.

> **Nota:** El espacio $L_2$ es un espacio de Hilbert, pero **no** un RKHS porque la función delta que tiene la propiedad reproductora
>
> $$f(x) = \int_S \delta(x - u)\, f(u)\, du$$
>
> no satisface la condición de integrabilidad cuadrática:
>
> $$\int_S \delta(u)^2\, du \not< \infty$$

---

> **📘 Definición.** $k: \mathcal{X} \times \mathcal{X} \to \mathbb{R}$ es un *núcleo* si:
> 1. $k$ es simétrico: $k(x,y) = k(y,x)$
> 2. $k$ es semidefinido positivo: $\forall\, x_1,\ldots,x_n \in \mathcal{X}$, la matriz de Gram $K_{ij} = k(x_i, x_j)$ es semidefinida positiva.

**Propiedades:**

1. $k(x,x) \geq 0$ (caso $n=1$ de la matriz de Gram)
2. $k(u,v) \leq \sqrt{k(u,u)\,k(v,v)}$ (desigualdad de Cauchy-Schwarz)

Para la segunda propiedad, con $n=2$:

$$a = \begin{bmatrix} k(v,v) \\ -k(u,v) \end{bmatrix}, \quad K = \begin{bmatrix} k(u,u) & k(u,v) \\ k(v,u) & k(v,v) \end{bmatrix} \succeq 0 \iff a^\top K a \geq 0$$

$$\iff \bigl[k(v,v)k(u,u) - k(u,v)^2\bigr]k(v,v) \geq 0$$

Como $k(v,v) \geq 0$, se concluye que $k(v,v)k(u,u) \geq k(u,v)^2$.

---

### Construir un RKHS

Dado un núcleo $k$, se define el mapeo de características del núcleo reproductor $\Phi: \mathcal{X} \to \mathbb{R}^{\mathcal{X}}$ como:

$$\Phi(x) = k(\cdot, x)$$

El espacio vectorial considerado es:


$$
\mathrm{span}\{\Phi(x) : x \in \mathcal{X}\} = \{ 
f(\cdot) = \sum_{i=1}^n \alpha_i k(\cdot, x_i) 
: n \in \mathbb{N},\ x_i \in \mathcal{X},\ \alpha_i \in \mathbb{R} 
\}
$$











Para $f = \sum_i \alpha_i k(\cdot, u_i)$ y $g = \sum_j \beta_j k(\cdot, v_j)$, se define:

$$\langle f, g \rangle = \sum_{i,j} \alpha_i \beta_j\, k(u_i, v_j)$$

Nótese que:

$$\langle f,\, k(\cdot,x) \rangle = \sum_i \alpha_i\, k(x, u_i) = f(x)$$

**Verificación del producto interno:**

1. **Simetría:** $\langle f,g \rangle = \sum_{ij} \alpha_i \beta_j k(u_i,v_j) = \sum_{ij} \beta_j \alpha_i k(v_j,u_i) = \langle g,f \rangle$
2. **Bilinealidad:** $\langle f,g \rangle = \sum_i \alpha_i g(u_i) = \sum_j \beta_j f(v_j)$
3. **Definición positiva:** $\langle f,f \rangle = \alpha^\top K \alpha \geq 0$ con igualdad ssi $f = 0$

Del punto 3 se deriva:

- $\langle f,g \rangle^2 \leq \langle f,f \rangle \langle g,g \rangle$
  - *Prueba:* $\forall a \in \mathbb{R}$, $\langle af+g,\, af+g \rangle = a^2\langle f,f\rangle + 2a\langle f,g\rangle + \langle g,g\rangle \geq 0$. El discriminante de esta cuadrática es no positivo, luego $\langle f,g\rangle^2 - \langle f,f\rangle\langle g,g\rangle \leq 0$.
- $|f(x)|^2 = \langle k(\cdot,x), f\rangle^2 \leq k(x,x)\langle f,f\rangle$, lo que implica que si $\langle f,f\rangle = 0$ entonces $f \equiv 0$.

---

> **📘 Definición.** Para un (compacto) $\mathcal{X} \subseteq \mathbb{R}^d$ y un espacio de Hilbert $\mathcal{H}$ de funciones $f: \mathcal{X} \to \mathbb{R}$, decimos que $\mathcal{H}$ es un *RKHS* si $\exists\, k: \mathcal{X} \to \mathbb{R}$ tal que:
> 1. $k$ tiene la propiedad reproductora: $f(x) = \langle f(\cdot),\, k(\cdot,x) \rangle$
> 2. $k$ genera $\mathcal{H} = \text{span}\{k(\cdot,x) : x \in \mathcal{X}\}$

---

### Teorema de Mercer

> **📐 Teorema 1.1 (Mercer).** Suponga que $k$ es un núcleo continuo semidefinido positivo en un conjunto compacto $\mathcal{X}$, y el operador integral $T_k: L_2(\mathcal{X}) \to L_2(\mathcal{X})$ definido por
>
> $$(T_k f)(\cdot) = \int_{\mathcal{X}} k(\cdot, x)\, f(x)\, dx$$
>
> es semidefinido positivo, es decir, $\forall f \in L_2(\mathcal{X})$:
>
> $$\int_{\mathcal{X}} k(u,v)\, f(u)\, f(v)\, du\, dv \geq 0$$
>
> Entonces existe una base ortonormal $\{\psi_i\}$ de $L_2(\mathcal{X})$ compuesta por autofunciones de $T_k$, con autovalores no negativos $\{\lambda_i\}$. Las autofunciones con autovalores no nulos son continuas en $\mathcal{X}$ y $k$ admite la representación:
>
> $$k(u,v) = \sum_{i=1}^{\infty} \lambda_i\, \psi_i(u)\, \psi_i(v)$$
>
> donde la convergencia es absoluta y uniforme:
>
> $$\lim_{n \to \infty} \sup_{u,v} \left| k(u,v) - \sum_{i=1}^{n} \lambda_i\, \psi_i(u)\, \psi_i(v) \right| = 0$$

---

**Analogía en el caso finito** ($\mathcal{X} = \{x_1, \ldots, x_n\}$):

Sea $K_{ij} = k(x_i, x_j)$, y $f: \mathcal{X} \to \mathbb{R}^n$ con $f_i = f(x_i)$. Entonces:

$$T_k f = \sum_{i=1}^n k(\cdot, x_i)\, f_i$$

$$\forall f,\; f^\top K f \geq 0 \implies K \succeq 0 \implies K = \sum_k \lambda_k v_k v_k^\top$$

Por lo tanto:

$$k(x_i, x_j) = K_{ij} = (V\Lambda V^\top)_{ij} = \sum_k \lambda_k (v_k)_i (v_k)_j = \sum_k \lambda_k \psi_k(x_i)\psi_k(x_j)$$

---

### Condiciones equivalentes sobre $k$

Para $k$ simétrico continuo en $\mathcal{X}$ compacto, las siguientes condiciones son equivalentes:

| # | Condición |
|---|-----------|
| 1 | Toda matriz de Gram es semidefinida positiva |
| 2 | $T_k$ es semidefinido positivo |
| 3 | $k(u,v) = \sum_i \lambda_i\, \psi_i(u)\, \psi_i(v)$ |
| 4 | $k$ es el núcleo reproductor de un RKHS de funciones en $\mathcal{X}$ |

---

<div align="center">

*Antono Lam García - Universidad Nacional de Ingeniería — Doctorado en Ciencias e Ingeniería Estadística*  
© 2023 Todos los derechos reservados

</div>
