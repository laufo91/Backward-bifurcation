# Backward-bifurcation
For epidemiological models, calculating the backward bifurcation represents having knowledge about disease control. In the [article](https://www.aimsciences.org/article/doi/10.3934/mbe.2004.1.361), they employ a technique that involves 
determining the values of $a$ and $b$ to know whether a backward bifurcation exists or not.

I use python, the libreries are
```python
from sympy.solvers import solve
from sympy import Symbol
import numpy as np
from sympy import sin, cos, Matrix
from sympy.matrices import hessian
from sympy import *
import sympy as sp
```

Consider the model described in the [Article](https://www.sciencedirect.com/science/article/abs/pii/S0025556400000031), then
1. Define the variaveis
```python
I = Symbol("I")
V = Symbol("V")
N = Symbol("N")
sigma=Symbol("sigma")
beta=Symbol("beta")
mu=Symbol("mu")
c=Symbol("c")
phi=Symbol("phi")
theta=Symbol("theta")
```
2. Define the system
```python
X = Matrix([(N-I-(1-sigma)*V)*(beta*I)/N-(mu+c)*I,phi*(N-I-V)-sigma*beta*V*I/N-(mu+theta)*V])
```

3. Determine the Jacobian with respect all variaveis
```python
Y = Matrix([I,V])
Z=X.jacobian(Y)
```
4. Subtitute the endemic point
```python
Z=Z.evalf(subs={I:0,V:phi*N/(mu+phi+theta)})
Z=simplify(Z)
```
5. In this model, $R_0=\frac{\beta}{\mu+c}\frac{\mu+θ+\sigma \phi}{\mu+\theta+\phi}$, equal $1$ and solve for any variable. In this case $\beta=(\mu+c) \frac{\mu+\theta+\phi}{\mu+θ+\sigma \phi}$, then substitue en last matrix.
```python
Z=Z.evalf(subs={beta:(mu+c)*(mu+theta+phi)/(mu+theta+sigma*phi)})
Z=simplify(Z)
```


