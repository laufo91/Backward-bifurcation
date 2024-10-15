# Backward-bifurcation
In epidemiological models, calculating the backward bifurcation provides valuable insights into disease control. In the [article](https://www.aimsciences.org/article/doi/10.3934/mbe.2004.1.361) the authors use a technique that involves determining the values of $a$ and $b$ to assess the existence of a backward bifurcation.

To describe the backward bifurcation, you need to know the disease-free equilibrium point and $R_0$. I used Python, and the libraries are:
```python
from sympy.solvers import solve
from sympy import Symbol
import numpy as np
from sympy import sin, cos, Matrix
from sympy.matrices import hessian
from sympy import *
import sympy as sp
```

Consider the model described in the [Article](https://www.sciencedirect.com/science/article/abs/pii/S0025556400000031), the ideas of programming next, but remember with each system de values changes.
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

3. Determine the Jacobian with respect to all variables.
```python
Y = Matrix([I,V])
Z=X.jacobian(Y)
```
4. Substitute the endemic point into the matrix.
```python
Z=Z.evalf(subs={I:0,V:phi*N/(mu+phi+theta)})
Z=simplify(Z)
```
5. In this model, $R_0=\frac{\beta}{\mu+c} \frac{\mu+\theta+\sigma \phi}{\mu+\theta+\phi}$ is set equal to 1, and solved for any variable. In this case, $\beta = (\mu+c) \frac{\mu+\theta+\phi}{\mu+\theta+\sigma \phi}$. Then, substitute this into the last matrix.
```python
Z=Z.evalf(subs={beta:(mu+c)*(mu+theta+phi)/(mu+theta+sigma*phi)})
Z=simplify(Z)
```
Rembember, if all calcultions are perfect, the determinat sholud be $0$.

6. Determine the left and right eigenvectors associated with the eigenvalue $0$.
```python
v=Z.left_eigenvects()[0][2][0]
w=Z.eigenvects()[0][2][0]
```
7. Verify that the dot product is $1$; if it is negative, apply the following:
```python
if v.dot(w)!=1:
  l=Symbol("l")
  escalar=solve((l*v).dot(w)-1, l)
  w=escalar[0]*w
v=simplify(v)
w=simplify(w)
```
8. To determine the value of $b$, proceed as follows:
```python
y = Matrix([beta])# beta was the value solved in the item 5. 
Matrix_b=X.jacobian(y)
Y = Matrix([I,V])
Matrix_b1=Matrix_b.jacobian(Y)
Matrix_b1=Matrix_b1.evalf(subs={I:0,V:phi*N/(mu+phi+theta),beta:(mu+c)*(mu+theta+phi)/(mu+theta+sigma*phi)})# Substitute beta as item 5 and endemic point
Matrix_b1
b=0
for j in range(2):#El 2 es el numero de ecuaciones
  for i in range(2):
    if Matrix_b1[j,i]!=0:
      b=b+Matrix_b1[j,i]*v[j]*w[i]
b=simplify(b)
print("b")
b
```
9. To determine the value of $a$, proceed as follows:
```python
aa=0
for pp in range(2):
  XX=Matrix([X[pp]])
  matrix_a=hessian(XX,(I,V))
  matrix_a=matrix_a.evalf(subs={I:0,V:phi*N/(mu+phi+theta),beta:(mu+c)*(mu+theta+phi)/(mu+theta+sigma*phi)}) # Substitute beta as item 5 and endemic point
  for i in range(0,2):
    for j in range(0,2):
      if matrix_a[i,j]!=0:
        aa=aa+matrix_a[i,j]*v[pp]*w[i]*w[j]
        aa=simplify(aa)
print("a")
aa=simplify(aa)
aa
```
If $a>0$ and $b>0$ the backward bifurcation exits.

