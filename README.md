# Backward-bifurcation
For epidemiological models, calculating the backward bifurcation represents having knowledge about disease control. In the article [GitHub Pages][(https://www.aimsciences.org/article/doi/10.3934/mbe.2004.1.361)], they employ a technique that involves 
determining the values of a and b to know whether a backward bifurcation exists or not.

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

see some examples:
1. Consider the model describe in the [Article](https://www.sciencedirect.com/science/article/abs/pii/S0025556400000031), luego
   ```python
a=5
b=10
```
    
