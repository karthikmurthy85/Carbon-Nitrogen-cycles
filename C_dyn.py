#%%
from numpy import negative
import sympy as sm
from sympy.matrices import expressions
from sympy.simplify.simplify import simplify
I, k_d, alpha, k_l, r_h, k_h,r_r,k_d = sm.symbols('I, k_d, alpha, k_l, r_h, k_h,r_r,k_d', positive=True)
C_b, C_l, C_h = sm.symbols('C_b, C_l, C_h', negative=False)
#%%
DCl = I + C_b*k_d - C_b*k_l*C_l
DCh = r_h*C_b*k_l*C_l - C_b*k_h*C_h
DCb = (1- r_h - r_r)*C_b*k_l*C_l + (1- r_r)*C_b*k_h*C_h -  C_b*k_d
DCb
#%%
eq1 = sm.Eq(DCl, 0)
equiCb = (sm.solve(eq1, C_b))
equiCb = (I/(C_l*k_l - k_d))
equiCb
#%%
eq2 = sm.Eq(DCh, 0)
equiCh = (sm.solve(eq2, C_h))
equiCh = C_l*k_l*r_h/k_h
equiCh
#%%
eq3 = sm.Eq(DCb, 0)
equiCl = (sm.solve(eq3, C_l))
equiCl = sm.simplify(-C_h*k_h*r_r + C_h*k_h - k_d)/(k_l*(r_h + r_r - 1))
equiCl
equiCl.subs(C_h, equiCh)
#%%
equilibria = (sm.solve( (eq1, eq2, eq3), C_l, C_h, C_b))
print(equilibria)
eqCl = equilibria[0][0]
eqCh = equilibria[0][1]
eqCb = equilibria[0][2]
eqCb
#%%
F = sm.Matrix([DCl,DCh,DCb])
F
#%%
F.jacobian([C_l,C_h,C_b])
#%%
J1 = F.jacobian([C_l,C_h,C_b]).subs([(C_l,equilibria[0][0]), (C_h,equilibria[0][1]), (C_b,equilibria[0][2])])
J1 = sm.simplify(J1)
#TrJ1 = J1[0,0] + J1[1,1] + J1[2,2]
#sm.simplify(TrJ1)
J1
#%%
a, b, c, d, e, f, g, lb = sm.symbols('a, b, c, d, e, f, g, lb')
F1 = sm.Matrix([[a-lb,0,b],[c,d-lb,0],[e,f,-lb]])
exp1 = F1.det()
exp11 = sm.Poly(exp1,lb)
cfs = (exp11.coeffs())
cfs[3]*-1
#cfs[3]
#cfs[1]
#cfs
#%%
from sympy.plotting import plot 
a3, a2, a1, a0, x = sm.symbols('a3, a2, a1, a0, x', real=True)
#fx = a3*x**3 + a2*x**2+a1*x+a0
#plot((1*x**3 + 1*x**2+1*x+100), (1*x**3 + 10*x**2+10*x+0))
plot(x, (1*x**3 + 25*x**2 + 40*x - 40,), (1*x**3 + 1*x**2 + 1*x - 1), (x, -50, 50))
# %%
#Understand the roots
A = cfs[1]*-1; B = cfs[2]*-1; C = cfs[3]*-1
p = B - (A**2/3)
q = (2*A**3/27) - (A*B/3) + C
del1 = (q**2/4) + (p**3/27)
#del = del1.subs([(a, J1[0]), (b, J1[2]), (c, J1[3]) ,  (e, J1[6]), (f, J1[7]), (g, J1[8])])
#del = del1.subs([(q, J1[0])])
dl = del1.subs([(a, J1[0]), (b, J1[2]), (c, J1[3]) , (d, J1[4]), (e, J1[6]), (f, J1[7])])
dl1 = sm.simplify(dl)
#dl2 = sm.diff(dl1, r_r)
dl1
#%%
B
# %%
intc = sm.simplify(J1[4]*J1[6] - (J1[3]*J1[7]))
intc