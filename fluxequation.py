import sympy
sympy.init_printing()

u_max, u_star, rho_max, rho_star, A, B = sympy.symbols('u_max u_star rho_max rho_star A B')

eq1 = sympy.Eq( 0,      u_max*rho_max*(1 - A*rho_max-B*rho_max**2) )
eq2 = sympy.Eq( 0,      u_max*(1 - 2*A*rho_star-3*B*rho_star**2) )
eq3 = sympy.Eq( u_star, u_max*(1 - A*rho_star - B*rho_star**2) )

## SUBTRACT: EQ2 - 3 EQ3
eq4 = sympy.Eq(eq2.lhs - 3*eq3.lhs, eq2.rhs - 3*eq3.rhs)

rho_sol = sympy.solve(eq4,rho_star)[0]
B_sol = sympy.solve(eq1,B)[0]

quadA = eq2.subs([(rho_star, rho_sol), (B,B_sol)])
A_sol = sympy.solve(quadA, A)


# NUMERIC SUBSTITUTIONS

aval1 = A_sol[0].evalf(subs={u_star: 1.5, u_max:2.0, rho_max:15} )
print("========================== A =====================")
print(aval1)


bval = B_sol.evalf(subs={rho_max:15, A:aval1} )
print("========================== B =====================")
print(bval)
