from dolfin import *

if not has_cgal():
    print "DOLFIN must be compiled with CGAL to run this demo."
    exit(0)

n1, n2 = 20, 20

# center of the ellipses
x1, y1 = 0.0, 0.0
x2, y2 = 0.0, 0.0
# major and minor axes
a1, b1 = 2.0, 1.0
a2, b2 = 1.0, 0.5

e1 = Ellipse(x1, y1, a1, b1, n1)
e2 = Ellipse(x2, y2, a2, b2, n2)
g2d = e1 - e2

mesh = Mesh(g2d, 20)
plot(mesh)

V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)

a = inner(grad(u), grad(v))*dx
L = v*dx

bc = DirichletBC(V, 0, "on_boundary")

u = Function(V)

# Functional
M = u*dx
tol = 1.0e-5

# Solve equation a = L with respect to u and the given boundary
# conditions, such that the estimated error (measured in M) is less
# than tol
problem = LinearVariationalProblem(a, L, u, bc)
solver = AdaptiveLinearVariationalSolver(problem, M)
solver.parameters["error_control"]["dual_variational_solver"]["linear_solver"] = "cg"
solver.solve(tol)

solver.summary()

# Plot solution(s)
plot(u.root_node(), title="Solution on initial mesh")
plot(u.leaf_node(), title="Solution on final mesh")
interactive()
