import z3

__all__ = [
       'getNewId',
       'getDifferentSolution',
       'getDifferentSolutionMatrix',
       'member_of',
        ]

def getNewId():
  return uuid.uuid4().int

def getDifferentSolution(sol,mod, *params):
  for t in params:
    sol.add(Or([t[i] != mod.eval(t[i]) for i in range(len(t))]))

# special case for a matrix; requires number of rows and columns
def getDifferentSolutionMatrix(sol,mod, x, rows, cols):
    sol.add(z3.Or([x[i,j] != mod.eval(x[i,j]) for i in range(rows) for j in range(cols)]))

def member_of(sol, e, v):
    sol.add(z3.Or([e == i for i in v]))
