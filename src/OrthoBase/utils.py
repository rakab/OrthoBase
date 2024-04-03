import z3

__all__ = [
       'getNewId',
       'getDifferentSolution',
       'getDifferentSolutionMatrix',
       'member_of',
       'generate_pairings',
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

def generate_pairings(perms, nconns):
    """
    Create all possible connections of indices.
    Primaraly designed to fill in the treapezoid in Eq. (4.13) of
    arXiv:1207.0609.

    Here is a recursive version, which is more readable, but less efficient:
    def generate_pairings(perms,res,i):
        if not len(perms):
            if i == 3:
                print(res)
            return
        for j,perm in enumerate(perms):
            remaining = [p for p in perms[j+1:] if ((p[0]!=perm[0] and p[0]!=perm[1]) and (p[1]!=perm[0] and p[1]!=perm[1]))]
            res[i] = perm
            generate_pairings(remaining,res,i+1)
    """
    stack = [(perms, [], 0)]
    result = list()

    while stack:
        current_perms, res, i = stack.pop()

        if not len(current_perms):
            if i == nconns:
                result.append(res)
            continue

        for j, perm in enumerate(current_perms):
            remaining = [p for p in current_perms[j+1:] if ((p[0]!=perm[0] and p[0]!=perm[1]) and (p[1]!=perm[0] and p[1]!=perm[1]))]
            new_res = res + [perm]
            stack.append((remaining, new_res, i+1))
    return result
