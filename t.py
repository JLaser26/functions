from sympy import symbols, sympify, diff, integrate

class function:
    def __init__(self, expr: str):
        """
        expr: string formula like "x - 2*w + 3*z"
        """
        self.expr_str = expr.replace("^", "**")     # allow ^ as power
        self.expr = sympify(self.expr_str)          # convert to a symbolic expression
        self.vars = sorted(list(self.expr.free_symbols), key=lambda s: s.name)

    def __str__(self):
        return str(self.expr)

    def tokens(self):
        """Return expression tokens split by space."""
        return self.expr_str.split()

    def evaluater(self, **values):
        """
        Evaluate the expression by substituting numerical values.
        Example: f.evaluater(x=2, w=1, z=3)
        """
        return float(self.expr.subs(values))

    def output(self):
        return self.expr

    def Domain(self):
        """
        Domain of a sympy expression. (General case is complicated.)
        For now: return all real numbers except where expression is undefined.
        """
        return f"All real numbers except points where the expression is undefined."

    def Range(self):
        """General range is very hard; placeholder."""
        return "Range depends on variable domain (not computed automatically)."

    def number_of_variables(self):
        return len(self.vars)

    def derivative(self, var=None):
        """
        Derivative w.r.t. one variable.
        If var is None and only one variable exists → use that.
        """
        if var is None:
            if len(self.vars) != 1:
                raise ValueError("Specify variable for partial derivative.")
            var = self.vars[0]
        else:
            var = symbols(var)

        return function(str(diff(self.expr, var)))

    def integral(self, var=None):
        """Indefinite integral."""
        if var is None:
            if len(self.vars) != 1:
                raise ValueError("Specify variable for integration.")
            var = self.vars[0]
        else:
            var = symbols(var)

        return function(str(integrate(self.expr, var)))

    # --- Operator overloading ---
    def __add__(self, other):
        return function(str(self.expr + other.expr))

    def __sub__(self, other):
        return function(str(self.expr - other.expr))

    def __mul__(self, other):
        return function(str(self.expr * other.expr))

    def __truediv__(self, other):
        return function(str(self.expr / other.expr))

    def __pow__(self, power):
        return function(str(self.expr ** power))

#testing
y = function("x - 2*w + 3*z")
q = function("x - 2*w + 3*z")

print(y, y.number_of_variables())
print("tokens:", y.tokens())

print("evaluate:", y.evaluater(x=1, w=2, z=3))

print("y + q =", y + q)
print("y - q =", y - q)

# Derivative
print("dy/dx =", y.derivative("x"))

# Integral
print("∫ y dx =", y.integral("x"))
