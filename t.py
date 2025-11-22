from sympy import symbols, sympify, diff, integrate, sin, cos, tan

# ------------------------------
# 🔵 SUPERCLASS: Generic Function
# ------------------------------
class Function:
    def __init__(self, expr: str):
        self.expr_str = expr.replace("^", "**")
        self.expr = sympify(self.expr_str)
        self.vars = sorted(self.expr.free_symbols, key=lambda s: s.name)

    def __str__(self):
        return str(self.expr)

    def evaluater(self, **values):
        return float(self.expr.subs(values))

    def output(self):
        return self.expr

    def number_of_variables(self):
        return len(self.vars)

    def derivative(self, var=None):
        if var is None:
            if len(self.vars) != 1:
                raise ValueError("Specify a variable for partial derivative.")
            var = self.vars[0]
        else:
            var = symbols(var)

        return self.__class__(str(diff(self.expr, var)))

    def integral(self, var=None):
        if var is None:
            if len(self.vars) != 1:
                raise ValueError("Specify variable for integration.")
            var = self.vars[0]
        else:
            var = symbols(var)

        return self.__class__(str(integrate(self.expr, var)))

    # Arithmetic operators return same class (important for polymorphism)
    def __add__(self, other):
        return self.__class__(str(self.expr + other.expr))

    def __sub__(self, other):
        return self.__class__(str(self.expr - other.expr))

    def __mul__(self, other):
        return self.__class__(str(self.expr * other.expr))

    def __truediv__(self, other):
        return self.__class__(str(self.expr / other.expr))

    def __pow__(self, power):
        return self.__class__(str(self.expr ** power))


# -----------------------------------
# 🔴 SUBCLASS: Polynomial Function
# -----------------------------------
class PolynomialFunction(Function):
    def degree(self):
        """Return polynomial degree (only if single-variable)."""
        if len(self.vars) != 1:
            raise ValueError("Degree only defined for single-variable polynomials.")

        return self.expr.as_poly().degree()


# -----------------------------------
# 🟢 SUBCLASS: Trigonometric Function
# -----------------------------------
class TrigFunction(Function):
    def simplify_trig(self):
        """Simplify the trigonometric expression."""
        return TrigFunction(str(self.expr.simplify()))

    def to_polynomial_if_possible(self):
        """Convert sin^2 + cos^2, etc., into simpler algebraic forms."""
        return TrigFunction(str(self.expr.simplify().rewrite('sqrt')))


# -----------------------------------
# 🟣 SUBCLASS: Rational Function
# -----------------------------------
class RationalFunction(Function):
    def numerator(self):
        return self.expr.as_numer_denom()[0]

    def denominator(self):
        return self.expr.as_numer_denom()[1]

    def vertical_asymptotes(self):
        """Points where denominator = 0."""
        den = self.denominator()
        return f"Solutions of: {den} = 0"

    def horizontal_asymptote(self):
        """Based on polynomial degree comparison."""
        num, den = self.expr.as_numer_denom()
        try:
            p_num = num.as_poly()
            p_den = den.as_poly()
        except:
            return "Not a rational polynomial function."

        if p_num.degree() < p_den.degree():
            return "y = 0"
        elif p_num.degree() == p_den.degree():
            return f"y = {p_num.LC()}/{p_den.LC()}"
        else:
            return "No horizontal asymptote (slant or higher)."



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
