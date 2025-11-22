from sympy import symbols, sympify, diff, integrate, sin, cos, tan
from sympy.parsing.sympy_parser import (
    parse_expr,
    standard_transformations,
    implicit_multiplication_application
)

# Transformations that allow implicit multiplication like "4x", "2sin(x)"
TRANSFORMS = standard_transformations + (implicit_multiplication_application,)


# ------------------------------
# 🔵 SUPERCLASS: Generic Function
# ------------------------------
class Function:
    def __init__(self, expr: str):
        # allow ^ for power
        self.expr_str = expr.replace("^", "**")

        # Parse expression with implicit multiplication enabled
        self.expr = parse_expr(self.expr_str, transformations=TRANSFORMS)

        # Detect variables
        self.vars = sorted(self.expr.free_symbols, key=lambda s: s.name)

    def __str__(self):
        return str(self.expr)

    def tokens(self):
        return self.expr_str.split()

    def evaluater(self, **values):
        evaluated = self.expr.subs(values)
        if evaluated.free_symbols:
            return evaluated
        return float(evaluated)

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
        if len(self.vars) != 1:
            raise ValueError("Degree only defined for single-variable polynomials.")
        poly = self.expr.as_poly()
        if poly is None:
            raise ValueError("Expression is not a polynomial.")
        return poly.degree()


# -----------------------------------
# 🟢 SUBCLASS: Trigonometric Function
# -----------------------------------
class TrigFunction(Function):
    def simplify_trig(self):
        return TrigFunction(str(self.expr.simplify()))

    def to_polynomial_if_possible(self):
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
        den = self.denominator()
        return f"Solutions of: {den} = 0"

    def horizontal_asymptote(self):
        num, den = self.expr.as_numer_denom()

        try:
            p_num = num.as_poly()
            p_den = den.as_poly()
        except:
            return "Not a rational polynomial function."

        if p_num is None or p_den is None:
            return "Not a rational polynomial function."

        deg_num = p_num.degree()
        deg_den = p_den.degree()

        if deg_num < deg_den:
            return "y = 0"
        elif deg_num == deg_den:
            return f"y = {p_num.LC()}/{p_den.LC()}"
        else:
            return "No horizontal asymptote (slant or higher)."


# ------------------------------
# TESTING
# ------------------------------

y = Function("x - 2w + 3z")
q = Function("x - 2w + 3z")

print(y, y.number_of_variables())
print("tokens:", y.tokens())
print("evaluate:", y.evaluater(x=1, w=2, z=3))

print("y + q =", y + q)
print("y - q =", y - q)

print("dy/dx =", y.derivative("x"))
print("∫ y dx =", y.integral("x"))

p = PolynomialFunction("x^3 - 4x + 1")
print("p:", p)
print("degree:", p.degree())
print("p' =", p.derivative())

t = TrigFunction("sin(x) + cos(x)")
print("t simplified:", t.simplify_trig())
print("t' =", t.derivative())

r = RationalFunction("(x^2 + 1)/(x - 2)")
print("vertical asymptotes:", r.vertical_asymptotes())
print("horizontal asymptote:", r.horizontal_asymptote())
