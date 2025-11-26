from sympy import (
    symbols, diff, integrate, sin, cos, tan,
    exp, log, Abs, S, Interval
)
from sympy.parsing.sympy_parser import (
    parse_expr,
    standard_transformations,
    implicit_multiplication_application
)
from sympy.calculus.util import continuous_domain, function_range

# Allow implicit multiplication like "4x", "2sin(x)"
TRANSFORMS = standard_transformations + (implicit_multiplication_application,)


# ------------------------------
# 🔵 SUPERCLASS: Generic Function
# ------------------------------
class Function:
    def __init__(self, expr: str):
        # allow "^" for power
        self.expr_str = expr.replace("^", "**")
        # parse with implicit multiplication enabled
        self.expr = parse_expr(self.expr_str, transformations=TRANSFORMS)
        # detect variables
        self.vars = sorted(self.expr.free_symbols, key=lambda s: s.name)

    def __str__(self):
        return str(self.expr)

    def tokens(self):
        """Return expression tokens split by whitespace from original string."""
        return self.expr_str.split()

    def evaluater(self, **values):
        """
        Evaluate expression by substituting values into variables.
        Example:
            f = Function("x - 2w + 3z")
            f.evaluater(x=1, w=2, z=3)
        """
        evaluated = self.expr.subs(values)
        if evaluated.free_symbols:
            return evaluated
        try:
            return float(evaluated.evalf())
        except TypeError:
            return evaluated

    def output(self):
        return self.expr

    def number_of_variables(self):
        return len(self.vars)

    # ---------- internal helper ----------
    def _ensure_symbol(self, var):
        """Accept either string or SymPy Symbol; always return Symbol."""
        if isinstance(var, str):
            return symbols(var)
        return var

    # ---------- Calculus ----------
    def derivative(self, var=None):
        if var is None:
            if len(self.vars) != 1:
                raise ValueError("Specify a variable for partial derivative.")
            var = self.vars[0]
        else:
            var = self._ensure_symbol(var)
        return self.__class__(str(diff(self.expr, var)))

    def integral(self, var=None):
        if var is None:
            if len(self.vars) != 1:
                raise ValueError("Specify variable for integration.")
            var = self.vars[0]
        else:
            var = self._ensure_symbol(var)
        return self.__class__(str(integrate(self.expr, var)))

    # ---------- Domain & Range ----------
    def domain(self, var=None):
        """
        Generic domain over ℝ using SymPy's continuous_domain.
        Works for most single-variable expressions.
        """
        if var is None:
            if len(self.vars) != 1:
                raise ValueError("Domain currently only for single-variable functions; pass var.")
            var = self.vars[0]
        else:
            var = self._ensure_symbol(var)
        return continuous_domain(self.expr, var, S.Reals)

    def range(self, var=None, domain=None):
        """
        Generic range using SymPy's function_range.
        For some functions (like Abs), SymPy may fail → we catch that.
        """
        if var is None:
            if len(self.vars) != 1:
                raise ValueError("Range currently only for single-variable functions; pass var.")
            var = self.vars[0]
        else:
            var = self._ensure_symbol(var)

        if domain is None:
            domain = self.domain(var)

        try:
            return function_range(self.expr, var, domain)
        except NotImplementedError:
            return "Range could not be computed symbolically."

    # You also wanted Domain/Range names:
    def Domain(self, var=None):
        return self.domain(var)

    def Range(self, var=None, domain=None):
        return self.range(var, domain)

    # ---------- Operator Overloading ----------
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
        poly = self.expr.as_poly()
        if poly is None:
            raise ValueError("Expression is not a polynomial.")
        return poly.degree()

    def domain(self, var=None):
        """
        For real polynomials in one variable: domain is ℝ.
        """
        return S.Reals

    def range(self, var=None, domain=None):
        """
        Use SymPy's function_range over ℝ.
        """
        if var is None:
            if len(self.vars) != 1:
                raise ValueError("Range only implemented for single-variable polynomials; pass var.")
            var = self.vars[0]
        else:
            var = self._ensure_symbol(var)

        if domain is None:
            domain = S.Reals

        return function_range(self.expr, var, domain)


# -----------------------------------
# 🟢 SUBCLASS: Trigonometric Function
# -----------------------------------
class TrigFunction(Function):
    def simplify_trig(self):
        return TrigFunction(str(self.expr.simplify()))

    def to_polynomial_if_possible(self):
        return TrigFunction(str(self.expr.simplify().rewrite('sqrt')))

    # domain() and range() are good using the base class implementations.


# -----------------------------------
# 🟣 SUBCLASS: Rational Function
# -----------------------------------
class RationalFunction(Function):
    def numerator(self):
        return self.expr.as_numer_denom()[0]

    def denominator(self):
        return self.expr.as_numer_denom()[1]

    def domain(self, var=None):
        """
        Domain: all real numbers where denominator ≠ 0.
        Base class continuous_domain already handles this.
        """
        return super().domain(var)

    def range(self, var=None, domain=None):
        """
        Range via function_range over that domain.
        """
        return super().range(var, domain)

    def vertical_asymptotes(self):
        den = self.denominator()
        return f"Solutions of: {den} = 0"

    def horizontal_asymptote(self):
        num, den = self.expr.as_numer_denom()
        try:
            p_num = num.as_poly()
            p_den = den.as_poly()
        except Exception:
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


# -----------------------------------
# 🔵 SUBCLASS: Exponential Function
# -----------------------------------
class ExponentialFunction(Function):
    # generic domain/range via base class works fine for exp-type expressions
    pass


# -----------------------------------
# 🟤 SUBCLASS: Logarithmic Function
# -----------------------------------
class LogFunction(Function):
    # continuous_domain does the x > 0 restriction automatically
    pass


# -----------------------------------
# 🟠 SUBCLASS: Absolute Value Function
# -----------------------------------
class AbsoluteFunction(Function):
    def range(self, var=None, domain=None):
        """
        Special-case Abs for a clean result.
        For Abs(linear), the range is [0, ∞).
        To keep things simple and robust, for Abs(something in ℝ) we return [0, ∞).
        """
        if var is None:
            if len(self.vars) != 1:
                raise ValueError("Range only implemented for single-variable |.|; pass var.")
            var = self.vars[0]
        else:
            var = self._ensure_symbol(var)

        if domain is None:
            domain = self.domain(var)

        # If it's literally Abs(something), we know output is ≥ 0
        if isinstance(self.expr, Abs):
            return Interval(0, S.Infinity)

        # Fallback to SymPy's general function_range
        try:
            return function_range(self.expr, var, domain)
        except NotImplementedError:
            return "Range could not be computed symbolically."


# ------------------------------
# TESTING
# ------------------------------
if __name__ == "__main__":
    # Generic Function
    y = Function("x - 2w + 3z")
    q = Function("x - 2w + 3z")

    print("y =", y, " | #vars =", y.number_of_variables())
    print("tokens:", y.tokens())
    print("evaluate y at (x=1, w=2, z=3):", y.evaluater(x=1, w=2, z=3))

    print("y + q =", y + q)
    print("y - q =", y - q)

    # Single-variable generic function: f(x) = x^2 + 1
    f = Function("x^2 + 1")
    print("\nf(x) =", f)
    print("Domain(f):", f.Domain())
    print("Range(f):", f.Range())

    print("f'(x):", f.derivative("x"))
    print("∫ f dx:", f.integral("x"))

    # Polynomial
    p = PolynomialFunction("x^3 - 4x + 1")
    print("\nPolynomial p(x) =", p)
    print("degree(p):", p.degree())
    print("Domain(p):", p.Domain())
    print("Range(p):", p.Range())
    print("p'(x):", p.derivative())

    # Trigonometric
    t = TrigFunction("sin(x) + cos(x)")
    print("\nTrig t(x) =", t)
    print("t simplified:", t.simplify_trig())
    print("t'(x):", t.derivative())
    print("Domain(t):", t.Domain())
    print("Range(t):", t.Range())

    # Rational
    r = RationalFunction("(x^2 + 1)/(x - 2)")
    print("\nRational r(x) =", r)
    print("vertical asymptotes:", r.vertical_asymptotes())
    print("horizontal asymptote:", r.horizontal_asymptote())
    print("Domain(r):", r.Domain())
    print("Range(r):", r.Range())

    # Exponential
    e = ExponentialFunction("exp(x)")
    print("\nExponential e(x) =", e)
    print("Domain(e):", e.Domain())
    print("Range(e):", e.Range())

    # Logarithmic
    g = LogFunction("log(x)")
    print("\nLog g(x) =", g)
    print("Domain(g):", g.Domain())
    print("Range(g):", g.Range())

    # Absolute
    h = AbsoluteFunction("Abs(x - 1)")
    print("\nAbsolute h(x) =", h)
    print("Domain(h):", h.Domain())
    print("Range(h):", h.Range())
