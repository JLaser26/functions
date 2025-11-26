from sympy import (
    symbols, diff, integrate,
    sin, cos, tan, exp, log, Abs,
    S
)
from sympy.parsing.sympy_parser import (
    parse_expr,
    standard_transformations,
    implicit_multiplication_application
)
from sympy.calculus.util import continuous_domain, function_range

# Transformations that allow implicit multiplication like "4x", "2sin(x)"
TRANSFORMS = standard_transformations + (implicit_multiplication_application,)


# ------------------------------
# 🔵 SUPERCLASS: Generic Function
# ------------------------------
class Function:
    def __init__(self, expr: str):
        # Allow '^' for powers
        self.expr_str = expr.replace("^", "**")

        # Parse expression with implicit multiplication enabled
        self.expr = parse_expr(self.expr_str, transformations=TRANSFORMS)

        # Detect variables (symbols)
        self.vars = sorted(self.expr.free_symbols, key=lambda s: s.name)

    def __str__(self):
        return str(self.expr)

    def tokens(self):
        """Return expression tokens split by whitespace from the original string."""
        return self.expr_str.split()

    def evaluater(self, **values):
        """
        Evaluate the expression by substituting numerical values.
        Example:
            f = Function("x - 2w + 3z")
            f.evaluater(x=1, w=2, z=3)
        """
        evaluated = self.expr.subs(values)
        if evaluated.free_symbols:
            # still symbolic
            return evaluated
        try:
            return float(evaluated.evalf())
        except TypeError:
            return evaluated

    def output(self):
        return self.expr

    def number_of_variables(self):
        return len(self.vars)

    # ---------- Calculus ----------
    def derivative(self, var=None):
        """
        Derivative w.r.t. specified variable.
        If var is None and there's only one variable, use that.
        """
        if var is None:
            if len(self.vars) != 1:
                raise ValueError("Specify a variable for partial derivative.")
            var = self.vars[0]
        else:
            var = symbols(var)
        return self.__class__(str(diff(self.expr, var)))

    def integral(self, var=None):
        """
        Indefinite integral w.r.t. specified variable.
        If var is None and there's only one variable, use that.
        """
        if var is None:
            if len(self.vars) != 1:
                raise ValueError("Specify variable for integration.")
            var = self.vars[0]
        else:
            var = symbols(var)
        return self.__class__(str(integrate(self.expr, var)))

    # ---------- Domain & Range ----------
    def domain(self, var=None):
        """
        Generic domain over ℝ using SymPy's continuous_domain.
        """
        if var is None:
            if len(self.vars) != 1:
                raise ValueError("Domain currently implemented only for single-variable functions; specify var.")
            var = self.vars[0]
        else:
            var = symbols(var)

        dom = continuous_domain(self.expr, var, S.Reals)
        return dom

    def range(self, var=None, domain=None):
        """
        Generic range using SymPy's function_range on a given domain.
        """
        if var is None:
            if len(self.vars) != 1:
                raise ValueError("Range currently implemented only for single-variable functions; specify var.")
            var = self.vars[0]
        else:
            var = symbols(var)

        if domain is None:
            domain = self.domain(var)

        rng = function_range(self.expr, var, domain)
        return rng

    # For backward compatibility with your earlier naming idea
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
        if var is None:
            if len(self.vars) != 1:
                raise ValueError("Domain for PolynomialFunction: specify var for multivariable case.")
        return S.Reals

    def range(self, var=None, domain=None):
        """
        Use SymPy's function_range for polynomials over ℝ.
        """
        if var is None:
            if len(self.vars) != 1:
                raise ValueError("Range currently implemented only for single-variable polynomials; specify var.")
            var = self.vars[0]
        else:
            var = symbols(var)

        if domain is None:
            domain = S.Reals

        return function_range(self.expr, var, domain)


# -----------------------------------
# 🟢 SUBCLASS: Trigonometric Function
# -----------------------------------
class TrigFunction(Function):
    def simplify_trig(self):
        """Simplify the trigonometric expression."""
        return TrigFunction(str(self.expr.simplify()))

    def to_polynomial_if_possible(self):
        """
        Try to rewrite using algebraic forms where possible.
        """
        return TrigFunction(str(self.expr.simplify().rewrite('sqrt')))

    def domain(self, var=None):
        """
        Domain of trig functions:
        - In general, use continuous_domain.
        (This will automatically exclude points where tan, sec etc. blow up.)
        """
        return super().domain(var)

    def range(self, var=None, domain=None):
        """
        Use function_range. For simple cases:
        - sin(x), cos(x) → [-1, 1]
        - tan(x) → ℝ
        SymPy handles many of these directly.
        """
        return super().range(var, domain)


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
        continuous_domain already does this, but we keep it explicit here.
        """
        return super().domain(var)

    def range(self, var=None, domain=None):
        """
        Range via function_range over rational function domain.
        """
        return super().range(var, domain)

    def vertical_asymptotes(self):
        """Points where denominator = 0 (symbolic equation)."""
        den = self.denominator()
        return f"Solutions of: {den} = 0"

    def horizontal_asymptote(self):
        """Determine horizontal asymptote based on polynomial degrees if possible."""
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
    def domain(self, var=None):
        """
        exp-based expressions are usually defined for all ℝ,
        but we still rely on continuous_domain in case of more complex expressions
        like exp(1/x), etc.
        """
        return super().domain(var)

    def range(self, var=None, domain=None):
        """
        For pure exp(x) or a*exp(x)+b, range is (0, ∞) or shifted,
        but in general we again rely on function_range.
        """
        return super().range(var, domain)


# -----------------------------------
# 🟤 SUBCLASS: Logarithmic Function
# -----------------------------------
class LogFunction(Function):
    def domain(self, var=None):
        """
        Logarithm is defined where its argument > 0.
        continuous_domain handles this (for real log).
        """
        return super().domain(var)

    def range(self, var=None, domain=None):
        """
        For log(x) (base e), range is all ℝ.
        But with more complex arguments, we use function_range.
        """
        return super().range(var, domain)


# -----------------------------------
# 🟠 SUBCLASS: Absolute Value Function
# -----------------------------------
class AbsoluteFunction(Function):
    def domain(self, var=None):
        """
        |x| is defined for all real x (normally),
        but we delegate to continuous_domain in case it's nested inside other ops.
        """
        return super().domain(var)

    def range(self, var=None, domain=None):
        """
        For |x|, range is [0, ∞), but more complex expressions may change that.
        So we use function_range generically.
        """
        return super().range(var, domain)


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

    # For domain & range: use a single-variable example
    f = Function("x^2 + 1")
    print("\nf(x) =", f)
    print("Domain(f):", f.Domain())
    print("Range(f):", f.Range())

    # Derivative & Integral
    print("dy/dx for f:", f.derivative("x"))
    print("∫ f dx:", f.integral("x"))

    # Polynomial
    p = PolynomialFunction("x^3 - 4x + 1")
    print("\nPolynomial p(x) =", p)
    print("degree(p):", p.degree())
    print("Domain(p):", p.Domain())
    print("Range(p):", p.Range())
    print("p' =", p.derivative())

    # Trigonometric
    t = TrigFunction("sin(x) + cos(x)")
    print("\nTrig t(x) =", t)
    print("t simplified:", t.simplify_trig())
    print("t' =", t.derivative())
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
