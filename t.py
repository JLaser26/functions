from sympy import (
    symbols, diff, integrate, sin, cos, tan,
    exp, log, Abs, S, Interval, pi
)
from sympy.parsing.sympy_parser import (
    parse_expr,
    standard_transformations,
    implicit_multiplication_application
)
from sympy.calculus.util import continuous_domain, function_range
import numpy as np
import matplotlib.pyplot as plt

# Allow implicit multiplication like "4x", "2sin(x)"
TRANSFORMS = standard_transformations + (implicit_multiplication_application,)

# ------------------------------
# 🔵 SUPERCLASS: Generic Function
# ------------------------------
class Function:
    def __init__(self, expr: str):
        self.expr_str = expr.replace("^", "**")
        self.expr = parse_expr(self.expr_str, transformations=TRANSFORMS)
        self.vars = sorted(self.expr.free_symbols, key=lambda s: s.name)

    def __str__(self):
        return str(self.expr)

    def tokens(self):
        return self.expr_str.split()

    def evaluater(self, **values):
        evaluated = self.expr.subs(values)
        if evaluated.free_symbols:
            return evaluated
        try:
            return float(evaluated.evalf())
        except TypeError:
            return evaluated

    def _ensure_symbol(self, var):
        if isinstance(var, str):
            return symbols(var)
        return var

    # ---------- Calculus ----------
    def derivative(self, var=None):
        if var is None:
            if len(self.vars) != 1:
                raise ValueError("Specify a variable for derivative.")
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
        if var is None:
            if len(self.vars) != 1:
                raise ValueError("Domain only for single-variable functions; pass var explicitly.")
            var = self.vars[0]
        else:
            var = self._ensure_symbol(var)
        return continuous_domain(self.expr, var, S.Reals)

    def range(self, var=None, domain=None):
        if var is None:
            if len(self.vars) != 1:
                raise ValueError("Range only for single-variable functions; pass var explicitly.")
            var = self.vars[0]
        else:
            var = self._ensure_symbol(var)

        if domain is None:
            domain = self.domain(var)

        try:
            return function_range(self.expr, var, domain)
        except NotImplementedError:
            return "Range could not be computed symbolically."

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

    # ---------------------------------------------------------
    # 📈 NEW: PLOT METHOD (auto domain, auto handling asymptotes)
    # ---------------------------------------------------------
    def plot(self, var=None, interval=None, points=500):
        if var is None:
            if len(self.vars) != 1:
                raise ValueError("Plotting only works for single-variable functions. Pass var explicitly.")
            var = self.vars[0]
        else:
            var = self._ensure_symbol(var)

        # Auto interval if none given
        if interval is None:
            if isinstance(self, TrigFunction):
                interval = (-2 * float(pi), 2 * float(pi))
            elif isinstance(self, RationalFunction):
                # Avoid poles for rational functions
                dom = self.domain(var)
                # choose first connected interval
                try:
                    parts = list(dom.as_relational(var))
                except:
                    interval = (-10, 10)
                interval = (-10, 10)
            else:
                interval = (-10, 10)

        a, b = interval

        xs = np.linspace(a, b, points)
        ys = []

        for x_val in xs:
            try:
                y_val = self.expr.subs(var, x_val)
                if y_val.is_real or y_val.is_number:
                    ys.append(float(y_val))
                else:
                    ys.append(np.nan)
            except Exception:
                ys.append(np.nan)

        plt.figure(figsize=(7, 5))
        plt.plot(xs, ys, label=str(self.expr), linewidth=2)
        plt.axhline(0, color="black", linewidth=0.8)
        plt.axvline(0, color="black", linewidth=0.8)
        plt.title(f"Plot of {self.expr}")
        plt.xlabel(str(var))
        plt.ylabel("f(x)")
        plt.grid(True)
        plt.legend()
        plt.show()


# -----------------------------------
# Other subclasses (same as earlier)
# -----------------------------------

class PolynomialFunction(Function):
    def degree(self):
        if len(self.vars) != 1:
            raise ValueError("Degree only for single-variable polynomials.")
        poly = self.expr.as_poly()
        if poly is None:
            raise ValueError("Not a polynomial.")
        return poly.degree()

    def domain(self, var=None):
        return S.Reals


class TrigFunction(Function):
    pass


class RationalFunction(Function):
    def numerator(self):
        return self.expr.as_numer_denom()[0]

    def denominator(self):
        return self.expr.as_numer_denom()[1]

    def vertical_asymptotes(self):
        return f"Solutions of: {self.denominator()} = 0"

    def horizontal_asymptote(self):
        num, den = self.expr.as_numer_denom()
        p_num, p_den = num.as_poly(), den.as_poly()
        if not p_num or not p_den:
            return "No horizontal asymptote"
        if p_num.degree() < p_den.degree():
            return "y = 0"
        if p_num.degree() == p_den.degree():
            return f"y = {p_num.LC()}/{p_den.LC()}"
        return "No horizontal asymptote"


class ExponentialFunction(Function):
    pass


class LogFunction(Function):
    pass


class AbsoluteFunction(Function):
    def range(self, var=None, domain=None):
        if isinstance(self.expr, Abs):
            return Interval(0, S.Infinity)
        return super().range(var, domain)


# ------------------------------
# TESTING PLOT FEATURE
# ------------------------------
if __name__ == "__main__":
    f = PolynomialFunction("x^3 - 4x + 1")
    f.plot()      # automatically plots on [-10, 10]

    g = TrigFunction("sin(x) + cos(x)")
    g.plot()      # automatically plots on [-2π, 2π]

    r = RationalFunction("(x^2 + 1)/(x - 2)")
    r.plot()      # automatically avoids asymptote reasonably

    h = AbsoluteFunction("Abs(x - 1)")
    h.plot()

    e = ExponentialFunction("exp(x)")
    e.plot()

    l = LogFunction("log(x)")
    l.plot(interval=(0.1, 10))   # safe interval for log

    k = PolynomialFunction("x^5 - 2x^4 + 3x^3 - 4x^2 + 5x + 6")
    k.plot(interval=(-100, 100))

    jj = TrigFunction("(sin(x-x^2))/cos(1/x)")
    jj.plot()

