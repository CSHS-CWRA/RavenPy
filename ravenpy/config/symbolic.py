from typing import Union

import pymbolic
from pymbolic.mapper.coefficient import CoefficientCollector
from pymbolic.mapper.evaluator import EvaluationMapper as EM
from pymbolic.mapper.stringifier import StringifyMapper
from pymbolic.primitives import Expression, Variable

# Type hint for symbolic expressions
RavenExp = Union[Variable, Expression, float, None]

# Global registry to store Ostrich TiedParams expressions:
#   * If None, it means that we are running in non-Ostrich mode, where
#     an unknown symbolic value should raise an error
#   * When running Ostrich, this registry is first set to an empty dict,
#     which signals that the unknown symbolic values should be kept in it
TIED_PARAMS_REGISTRY = None


class TiedParamsMapper(pymbolic.mapper.stringifier.StringifyMapper):
    """Return a string representation of an expression.

    This is used to identify expressions for Ostrich.

    Notes
    -----
    No guarantee that this returns unique identifiers.
    """

    def format(self, s, *args):
        return (s % args).replace(".", "d")

    def map_sum(self, expr, enclosing_prec, *args, **kwargs):
        return self.join_rec("_add_", expr.children, 1, *args, **kwargs)

    def map_product(self, expr, enclosing_prec, *args, **kwargs):
        return self.join_rec("_mul_", expr.children, 1, *args, **kwargs)

    def map_quotient(self, expr, enclosing_prec, *args, **kwargs):
        return self.join_rec(
            "%s_over_%s", [expr.numerator, expr.denominator], 1, *args, **kwargs
        )


def parse_symbolic(value, **kwds):
    if isinstance(value, dict):
        return {k: parse_symbolic(v, **kwds) for k, v in value.items()}

    elif isinstance(value, (list, tuple)):
        return [parse_symbolic(v, **kwds) for v in value]

    elif isinstance(value, (Variable, Expression)):
        try:
            # Convert to numerical value
            return EM(context=kwds)(value)

        except pymbolic.mapper.evaluator.UnknownVariableError:
            # Convert to string, an identifier for expressions, and the variable name for variables.
            if isinstance(value, Expression) and not isinstance(value, Variable):

                # If this registry is not defined, it means that we are
                # running in a non-Ostrich context, where an unknown symbolic
                # value should be an error
                if TIED_PARAMS_REGISTRY is None:
                    raise ValueError(f"{value} is not defined")
                else:
                    # We are running Ostrich, with which we want to collect
                    # tied params in the registry
                    key = "par_" + TiedParamsMapper()(value).lower()
                    TIED_PARAMS_REGISTRY[key] = value
                    return key

            # Convert to expression string
            return StringifyMapper()(value)
        return value

    else:
        return value
