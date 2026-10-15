from .expression import expression_references


def parse_expression(expr):
    """Return the names and dotted probe-part names an expression reads.

    Uses the same parser as the expression evaluator, so invalid or unsafe
    expressions raise ExpressionError here as well.
    """
    return expression_references(expr)
