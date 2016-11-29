import re
from enum import Enum
from sys import argv

### PARSING

def parse_expr(str):
    """
    Parse a string consisting of +-separated terms.

    >>> parse_expr("c[k1].h[p1mu] + a[k1].ad2f[p2]")
    [('+', [('C', 'c', 'k1'), ('A', 'h', 'p1mu')]), ('+', [('A', 'a', 'k1'), ('A', 'ad2f', 'p2')])]
    """
    return [parse_term(t.strip()) for t in str.split('+')]

def parse_term(str):
    """
    Parse a string consisting of dot-separated factors. The result may be
    negated or not.

    >>> parse_term("c[k1]")
    ('+', [('C', 'c', 'k1')])

    >>> parse_term("-c[k1]")
    ('-', [('C', 'c', 'k1')])

    >>> parse_term("c[k1].h[p1mu]")
    ('+', [('C', 'c', 'k1'), ('A', 'h', 'p1mu')])

    """
    if str[0] == '-':
        t = parse_term(str[1:])
        factors = t[1]
        return ('-', factors)
        return t
    else:
        return ('+', [parse_factor(part) for part in str.split('.')])

def parse_factor(input):
    """
    Given a string that looks like "c[k1]" or "ad2f[p1mu]", return
    a triple with the
        * operation type ('C' for constructor, 'A' for annihilator),
        * operation ("c" and "ad2f" in the examples), and
        * operand ("k1" and "p1mu" in the examples).

    >>> parse_factor("c[k1]")
    ('C', 'c', 'k1')

    >>> parse_factor('ad2f[p1mu]')
    ('A', 'ad2f', 'p1mu')

    >>> parse_factor('h[k1]')
    ('A', 'h', 'k1')

    """

    match = re.search('(\w+)\[(\w+)\]', input)
    operation = match.group(1)
    operand = match.group(2)
    op_type = classify(operation)

    return (op_type, operation, operand)

def classify(op):
    """
    Given an operation name, decide whether it is a constructor or
    an annihilator.  The convention is that constructors are operations
    starting with 'c', and all other operations are annihilators.

    >>> classify("c")
    'C'

    >>> classify("c2df")
    'C'

    >>> classify("a")
    'A'

    >>> classify("h")
    'A'

    """

    if op[0] is 'c':
        return 'C'
    else:
        return 'A'


### UNPARSING

def show_expr(expr):
    """
    Convert an expression back into a string.

    >>> show_expr([])
    ''

    >>> show_expr(parse_expr('c1ds[p1m].a1ds[k1]+cds[p1m].ads[k2]'))
    'c1ds[p1m].a1ds[k1] + cds[p1m].ads[k2]'
    """
    return " + ".join([ show_term(term) for term in expr ])

def print_expr(expr):
    print(" +\n".join([ show_term(term) for term in expr ]))


def show_term(term):
    """
    Convert a term back into a string.

    >>> show_term(parse_term('c1ds[p1m].a1ds[k1]'))
    'c1ds[p1m].a1ds[k1]'

    >>> show_term(parse_term('-h[b].h[c]'))
    '-h[b].h[c]'

    """
    if term[0] == '-':
        positive_term = ('+', term[1])
        str = show_term(positive_term)
        return '-' + str
    else:
        return ".".join([ show_factor(factor) for factor in term[1]])

def show_factor(factor):
    """
    Convert a factor back into a string.

    >>> show_factor(parse_factor('c1ds[p1m]'))
    'c1ds[p1m]'

    """
    return "{1}[{2}]".format(*factor)


### CALCULATIONS

def related_factors(f1, f2):
    """
    Given two parsed factors, decide whether they are the creator and
    annihilator (in either order) of the same type. The operators do
    not have to be acting on the same operand.

    >>> related_factors(parse_factor("cd2f[p2]"), parse_factor("ad2f[p1]"))
    True

    >>> related_factors(parse_factor("ad2f[p1]"), parse_factor("cd2f[p1]"))
    True

    >>> related_factors(parse_factor("ad2f[p1]"), parse_factor("ad2f[p1]"))
    False

    >>> related_factors(parse_factor("h[p1]"), parse_factor("a[p1]"))
    False

    """

    if op_type(f1) == 'C' and op_type(f2) == 'A':
        f1, f2 = f2, f1

    if op_type(f1) == 'A' and op_type(f2) == 'C':
        op1 = op(f1)
        op2 = op(f2)
        return op1[1:] == op2[1:]
    else:
        return False

def is_spinor_factor(factor):
    """
    A parsed factor represents a spinor if its operation name ends in f.

    >>> is_spinor_factor(parse_factor("cd2f[p2]"))
    True

    >>> is_spinor_factor(parse_factor("ad1s[k1]"))
    False
    """
    return op(factor)[-1] == 'f'


### NORMAL ORDERING ALGORITHM

def normal_order(expr):
    """
    Return the normal ordering on a parsed expression expr.

    >>> print_expr(normal_order(parse_expr('ads[k2].cds[p1m]')))
    DiracDelta[-k2+p1m] +
    cds[p1m].ads[k2]

    >>> print_expr(normal_order(parse_expr('a1ds[k1].c2ds[p1m]')))
    c2ds[p1m].a1ds[k1]

    """

    while True:
        next_expr = step_expr(expr)
        if next_expr == None:
            return expr
        expr = next_expr

def step_expr(expr):
    """
    Return the next step in the normal ordering for the given expression. If
    no next step exists, return None. If there are multiple options, pick the
    leftmost option.

    >>> show_expr(step_expr(parse_expr("a1ds[k1].c2ds[p1m]")))
    'c2ds[p1m].a1ds[k1]'

    >>> show_expr(step_expr(parse_expr("a1ds[k1].c1ds[p1m]")))
    'DiracDelta[-k1+p1m] + c1ds[p1m].a1ds[k1]'

    >>> step_expr(parse_expr('h[a].h[b]'))
    """
    transform = find_transformation(expr)
    if transform:
        term_idx, factor_idx, swap_type = transform
        pre, term, post = expr[0:term_idx], expr[term_idx], expr[term_idx+1:]
        new_terms = swap_factors(term, factor_idx, swap_type)
        return pre + new_terms + post


def swap_factors(term, idx, swap_type):
    if swap_type == 'commutative':
        return commutative_swap(idx, term)
    else:
        return delta_swap(idx, term)


def find_transformation(expr):
    """
    Return the type and position of the next normal ordering step
    for the given expression. If no next step exists, return None.
    If there are multiple options, pick the leftmost option. The
    tuple returned, if any, indicates respectively the index of the
    term in the expression, the index of the factor in the term, and
    the type of swap ('delta' or 'commutative').

    >>> find_transformation(parse_expr('h[a].h[b] + a1ds[k1].c2ds[p1m]'))
    (1, 0, 'commutative')

    >>> find_transformation(parse_expr('a1ds[k1].c1ds[p1m] + h[a].h[b]'))
    (0, 0, 'delta')

    >>> find_transformation(parse_expr('h[a].h[b]'))

    """
    for term_idx, term in enumerate(expr):
        pos = find_swap(term)
        if pos >= 0:
            factors = term[1]
            f1, f2 = factors[pos], factors[pos+1]
            swap_type = 'commutative' if commutative(f1, f2) else 'delta'
            return (term_idx, pos, swap_type)


def commutative(f1, f2):
    """
    Given two factors representing an annihilator and a creator,
    determine whether they are commutative or not.  The rule is
    that they commute iff their names differ after the first
    letter.

    >>> commutative(parse_factor('a1ds[p1m]'), parse_factor('c1ds[k]'))
    False

    >>> commutative(parse_factor('a2ds[p1m]'), parse_factor('c1ds[k]'))
    True
    """

    return op(f1)[1:] != op(f2)[1:]



def find_swap(term):
    """
    If a term contains a creator immediately to the right of an annihilator,
    then the creator can be moved left with either a Dirac Delta or a
    Commutative Swap.

    If such a pair exists in the term, find_transformation() returns the
    index of the first factor in the term. If there are multiple options,
    it returns the leftmost one.  If no such pair exists, the function returns
    -1.

    >>> find_swap(parse_term("a1ds[k1].c2ds[p1m]"))
    0

    >>> find_swap(parse_term("a1ds[k1].c2ds[p1m].a1ds[k1].c2ds[p1m]"))
    0

    >>> find_swap(parse_term("a1df[k2].a1ds[k1].c2ds[p1m]"))
    1

    >>> find_swap(parse_term("c1ds[k1].a2ds[p1m]"))
    -1

    """

    factors = term[1]
    for i in range(len(factors) - 1):
        if op_type(factors[i]) == 'A' and op_type(factors[i+1]) == 'C':
            return i

    return -1

def commutative_swap(pos, term):
    """
    Swap two factors in a term when they commute. Returns a list
    containing the new term, ready to be spliced.

    >>> show_expr(commutative_swap(1, parse_term("h[p].a[q].c[r].h[s]")))
    'h[p].c[r].a[q].h[s]'

    """

    sign, factors = term
    pre, f1, f2, post = factors[0:pos], factors[pos], factors[pos+1], factors[pos+2:]

    return [ (sign, pre + [f2, f1] + post) ]

def delta_swap(pos, term):
    """
    Swap two related factors in a term with a Dirac Delta. Returns the
    terms the given term is to be replaced by.

    >>> show_expr(delta_swap(0, parse_term("a[k].c[m]")))
    'DiracDelta[-k+m] + c[m].a[k]'

    >>> show_expr(delta_swap(1, parse_term("a[k2].a[k1].c[p1].c[p2]")))
    'DiracDelta[-k1+p1].a[k2].c[p2] + a[k2].c[p1].a[k1].c[p2]'

    >>> show_expr(delta_swap(1, parse_term("h[a].h[b].h[c].h[d]")))
    'DiracDelta[-b+c].h[a].h[d] + h[a].h[c].h[b].h[d]'

    """

    sign, factors = term
    pre, f1, f2, post = factors[0:pos], factors[pos], factors[pos+1], factors[pos+2:]

    sign1 = '+'
    sign2 = '-' if is_spinor_factor(f2) else '+'

    delta = ('', 'DiracDelta', "-{0}+{1}".format(arg(f1), arg(f2)))
    term1 = (sign1, [delta] + pre + post)
    term2 = (sign2, pre + [f2, f1] + post)

    return [term1, term2]

def consecutive_pairs(list):
    """
    Utility function to return each consecutive pair from a list

    >>> consecutive_pairs([1,2,3,4,5])
    [(1, 2), (2, 3), (3, 4), (4, 5)]

    >>> consecutive_pairs(['a', 'b', 'c'])
    [('a', 'b'), ('b', 'c')]

    >>> consecutive_pairs([2])
    []

    """

    return [(list[i], list[i+1]) for i in range(len(list) - 1)]


### ACCESSORS

def op_type(factor):
    return factor[0]

def op(factor):
    return factor[1]

def arg(factor):
    return factor[2]


if __name__ == "__main__":
    if len(argv) != 2:
        print("Usage: {0} expr".format(argv[0]))
        print("Example:")
        print("{0} ad1s[k1].ad1s[k2].ah1s[p1mu].cd1s[p1m].ad1s[p2m].cd1s[q1].cd1s[q2].ch1s[q3]".format(argv[0]))
    elif argv[1] == '--selftest':
        import doctest
        doctest.testmod()
    else:
        expr = parse_expr(argv[1])
        result = normal_order(expr)
        print_expr(result)
