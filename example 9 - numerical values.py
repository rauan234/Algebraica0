from Algebraica import cp, function, numerical
import Algebraica
# with Algebraica, you can calculate numerical values of certain expressions.
# for that, use C.numerical() method, that returns the numerical value of C.
# the method returns the numerical value in class Numerical format.
# to learn about class Numerical, visit example 2.


print('sqrt(2) =', cp( (1, ('sqrt', cp(2))) ).numerical())
print('(2 + 6 i) ** (3 + 1 i) =', cp( (1, ('pow', cp(2, 6), cp(3, 1))) ).numerical())
print()
# this is how to calculate power functions of type pow(u, v).
# if u is real and positive and v is real, pow(u, v) is calculated just as numpy.pow(u.rl, v.rl).
# otherwise, the answer is calculated by representing u as a complex exponential.

print('ln(3) =', cp( (1, ('ln', cp(3))) ).numerical())
print('ln(3 + 4 i) =', cp( (1, ('ln', cp(3, 4))) ).numerical())
print()
# natural logs ln(x) are calculated as numpy.log(x), if x is real and positive.
# otherwise, ln(x) = numpy.log(|x|) + i x.arg().

print('Rl(5 + 3 i) =', cp( (1, ('Rl', cp(5, 3))) ).numerical())
print('Im(5 + 3 i) =', cp( (1, ('Im', cp(5, 3))) ).numerical())
print('abs(5 + 3 i) =', cp( (1, ('abs', cp(5, 3))) ).numerical())
print('arg(5 + 3 i) =', cp( (1, ('arg', cp(5, 3))) ).numerical())
print('sign(5 + 3 i) =', cp( (1, ('sign', cp(5, 3))) ).numerical())
print()
# here is how to get the real and the imaginary parts, and also arg and abs ans sign

print('sign(1) =', cp( (1, ('sign', cp(1))) ).numerical())
print('sign(-4 + 2 i) =', cp( (1, ('sign', cp(-4, 2))) ).numerical())
print()
# in Algebraica, sign(x) is defined as sign(x) = x / abs(x)

print('sin(2) =', cp( (1, ('sin', cp(2))) ).numerical())
print('cos(2) =', cp( (1, ('cos', cp(2))) ).numerical())
print('sec(2) =', cp( (1, ('sec', cp(2))) ).numerical())
print('csc(2) =', cp( (1, ('csc', cp(2))) ).numerical())
print('tan(2) =', cp( (1, ('tan', cp(2))) ).numerical())
print('cot(2) =', cp( (1, ('cot', cp(2))) ).numerical())
print()
# trigonometric functions of real arguments sin(x) are calculated as numpy.sin(x)

print('sin(2 + 1 i) =', cp( (1, ('sin', cp(2, 1))) ).numerical())
print('cos(2 + 1 i) =', cp( (1, ('cos', cp(2, 1))) ).numerical())
print('sec(2 + 1 i) =', cp( (1, ('sec', cp(2, 1))) ).numerical())
print('csc(2 + 1 i) =', cp( (1, ('csc', cp(2, 1))) ).numerical())
print('tan(2 + 1 i) =', cp( (1, ('tan', cp(2, 1))) ).numerical())
print('cot(2 + 1 i) =', cp( (1, ('cot', cp(2, 1))) ).numerical())
print()
# if x is complex, sin(x) = numpy.sin(Rl(x)) * numpy.cosh(Im(x)) + i numpy.cos(Rl(x)) * numpy.sinh(Im(x))


print('sinh(2) =', cp( (1, ('sinh', cp(2))) ).numerical())
print('cosh(2) =', cp( (1, ('cosh', cp(2))) ).numerical())
print('sech(2) =', cp( (1, ('sech', cp(2))) ).numerical())
print('csch(2) =', cp( (1, ('csch', cp(2))) ).numerical())
print('tanh(2) =', cp( (1, ('tanh', cp(2))) ).numerical())
print('coth(2) =', cp( (1, ('coth', cp(2))) ).numerical())
print()

print('sinh(2 + 1 i) =', cp( (1, ('sinh', cp(2, 1))) ).numerical())
print('cosh(2 + 1 i) =', cp( (1, ('cosh', cp(2, 1))) ).numerical())
print('sech(2 + 1 i) =', cp( (1, ('sech', cp(2, 1))) ).numerical())
print('csch(2 + 1 i) =', cp( (1, ('csch', cp(2, 1))) ).numerical())
print('tanh(2 + 1 i) =', cp( (1, ('tanh', cp(2, 1))) ).numerical())
print('coth(2 + 1 i) =', cp( (1, ('coth', cp(2, 1))) ).numerical())
print()
# with hyperbolic trig functions, the process is similar


print('arcsin(1/2) =', cp( (1, ('arcsin', cp((1, 2)))) ).numerical())
print('arccos(1/2) =', cp( (1, ('arccos', cp((1, 2)))) ).numerical())
print('arcsec(2) =', cp( (1, ('arcsec', cp((2, 1)))) ).numerical())
print('arccsc(2) =', cp( (1, ('arccsc', cp((2, 1)))) ).numerical())
print('arctan(1/2) =', cp( (1, ('arctan', cp((1, 2)))) ).numerical())
print('arccot(1/2) =', cp( (1, ('arccot', cp((1, 2)))) ).numerical())
print()

print('arcsin(2 + 1 i) =', cp( (1, ('arcsin', cp(2, 1))) ).numerical())
print('arccos(2 + 1 i) =', cp( (1, ('arccos', cp(2, 1))) ).numerical())
print('arcsec(2 + 1 i) =', cp( (1, ('arcsec', cp(2, 1))) ).numerical())
print('arccsc(2 + 1 i) =', cp( (1, ('arccsc', cp(2, 1))) ).numerical())
print('arctan(2 + 1 i) =', cp( (1, ('arctan', cp(2, 1))) ).numerical())
print('arccot(2 + 1 i) =', cp( (1, ('arccot', cp(2, 1))) ).numerical())
print()
# and with inverse trig it`s the same


print('arcsinh(1/2) =', cp( (1, ('arcsinh', cp((1, 2)))) ).numerical())
print('arccosh(1/2) =', cp( (1, ('arccosh', cp((1, 2)))) ).numerical())
print('arcsech(1/2) =', cp( (1, ('arcsech', cp((1, 2)))) ).numerical())
print('arccsch(1/2) =', cp( (1, ('arccsch', cp((1, 2)))) ).numerical())
print('arctanh(1/2) =', cp( (1, ('arctanh', cp((1, 2)))) ).numerical())
print('arccoth(1/2) =', cp( (1, ('arccoth', cp((1, 2)))) ).numerical())
print()

print('arcsinh(2 + 1 i) =', cp( (1, ('arcsinh', cp(2, 1))) ).numerical())
print('arccosh(2 + 1 i) =', cp( (1, ('arccosh', cp(2, 1))) ).numerical())
print('arcsech(2 + 1 i) =', cp( (1, ('arcsech', cp(2, 1))) ).numerical())
print('arccsch(2 + 1 i) =', cp( (1, ('arccsch', cp(2, 1))) ).numerical())
print('arctanh(2 + 1 i) =', cp( (1, ('arctanh', cp(2, 1))) ).numerical())
print('arccoth(2 + 1 i) =', cp( (1, ('arccoth', cp(2, 1))) ).numerical())
print()
# as well as with inverse hyperbolic trig


X = cp( (1, ('sin', cp('x'))) )
# let`s say you want to compute the numerical value of sin(x)

try:
    print(X.numerical())
except Exception as exc:
    print(exc.get_function())
    print(exc)
print()
# but since the numerical value for x is not known, Algebraica will throw an exception.
# if you catch it and print, you can see what exactly went wrong by printing the exception.
# also, you can access the problematic function directly with exc.get_function().

# you can rectify this problem by calling X.numerical() with additional arguments.
print(X.numerical( function('x'), numerical(1.57) ))
print()
# in this case, X.numerical is called with two arguments, the first
# argument representing the function whose numerical value is known,
# and the second argument representing its` numerical value.
# function(*args) and numerical(rl, im=0) are just shortcuts for __init__.

A = cp( (1, ('pow', cp('x'), cp('y'))) )
print('A =', A)
print(A.numerical([
    (function('x'), numerical(2)),
    (function('y'), numerical(5, 3))
]))
print()
# if you have multiple variables to substitute for, you can call numerical with a list
# of tuples, where each tuple containt a function to be substituted for and the value
# that should be substituted.

B = cp( (1, ('f', cp(4))) )
print('B =', B)
print(B.numerical( function('f', cp(4)), numerical(1) ))
print()
# you can not only substitute for variables, but also for functions

# and also derivatives of functions
C = cp(function('g', cp('z'), partials=[2]))
print('C =', C)
print(C.numerical( function('g', cp('z'), partials=[2]), numerical(7) ))
print()

try:
    C.numerical(1)
except Exception as exc:
    print(exc)
print()
# the arguments for Complex.numerical must always be in either one of
# the two acceptable forms given above.
# otherwise, an exception is thrown.


input()
