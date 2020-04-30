from Algebraica import Rational, Numerical, Function, Term, Irrational, Complex, Numerical
from Algebraica import cp, function, numerical
import Algebraica
# this example explains the general structure of Algebraica.
# class Complex is defined using a few slave classes, namely:
# Irrational, Term, Function, Rational


# the most basic subsidiary class in Algebraica is Rational
# it represents a fraction with integer numerator and an
# integer denominator

A = Rational(1, 3)
B = Rational(7, 4)
print('A =', A)
print('B =', B)
print()
# you can create a class Rational object by writing an expression of form
# Rational( N, D ), where N and D are the numerator and the denominator

print('A + B =', A + B)
print('A - B =', A - B)
print('A * B =', A * B)
print('A / B =', A / B)
print()
# here are the basic arithmetic operations between two Rationals

print('A + 3 =', A + 3)
print('A + 4.2 =', A + 4.2)
print('A - 3 =', A - 3)
print('A - 4.2 =', A - 4.2)
print()

print('A * 3 =', A * 3)
print('A * 4.2 =', A * 4.2)
print('A / 3 =', A / 3)
print('A / 4.2 =', A / 4.2)
print()
# you can also do arithmetics with Rationals and integers or floats

C = Rational(7)
print('C =', C)
print()
# of you want, you can create Rational class object by putting
# in just one number instead of two.
# in this case, the number you entered will be the numerator, and the
# denominator will be set equal to one.

D = Rational(6.4)
print('D =', D)
print()
# you can declare a Rational using a float number.
# then, the float will be first converted to a decimal fraction and then simplified.
# for example, 6.4 -> 64/10 -> 32/5

print('A =', A)
A *= B
print('A =', A)
A /= B
print('A =', A)
A += B
print('A =', A)
A -= B
print('A =', A)
print()
# here are some assignment operators

C += 1.2
D += 4
print('C =', C)
print('D =', D)
print()

C -= 3
D -= 7.4
print('C =', C)
print('D =', D)
print()

C *= 5
D *= 9.7
print('C =', C)
print('D =', D)
print()

C /= 2.6
D /= 9
print('C =', C)
print('D =', D)
print()
# there are also assignment operators for class Rational that
# take floats or integers as an argument.

print('(A == B) =', A == B)
print('(A != B) =', A != B)
print('(A >= B) =', A >= B)
print('(A <= B) =', A <= B)
print('(A > B) =', A > B)
print('(A < B) =', A < B)
print()
# class Rational objects can be compared with each other

print('(A == 1) =', A == 1)
print('(A != 2.5) =', A != 2.5)
print('(A >= 3.2) =', A >= 3.2)
print('(A <= 4.9) =', A <= 4.9)
print('(A > 7) =', A > 7)
print('(A < 9.0) =', A < 9.0)
print()
# class Rational objects can also be compared with integers and floats

print('-C =', -C)
print('-D =', -D)
print()
# by putting a minus sign in front of the Rational, you can get its` negative

print('~C =', ~C)
print('~D =', ~D)
print()
# by putting a tilda, you can reverse the fraction, e.g. switch the numerator
# and the denominator, so that, for example, 3/4 becomes 4/3

print(Rational(12, 245).prime_factors())
print()
# this is how Rationals can be decomposed into primes

print('2/4 =', Rational(2, 4))
print()
# note that the fractions are automatically simplified

print('1/(-2) =', Rational(1, -2))
print()
# also, the minus sign is automatically transfered from
# the numerator to the denominator

del A
del B
del C
del D


# the next subsidiary class is called Function.
# it carries information about one specific function.
# more precisely, it tells us the name of the function and what
# arguments were given to it, such as sqrt[ 2 ].

A = Function( ( 'sqrt', cp(2) ) )
print('A =', A)
print()
# in order to initialize a class Function object, we should call __init__
# with argument being a tuple, whose 0th entry tells the name of the
# function and whose other entries are the arguments of the function.

B = Function( ( 'pi' ) )
C = Function(   'pi'   )
print('B =', B)
print('C =', C)
print()
# sometimes you want you function to just be constant, so that there is no arguments.
# you can do it by either giving a tuple with only one element or by entering the name of the constant.
# there is no difference, because __init__ treats both equally.

# instead of Function.__init__ you might want to use function, because it`s
# a bit more convenient.
# the process of using is almost the same as with __init__, except that the
# arguments are not inserted into a tuple.

print('(A == B) =', A == B)
print('(B == C) =', B == C)
print('(A == C) =', A == C)
print('(sqrt(2) == sqrt(3)) =', function( 'sqrt', cp(2) ) == 
      function( 'sqrt', cp(3) ) )
print('exp(1) == exp(1) =', function( 'exp', cp(1) ) == 
      function( 'exp', cp(1) ) )
print()
# this is how Functions are compared.
# in order for two Functions to be equal, they must have equal names,
# equal number of arguments, and the arguments at each position must be the same.

# class Function is extremely important for Algebraica. It contains
# the most interesting part of the program, as you will see in further examples.

del A
del B
del C


# by combining Rational and Function we get Term.
# a Term is basically a Rational coefficient multiplied by
# a series of Fuctions.
# for example, 4/7 * sqrt(3) * ln(2) is a Term

R = Rational( 7, 2 )
functions = [
    function('sqrt', cp(3)),
    function('ln', cp(10)),
    function( 'pow', cp(2), cp( (2, 3) ) ),
    Function( 'pi' )
]
A = Term( R, functions )
del R
del functions
print('A =', A)
print()
# in order to declare a Term, write an expression of a form
# Term( R, functions ), where R is the rational coefficient and
# functions is a list of all functions that are multiplied together.
# the end result is R * functions[0] * functions[1]...

R = Rational( 3, 8 )
functions = [
    function( 'sqrt', cp('pi') )
]
B = Term( R, functions )
print('B =', B)
print()
# here is one more Term to play with

# ! print(A + B) !
# ! print(A - B) !
# ! print(A / B) !
# Term objects cannot be added, subtracted or divided

print('A * B =', A * B)
print()
# they can only be multiplied.
# during multiplication, the Rational coefficients of the Term objects
# are multiplied and their function lists are added together.
# the reason why Terms can only be multiplied and not anything else is because that`s the only
# way they are utilized.
# in other words, we simply don`t need anything except multiplication.

del A
del B


# the reason why Terms only need to be multiplied is because they are just
# slaves of class Irrational, which in turn takes care of addition, subtraction,
# multiplication, and lots of other things.
# a class Irrational object is simply a sum of many class Term objects.
# and yes, name Irrational might be rather misleading in this case.

t1a = Term(
    Rational(1, 2),
    [ Function('pi') ]
)
t2a = Term(
    Rational(1),
    []
)
t3a = Term(
    Rational(3),
    [ function('ln', cp(3)) ]
)
A = Irrational( t1a, t2a, t3a )
print('A =', A)
print()
# in order to declare a class Irrational object, enter a few
# class Term objects in __init__ method and that`s it.
# the result will be t1a + t2a + t3a.
# of course, you can have not just three terms, but as many of them as you want.

B = Irrational(
    Term(
        Rational(2),
        [
            function('sqrt', cp(2)),
            function('exp', cp(0, 1))
        ]
    )
)
print('B =', B)
print()
# let`s create one more Irrational so that we could perform operations with A and B

print('A + B =', A + B)
print('A - B =', A - B)
print('A * B =', A * B)
print()
# here is some basic arithmetics.
# when adding two Irrationals, the result is simply a class Irrational
# object whose terms are the sum of the terms of the first two Irrationals.
# when subtracting, the procedure is the same, except that all terms of B
# undergo sign change before being added to the terms of A.
# multiplication of two class Irrational object happens just by
# expanding the bracket, e.g. (t1a + t2a + t3a) * (t1b + t2b) =
# = t1a * t1b + t1a * t2b + t2a * t1b + t2a * t2b + t3a * t1b + t3a * t2b.
# class Irrational objects cannot be divided.
# Complex will take care of it later.

print('(A == B) =', A == B)
print('(A != B) =', A != B)
print()
# two Irrationals can be compared to each other

print('-A =', -A)
print()
# that`s how you can change the sign of a class Irrational object

print('A =', A)
A *= B
print('A =', A)
A += B
print('A =', A)
A -= B
print('A =', A)
print()
# here are some assignment operators

print('A + 2 =', A + 2)
print('A + 1.6 =', A + 1.6)
print('A - 2 =', A - 2)
print('A - 1.6 =', A - 1.6)
print()
# you can add or subtract ints and floats from Irrationals

print('A * 2 =', A * 2)
print('A * 1.6 =', A * 1.6)
print('A / 2 =', A / 2)
print('A / 1.6 =', A / 1.6)
print()
# Irrationals can be multiplied and divided by int and float numbers

print('A + pi =', A + 'pi')
print('A - pi =', A - 'pi')
print('A * pi =', A * 'pi')
print()
# you can also do addition, subtraction and multiplication
# between an Irrational and a string

del A
del B


# now let`s look at the main class, for which everything is built for - Complex
rt1c = Rational( 7, 2 )
rt2c = Rational( 4, 3 )
ti1c = Rational( 9, 13 )
tr2f1 = function('pi',)
ti1f1 = function('sqrt', cp(2))
ti1f2 = function('ln', cp(3))
tr1 = Term( rt1c, [] )
tr2 = Term( rt2c, [tr2f1] )
ti1 = Term( ti1c, [ti1f1, ti1f2 ] )
real = Irrational( tr1, tr2 )
imag = Irrational (ti1 )
A = Complex( real, imag)
print('A =', A)
print()
# in order to create a Complex, call __init__ with two Irrational arguments,
# where the first argument is the real part of the number and the second argument
# is the imaginary part.
# as you can see, in order to construct one single Complex we need all that arsenal
# we`ve prepared.

B = Complex(
    Irrational(
        Term(
            Rational( 1, 2 ),
            [
                function('a', cp(3)),
                function('b', cp(4))
            ]
            ),
        Term(
            Rational( 5, 6 ),
            [
                function('c', cp(7)),
                function('d', cp(8))
            ]
            )
        ),
    Irrational(
        Term(
            Rational( 9, 10 ),
            [
                function('f', cp(11)),
                function('g', cp(12))
            ]
            ),
        Term(
            Rational( 13, 14 ),
            [
                function('h', cp(15)),
                function('j', cp(16))
            ]
            ),
        Term(
            Rational( 17, 18 ),
            [
                function('k', cp(19)),
                function('l', cp(20))
            ]
            )
        )
)
print('B =', B)
print()
# here is how you would create a Complex by using __init__.
# but if you don`t need to build that complicated structures, you might use a little
# shortcut prepared to make life easier, namely - the cp(*args) function.

C = cp(3)
print('C =', C)
print()
# you can look up how to use cp in example 0

D = Rational(1, 2)
print('D =', cp(D))
print()

E = Function(('sqrt', cp(3)))
print('E =', cp(E))
print()

F = Term(D, [E])
print('F =', cp(F))
print()

G = Irrational(F)
print('G =', cp(G))
print()
# you can create Complex by calling the cp(*args) function
# with argument being an object of one of the subsidiary classes

H = cp((8, 3), (1, ('sqrt', cp(2))))
J = Rational(4, 7)
print('H * J =', H * J)
print('H / J =', H / J)
print()
# class Complex objects can be multiplied and divided by Rationals.
# in order to see more examples of operators and other things, visit example 1.


# in addition to class Complex, Algebraica also feaures one more independent class called Numerical.
# class Numerical might be useful, if instead of formulas you want numerical values of certian expressions.

print(numerical(2))
print()
# numerical(r) returns a class Numerical object equal to r.
# here r must be either int or float.

A = numerical(1.4, 6.2)
B = Numerical(-0.3, 2)
print('A =', A)
print('B =', B)
print()
# to create Numericals with a non-zero imaginary part, you can
# use both Numerical.__init__ or numerical - there is no difference.

print('A + B =', A + B)
print('A - B =', A - B)
print('A * B =', A * B)
print('A / B =', A / B)
print('A ^ B =', A ** B)
print()
# here are the binary operations between two Numericals

print('A + 3 =', A + 3)
print('A - 3 =', A - 3)
print('A * 3 =', A * 3)
print('A / 3 =', A / 3)
print('A ^ 3 =', A ** 3)
print()
# and between Numericals and ints

print('A + 4.2 =', A + 4.2)
print('A - 4.2 =', A - 4.2)
print('A * 4.2 =', A * 4.2)
print('A / 4.2 =', A / 4.2)
print('A ^ 4.2 =', A ** 4.2)
print()
# and floats

print('-A =', -A)
print('~A =', ~A)
print()
# -(a + b i) = -a - b i
# ~(a + b i) = a - b i


print(A.wolfram_lang())
print(A.compileable())
print()
# class Numerical objects can be printed in a form suitble for Wolfram Mathematica,
# or in a form Python 3 can understand


# Algebraica also includes functions that take Numerical as an argument
# and return another Numerical
print('(6 + 7.2 i) ** (2 + 3 i) =', Algebraica.cpower(Numerical(6, 7.2), Numerical(2, 3)))
print('sqrt(6 + 7.2 i) =', Algebraica.csqrt(Numerical(6, 7.2)))
print()
# this is the power function.
# the "c" in "cpower" stands for complex.
# it`s added so that the complex power function wouldn`t be confused
# with the real-domain numpy.power.

print('Rl(6 + 7.2 i) =', Numerical(6, 7.2).rl)
print('Im(6 + 7.2 i) =', Numerical(6, 7.2).im)
print('abs(6 + 7.2 i) =', abs(Numerical(6, 7.2)))
print('arg(6 + 7.2 i) =', Numerical(6, 7.2).arg())
print()
# here is how to find real and imaginary part, and also
# the absoute value and the complex argument 

print('sign(6 + 7.2 i) =', Numerical(6, 7.2).get_sign())
print('sign(-6 + 7.2 i) =', Numerical(-6, 7.2).get_sign())
print('sign(6 - 7.2 i) =', Numerical(6, -7.2).get_sign())
print('sign(-6 - 7.2 i) =', Numerical(-6, -7.2).get_sign())
print()
# sign(z) is defined to be z / abs(z)

print('exp(6 + 7.2 i) =', Algebraica.cexp(Numerical(6, 7.2)))
print('ln(6 + 7.2 i) =', Algebraica.cln(Numerical(6, 7.2)))
print('log(6 + 7.2 i, 2 + i) =', Algebraica.clog(Numerical(6, 7.2), Numerical(2, 1)))
print()
# here are the exponential function and the natural log, and also log with arbitrary base.
# log of a with base b is written as log(a, b) and equals Ln(a) / Ln(b).

print('sin(6 + 7.2 i) =', Algebraica.csin(Numerical(6, 7.2)))
print('cos(6 + 7.2 i) =', Algebraica.ccos(Numerical(6, 7.2)))
print('sec(6 + 7.2 i) =', Algebraica.csec(Numerical(6, 7.2)))
print('csc(6 + 7.2 i) =', Algebraica.ccsc(Numerical(6, 7.2)))
print('tan(6 + 7.2 i) =', Algebraica.ctan(Numerical(6, 7.2)))
print('cot(6 + 7.2 i) =', Algebraica.ccot(Numerical(6, 7.2)))
print()
# and trig functions

print('arcsin(6 + 7.2 i) =', Algebraica.carcsin(Numerical(6, 7.2)))
print('arccos(6 + 7.2 i) =', Algebraica.carccos(Numerical(6, 7.2)))
print('arcsec(6 + 7.2 i) =', Algebraica.carcsec(Numerical(6, 7.2)))
print('arccsc(6 + 7.2 i) =', Algebraica.carccsc(Numerical(6, 7.2)))
print('arctan(6 + 7.2 i) =', Algebraica.carctan(Numerical(6, 7.2)))
print('arccot(6 + 7.2 i) =', Algebraica.carccot(Numerical(6, 7.2)))
print()
# and inverse trig

print('sinh(6 + 7.2 i) =', Algebraica.csinh(Numerical(6, 7.2)))
print('cosh(6 + 7.2 i) =', Algebraica.ccosh(Numerical(6, 7.2)))
print('sech(6 + 7.2 i) =', Algebraica.csech(Numerical(6, 7.2)))
print('csch(6 + 7.2 i) =', Algebraica.ccsch(Numerical(6, 7.2)))
print('tanh(6 + 7.2 i) =', Algebraica.ctanh(Numerical(6, 7.2)))
print('coth(6 + 7.2 i) =', Algebraica.ccoth(Numerical(6, 7.2)))
print()
# and hyperbolic trig

print('arcsinh(6 + 7.2 i) =', Algebraica.carcsinh(Numerical(6, 7.2)))
print('arccosh(6 + 7.2 i) =', Algebraica.carccosh(Numerical(6, 7.2)))
print('arcsech(6 + 7.2 i) =', Algebraica.carcsech(Numerical(6, 7.2)))
print('arccsch(6 + 7.2 i) =', Algebraica.carccsch(Numerical(6, 7.2)))
print('arctanh(6 + 7.2 i) =', Algebraica.carctanh(Numerical(6, 7.2)))
print('arccoth(6 + 7.2 i) =', Algebraica.carccoth(Numerical(6, 7.2)))
print()
# and inverse hyperbolic trig

# all functions return the principal value

# you can also convert a class Complex object into a class Numerical object.
# for more information, visit example 9.

input()
