from Algebraica import cp
from Algebraica import Rational, Function, Term, Irrational, Complex


# the mosr basic subsidiary class in Algebraica is Rational
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
# here are the basic arithmetic operations between two Rational

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
# you can also do arithmetics with integers and floats

C = Rational(7)
print('C =', C)
print()
# of you want, you can create Rational class object by putting in one number
# in this case, the number you entered will be the numerator, and the
# denominator will be set equal to one

D = Rational(6.3)
print('D =', D)
print()
# in a similar way, you can declare a Rational using a float number
# then, the float will be automatically converted into a fraction

C += 1.2
D += 4
print('C =', C)
print('D =', D)
print()

C -= 3
D -= 7.4
print('C =', C)
print('D =', D)

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
# class Rational objects can be compared with integers and floats

print('-C =', -C)
print('-D =', -D)
print()
# by putting a minus sign in front of the Rational, you can change its` sign

print('~C =', ~C)
print('~D =', ~D)
print()
# by putting a tilda, you can reverse the fraction, e.g. switch the numerator
# and the denominator, so that, for example, 3/4 becomes 4/3

print('2/4 =', Rational(2, 4))
# note that the fractions are automatically simplified

print('1/(-2) =', Rational(1, -2))
print()
# also, the minus sign is automatically transfered from
# the numerator to the denominator

del A
del B
del C
del D


# the next subsidiary class is called Function
# it carries information about one specific function
# more precisely, it tells us the name of the function and what
# arguments were given to it, such as sqrt[ 2 ]

A = Function( ( 'sqrt', cp(2) ) )
print('A =', A)
print()
# in order to declare a Function, we need to enter a tuple as an argument for __init__
# the 0th element of the tuple is the name of the function, always
# the first and further elements of the tuple are the arguments of the function
# they must all be Complex

B = Function( ( 'pi' ) )
C = Function(   'pi'   )
print('B =', B)
print('C =', C)
# sometimes you want you function to just be constant, so that there is no arguments
# you can do it by either giving a tuple with only one element or by entering the name of the constant
# there is no difference, because __init__ treats both equally

print('(A == B) =', A == B)
print('(B == C) =', B == C)
print('(A == C) =', A == C)
print('(sqrt(2) == sqrt(3)) =', Function( ( 'sqrt', cp(2) ) ) == 
      Function( ( 'sqrt', cp(3) ) ))
print('exp(1) == exp(1) =', Function( ( 'exp', cp(1) ) ) == 
      Function( ( 'exp', cp(1) ) ))
print()
# function class does not have much methods
# all it can do is to be created, print itself and compare itself to others
# the comparison is conducted by comparing the names and the arguments of the two functions
# if their names or any of the arguments are different, the result is False
# otherwise it`s true

del A
del B
del C


# by combining Rational and Function we get Term
# Term is basically a Fational coefficient multiplied by
# a series of Fuction

R = Rational( 7, 2 )
functions = [
    Function( ('sqrt', cp(3)) ),
    Function( ('ln', cp(10)) ),
    Function( ( 'pow', cp(2), cp( (2, 3) ) ) ),
    Function( 'pi' )
]
A = Term( R, functions )
del R
del functions
print('A =', A)
print()
# in order to declare a Term, write an expression of form
# Term( R, functions ), where R is the rational coefficient and
# functions is a list of all functions that are multiplied together
# the end result is R * functions[0] * functions[1]...

R = Rational( 3, 8 )
functions = [
    Function( ( 'sqrt', cp('pi') ) )
]
B = Term( R, functions )
print('B =', B)
print()

# ! print(A + B) !
# ! print(A - B) !
# ! print(A / B) !
# Term objects cannot be added, subtracted or divided

print('A * B =', A * B)
print()
# they can only be multiplied
# during multiplication, the Rational coefficients of the Term objects
# are multiplied and their functions lists are added together
# the reason why Terms can only be multiplied is because that`s the only
# way they can be utilized. in other words, we simply don`t need
# anything except multiplication

del A
del B


# the reason why Terms only need to be multiplied is because they are just
# slaves of class Irrational, which in turn takes care of addition, subtraction,
# multiplication, division and lots of other things
# a class Irrational object is simply a sum of many class Term objects

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
    [ Function( ('ln', cp(3)) ) ]
)
A = Irrational( t1a, t2a, t3a )
print('A =', A)
print()
# in order to declare a class Irrational object, enter a few
# class Term objects in __init__ method and that`s it
# the result will be t1a + t2a + t3a
# of course, you can have not three but any natural number of terms

B = Irrational(
    Term(
        Rational(2),
        [
            Function( ('sqrt', cp(2) ) ),
            Function( ('exp', cp(0, 1) ) )
        ]
    )
)
print('B =', B)
print()

print('A + B =', A + B)
print('A - B =', A - B)
print('A * B =', A * B)
print()
# here is some basic arithmetics
# when adding two Irrationals, the result is simply a class Irrational
# object whose terms are the sum of the terms of the first two Irrationals
# when subtracting, the procedure is the same, except that all terms of B
# undergo sign change before being added to the terms of A
# multiplication of two class Irrational object happens just by
# expanding the bracket, e.g. (t1a + t2a + t3a) * (t1b + t2b) =
# = t1a * t1b + t1a * t2b + t2a * t1b + t2a * t2b + t3a * t1b + t3a * t2b
# as you can see, Irrational class objects cannot be divided
# class Complex will take care of it later

print('(A == B) =', A == B)
print('(A != B) =', A != B)
print()
# class Irrational objects can be compared

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
# ! print(A / 'pi')
# Irrationals can even be multiplied by strings, which are automaticaly
# turned into class Function objects
# though dividing by string is imposible

del A
del B


# now let`s look at the main class, for which everything is built for - Complex
rt1c = Rational( 7, 2 )
rt2c = Rational( 4, 3 )
ti1c = Rational( 9, 13 )
tr2f1 = Function( ('pi',) )
ti1f1 = Function( ('sqrt', cp(2)) )
ti1f2 = Function( ('ln', cp(3)) )
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
# is the imaginary part
# as you can see, in order to construct one single Complex we need all that arsenal
# we`ve prepared

B = Complex(
    Irrational(
        Term(
            Rational( 1, 2 ),
            [
                Function( ('a', cp(3)) ),
                Function( ('b', cp(4)) )
            ]
            ),
        Term(
            Rational( 5, 6 ),
            [
                Function( ('c', cp(7)) ),
                Function( ('d', cp(8)) )
            ]
            )
        ),
    Irrational(
        Term(
            Rational( 9, 10 ),
            [
                Function( ('f', cp(11)) ),
                Function( ('g', cp(12)) )
            ]
            ),
        Term(
            Rational( 13, 14 ),
            [
                Function( ('h', cp(15)) ),
                Function( ('j', cp(16)) )
            ]
            ),
        Term(
            Rational( 17, 18 ),
            [
                Function( ('k', cp(19)) ),
                Function( ('l', cp(20)) )
            ]
            )
        )
)
print('B =', B)
print()
# here is how you would create a Complex by using __init__
# but if you don`t need to build that complicated numbers, you might use a little
# shortcut prepared to make life easier, namely - the cp(*args) function

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
# you can also create Complex by calling the cp(inp) function
# with argument being an object of one of the subsidiary classes

H = cp((8, 3), ((1, 1), ('sqrt', cp(2))))
J = Rational(4, 7)
print('H * J =', H * J)
print('H / J =', H / J)
print()
# you can perform multiplication and division with Complex and Rational

# so that`s how it is
# this example demonstrates the structre of Complex
# in order to see the examples of operators and other things, visit example 1

input()
