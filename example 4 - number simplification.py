from Algebraica import cp
import Algebraica


A = cp( ((1, 1), ('pow', cp(36), cp( (9, 5) ))) )
print('A =', A)
print(A.simplified())
print()
# simplify() method, as you can judje by its` name simplifies a Complex
# one of the examples of such simplification is power function

B = cp( ((1, 1), ('ln', cp(1))) )
print('B =', B)
print(B.simplified())
print()
# ln of 1 is automatically converted to 0

C = cp( ((1, 1), ('ln', cp('e'))) )
print('C = ', end='')
print(C)
C.simplify()
print(C)
print()
# and ln of e is converted to 1

print('sin   0 =', cp( ((1, 1), ('sin', cp( ((0, 1), 'pi') ))) ).simplified())
print('sin 180 =', cp( ((1, 1), ('sin', cp( ((1, 1), 'pi') ))) ).simplified())
print('sin 120 =', cp( ((1, 1), ('sin', cp( ((2, 3), 'pi') ))) ).simplified())
print('sin  90 =', cp( ((1, 1), ('sin', cp( ((1, 2), 'pi') ))) ).simplified())
print('sin  60 =', cp( ((1, 1), ('sin', cp( ((1, 3), 'pi') ))) ).simplified())
print('sin  45 =', cp( ((1, 1), ('sin', cp( ((1, 4), 'pi') ))) ).simplified())
print('sin  30 =', cp( ((1, 1), ('sin', cp( ((1, 6), 'pi') ))) ).simplified())
print()
# here are some sin values

print('cos   0 =', cp( ((1, 1), ('cos', cp( ((0, 1), 'pi') ))) ).simplified())
print('cos 180 =', cp( ((1, 1), ('cos', cp( ((1, 1), 'pi') ))) ).simplified())
print('cos 120 =', cp( ((1, 1), ('cos', cp( ((2, 3), 'pi') ))) ).simplified())
print('cos  90 =', cp( ((1, 1), ('cos', cp( ((1, 2), 'pi') ))) ).simplified())
print('cos  60 =', cp( ((1, 1), ('cos', cp( ((1, 3), 'pi') ))) ).simplified())
print('cos  45 =', cp( ((1, 1), ('cos', cp( ((1, 4), 'pi') ))) ).simplified())
print('cos  30 =', cp( ((1, 1), ('cos', cp( ((1, 6), 'pi') ))) ).simplified())
print()
# here are cosines

print('tan   0 =', cp( ((1, 1), ('tan', cp( ((0, 1), 'pi') ))) ).simplified())
print('tan 180 =', cp( ((1, 1), ('tan', cp( ((1, 1), 'pi') ))) ).simplified())
print('tan 120 =', cp( ((1, 1), ('tan', cp( ((2, 3), 'pi') ))) ).simplified())
print('tan  60 =', cp( ((1, 1), ('tan', cp( ((1, 3), 'pi') ))) ).simplified())
print('tan  45 =', cp( ((1, 1), ('tan', cp( ((1, 4), 'pi') ))) ).simplified())
print('tan  30 =', cp( ((1, 1), ('tan', cp( ((1, 6), 'pi') ))) ).simplified())
print()
# and tangents

print('ctg 120 =', cp( ((1, 1), ('ctg', cp( ((2, 3), 'pi') ))) ).simplified())
print('ctg  60 =', cp( ((1, 1), ('ctg', cp( ((1, 3), 'pi') ))) ).simplified())
print('ctg  45 =', cp( ((1, 1), ('ctg', cp( ((1, 4), 'pi') ))) ).simplified())
print('ctg  30 =', cp( ((1, 1), ('ctg', cp( ((1, 6), 'pi') ))) ).simplified())
print()
# and cotangents

print('arcsin 0 =', cp( ((1, 1), ('arcsin', cp( ((0, 1)) ))) ).simplified())
print('arcsin 1/2 =', cp( ((1, 1), ('arcsin', cp( ((1, 2)) ))) ).simplified())
print('arcsin 1/sqrt(2) =', cp( ((1, 1), ('arcsin', cp( ((1, 2), ('sqrt', cp(2))) ))) ).simplified())
print('arcsin sqrt(3)/2 =', cp( ((1, 1), ('arcsin', cp( ((1, 2), ('sqrt', cp(3))) ))) ).simplified())
print('arcsin 1 =', cp( ((1, 1), ('arcsin', cp( ((1, 1)) ))) ).simplified())
print()
# here are some arcsin values

print('arccos 0 =', cp( ((1, 1), ('arccos', cp( ((0, 1)) ))) ).simplified())
print('arccos 1/2 =', cp( ((1, 1), ('arccos', cp( ((1, 2)) ))) ).simplified())
print('arccos 1/sqrt(2) =', cp( ((1, 1), ('arccos', cp( ((1, 2), ('sqrt', cp(2))) ))) ).simplified())
print('arccos sqrt(3)/2 =', cp( ((1, 1), ('arccos', cp( ((1, 2), ('sqrt', cp(3))) ))) ).simplified())
print('arccos 1 =', cp( ((1, 1), ('arccos', cp( ((1, 1)) ))) ).simplified())
print()
# here are arccosines

print('arctan 0 =', cp( ((1, 1), ('arctan', cp( ((0, 1)) ))) ).simplified())
print('arctan 1/sqrt(3) =', cp( ((1, 1), ('arctan', cp( ((1, 3), ('sqrt', cp(3))) ))) ).simplified())
print('arctan 1 =', cp( ((1, 1), ('arctan', cp( ((1, 1)) ))) ).simplified())
print('arctan sqrt(3) =', cp( ((1, 1), ('arctan', cp( ((1, 1), ('sqrt', cp(3))) ))) ).simplified())
print()
# and arctangents

print('sinh 0 =', cp( ((1, 1), ('sinh', cp(0)))).simplified())
print('cosh 0 =', cp( ((1, 1), ('cosh', cp(0)))).simplified())
print('tanh 0 =', cp( ((1, 1), ('tanh', cp(0)))).simplified())
print()
# some hyperbolic trig function values

print('sinh -f =', cp( ((1, 1), ('sinh', -cp('f'))) ).simplified())
print('cosh -f =', cp( ((1, 1), ('cosh', -cp('f'))) ).simplified())
print('tanh -f =', cp( ((1, 1), ('tanh', -cp('f'))) ).simplified())
print()
# hyperbolic trig functions of negative argument

print(cp( ((1, 1), ('sin', cp( ((-1, 7), 'pi') ))) ), '= ', end='')
print(cp( ((1, 1), ('sin', cp( ((-1, 7), 'pi') ))) ).simplified())
print()
# minus signs in sin and tan are taken out

print(cp( ((1, 1), ('cos', cp((1, 4))))).simplified())
print(cp( ((1, 1), ('ctg', cp((1, 4))))).simplified())
print()
# cos and ctg are replaced with sin and tan
# if you want to turn it off, switch ReplaceTrig to False
# for more info, look up example 6

print('sin(11/3 pi) =', cp( ((1, 1), ('sin', cp( ((11, 3), 'pi') ))) ).simplified())
print('cos(11/3 pi) =', cp( ((1, 1), ('cos', cp( ((11, 3), 'pi') ))) ).simplified())
print('tan(11/3 pi) =', cp( ((1, 1), ('tan', cp( ((11, 3), 'pi') ))) ).simplified())
print('ctg(11/3 pi) =', cp( ((1, 1), ('ctg', cp( ((11, 3), 'pi') ))) ).simplified())
print()
# trig functions of type sin(2 * pi + f) are simplified to sin(f)

D = cp( ((1, 1), ('pow', cp( ((1, 1), [
    ('sqrt', cp(3)),
    ('ln', cp(2))
    ]) ), cp( (9, 4) )) ) )
print('D =', D)
print(D.simplified())
print()
# factors inside a power function are taken out into endarate
# power functions

E = (cp(0, 1) * D) ** cp((7, 3))
print('E =', E)
print(E.simplified())
print()
# arguments of functions are also automatically simplified

print('(pi ** 1/2) ** 2 =', cp( ((1, 1), ('pow', cp( ((1, 1), ('pow', cp('pi'), cp((1, 2)))) ), cp(2)) )).simplified())
print('(a ** 1/2) ** 2 =',cp( ((1, 1), ('pow', cp( ((1, 1), ('pow', cp('a'), cp((1, 2)))) ), cp(2)) )).simplified())
print()
# pow( pow(a, b), c ) is simplified to pow(a, b * c) if a is always positive

print('(-f) ** 12 =', cp( ((1, 1), ('pow', -cp('f'), cp(12))) ).simplified())
print('(-f) ** 11 =', cp( ((1, 1), ('pow', -cp('f'), cp(11))) ).simplified())
print('(if) ** 12 =', cp( ((1, 1), ('pow', cp(0, 'f'), cp(12))) ).simplified())
print('(if) ** 11 =', cp( ((1, 1), ('pow', cp(0, 'f'), cp(11))) ).simplified())
print()
# that`s how minus signs and imaginary units are taken out of the power functions

F = cp( ((1, 1), ('pow', D, E**2)) )
print('F =', F)
print(F.simplified())
print()

G = cp( ((1, 1), ('pow', cp( ((1, 1), ('pow', cp('f'), cp((2)))) ), cp((1, 2)))) )
print('G =', G)
print(G.simplified())
print()
# like that

H = cp( ((1, 1), ('ln', cp( ((7, 12), [
    ('sqrt', cp(3)),
    ('sin', cp(2))
])) )) )
print('H =', H)
print(H.simplified())
print()
# factors inside logarithms are taken out

I = cp( ((1, 1), ('ln', cp((2, 13)) )) )
print('I =', I)
print(I.simplified())
print()
# just like that

J = cp( ((1, 1), ('exp', G)) )
print('J =', J)
print(J.simplified())
print()
# if there is a ln inside an exp, the thing is simplified

K = cp( ((1, 1), ('pow', cp(6), cp((1, 3)))) )
print('K =', K)
print(K.simplified())
print()
# numbers inside pow functions are desomposed into primes

L = cp( ((1, 1), ('sqrt', cp(2))) )
print('L =', L)
print(L.simplified())
print()
# sqrt(x) is automatically changed to pow(x, 1/2)
# during displaying, though, it`s shown as sqrt

M = cp( ((1, 1), ('log', cp(3), cp(10))) )
print('M =', M)
print(M.simplified())
print()
# log of arbitrary base is expressed through ln

N = cp( ((3, 1), ('pow', cp((2, 7)), cp(-2))) )
print('N =', N)
print(N.simplified())
print()
# power functions with integer exponents are calculated

O = cp( ((1, 1), ('pow', cp(36) * cp('e') * cp('pi'), cp((3, 2)))) )
print('O =', O)
print(O.simplified())
print()
# exponents with base e are written as exp

print('sign( 9/8 ) =', cp( ((1, 1), ('sign', cp((9, 8)))) ).simplified())
print('sign( -1/5 ) =', cp( ((1, 1), ('sign', cp((-1, 5)))) ).simplified())
print()
# signs of some numbers are calculated on spot

print('sign( pi ) =', cp( ((1, 1), ('sign', cp('pi'))) ).simplified())
print('sign( e ) =', cp( ((1, 1), ('sign', cp('e'))) ).simplified())
print()
# some functions are known to always be positive

print('sign(exp(3)) =', cp( ((1, 1), ('sign', cp( ((1, 1), ('exp', cp(3))) ))) ).simplified())
print('sign(cosh(-2)) =', cp( ((1, 1), ('sign', cp( ((1, 1), ('cosh', cp(-2))) ))) ).simplified())
print()
# some functions are known to be positive if their argument is real

print('sign(ln(2)) =', cp( ((1, 1), ('sign', cp( ((1, 1), ('ln', cp(2))) ))) ).simplified())
print('sign(sinh(3)) =', cp( ((1, 1), ('sign', cp( ((1, 1), ('sinh', cp(3))) ))) ).simplified())
print('sign(arctan(2)) =', cp( ((1, 1), ('sign', cp( ((1, 1), ('arctan', cp(2))) ))) ).simplified())
print()
# some functins are known to be positive/negative for certain values of their arguments

print('sign(pow(2, 13/4)) =', cp( ((1, 1), ('sign', cp( ((1, 1), ('pow', cp(2), cp((13, 4)))) ))) ).simplified())
print('sign(pow(C, 13/4)) =', cp( ((1, 1), ('sign', cp( ((1, 1), ('pow', cp('C'), cp((13, 4)))) ))) ).simplified())
print('sign(pow(2, D)) = ', cp( ((1, 1), ('sign', cp( ((1, 1), ('pow', cp(2), cp('D'))) ))) ).simplified())
print()
# power function is positive if the power is real and the base is real and positive
# if there is no guarantee if both of these consitions hold, the sign remains unknown

print('sign( -e * pi * ln(4) ) =', cp( ((1, 1), ('sign', cp( ((-1, 1), ['e', 'pi', ('ln', cp(4))]) ))) ).simplified())
print()
# if a sign function takes product of several real functions as an argument,
# the result is a product of signs of the functions

print('abs(-n) =', cp( ((1, 1), ('abs', -cp('n'))) ).simplified())
print('abs(n i) =', cp( ((1, 1), ('abs', cp(0, 'n'))) ).simplified())
print()
# abs(-n) = abs(n)
# abs(n i) = abs(n)

print('abs(a * b * c) =', cp( ((1, 1), ('abs', cp( ((1, 1), ['a', 'b', 'c']) ))) ).simplified())
print()
# abs(a * b * c...) = abs(a) * abs(b) * abs(c)...

print('abs(3) =', cp( ((1, 1), ('abs', cp(3))) ).simplified())
print('abs(i) =', cp( ((1, 1), ('abs', cp(0, 1))) ).simplified())
print()

print('abs(-pi) =', cp( ((1, 1), ('abs', -cp('pi'))) ).simplified())
print('abs(f) =', cp( ((1, 1), ('abs', cp('f'))) ).simplified())
print()
# abs(r) = +-r, where r is a real number
# if there is no guarantee r is real, simplification is not done

print('abs(e + 3 i) =', cp( ((1, 1), ('abs', cp('e', 3))) ).simplified())
print('abs(2 + 9 i) =', cp( ((1, 1), ('abs', cp(2, 9))) ).simplified())
print()
# absolute values of type( a + b i ), where a and b are real numbers
# are simplified as sqrt( a**2 + b**2 )

print('Rl( 3 + 8 i) =', cp( ((1, 1), ('Rl', cp(3, 8) )) ).simplified())
print('Im( 3 + 8 i) =', cp( ((1, 1), ('Im', cp(3, 8) )) ).simplified())
print('Rl( -pi ) =', cp( ((1, 1), ('Rl', cp( ((-1, 1), 'pi') ) )) ).simplified())
print('Im( -pi ) =', cp( ((1, 1), ('Im', cp( ((-1, 1), 'pi') ) )) ).simplified())
print('Rl( 2/3 e ) =', cp( ((1, 1), ('Rl', cp( ((2, 3), 'e') ) )) ).simplified())
print('Im( 2/3 e ) =', cp( ((1, 1), ('Im', cp( ((2, 3), 'e') ) )) ).simplified())
print()
# Rl and Im functions are simplified if their argument is well-behaved
# explicitly 'well-behaveness' means that the real part of the Complex
# is guaranteed to be real, and the same applies to the imaginary part

print('Rl( ln(2) + 9 + 4 i ) =', cp( ((1, 1), ('Rl', cp([ ('ln', cp(2)), 9 ], 4) )) ).simplified())
print('Rl( ln(2) + 9 + 4 i ) =', cp( ((1, 1), ('Im', cp([ ('ln', cp(2)), 9 ], 4) )) ).simplified())
print()
# if Rl or Im functions take sum as an argument, they are simplified in the above way

print('Rl(b i) =', cp( ((1, 1), ('Rl', cp(0, 'b')))).simplified())
print('Im(b i) =', cp( ((1, 1), ('Im', cp(0, 'b')))).simplified())
# Rl and Im are also simplified like that

print('pow(3 + 7 i, 8) =', cp( ((1, 1), ('pow', cp(3, 7), cp(8))) ).simplified())
print()
# expressions of type pow( a + b i, n ), where a and b are rational and n is real,
# are calculated

print('exp(0) =', cp( ((1, 1), ('exp', cp(0))) ).simplified())
print('exp(1) =', cp( ((1, 1), ('exp', cp(1))) ).simplified())
print()
# exp functions are calculated if the argument is either 0 or 1

print('exp(39/11 i) =', cp( ((1, 1), ('exp', cp(0, (39, 11)))) ).simplified())
print('exp(a i) =', cp( ((1, 1), ('exp', cp(0, 'a'))) ).simplified())
print('exp(pi i) =', cp( ((1, 1), ('exp', cp(0, 'pi'))) ).simplified())
print()
# expressions of type exp(f i) are re-written as cos(f) + i sin(f),
# if f is guaranteed to be a real number
# you can turn this off using ExpandComplexExponents setting (view example 6)

print('pow(exp(a), b) =', cp( ((1, 1), ('pow', cp( ((1, 1), ('exp', cp('a'))) ), cp('b'))) ).simplified())
print()
# pow(exp(a), b) = exp(a * b)

print('ln(24) =', cp( ((1, 1), ('ln', cp(24))) ).simplified())
# numbers inside logarithms are decomposed into primes

print('ln(2/3 pi) =', cp( ((1, 1), ('ln', cp((2, 3)) * cp('pi'))) ).simplified())
print('ln(pow(pi, 3)) =', cp( ((1, 1), ('ln', cp( ((1, 1), ('pow', cp('pi'), cp(3))) ))) ).simplified())
print()
# logarithms are simplified in the following way

print('arg(2) =', cp( ((1, 1), ('arg', cp(2))) ).simplified())
print('arg(3 i) =', cp( ((1, 1), ('arg', cp(0, 3))) ).simplified())
print('arg(-3 i) =', cp( ((1, 1), ('arg', cp(0, -3))) ).simplified())
print('arg(2 + 3 i) =', cp( ((1, 1), ('arg', cp(2, 3))) ).simplified())
print('arg(-2 + 3 i) =', cp( ((1, 1), ('arg', cp(-2, 3))) ).simplified())
print('arg(2 + -3 i) =', cp( ((1, 1), ('arg', cp(2, -3))) ).simplified())
print('arg(-2 + -3 i) =', cp( ((1, 1), ('arg', cp(-2, -3))) ).simplified())
print()
# arg(a + b i) is simplified to arctan(b/a) if a an b are real and rational and a != 0
# if a = 0, arg(bi) = pi/2 * sign(b)

print('arcsin(-2/3 f) =', cp( ((1, 1), ('arcsin', cp(((-2, 3), 'f')))) ).simplified())
print('arccos(-2/3 f) =', cp( ((1, 1), ('arccos', cp(((-2, 3), 'f')))) ).simplified())
print('arctan(-2/3 f) =', cp( ((1, 1), ('arctan', cp(((-2, 3), 'f')))) ).simplified())
# that`s how arc-trig functions of negative arguments are simplified

P = cp(
    [(2, 3), ((7, 5), 'pi'), ((4, 3), ('exp', cp(3))),
         ((1, 2), [
             'pi',
             ('exp', cp(4)),
             ('sqrt', cp((1, 1), [
                 1,
                 (3, 7),
                 ((9, 1), ('log', cp(3), cp(11))),
                 ((5, 2), ('pow', cp(1, 2), cp(3, 4)))
             ]))  # these are added
         ])  # these things are multiplied
     ],  # these things are added
    0
)
Q = cp( ((1, 1), ('exp', P)) )
print('Q =', Q)
print(Q.simplified())
print()
# this is an example of a complicated thing being simplified

print((cp(0, 1) ** cp(0, 1)).simplified())
print()
# and, of course, we can also do i ** i

R = cp( ((1, 1), ('log', cp(3), cp('e'))) )
print('R =', R)
R.simplify()
print(R)
print()
# you might want to use .simplify() instead of .simplified()
# the difference between them is that .simplified() does not modify the number
# and returns its` simplified version
# .simplify(), on the other hand, simplifies the initial number
# and does not return anything

'''
S = cp(6, -2)
T = cp(19, 3)
U = S ** 3 + T ** 4 + T * S
V = (S / T**2) ** U
U += V
W = cp(('exp', U))
print('W =', W)
print(W.simplified())
print()
# this is how complicated it might get
# be patient with this one as it takes quite a while to print
'''

# as you might have noticed, the .simplify() function is very primitive
# and clearly needs modernization
# for now, it is the Ahilles heel of Algebraica as it has a lot of bugs in it, so it`s unreliable
# also, the .simplify() method is not very effective as often it does not help
# to make thins look simpler
input()
