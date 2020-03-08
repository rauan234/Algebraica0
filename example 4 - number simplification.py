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

print(cp( ((1, 1), ('sin', cp( ((-1, 7), 'pi') ))) ), '= ', end='')
print(cp( ((1, 1), ('sin', cp( ((-1, 7), 'pi') ))) ).simplified())
print()
# minus signs in sin and tan are taken out

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

Algebraica.MultiplyPowers = True
F = (cp(0, 1) * D) ** cp((7, 3))
print('F =', F)
print(F.simplified())
print()
# if you want, you can turn on MultiplyPowers
# it allows to simplify things even further
# for example, (x ** 2) ** 1/2 becomes x

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
# sqrt is automatically changed to sqrt

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

O = cp( ((1, 1), ('pow', cp('e'), cp((3, 2)))) )
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

print('sign( e * pi * ln(4) ) =', cp( ((1, 1), ('sign', cp( ((1, 1), ['e', 'pi', ('ln', cp(4))]) ))) ).simplified())
print()
# if a sign function takes product of several real functions as an argument,
# the result is a product of signs of the functions

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

print('exp(0) =', cp( ((1, 1), ('exp', cp(0))) ).simplified())
print('exp(1) =', cp( ((1, 1), ('exp', cp(1))) ).simplified())
print()
# exp functions are calculated if the argument is either 0 or 1

print('ln(2/3 pi) =', cp( ((1, 1), ('ln', cp((2, 3)) * cp('pi'))) ).simplified())
print('ln(pow(pi, 3)) =', cp( ((1, 1), ('ln', cp( ((1, 1), ('pow', cp('pi'), cp(3))) ))) ).simplified())
print()
# logarithms are simplified in the following way

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

# as you might have noticed, the .simplify() function is very primitive
# and clearly needs modernization
# for now, it is the Ahilles heel of Algebraica as it has a lot of bugs in it, so it`s unreliable
# also, the .simplify() method is not very effective as often it does not help
# to make thins look simpler
input()
