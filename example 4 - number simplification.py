from Algebraica import cp
import Algebraica

# Algebraica automatically simplifies the answers.
Algebraica.AutoSimplify = False
# in order to explore how simplification is done,
# let`s turn off the auto-simplification function.
# to learn more, visit example 6.


print('36 ** (9/5) =', cp( (1, ('pow', cp(36), cp( (9, 5) ))) ).simplified())
print()
# this is how 36 ** (9/5) is simplified.
# this is done by first decomposing 36 into its` prime components and then using
# the exponential law, pow(a * b, p) = pow(a, p) * pow(b, p)

print('pow(a * b * c, P) =', cp( (1, ('pow', cp( (1, ['a', 'b', 'c']) ), cp('P'))) ).simplified())
print()
# that`s one more example of the exponential law, this time with letters

print('pow(a, b) * pow(a, c) =', (cp( (1, ('pow', cp('a'), cp('b'))) ) *
                                  cp( (1, ('pow', cp('a'), cp('c'))) )).simplified())
print()
# pow(a, b) * pow(a, c) = pow(a, b + c)

print('a * a =', (cp('a') * cp('a')).simplified())
print()
# a * a = a ** 2

print('pow(a, b) * a =', (cp( (1, ('pow', cp('a'), cp('b'))) ) * cp('a')).simplified())
print()
# pow(a, b) * a = pow(a, b + 1)

print('pow( 8 a + 42 b, 1/3) =', cp( (1, ('pow', cp((8, 'a')) + cp((42, 'b')), cp((1, 3)) )) ).simplified())
print()
# common rational coefficients are taken out

print('(pi ** 1/2) ** 2 =', cp( (1, ('pow', cp( (1, ('pow', cp('pi'), cp((1, 2)))) ), cp(2)) )).simplified())
print('(a ** 1/2) ** 2 =',cp( (1, ('pow', cp( (1, ('pow', cp('a'), cp((1, 2)))) ), cp(2)) )).simplified())
print()
# pow( pow(a, b), c ) is simplified to pow(a, b * c) if a is always positive

print('(-f) ** 12 =', cp( (1, ('pow', -cp('f'), cp(12))) ).simplified())
print('(-f) ** 11 =', cp( (1, ('pow', -cp('f'), cp(11))) ).simplified())
print('(if) ** 12 =', cp( (1, ('pow', cp(0, 'f'), cp(12))) ).simplified())
print('(if) ** 11 =', cp( (1, ('pow', cp(0, 'f'), cp(11))) ).simplified())
print()
# that`s how minus signs and imaginary units are taken out of the power functions

print('pow(3 + 5 i, 2) =', cp( ('pow', cp(3, 5), cp(2)) ).simplified())
print('pow(3 + 5 i, -2) =', cp( ('pow', cp(3, 5), cp(-2)) ).simplified())
print()
# power functions with integer powers and bases of type a + bi,
# where a and b are rational, are calculated

print('pow(e, P) =', cp( ('pow', cp('e'), cp('P')) ).simplified())
print()
# pow(e, P) = exp(P)
# (by the definition of an exponential function)


print('exp(0) =', cp( (1, ('exp', cp(0))) ).simplified())
print('exp(1) =', cp( (1, ('exp', cp(1))) ).simplified())
print()
# exp functions are calculated if the argument is either 0 or 1

print('exp(ln(Z)) =', cp( (1, ('exp', cp( (1, ('ln', cp('Z'))) ))) ).simplified())
print()
# exp(ln(Z)) = Z

print('exp(2 + 4 ln(3) + 8 i) =', cp( (1, ('exp', cp(2) + cp( ((4, 1), ('ln', cp(3))) ) + cp(0, 8)))).simplified())
print()
# expressions like exp(P + ln(Z)) are simplified to exp(P) * Z

print('exp(2 + 39/11 i) =', cp( (1, ('exp', cp(2, (39, 11)))) ).simplified())
print('exp(3 + a i) =', cp( (1, ('exp', cp(3, 'a'))) ).simplified())
print('exp(e + pi i) =', cp( (1, ('exp', cp('e', 'pi'))) ).simplified())
print()
# if ExpandComplexExponentials setting is active (and it, by default, is),
# expressions of type exp(A + f i) are re-written as exp(A) * ( cos(f) + i sin(f) ),
# if f and A are guaranteed to be real numbers.
# more information in example 6.

print('exp(a) * exp(b) =', (cp( (1, ('exp', cp('a'))) ) * cp( (1, ('exp', cp('b'))) )).simplified())
print()
# exp(a) * exp(b) is simplified to exp(a + b).
# it is done not by Function.alternate_form() but by Term.collect_exponentials().

print('pow(exp(a), 2) =', cp( (1, ('pow', cp( (1, ('exp', cp('a'))) ), cp(2))) ).simplified())
print()
# pow(exp(a), 2) = exp(2 * a)


print('ln(1) =', cp( (1, ('ln', cp(1))) ).simplified())
print('ln(-1) =', cp( (1, ('ln', cp(-1))) ).simplified())
print('ln(i) =', cp( (1, ('ln', cp(0, 1))) ).simplified())
print('ln(e) =', cp( (1, ('ln', cp('e'))) ).simplified())
print()
# natural logarithms of some arguments are calculated

print('ln(24) =', cp( (1, ('ln', cp(24))) ).simplified())
print('ln(2/3 pi) =', cp( (1, ('ln', cp((2, 3)) * cp('pi'))) ).simplified())
print()
# ln(a * b * c ...) = ln(a) + ln(b) + ln(c) ...

print('ln(pow(pi, 3)) =', cp( (1, ('ln', cp( (1, ('pow', cp('pi'), cp(3))) ))) ).simplified())
print()
# ln(pow(a, b)) = b * ln(a)

print('ln(exp(Z)) =', cp( (1, ('ln', cp( (1, ('exp', cp('Z'))) ))) ).simplified())
print()
# ln(exp(Z)) = Z


print('log10(2) =', cp( (1, ('log', cp(2), cp(10))) ).simplified())
print()
# log_B(N) = ln(N) / ln(B).
# in order to enter expression line log10(2), use cp( (1, ('log', cp(2), cp(10))) ).
# note that the base of the logarithm comes second.


print('Rl( 3 + 8 i) =', cp( (1, ('Rl', cp(3, 8) )) ).simplified())
print('Im( 3 + 8 i) =', cp( (1, ('Im', cp(3, 8) )) ).simplified())
print('Rl( -pi ) =', cp( (1, ('Rl', cp( ((-1, 1), 'pi') ) )) ).simplified())
print('Im( -pi ) =', cp( (1, ('Im', cp( ((-1, 1), 'pi') ) )) ).simplified())
print('Rl( 2/3 i e ) =', cp( (1, ('Rl', cp( 0, ((2, 3), 'e') ) )) ).simplified())
print('Im( 2/3 i e ) =', cp( (1, ('Im', cp( 0, ((2, 3), 'e') ) )) ).simplified())
print()
# Rl and Im functions are simplified if their argument has a form a + b i,
# where a and b are both guaranteed to be real

print('Rl( ln(2) + 9 + 4 i ) =', cp( (1, ('Rl', cp([ ('ln', cp(2)), 9 ], 4) )) ).simplified())
print('Im( ln(2) + 9 + 4 i ) =', cp( (1, ('Im', cp([ ('ln', cp(2)), 9 ], 4) )) ).simplified())
print()
# if Rl or Im functions take sum as an argument, they are simplified in the above shown way

print('Rl(b i) =', cp( (1, ('Rl', cp(0, 'b')))).simplified())
print('Im(b i) =', cp( (1, ('Im', cp(0, 'b')))).simplified())
print()
# Rl(b i) = -Im(b) and Im(b i) = Rl(b).
# this is done for any b, real or complex.


print('abs(-n) =', cp( (1, ('abs', -cp('n'))) ).simplified())
print('abs(n i) =', cp( (1, ('abs', cp(0, 'n'))) ).simplified())
print()
# abs(-n) = abs(n)
# abs(n i) = abs(n)

print('abs(a * b * c) =', cp( (1, ('abs', cp( (1, ['a', 'b', 'c']) ))) ).simplified())
print()
# abs(a * b * c...) = abs(a) * abs(b) * abs(c)...

print('abs(-pi) =', cp( (1, ('abs', -cp('pi'))) ).simplified())
print('abs(-f) =', cp( (1, ('abs', -cp('f'))) ).simplified())
print()
# abs(-r) = r, where r is a real number.
# if there is no guarantee r is real, simplification is not done.

print('abs(e + 3 i) =', cp( (1, ('abs', cp('e', 3))) ).simplified())
print('abs(2 + 9 i) =', cp( (1, ('abs', cp(2, 9))) ).simplified())
print()
# absolute values of type( a + b i ), where a and b are real numbers,
# are simplified as sqrt( a**2 + b**2 )


print('arg(2) =', cp( (1, ('arg', cp(2))) ).simplified())
print('arg(3 i) =', cp( (1, ('arg', cp(0, 3))) ).simplified())
print('arg(-2) =', cp( (1, ('arg', cp(-2))) ).simplified())
print('arg(-3 i) =', cp( (1, ('arg', cp(0, -3))) ).simplified())
print()
# arg(a) = pi/2 * (1 - sign(a)), if a is real.
# arg(b i) = pi/2 * sign(b), if b is real.

print('arg(-3 i) =', cp( (1, ('arg', cp(0, -3))) ).simplified())
print('arg(2 + 3 i) =', cp( (1, ('arg', cp(2, 3))) ).simplified())
print('arg(-2 + 3 i) =', cp( (1, ('arg', cp(-2, 3))) ).simplified())
print('arg(2 + -3 i) =', cp( (1, ('arg', cp(2, -3))) ).simplified())
print('arg(-2 + -3 i) =', cp( (1, ('arg', cp(-2, -3))) ).simplified())
print()
# arg(a + bi) is simplified using arctan, if a and b are both real and rational


print('sign( 9/8 ) =', cp( (1, ('sign', cp((9, 8)))) ).simplified())
print('sign( -1/5 ) =', cp( (1, ('sign', cp((-1, 5)))) ).simplified())
print('sign( 9/8 i ) =', cp( (1, ('sign', cp(0, (9, 8)))) ).simplified())
print('sign( -1/5 i ) =', cp( (1, ('sign', cp(0, (-1, 5)))) ).simplified())
print('sign( 9/8 - 1/5 i) =', cp( (1, ('sign', cp( ((9, 8)), ((-1, 5)) ))) ).simplified())
print()
# signs of some numbers are calculated on spot

print('sign( pi ) =', cp( (1, ('sign', cp('pi'))) ).simplified())
print('sign( e ) =', cp( (1, ('sign', cp('e'))) ).simplified())
print()
# some functions are known to always be positive

print('sign(exp(3)) =', cp( (1, ('sign', cp( (1, ('exp', cp(3))) ))) ).simplified())
print('sign(cosh(-2)) =', cp( (1, ('sign', cp( (1, ('cosh', cp(-2))) ))) ).simplified())
print()
# some functions are known to be positive if their argument is real

print('sign(ln(2)) =', cp( (1, ('sign', cp( (1, ('ln', cp(2))) ))) ).simplified())
print('sign(sinh(3)) =', cp( (1, ('sign', cp( (1, ('sinh', cp(3))) ))) ).simplified())
print('sign(arctan(2)) =', cp( (1, ('sign', cp( (1, ('arctan', cp(2))) ))) ).simplified())
print()
# some functins are known to be positive/negative for certain values of their arguments

print('sign(pow(2, 13/4)) =', cp( (1, ('sign', cp( (1, ('pow', cp(2), cp((13, 4)))) ))) ).simplified())
print('sign(pow(C, 13/4)) =', cp( (1, ('sign', cp( (1, ('pow', cp('C'), cp((13, 4)))) ))) ).simplified())
print('sign(pow(2, D)) = ', cp( (1, ('sign', cp( (1, ('pow', cp(2), cp('D'))) ))) ).simplified())
print()
# power function is positive if the power is real and the base is real and positive.
# if there is no guarantee if both of these consitions hold, the sign remains unknown.

print('sign( -e * pi * ln(4) ) =', cp( (1, ('sign', cp( ((-1, 1), ['e', 'pi', ('ln', cp(4))]) ))) ).simplified())
print()
# if a sign function takes product of several real functions as an argument,
# the result then is a product of signs of the functions


print('sin   0 =', cp( (1, ('sin', cp( ((0, 1), 'pi') ))) ).simplified())
print('sin 180 =', cp( (1, ('sin', cp( (1, 'pi') ))) ).simplified())
print('sin 120 =', cp( (1, ('sin', cp( ((2, 3), 'pi') ))) ).simplified())
print('sin  90 =', cp( (1, ('sin', cp( ((1, 2), 'pi') ))) ).simplified())
print('sin  60 =', cp( (1, ('sin', cp( ((1, 3), 'pi') ))) ).simplified())
print('sin  45 =', cp( (1, ('sin', cp( ((1, 4), 'pi') ))) ).simplified())
print('sin  30 =', cp( (1, ('sin', cp( ((1, 6), 'pi') ))) ).simplified())
print()
# here are some sin values

print('cos   0 =', cp( (1, ('cos', cp( ((0, 1), 'pi') ))) ).simplified())
print('cos 180 =', cp( (1, ('cos', cp( (1, 'pi') ))) ).simplified())
print('cos 120 =', cp( (1, ('cos', cp( ((2, 3), 'pi') ))) ).simplified())
print('cos  90 =', cp( (1, ('cos', cp( ((1, 2), 'pi') ))) ).simplified())
print('cos  60 =', cp( (1, ('cos', cp( ((1, 3), 'pi') ))) ).simplified())
print('cos  45 =', cp( (1, ('cos', cp( ((1, 4), 'pi') ))) ).simplified())
print('cos  30 =', cp( (1, ('cos', cp( ((1, 6), 'pi') ))) ).simplified())
print()
# here are cosines

print('sec   0 =', cp( (1, ('sec', cp( ((0, 1), 'pi') ))) ).simplified())
print('sec 180 =', cp( (1, ('sec', cp( (1, 'pi') ))) ).simplified())
print('sec 120 =', cp( (1, ('sec', cp( ((2, 3), 'pi') ))) ).simplified())
print('sec  60 =', cp( (1, ('sec', cp( ((1, 3), 'pi') ))) ).simplified())
print('sec  45 =', cp( (1, ('sec', cp( ((1, 4), 'pi') ))) ).simplified())
print('sec  30 =', cp( (1, ('sec', cp( ((1, 6), 'pi') ))) ).simplified())
print()
# and secants

print('csc 120 =', cp( (1, ('csc', cp( ((2, 3), 'pi') ))) ).simplified())
print('csc  90 =', cp( (1, ('csc', cp( ((1, 2), 'pi') ))) ).simplified())
print('csc  60 =', cp( (1, ('csc', cp( ((1, 3), 'pi') ))) ).simplified())
print('csc  45 =', cp( (1, ('csc', cp( ((1, 4), 'pi') ))) ).simplified())
print('csc  30 =', cp( (1, ('csc', cp( ((1, 6), 'pi') ))) ).simplified())
print()
# and cosecants

print('tan   0 =', cp( (1, ('tan', cp( ((0, 1), 'pi') ))) ).simplified())
print('tan 180 =', cp( (1, ('tan', cp( (1, 'pi') ))) ).simplified())
print('tan 120 =', cp( (1, ('tan', cp( ((2, 3), 'pi') ))) ).simplified())
print('tan  60 =', cp( (1, ('tan', cp( ((1, 3), 'pi') ))) ).simplified())
print('tan  45 =', cp( (1, ('tan', cp( ((1, 4), 'pi') ))) ).simplified())
print('tan  30 =', cp( (1, ('tan', cp( ((1, 6), 'pi') ))) ).simplified())
print()
# and tangents

print('cot 120 =', cp( (1, ('cot', cp( ((2, 3), 'pi') ))) ).simplified())
print('cot  90 =', cp( (1, ('cot', cp( ((1, 2), 'pi') ))) ).simplified())
print('cot  60 =', cp( (1, ('cot', cp( ((1, 3), 'pi') ))) ).simplified())
print('cot  45 =', cp( (1, ('cot', cp( ((1, 4), 'pi') ))) ).simplified())
print('cot  30 =', cp( (1, ('cot', cp( ((1, 6), 'pi') ))) ).simplified())
print()
# and cotangents


print('sin(-1/7 pi) =', cp( (1, ('sin', cp( ((-1, 7), 'pi') ))) ).simplified())
print('cos(-1/7 pi) =', cp( (1, ('cos', cp( ((-1, 7), 'pi') ))) ).simplified())
print('sec(-1/7 pi) =', cp( (1, ('sec', cp( ((-1, 7), 'pi') ))) ).simplified())
print('csc(-1/7 pi) =', cp( (1, ('csc', cp( ((-1, 7), 'pi') ))) ).simplified())
print('tan(-1/7 pi) =', cp( (1, ('tan', cp( ((-1, 7), 'pi') ))) ).simplified())
print('cot(-1/7 pi) =', cp( (1, ('cot', cp( ((-1, 7), 'pi') ))) ).simplified())
print()
# this is how the minus sign in a trig function argument is dealt with


print('sin(11/3 pi) =', cp( (1, ('sin', cp( ((11, 3), 'pi') ))) ).simplified())
print('cos(11/3 pi) =', cp( (1, ('cos', cp( ((11, 3), 'pi') ))) ).simplified())
print('sec(11/3 pi) =', cp( (1, ('sec', cp( ((11, 3), 'pi') ))) ).simplified())
print('csc(11/3 pi) =', cp( (1, ('csc', cp( ((11, 3), 'pi') ))) ).simplified())
print('tan(11/3 pi) =', cp( (1, ('tan', cp( ((11, 3), 'pi') ))) ).simplified())
print('cot(11/3 pi) =', cp( (1, ('cot', cp( ((11, 3), 'pi') ))) ).simplified())
print()
# trig functions of type sin(2 * pi + f) are simplified to sin(f)


print('tan(x)^2 * csc(x) =', cp( (1, [('pow', cp( (1, ('tan', cp('x'))) ), cp(2)), ('csc', cp('x'))]) ).simplified())
print()
# products of trigonometric functions are re-written

print('tanh(x)^2 * csch(x) =', cp( (1, [('pow', cp( (1, ('tanh', cp('x'))) ), cp(2)), ('csch', cp('x'))]) ).simplified())
print()
# as well as products of hyperbolic trig functions


print('arcsin 0 =', cp( (1, ('arcsin', cp( ((0, 1)) ))) ).simplified())
print('arcsin 1/2 =', cp( (1, ('arcsin', cp( ((1, 2)) ))) ).simplified())
print('arcsin 1/sqrt(2) =', cp( (1, ('arcsin', cp( ((1, 2), ('sqrt', cp(2))) ))) ).simplified())
print('arcsin sqrt(3)/2 =', cp( (1, ('arcsin', cp( ((1, 2), ('sqrt', cp(3))) ))) ).simplified())
print('arcsin 1 =', cp( (1, ('arcsin', cp( (1) ))) ).simplified())
print()
# here are some arcsin values

print('arccos 0 =', cp( (1, ('arccos', cp( ((0, 1)) ))) ).simplified())
print('arccos 1/2 =', cp( (1, ('arccos', cp( ((1, 2)) ))) ).simplified())
print('arccos 1/sqrt(2) =', cp( (1, ('arccos', cp( ((1, 2), ('sqrt', cp(2))) ))) ).simplified())
print('arccos sqrt(3)/2 =', cp( (1, ('arccos', cp( ((1, 2), ('sqrt', cp(3))) ))) ).simplified())
print('arccos 1 =', cp( (1, ('arccos', cp( (1) ))) ).simplified())
print()
# here are arccosines

print('arccsc 1 =', cp( (1, ('arccsc', cp(1))) ).simplified())
print('arccsc 2/3 sqrt(3) =', cp( (1, ('arccsc', cp( ((2, 3), ('sqrt', cp(3))) ))) ).simplified())
print('arccsc sqrt(2) =', cp( (1, ('arccsc', cp( (1, ('sqrt', cp(2))) ))) ).simplified())
print('arccsc 2 =', cp( (1, ('arccsc', cp(2))) ).simplified())
print('arccsc -2 =', cp( (1, ('arccsc', cp(-2))) ).simplified())
print()
# and arccosecants

print('arctan 0 =', cp( (1, ('arctan', cp( ((0, 1)) ))) ).simplified())
print('arctan 1/sqrt(3) =', cp( (1, ('arctan', cp( ((1, 3), ('sqrt', cp(3))) ))) ).simplified())
print('arctan 1 =', cp( (1, ('arctan', cp( (1) ))) ).simplified())
print('arctan sqrt(3) =', cp( (1, ('arctan', cp( (1, ('sqrt', cp(3))) ))) ).simplified())
print()
# and arctangents


print('arcsin(-2/3 f) =', cp( (1, ('arcsin', cp(((-2, 3), 'f')))) ).simplified())
print('arccos(-2/3 f) =', cp( (1, ('arccos', cp(((-2, 3), 'f')))) ).simplified())
print('arcsec(-2/3 f) =', cp( (1, ('arcsec', cp(((-2, 3), 'f')))) ).simplified())
print('arccsc(-2/3 f) =', cp( (1, ('arccsc', cp(((-2, 3), 'f')))) ).simplified())
print('arctan(-2/3 f) =', cp( (1, ('arctan', cp(((-2, 3), 'f')))) ).simplified())
print('arccot(-2/3 f) =', cp( (1, ('arccot', cp(((-2, 3), 'f')))) ).simplified())
print()
# that`s how arc-trig functions of negative arguments are simplified


print('sinh 0 =', cp( (1, ('sinh', cp(0)))).simplified())
print('cosh 0 =', cp( (1, ('cosh', cp(0)))).simplified())
print('sech 0 =', cp( (1, ('sech', cp(0)))).simplified())
print('tanh 0 =', cp( (1, ('tanh', cp(0)))).simplified())
print()
# some hyperbolic trig function values


print('sinh -f =', cp( (1, ('sinh', -cp('f'))) ).simplified())
print('cosh -f =', cp( (1, ('cosh', -cp('f'))) ).simplified())
print('sech -f =', cp( (1, ('sech', -cp('f'))) ).simplified())
print('csch -f =', cp( (1, ('csch', -cp('f'))) ).simplified())
print('tanh -f =', cp( (1, ('tanh', -cp('f'))) ).simplified())
print('coth -f =', cp( (1, ('coth', -cp('f'))) ).simplified())
print()
# that`s how hyperbolic trig functions of negative argument are simplified


print('arcsinh(0) =', cp( (1, ('arcsinh', cp(0))) ).simplified())
print('arccosh(0) =', cp( (1, ('arccosh', cp(0))) ).simplified())
print('arccosh(1) =', cp( (1, ('arccosh', cp(1))) ).simplified())
print('arcsech(1) =', cp( (1, ('arcsech', cp(1))) ).simplified())
print('arctanh(0) =', cp( (1, ('arctanh', cp(0))) ).simplified())
print('arccoth(0) =', cp( (1, ('arccoth', cp(0))) ).simplified())
print()
# here are some values of inverse hyperbolic trig functions


print('arcsinh(-f) =', cp( (1, ('arcsinh', -cp('f'))) ).simplified())
print('arccosh(-f) =', cp( (1, ('arccosh', -cp('f'))) ).simplified())
print('arcsech(-f) =', cp( (1, ('arcsech', -cp('f'))) ).simplified())
print('arccsch(-f) =', cp( (1, ('arccsch', -cp('f'))) ).simplified())
print('arctanh(-f) =', cp( (1, ('arctanh', -cp('f'))) ).simplified())
print('arccoth(-f) =', cp( (1, ('arccoth', -cp('f'))) ).simplified())
print()
# arcsinh(-f) = -arcsinh(f)


print('pow(N, 1/2) =', cp( (1, ('pow', cp('N'), cp((1, 2)))) ).simplified())
print()
# pow(N, 1/2) is displayed as sqrt(N).
# the important word here is "displayed", because the function pow(N, 1/2) is
# stored as Function( ('pow', cp('N'), cp((1, 2))) ) and not as Function( ('sqrt', cp('N')) ).

S = cp( (1, ('sqrt', cp(2))) ).simplified()
print('S =', S)
print(S.rl.terms[0].irt[0].name)
print()
# Function( ('sqrt', N) ) is converted to Function( ('pow', N, cp((1, 2))) ), but displayed as sqrt[N].
# you can see the true name of the function by accessing it directly.


A = cp( (1, ('ln',
             cp( (1, ('sin', cp('pi')/3)) ))
         ) )
B = cp( (1, ('sqrt', A)) )
print('B =', B)
print(B.simplified())
print()
# when simplifying a certain function, its` arguments are simplified as well


R = cp( (1, ('log', cp(3), cp('e'))) )
print('R =', R)
R.simplify()
print(R)
print()
# you might want to use .simplify() instead of .simplified().
# the difference between them is that .simplified() does not modify the number
# and returns its` simplified version.
# .simplify(), on the other hand, simplifies the initial number
# and does not return anything.


M = cp(4, -2)
T = cp(7, 3)
U = M ** 3 + T ** 4 + T * M
V = (M / T**2) ** U
U += V
W = cp(('exp', U))
print('W =', W)
print(W.simplified())
print()
# this is an example of how useful simplification might be.
# be careful with simplification, though, because it might have some bugs in it.


# in the next example, №5, you can learn how Algebraica lets you know if you
# have entered any invalid expressions.
input()
