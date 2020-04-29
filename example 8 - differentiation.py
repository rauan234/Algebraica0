from Algebraica import cp, function
import Algebraica
# Algebraica allows to automatically differentiate expressions.
# below you can see some examples of that.


print('d/dx[x] =', cp( (1, 'x') ).differentiated('x'))
print('d/dx[2x] =', cp( (2, 'x') ).differentiated('x'))
print('d/dx[ax] =', cp( (1, ['a', 'x']) ).differentiated('x'))
print()
# the derivative of constant times x is the constant itself

print('d/dx[x ** p] =', cp( (1, ('pow', cp('x'), cp('p'))) ).differentiated('x'))
print('d/dx[b ** x] =', cp( (1, ('pow', cp('b'), cp('x'))) ).differentiated('x'))
print()
# that`s how power functions are differentiated

print('d/dx[exp(x)] =', cp( (1, ('exp', cp('x'))) ).differentiated('x'))
print()
# d/dx[exp x] = exp x

print('d/dx[ln(x)] =', cp( (1, ('ln', cp('x'))) ).differentiated('x'))
print('d/dx[log2(x)] =', cp( (1, ('log', cp('x'), cp(2))) ).differentiated('x'))
print('d/dx[logx(b)] =', cp( (1, ('log', cp('b'), cp('x'))) ).differentiated('x'))
print()
# d/dx[ln x] = 1 / x
# d/dx[log2 x] = 1 / (ln(2) x)
# d/dx[logx(b)] = d/dx[ln(b) / ln(x)] = -ln(b) / (x ln(x)^2)

print('d/dx[Rl(f(x))] =', cp( (1, ('Rl', cp( (1, ('f', cp('x'))) ))) ).differentiated('x'))
print('d/dx[Im(f(x))] =', cp( (1, ('Im', cp( (1, ('f', cp('x'))) ))) ).differentiated('x'))
print()
# d/dx[Rl(f(x))] = Rl(d/dx[f(x)])
# d/dx[Im(f(x))] = Im(d/dx[f(x)])


print('d/dx[sin(x)] =', cp( (1, ('sin', cp('x'))) ).differentiated('x'))
print('d/dx[cos(x)] =', cp( (1, ('cos', cp('x'))) ).differentiated('x'))
print('d/dx[sec(x)] =', cp( (1, ('sec', cp('x'))) ).differentiated('x'))
print('d/dx[csc(x)] =', cp( (1, ('csc', cp('x'))) ).differentiated('x'))
print('d/dx[tan(x)] =', cp( (1, ('tan', cp('x'))) ).differentiated('x'))
print('d/dx[cot(x)] =', cp( (1, ('cot', cp('x'))) ).differentiated('x'))
print()
# here are the derivatives of trig functions

print('d/dx[arcsin(x)] =', cp( (1, ('arcsin', cp('x'))) ).differentiated('x'))
print('d/dx[arccos(x)] =', cp( (1, ('arccos', cp('x'))) ).differentiated('x'))
print('d/dx[arcsec(x)] =', cp( (1, ('arcsec', cp('x'))) ).differentiated('x'))
print('d/dx[arccot(x)] =', cp( (1, ('arccsc', cp('x'))) ).differentiated('x'))
print('d/dx[arctan(x)] =', cp( (1, ('arctan', cp('x'))) ).differentiated('x'))
print('d/dx[arccot(x)] =', cp( (1, ('arccot', cp('x'))) ).differentiated('x'))
print()
# and inverse trig

print('d/dx[sinh(x)] =', cp( (1, ('sinh', cp('x'))) ).differentiated('x'))
print('d/dx[cosh(x)] =', cp( (1, ('cosh', cp('x'))) ).differentiated('x'))
print('d/dx[sech(x)] =', cp( (1, ('sech', cp('x'))) ).differentiated('x'))
print('d/dx[csch(x)] =', cp( (1, ('csch', cp('x'))) ).differentiated('x'))
print('d/dx[tanh(x)] =', cp( (1, ('tanh', cp('x'))) ).differentiated('x'))
print('d/dx[coth(x)] =', cp( (1, ('coth', cp('x'))) ).differentiated('x'))
print()
# and hyperbolic trig

print('d/dx[arcsinh(x)] =', cp( (1, ('arcsinh', cp('x'))) ).differentiated('x'))
print('d/dx[arccosh(x)] =', cp( (1, ('arccosh', cp('x'))) ).differentiated('x'))
print('d/dx[arcsech(x)] =', cp( (1, ('arcsech', cp('x'))) ).differentiated('x'))
print('d/dx[arccsch(x)] =', cp( (1, ('arccsch', cp('x'))) ).differentiated('x'))
print('d/dx[arctanh(x)] =', cp( (1, ('arctanh', cp('x'))) ).differentiated('x'))
print('d/dx[arccoth(x)] =', cp( (1, ('arccoth', cp('x'))) ).differentiated('x'))
print()
# and hyperbolic inverse trig


F = cp( (1, ('f', cp('x'))) )
print('F =', F)
F.differentiate('x')
print(F)
F.differentiate('x')
print(F)
F.differentiate('x')
print(F)
print()
# if the derivative of a certain function is unknown, it is written
# in a form D< f(x) >, which is equivalent to traditional f`(x).
# any higher derivatives are written as Dn< f(x) >, where
# n is the order of the derivative.

F_alt = cp(function('f', cp('x'), partials=[3]))
print('F =', F_alt)
print()
# alternatively, you can initialize derivatives of f of any order
# using the above shown syntax.

G = cp( (1, ('g',
             cp( (1, ('u', cp('x'))) ),
             cp( (1, ('v', cp('x'))) ),
             cp( (1, ('w', cp('x'))) ),
             )) )
print('G =', G)
print(G.differentiated('x'))
print()
# if the differentiated function depends on multiple variables,
# say G(u, v, w), partial of G with respect to v is written as D|0 1 0|< G(u, v, w) >.
# here the sequence 0 1 0 denoted that G was differentiated once with respect to
# argument 1 and zero times with respect to arguments 0 and 2.

H = cp( (1, ('ln', cp('y'))) )
print('H =', H)
print(H.differentiated('x'))
print(H.differentiated('y'))
print()
# if a certain function and all its` arguments do not have an explicit
# dependency on the variable it is being differentiated with respect to,
# the derivative then is zero.

I = cp( (1, ('ln', cp( (1, ('y', cp('x'))) ))) )
print('I =', I)
print(I.differentiated('x'))
print(I.differentiated('y'))
print()
# however, if y is declared not as an independent variable but as a function
# of x, the situation changes

J = cp( (1, ('J', cp('x'), cp('y'))) )
print('J =', J)
J.differentiate('x')
J.differentiate('x')
J.differentiate('y')
J.differentiate('y')
J.differentiate('y')
print(J)
print()
# that`s how higher order derivatives are displayed

J_alt = cp(function('J', cp('x'), cp('y'), partials=[2, 3]))
print('J = ', J_alt)
print()
# that`s another way of initializing J


print(J.wolfram_lang())
print(J.compileable())
print()
# derivatives of unknown functions can also be printed
# in wolfram language or in a compileable form.


input()
