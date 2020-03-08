from Algebraica import cp
# cp() is a function, that allows to create
# objects of class Algebraica.Complex


A = cp(1)
print('A =', A)
print()
# the above function makes A equal 1 + 0 i
# if you want to initialize a Complex with
# only the real part, simply use cp(R), and
# it will create a comlex number R + 0 i

B = cp(2, 3)
print('B =', B)
print()
# if you need a complex number that actually has
# an imaginary component, call cp(R, I) and
# you will get R + I i

C = cp(4.5, 7)
print('C =', C)
print()
# you can create complex numbers using floats

D = cp(
    (1, 2),
    (2, 3)
)
print('D =', D)
print()
# or, instead of floats, you might want to use integer fractions
# the avobe funcion creates a number 1/2 + 2/3 i

E = cp('pi')
print('E =', E)
print()
# you can create irrational numbers just like that
# note, that there is no restrictions on the names

F = cp('Muad`dib')
print('F =', F)
print()
# you can choose arbitrary names

G = cp(
    ((3, 2), 'pi'),
    (1, 5)
)
print('G =', G)
print()
# the irrational might also contain a fraction
# in order to do such thing, write the expression in a form
# cp( (N, D), name )  where N, D are the numerator and denominator
# and name is the name of your constant

H = cp(
    ((4, 9), ('ln', cp(7))),
    1
)
print('H =', H)
print()
# besides constants, you might want to add functions, such as ln
# for that, write the expression in a form cp( (N, D), (name, *args) )
# there can be multiple arguments, so don`t hesitate to add them if needed
# note that args must be objects of class Complex

I = cp(
    ((4, 5), ('log', cp('e'), cp(2))),
    ((9, 7), ('sqrt', cp(2)))
)
print('I =', I)
print()
# that`s one more example of how you can create a complex number

J = cp(
    ((2, 1), [('sqrt', cp(3))]),
    1
)
print('J =', J)
print()
# if you love square brackets, you can write your expressions like
# cp( (N, D), [ (name, *args) ] )
# but really, why bother with square brackets?

K = cp(
    ((1, 1), [
               'pi',
              ('sqrt', cp(0, 2)),
              ('sin', cp(((1, 2), ['pi'])))
              ]),
    (2, 9)
)
print('K =', K)
print()
# well, with them you can have multiple irrational numbers multiplied togethet
# for example, the above expression gives out a number
# 1/1 * pi * sqrt(2 i) * sin(1/2 * pi) + 2/9 i

L = cp(
    ['pi', ((1, 2), [('ln', cp(11)), ('sqrt', cp(2))]), ('sqrt', cp(5))],
    [(1, 2), 'e']
)
print('L =', L)
print()
# if you really like lists, you can create numbers with multiple terms, like this one
# to do it, use cp( [*real_terms], [*comp_terms] ), where real_terms and complex_terms
# are things that will be added together eventually
# each of there terms is to be declared just as terms used to be declared before, aka
# in a form cp( (N, D), [ (name1, *args1), (name2, *args2) ], ... ) or in any of the shorter ways

M = cp(
    [(2, 3), ((7, 5), 'pi'), ((4, 3), ('exp', cp(3))),
         ((1, 2), [
             'pi',
             ('exp', cp(4)),
             ('sqrt', cp((1, 1), [
                 1,
                 (3, 7),
                 ((9, 1), ('ln', cp(3))),
                 ((5, 2), ('pow', cp(1, 2), cp(3, 4)))
             ]))  # these are added
         ])  # these things are multiplied
     ],  # these things are added
    0
)
print('M =', M)
print()
# this is how complicated it can get

N = cp(
    ((7, 5), [('exp', M)])
)
print('N =', N)
print()
# if you want, you can use M to define N
# this way Complex can be 'embeded' into one another

input()
