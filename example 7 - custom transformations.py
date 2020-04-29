import Algebraica
from Algebraica import cp


# let`s say you want to simplify a certain expression, but
# it turns out that Algebraica`s innate simplification mechaism is not sufficient.
# if that`s the case, you might want to try to use custom transformations.

# transformations work exactly the same way as simplification.
# if, for example, the transformation changes cosh(x) to (exp(x) + exp(-x)) / 2,
# Algebraica simply substitutes (exp(x) + exp(-x)) / 2 everywhere it sees a cosh.

# more concretely, a transformation function (let`s call it tr(func)) takes in a
# function and gives out a class Function object and returns a class Complex object
# that is supposed to have the same value as the initial function.

# for example, this is how a transformation that changes cosh to exp looks.
def tr(func):
    if(func.name == 'cosh'):
        return (cp( ((1, 2), ('exp', func.args[0])) ) + 
                cp( ((1, 2), ('exp', -func.args[0])) ))
    
    else:
        return cp(func)
# if the function is a hyperbolic cosine, it substitutes the formula.
# otherwise, it simply returns what it took in.

# a class Complex object consists of two class Irrational objects, which, in
# turn, consist of multiple class Term objects. C = A + Bi, A = t1 + t2 + t3..., B = t4  + t5 + t6...
# a class Term object is a rational coefficient multiplied by several functions.
# for example, 3/11 * ln(2) * cosh(1) * sqrt(5) is a term. let`s call it T for future purposes.
# a transformation tr takes T and transforms it into a class Complex object N.
# by the definition of a transformation, N = 3/11 * tr( ln(2) ) * tr( cosh(0) ) * tr( sqrt(5) ).
# after computing, we get N = 3/22 * e * ln(2) * sqrt(5) + 3/22 * exp(-1) * ln(2) * sqrt(5).
# the same procedure is done with every term constituting a class Complex object.

# this is how to apply a transformation
C = cp( ((3, 11), [('ln', cp(2)), ('cosh', cp(1)), ('sqrt', cp(5))]) )
print('C =', C)
print(C.transformed(tr))
print()
# C.transformed(tr) returns a transformed version of C and
# does not modify the initial number.


def tr1(func):
    if((func.name == 'a') and (func.args[0] == cp('b'))):
        return cp('b')
    
    else:
        return cp(func)
D = cp( (1, ('a', cp('b'))) )
print('D =', D)
D.transform(tr1)
print(D)
print()
# D.transform(tr1) does not return anything but
# changes D to D.transformed(tr1).


E = cp( (1, ('a', cp( (1, ('a', cp( (1, ('a', cp('b'))) ))) ))) )
print('E =', E)
print(E.transformed(tr1))
print()
# transformation is a recursive process.
# Algebraica will continue applying the transformation
# until it does not yield any changes.


def forbidden_transformation(func):
    if(func.name == 'f'):
        return cp('g')
    elif(func.name == 'g'):
        return cp('f')
    
    else:
        return cp(self)
    
F = cp('f')
'''
# if run this code, it will enter an infinite cycle, because forbidden_transformation
# first changes f to g, then g to f, and then again to g, and this process never stops.
F.transform(forbidden_transformation)
'''
# if the transformation applied is self-contradictive,
# the program enters an infinite cycle.
# be careful with transformations.


input()
