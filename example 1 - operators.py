from Algebraica import cp
# if creating complex numbers was all Algebraica was usable for, it would be pretty useless.
# luckily, Algebraica allows to perform a variety of operations on complex numbers.


A = cp(2, 3)
B = cp(13, 7)
# first, here are two class Complex objects

print('A =', A)
print('B =', B)
print()
# here are the values of A and B

print('-A =', -A)
print()
# this operator turns a + bi to (-a) + (-b) i

print('~A =', ~A)
print()
# this one returns a complex conjugate of A.
# if A = a + b i, then ~A = a + (-b) i.

print('|A| =', abs(A))
print()
# this is the absolute value function

print('arg(A) =', A.arg())
print()
# this is how you find the complex argument


print('A + 2 =', A + 2)
print('A - 2 =', A - 2)
print('A * 2 =', A * 2)
print('A / 2 =', A / 2)
print()
# arithmerics can be done between Complex and ints

print('A + 1.6 =', A + 1.6)
print('A - 1.6 =', A - 1.6)
print('A * 1.6 =', A * 1.6)
print('A / 1.6 =', A / 1.6)
print()
# and floats

print('A ** +3 =', A ** 3)
print('A ** 0 =', A ** 0)
print('A ** -3 =', A ** -3)
print('A ** B =', A ** B)
print()
# that`s how class Complex objects can be raised to an integer power


print('(A == B) =', A == B)
print('(A != B) =', A != B)
print()
# class Complex objects can be compared for equality.
# there is no operators like (A > B), because it is not
# clear what it should return when given a pair of complex numbers.

print('A + B =', A + B)
print('A - B =', A - B)
print('A * B =', A * B)
print('A / B =', A / B)
print()
# here are the basic arithmetic operations involving two class Complex objects

print(A ** B)
print(B ** A)
print()
# class Complex objects can also be raised to a complex power

i = cp(0, 1)  # the imaginary unit
print(i ** i)
print()
# this, for example, is i ** i


print('A * pi =', A * 'pi')
print()
# you can even multiply Complex by string!
# the result of this operation is the same as A * cp('pi').

print('A * stuff =', A * (['pi', [((1, 2), ('ln', cp(11))), ('sqrt', cp(2))], ('sqrt', cp(5))], [(1, 2), 'e']))
print()
# actually, you can write down any expression,
# just as if you were creating a class Complex object.
# in other words, A * (something) = A * cp(something).


C = A + B ** 3
D = (A / B) ** C
C -= D
E = cp(('exp', C))
print('C =', C)
print('D =', D)
print('E =', E)
print()
# this is how you might want to use class Complex objects for your calculations.
# notice, that the answer is automatically simplified.
# view example 4 to learn more about simplification.

input()
