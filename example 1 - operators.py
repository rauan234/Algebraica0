from Algebraica import cp
# if creating complex numbers was all Algebraica was usable for, it would be pretty useless
# luckily, Algebraica allows to perform a variety of operations with complex numbers


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
# if A = a + b i, ~A = a + (-b) i

print('A + B =', A + B)
print()
# that`s how two class Complex objects are added together

print('A - B =', A - B)
print()
# that`s how two class Complex objects can be subtracted from each other

print('A * B =', A * B)
print()
# that`s multiplication of two complex number

print('A * 2 =', A * 2)
print('A / 2 =', A / 2)
print()
# complex numbers can be multiplied and divided by integers

print('A * 1.6 =', A * 1.6)
print('A / 1.6 =', A / 1.6)
print()
# and floats

print('A * pi =', A * 'pi')
print()
# you can even multiply Complex by string!
# the result of this operation is A * cp('pi').

print('A * stuff =', A * (['pi', [((1, 2), ('ln', cp(11))), ('sqrt', cp(2))], ('sqrt', cp(5))], [(1, 2), 'e']))
print()
# actually, you can write down any expression,
# just as if you were creating a class Complex object.
# in other words, A * (something) = A * cp(something)

print('A / B =', A / B)
print()
# class Complex objects can be divided by other class Complex objects

print('A ** +3 =', A ** 3)
print('A ** 0 =', A ** 0)
print('A ** -3 =', A ** -3)
print()
# they can be raised to integer powers

print('A ** B =', A ** B)
print()
# they can also be raised to a complex power

print('|A| =', abs(A))
print()
# this is the absolute value function

print('arg(A) =', A.arg())
print()
# this is how you find the complex argument

print('(A == B) =', A == B)
print()
# comparing two numbers

print('(A != B) =', A != B)
print()
# checking if the two numbers are different

C = A + B ** 7
D = (A / B) ** C
C -= D
E = cp(('exp', C))
print('C =', C)
print('D =', D)
print('E =', E)
print()
# doing operations with variables
# as you might have noticed, the result sometimes is a real mess
# so after getting the answer it`s useful to re-write it on paper and try to simplify
# or try auto-simplification described in example 4

input()
