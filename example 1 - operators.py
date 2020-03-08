from Algebraica import cp


A = cp(2, 3)
B = cp(13, 7)

print('A =', A)
print('A =', B)
print()
# just printing the initial numbers

print('-A =', -A)
print()
# multiplying number by -1

print('~A =', ~A)
print()
# producing a complex conjugate

print('A + B =', A + B)
print()
# adding two numbers together

print('A - B =', A - B)
print()
# subtracting two numbers from each other

print('A * B =', A * B)
print()
# multiplying two numbers

print('A * 2 =', A * 2)
print()
# multiplying complex by int

print('A * 1.6 =', A * 1.6)
print()
# multiplying complex by float

print('A * pi =', A * 'pi')
print()
# you can even multiply Complex by string!

print('A * stuff =', A * (['pi', [((1, 2), ('ln', cp(11))), ('sqrt', cp(2))], ('sqrt', cp(5))], [(1, 2), 'e']))
print()
# actually, you can write down any expression,
# just as if you were creating a class Complex object

print('A / B =', A / B)
print()
# dividing two numbers

print('A ** +3 =', A ** 3)
print('A ** -3 =', A ** -3)
print()
# raising to integer power

print('A ** B =', A ** B)
print()
# raising complex to complex power

print('Rl(A) =', A.real())
print('Im(A) =', A.imag())
print()
# separating real and imaginary parts

print('|A| =', abs(A))
print()
# finding the absolute value

print('arg(A) =', A.arg())
print()
# finding the complex argument

print('(A == B) =', A == B)
print()
# comparing two numbers

print('(A != B) =', A != B)
print()
# checking if the two numbers are different

print('(A.rl == B.rl) =', A.rl == B.rl)
print('(A.rl != B.rl) =', A.rl != B.rl)
print('(A.im == B.im) =', A.im == B.im)
print('(A.im != B.im) =', A.im != B.im)
print()
# comparing real/imaginary parts

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
