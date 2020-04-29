from Algebraica import cp, make_irrational
from Algebraica import Rational, Function, Term, Irrational, Complex
# in this example, I will tell some details that didn`t quite fit to
# the last two examples but needed to be told.


A = cp( [ 1, (3, 7), ((1, 2), 'pi'), 0, ((9, 11), 'pi') ] )
print('A =', A)
print()
# notice the auto-simplification function that removes the zero-terms.
# also, note that terms are automatically grouped together, so that 1/2 pi + 9/11 pi becomes 29/22 pi.


B = Term(
    Rational(1, 1),
    [
        Function(( 'a', cp((2, 3)) )),
        Function(( 'd', cp(6, 2) )),
        Function(( 'b', cp(9, 4) )),
        Function(( 'c', cp( ((4, 7), 'pi'), 2 ) ))
    ]
)
print('B =', B)
print()
# in order for comparison operation to work, e.g. to make
# pi * sqrt(2) equal to sqrt(2) * pi,
# functions inside Terms are sorted in alphabet order


C = make_irrational([
    ( (1, 2), [
        'a',
        'e',
        'f'
        ] ),
    ( (4, 3), [
        'g',
        'b'
        ] ),
    ( (7, 4), [
        'h',
        'j',
        'c',
        'd'
    ] )
    ])
print('C =', C)
print()
# Terms inside Irrationals are sorted in alphabet order, too


D = Term(
    Rational(1, 1),
    [
        Function(( 'f', cp((1, 7)) )),
        Function(( 'f', cp((1, 7)) )),
        Function(( 'f', cp((1, 7)) )),
        Function(( 'f', cp((1, 7)) )),        
        Function(( 'g', cp((5, 2)) ))
    ]
)
print('D =', D)
print()
# identical functions are stacked together


E = make_irrational([
    ( (1, 2), [
        ( 'ln', cp(10) ),
        ( 'sqrt', cp(2) ),
        'pi'
    ]),
    ( (9, 4), [
        ( 'sqrt', cp(2) ),
        'pi',
        ( 'ln', cp(10) )
    ] ),
    ( (11, 13), [
        ( 'sqrt', cp(2) ),
        ( 'ln', cp(10) ),
        'pi'
    ] )    
])
print('E =', E)
print()
# similarly, terms with equal irrational parts are added together


F = Irrational()
print('F =', F)
print()
# if an Irrational has no terms, it is shown as 0


print('log3(7) =', cp( (1, ('log', cp(7), cp(3))) ).simplified())
print()
# log3(7) should be written as cp( (1, ('log', cp(7), cp(3))) ),
#                       not as cp( (1, ('log', cp(3), cp(7))) )


print('x ** 6 =', cp( (1, ('pow', cp('x'), cp(6))) ))
print('x ** -6 =', cp( (1, ('pow', cp('x'), cp(-6))) ))
print('x ** 1/2 =', cp( (1, ('pow', cp('x'), cp((1, 2)))) ))
print('x ** -1/2 =', cp( (1, ('pow', cp('x'), cp((-1, 2)))) ))
print('x ** -1 =', cp( (1, ('pow', cp('x'), cp(-1))) ))
print()
# power functions pow(base, pow) with pow being equal to -1, -1/2, 1/2 or being an
# integer are shown in a special way.


input()
