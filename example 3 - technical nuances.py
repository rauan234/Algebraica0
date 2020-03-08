from Algebraica import cp, make_irrational
from Algebraica import Rational, Function, Term, Irrational, Complex


# cp(*args) function relies completely on make_irrational(inp)
'''
def cp(*args):
    if(len(args) == 1):
        rl = args[0]
        
        return Complex(
            make_irrational(rl),
            make_irrational(0)
        )
        
    elif(len(args) == 2):
        rl = args[0]
        im = args[1]
        
        c = Complex(
            make_irrational(rl),
            make_irrational(im)
        )
        
    else:
        raise RuntimeError
    
    return c
'''
# as you can see, it creates a Complex by simply converting its` inputs
# into Irrational and then calling Complex.__init__


# make_irrational(inp) is a great pile of if-statements
# it was precariously designed to be able to sort out
# any trash you might want to throw in it, as long as it makes
# just a little bit of sense
# it takes care of 13 input formats occupies 100 lines of code
# the function itself is intimidatingly messy, so I decided not to
# include it in this example
# if you want, you can find the function in Algebraica.py


A = make_irrational( [ 1, (3, 7), ((1, 2), 'pi'), 0, ((9, 11), 'pi') ] )
print('A =', A)
print()
# note the auto-simplification function that removes the zero-terms from Irrationals
# also, terms are automatically grouped together, so that 1/2 pi + 9/11 pi becomes 29/22 pi


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
# the functions that are multiplied together to form a term are
# sorted in alphabet order

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
# likewise, terms are sorted in alphabet order

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

G = cp(0, 'f')
print('G =', G)
print('~G =', ~G)
print('arg(G) =', G.arg())
print('arg(~G) =', (~G).arg())
print()
# if the real part of the input of the arg function is zero,
# the result is pi/2 * sign( im )


input()
