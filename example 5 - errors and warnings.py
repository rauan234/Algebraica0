from Algebraica import cp, make_irrational
from Algebraica import Rational, Function, Term, Irrational, Complex
import Algebraica


# Algebraica lets you know if you did anything wrong
# if Algebraica finds whatever you did illegal, it raises an exception

try:
    Rational('h')
except Exception as exc:
    print(exc)

try:
    Rational(3, ())
except Exception as exc:
    print(exc)

try:
    Rational(1, 2, 3)
except Exception as exc:
    print(exc)
    
try:
    Function( [] )
except Exception as exc:
    print(exc)
    
try:
    Function( ('sqrt', 3) )
except Exception as exc:
    print(exc)

try:
    Function( (1, cp(3)) )
except Exception as exc:
    print(exc)
    
try:
    Function( ('') )
except Exception as exc:
    print(exc)
    
try:
    Term( [], {} )
except Exception as exc:
    print(exc)
    
try:
    Irrational('Cause I`ve been resting for this testing',
               'digesting',
               'every word that experts say')
except Exception as exc:
    print(exc)
    
try:
    Complex('Breathe', 'Breathe in the air')
except Exception as exc:
    print(exc)
    
# the most obvious way to initiate an exception is to use __init__ with invalid arguments


try:
    make_irrational(1, 2)
except Exception as exc:
    print(exc)
    
try:
    make_irrational([range])
except Exception as exc:
    print(exc)
    
try:
    make_irrational((1, 2, 3))
except Exception as exc:
    print(exc)
    
try:
    make_irrational((6, {'Capital of Kazakhstan': 'Astana'}))
except Exception as exc:
    print(exc)
    
try:
    make_irrational( ('ln', set) )
except Exception as exc:
    print(exc)

try:
    make_irrational( ((1, 1), ('sqrt', 2)) )
except Exception as exc:
    print(exc)

try:
    make_irrational( ((1, 1), (2, cp(2))) )
except Exception as exc:
    print(exc)
    
try:
    cp('I talk to the wind', 'My words followed, carried away', 'The wind does not hear', 'The wind cannot hear')
except Exception as exc:
    print(exc)
    
# cp and make_irrational functions must be fed a very specific kind of arguments


R = Rational(1, 2)

try:
    R == 'f'
except Exception as exc:
    print(exc)
    
try:
    R != [4, 2, 91]
except Exception as exc:
    print(exc)
    
try:
    R + range(3)
except Exception as exc:
    print(exc)
    
try:
    R - (3.5)
except Exception as exc:
    print(exc)
    
try:
    R * {'My favorite album': 'Pink Floyd: Wish you were here'}
except Exception as exc:
    print(exc)
    
try:
    R / int
except Exception as exc:
    print(exc)
    
try:
    R == list
except Exception as exc:
    print(exc)    
    
try:
    R != ()
except Exception as exc:
    print(exc)
    
try:
    R > Algebraica
except Exception as exc:
    print(exc)
    
try:
    R < sorted
except Exception as exc:
    print(exc)
    
try:
    R >= ['e']
except Exception as exc:
    print(exc)

try:
    R <= '\0'
except Exception as exc:
    print(exc)
    
# here is how illegal operations with class Rational objects cause exceptions


try:
    Function( ('pi') ) == 3.14
except Exception as exc:
    print(exc)
    
# comparing a class Function object to an object of any other class is illegal


T = Term(
    Rational(4, 11),
    [
        Function( ('sqrt', cp(2)) ),
        Function( ('ln', cp(5)) )
    ]
)

try:
    T * ['And the lamb', 'lies down', 'on Broadway']
    
except Exception as exc:
    print(exc)
    
try:
    T / float
except Exception as exc:
    print(exc)

try:
    T == open

except Exception as exc:
    print(exc)
    
# operators of class Term can only accept arguments of certain types


I = make_irrational((4, 7))

try:
    I / 'Mother, do you think they`ll drop the bomb'
except Exception as exc:
    print(exc)
    
try:
    I + ('Who, do you think, was a better composer?', 'Chopin', 'or', 'Liszt', '?')
except Exception as exc:
    print(exc)

try:
    I - {'Riders on the storm': 'like a dog without a bone',
         'Or an acter out on loan': 'riders on the storm'}
except Exception as exc:
    print(exc)
    
try:
    I * ('Come on baby light my fire', 'Try to set the night on fire')
except Exception as exc:
    print(exc)
    
try:
    I == ['My spirits are low in the depths of despair', 'My life blood', 'Spills over...']
except Exception:
    pass

# exceptions are raised when you try to perform illegal operations with class Irrational objects


try:
    cp(2) ** 't'
except Exception as exc:
    print(exc)

try:
    cp('I have become comfortably numb') / (3,)
except Exception as exc:
    print(exc)
    
try:
    cp(6) == 'And the lamb lies down on Broadway'
except Exception as exc:
    print(exc)
    
try:
    cp('This is a roundabout') - 'it`s gonna spin you out and out'
except Exception as exc:
    print(exc)
    
try:
    cp('The moments seemed lost in all the noise') + 'A snow storm a stimulating voice'
except Exception as exc:
    print(exc)

# operations involving Complex are picky when it comes to the type of arguments to give them


try:
    Rational(1, 0)
except Exception as exc:
    print(exc)
    
# creating a Rational with zero denominator is prohibited


try:
    Function( ('pi', cp(2)) )
except Exception as exc:
    print(exc)
    
try:
    Function( ('cosh', cp(6), cp(4)) )
except Exception as exc:
    print(exc)
    
try:
    Function( ('pow') )
except Exception as exc:
    print(exc)
    
# if you create a function and its` name and the number of arguments are inconsistent,
# an exception is thrown

# if you want to use pi as a function of one variable, you can change settings
# so that pi won`t be recognized as 3.1415... any longer
# look up example 6 to learn more
