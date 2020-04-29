from Algebraica import cp
import Algebraica
# Algebraica library is tunable.
# in this example, we will explore this possibility.


# Algebraica automatically simplifies certain expressions
print('ln(e) =', cp( (1, ('ln', cp('e'))) ))

# if you don`t want this, turn AutoSimpify off
Algebraica.AutoSimplify = False

print('ln(e) =', cp( (1, ('ln', cp('e'))) ))
print()


# let`s say you want to create a variable that would tell
# how many kilograms of silver a bakn has in its` reserves
# and how would you call it? arg, of course!
# but there is a problem with that name

try:
    cp('arg')
except Exception as exc:
    print(exc)
print()

# namely, Algebraica recognizes arg as the complex argument function

# in order to fix that issue, we can change the name which
# Algebraica assigns to the argument function, say to __arg__
Algebraica.arg_function_name = '__arg__'

# now everything`s okay
try:
    cp('arg')
except Exception as exc:
    print(exc)
print()

# you can change names of all constants and functions in a similar way by using a special function
Algebraica.set_names_to_secure()

# in order to see what you have done, call the show_names() function
Algebraica.show_names()
print()

# now you can go back by usint set_names_to_default() function
Algebraica.set_names_to_default()
Algebraica.show_names()
print()

# there is no restrictions on functions` names
Algebraica.log_function_name = 'dog'
Algebraica.ln_function_name = 'In'

print(cp( (1, ('dog', cp(8), cp(2))) ).simplified())
print()
# everything still works like it did before

Algebraica.set_names_to_default()

Algebraica.sin_function_name = 'I_can`t_spare_a_buck_but_that`s'
Algebraica.pi_constant_name = 'a_rock'
print(cp( (1, ('I_can`t_spare_a_buck_but_that`s', cp( ((1, 2), 'a_rock'))))))
print(cp( (1, ('I_can`t_spare_a_buck_but_that`s', cp( ((1, 2), 'a_rock'))))).simplified())
print()
Algebraica.set_names_to_default()
# just one more example


# Also, you can adjust the indentation step
Algebraica.TabulationStep = 6
print(cp( (1, ('ln', cp(1, ((2, 3), ('sqrt', cp(3)))))) ))

Algebraica.TabulationStep = 2
print(cp( (1, ('ln', cp(1, ((2, 3), ('sqrt', cp(3)))))) ))

Algebraica.TabulationStep = 4
print(cp( (1, ('ln', cp(1, ((2, 3), ('sqrt', cp(3)))))) ))
print()


# if you, like me, prefer tau over pi, you can turn on the TauMode
Algebraica.TauMode = True
print(cp('pi').simplified())
print(cp( (1, ('arcsin', cp(1)))).simplified())
print()
# if TauMode is on, pi is expressed as tau/2

# if TauMode is off, tau is expressed as 2 pi
Algebraica.TauMode = False
print(cp('tau').simplified())
print()


# you can also change the precision with which class Numerical objects are displayed
print(Algebraica.Numerical(1/3, 0))
Algebraica.NumericalNumbersPrintLength = 7
print(Algebraica.Numerical(1/3, 0))
print()


# if you want to see stuff in standard blackwhite, switch ColorfulPrint to False
Algebraica.ColorfulPrint = False
print( cp( (1, ('ln', cp(2))) ))
print()


# if you want, you can turn on the use of Euler`s formula for complex exponentials
Algebraica.ExpandComplexExponentials = True
print('exp(i) =', cp( (1, ('exp', cp(0, 1))) ).simplified())
print()


# now that we have done so much mess it makes sense to restore order
Algebraica.return_to_default()

# this function returns all settings to what they were in the beginning


input()
