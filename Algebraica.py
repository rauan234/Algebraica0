import unicodedata
import colorama
import numpy
import math

colorama.init(convert=True)



global TabulationStep
TabulationStep = 4
# this refers to Complex.__str__ function

global TauMode
TauMode = False
# if true, tau = 2 pi will be used instead of pi

global AutoSimplify
AutoSimplify = True
# if true, all class Complex are simplified before each print
# and numerical value calculation

global ExpandComplexExponentials
ExpandComplexExponentials = False
# if true, expressions of type exp(f i), where f is a real number, will
# be re-written as cos(f) + i sin(f)

global ExpressArg
ExpressArg = True
# if true, arg(z) will be represented through arctan(z)
# if z is of form (a/b + c/d i), where a, b, c, and d are integers

global NumericalNumbersPrintLength
NumericalNumbersPrintLength = 12
# determines the resolution with which Numericals will be printed

global ColorfulPrint
ColorfulPrint = True
# if true, numbers will be printed in color

global e_constant_name
global pi_constant_name
global tau_constant_name
global sin_function_name
global cos_function_name
global sec_function_name
global csc_function_name
global tan_function_name
global cot_function_name
global sinh_function_name
global cosh_function_name
global sech_function_name
global csch_function_name
global tanh_function_name
global coth_function_name
global arcsin_function_name
global arccos_function_name
global arcsec_function_name
global arccsc_function_name
global arctan_function_name
global arccot_function_name
global arcsinh_function_name
global arccosh_function_name
global arcsech_function_name
global arccsch_function_name
global arctanh_function_name
global arccoth_function_name
global ln_function_name
global log_function_name
global pow_function_name
global sqrt_function_name
global exp_function_name
global rl_function_name
global im_function_name
global arg_function_name
global sign_function_name

e_constant_name =         'e'
pi_constant_name =        'pi'
tau_constant_name =       'tau'
sin_function_name =       'sin'
cos_function_name =       'cos'
csc_function_name =       'csc'
sec_function_name =       'sec'
tan_function_name =       'tan'
cot_function_name =       'cot'
sinh_function_name =      'sinh'
cosh_function_name =      'cosh'
sech_function_name =      'sech'
csch_function_name =      'csch'
tanh_function_name =      'tanh'
coth_function_name =      'coth'
arcsin_function_name =    'arcsin'
arccos_function_name =    'arccos'
arcsec_function_name =    'arcsec'
arccsc_function_name =    'arccsc'
arctan_function_name =    'arctan'
arccot_function_name =    'arccot'
arcsinh_function_name =   'arcsinh'
arccosh_function_name =   'arccosh'
arcsech_function_name =   'arcsech'
arccsch_function_name =   'arccsch'
arctanh_function_name =   'arctanh'
arccoth_function_name =   'arccoth'
ln_function_name =        'ln'
log_function_name =       'log'
pow_function_name =       'pow'
sqrt_function_name =      'sqrt'
exp_function_name =       'exp'
rl_function_name =        'Rl'
im_function_name =        'Im'
arg_function_name =       'arg'
abs_function_name =       'abs'
sign_function_name =      'sign'
# you can change names of the constants if you want



global InitialMaxPrime
InitialMaxPrime = 10000
# when decomposing numbers into primes, Algebraica first generates
# a list of primes using Eratosphenes` method.
# InitialMaxPrime determines the length of that list.

global PrimeNumberSearchStep
PrimeNumberSearchStep = 10000
# if the list of prime numbers is not sufficient, it is
# expanded when needed.
# PrimeNumberSearchStep determines by how much it will be expanded.

global HighestAllowedPrime
HighestAllowedPrime = 100000
# decomposing large numbers into primes can be a very hard task.
# that`s why there is a limit to how far Algebraica will go.
# namely, it won`t search for primes that are higher than HighestAllowedPrime.

global NumericalComparisonMargin
NumericalComparisonMargin = 10**-12
# if two class Numerical objects are different by less than NumericalComparisonMargin,
# they are considered equal



class Rational:
    def __init__(self, *args):  # object of class Rational represents a fraction
                                # it has an integer numerator and an integer denominator
        if(len(args) == 1):
            f = args[0]
            
            if isinstance(f, int):
                self.num = f
                self.den =  1
                
                self.reduce()
            
            elif isinstance(f, float):
                self.num, self.den = float_to_frac(f)
                
                self.reduce()
            
            else:
                raise TypeError('Wrong argument type: expected int or float, got ' + str(type(f)))
        
        elif(len(args) == 2):
            numerator = args[0]
            denominator = args[1]
            
            if (isinstance(numerator, int) and isinstance(denominator, int)):
                if(denominator != 0):
                    self.num = numerator
                    self.den = denominator
                    
                    self.reduce()
                
                else:
                    raise ValueError('Cannot create a Rational with a 0 denominator')
            
            else:
                raise TypeError('Wrong argument type: expected (int, int), got ' + str(type(numerator)) + ', ' + str(type(denominator)))
        
        else:
            raise ValueError('Wrong number of arguments: expected 1 or 2, got ' + str(len(args)))
        
    def __str__(self):
        if ColorfulPrint:
            if(self.den == 1):
                return colorama.Fore.WHITE + str(self.num) + colorama.Style.RESET_ALL
            
            elif(self.num == 0):
                return colorama.Fore.WHITE + '0' + colorama.Style.RESET_ALL
            
            else:
                return colorama.Fore.WHITE + str(self.num) + colorama.Style.NORMAL +\
                       '/' + colorama.Style.BRIGHT + str(self.den) + colorama.Style.RESET_ALL
        
        else:
            if(self.den == 1):
                return str(self.num)
            
            elif(self.num == 0):
                return '0'
            
            else:
                return str(self.num) + '/' + str(self.den)      
            
    def compileable(self):
        if(self == 1):
            return '1'
        
        else:
            return '(' + str(self.num) + ', ' + str(self.den) + ')'
        
    def __abs__(self):
        # abs( N / D ) = abs(N) / abs(D)
        return cp( (abs(self.num), abs(self.den)) )
    
    def __add__(first, second):
        if(isinstance(second, Rational)):
            return Rational(first.num * second.den + first.den * second.num,
                            first.den * second.den)
        
        elif(isinstance(second, int)):
            return Rational(
                first.num + second * first.den,
                first.den
            )
        
        elif(isinstance(second, float)):
            return first + Rational(second)
        
        else:
            raise TypeError('Wrong argument type: ' + str(type(second)))
    def __iadd__(first, second):
        return first + second
    
    def __mul__(first, second):
        if(isinstance(first, Rational) and
           isinstance(second, Rational)):
            
            return Rational(first.num * second.num,
                            first.den * second.den)
        
        elif(isinstance(first, Rational) and
             isinstance(second, int)):
            
            return Rational(first.num * second, first.den)
        
        elif(isinstance(first, Rational) and
             isinstance(second, float)):
            N, D = float_to_frac(second)
            
            return Rational(first.num * N, first.den * D)
        
        else:
            raise TypeError('Wrong argument type: ' + str(type(second)))
    def __imul__(first, second):
        return first * second
    
    def __pow__(first, second):
        if not isinstance(second, int):
            raise TypeError('Wrong argument type: ' + str(type(second)))
        
        if(second == 0):
            return Rational(1, 1)
        
        elif(second < 0):
            return Rational(
                first.den ** (-second),
                first.num ** (-second)
            )
        
        elif(second > 0):
            return Rational(
                first.num ** second,
                first.den ** second
            )
    
    def __mod__(first, second):
        if isinstance(second, int):
            return Rational(first.num % (first.den * second), first.den)
        
        else:
            raise TypeError('Wrong argument type: expected int, got ' + str(type(second)))
        
    def __sub__(first, second):
        if(isinstance(second, Rational)):
            return Rational(first.num * second.den - first.den * second.num,
                            first.den * second.den)
        
        elif(isinstance(second, int)):
            return Rational(
                first.num - second * first.den,
                first.den
            )
        
        elif(isinstance(second, float)):
            return first - Rational(second)
        
        else:
            raise TypeError('Wrong argument type: ' + str(type(second)))
    def __isub__(first, second):
        return first - second
    
    def __truediv__(first, second):
        if(isinstance(second, Rational)):
            return Rational(first.num * second.den,
                            first.den * second.num)
        
        elif isinstance(second, int):
            return Rational(first.num,
                            first.den * second)
        
        elif isinstance(second, float):
            N, D = float_to_frac(second)
            
            return Rational(first.num * D,
                            first.den * N)
        
        else:
            raise TypeError('Wrong argument type: ' + str(type(second)))
    def __idiv__(first, second):
        return first / second
    
    def __eq__(first, second):
        if isinstance(second, Rational):
            return ((first.num == second.num) and (first.den == second.den))
        
        elif isinstance(second, int):
            return ((first.num == first.den * second))
        
        elif isinstance(second, float):
            N, D = float_to_frac(second)
            return ((first.num * D == first.den * N))        
        
        else:
            raise TypeError('Wrong argument type: ' + str(type(second)))
    def __ne__(first, second):
        return not (first == second)
    
    def __gt__(first, second):
        if isinstance(second, Rational):
            return (first.num * second.den > second.num * first.den)
        
        elif isinstance(second, int):
            return (first.num > second * first.den)
        
        elif isinstance(second, float):
            N, D = float_to_frac(second)
            
            return (first.num * D > N * first.den)    
        
        else:
            raise TypeError('Wrong argument type: ' + str(type(second)))
        
    def __lt__(first, second):
        if isinstance(second, Rational):
            return (first.num * second.den < second.num * first.den)
        
        elif isinstance(second, int):
            return (first.num < second * first.den)
        
        elif isinstance(second, float):
            N, D = float_to_frac(second)
            
            return (first.num * D < N * first.den)    
        
        else:
            raise TypeError('Wrong argument type: ' + str(type(second)))
        
    def __ge__(first, second):
        if isinstance(second, Rational):
            return (first.num * second.den >= second.num * first.den)
        
        elif isinstance(second, int):
            return (first.num >= second * first.den)
        
        elif isinstance(second, float):
            N, D = float_to_frac(second)
            
            return (first.num * D >= N * first.den)    
        
        else:
            raise TypeError('Wrong argument type: ' + str(type(second)))
        
    def __le__(first, second):
        if isinstance(second, Rational):
            return (first.num * second.den <= second.num * first.den)
        
        elif isinstance(second, int):
            return (first.num <= second * first.den)
        
        elif isinstance(second, float):
            N, D = float_to_frac(second)
            
            return (first.num * D <= N * first.den)    
        
        else:
            raise TypeError('Wrong argument type: ' + str(type(second)))
    
    def __neg__(self):
        return Rational(-self.num, self.den)
    
    def __invert__(self):
        return Rational(self.den, self.num)
    
    def reduce(self):
        if(self.den < 0):
            self.den *= -1
            self.num *= -1
            
        if((self.den == 1) or (self.num == 1)):  # in this case the fraction
                                                     # cannot be reduced further
            return
        
        gcd = math.gcd(abs(self.num), abs(self.den))
        if(gcd != 1):
            self.num = self.num // gcd
            self.den = self.den // gcd
    
    def change_sign(self):
        self.num *= -1
        
    def get_sign(self):
        return sign(self.num) * sign(self.den)
        
    def copy(self):
        return Rational(self.num, self.den)
    
    def might_be_nonzero(self):
        return self.num != 0
    
    def proper_fraction(self):
        # turns a/b to n + c/d
        # where a, b, c, d, n are integers
        # and c/d < 1
        # for example, 5/3 = 1 + 2/3
        
        n = self.num // self.den
        if(n > 0):
            return n, self - n
        else:
            return 0, self
    
    def prime_factors(self):
        pnum = prime_factors(self.num)
        pden = prime_factors(self.den)
        
        primes = pnum
        for p in pden:
            primes.append( (p[0], -p[1]) )
        
        return primes
    
class Numerical:
    def __init__(self, new_rl, new_im):
        if( (isinstance(new_rl, float) or isinstance(new_rl, int)) and
            (isinstance(new_im, float) or isinstance(new_im, int)) ):
            self.rl = float(new_rl)
            self.im = float(new_im)
        
        else:
            raise TypeError('Wrong arguments type: expected float, got ' + str(type(new_rl)) + ' and ' + str(type(new_im)))
    
    def __str__(self):
        out = ''
        
        out += Str(self.rl, NumericalNumbersPrintLength)
        
        if(self.im > 0):
            if ColorfulPrint:
                out += colorama.Fore.WHITE + ' + ' + colorama.Style.RESET_ALL
            else:
                out += ' + '
                
            out += Str(self.im, NumericalNumbersPrintLength)
            
            if ColorfulPrint:
                out += colorama.Fore.MAGENTA + ' i' + colorama.Style.RESET_ALL
            else:
                out += ' i'
        
        elif(self.im < 0):
            if ColorfulPrint:
                out += colorama.Fore.WHITE + ' - ' + colorama.Style.RESET_ALL
            else:
                out += ' - '
                
            out += Str(-self.im, NumericalNumbersPrintLength)
            
            if ColorfulPrint:
                out += colorama.Fore.MAGENTA + ' i' + colorama.Style.RESET_ALL
            else:
                out += ' i'
        
        else:
            pass
        
        return out
    
    def __eq__(first, second):
        if isinstance(second, Numerical):
            return ((abs(first.rl - second.rl) < NumericalComparisonMargin) and
                    (abs(first.im - second.im) < NumericalComparisonMargin))
        
        elif(isinstance(second, int) or
             isinstance(second, float)):
            return ((abs(first.rl - second) < NumericalComparisonMargin) and
                    (abs(first.im) < NumericalComparisonMargin))
        
        else:
            raise TypeError('Can`t compare class Numerical object to ' + str(type(second)))
    
    def __ne__(first, second):
        return not (first == second)
    
    def __neg__(self):
        return Numerical(-self.rl, -self.im)
    
    def __invert__(self):
        return Numerical(self.rl, -self.im)
    
    def reciprocal(self):  # returns 1 / self
        abs_2 = self.rl**2 + self.im**2
        
        return Numerical( self.rl / abs_2, -self.im / abs_2 )
    
    def get_sign(self):
        return self / abs(self)
    
    def __add__(first, second):
        if isinstance(second, Numerical):
            return Numerical(first.rl + second.rl,
                             first.im + second.im)
        
        elif(isinstance(second, int) or isinstance(second, float)):
            return Numerical(
                first.rl + second,
                first.im
            )
        
        else:
            raise TypeError('Wrong argument type: ' + str(type(second)))
    
    def __sub__(first, second):
        if isinstance(second, Numerical):
            return Numerical(first.rl - second.rl,
                             first.im - second.im)
        
        elif(isinstance(second, int) or isinstance(second, float)):
            return Numerical(first.rl - second,
                             first.im)
        
        else:
            raise TypeError('Wrong argument type: ' + str(type(second)))
    
    def __mul__(first, second):
        if isinstance(second, Numerical):
            return Numerical(first.rl * second.rl - first.im * second.im,
                             first.rl * second.im + first.im * second.rl)
        
        elif(isinstance(second, int) or isinstance(second, float)):
            return Numerical(first.rl * second,
                             first.im * second)
        
        else:
            raise TypeError('Wrong argument type: ' + str(type(second)))
    
    def __imul__(first, second):
        return first * second
    
    def __truediv__(first, second):
        if isinstance(second, Numerical):
            abs_2 = second.rl ** 2 + second.im ** 2
            
            return Numerical( (first.rl * second.rl + first.im * second.im) / abs_2,
                              (-first.rl * second.im + first.im * second.rl) / abs_2 )
        
        elif(isinstance(second, int) or isinstance(second, float)):
            return Numerical(first.rl / second,
                             first.im / second)
        
        else:
            raise TypeError('Wrong argument type: ' + str(type(second)))
    
    def __pow__(first, second):
        if(isinstance(second, Numerical)):
            if( (first.rl > 0) and (first.im == 0) and
                (second.im == 0) ):
                return Numerical(numpy.power(first.rl, second.rl), 0.0)
            
            else:
                log_abs = numpy.log(first.rl**2 + first.im**2) / 2
                arg = first.arg()
                
                r = numpy.exp(log_abs * second.rl - arg * second.im)
                f = log_abs * second.im + arg * second.rl
                return Numerical(
                    r * numpy.cos(f),
                    r * numpy.sin(f)
                )
        
        elif(isinstance(second, float) or isinstance(second, int)):
            if( (first.rl > 0) and (first.im == 0) ):
                return Numerical(numpy.power(first.rl, second), 0.0)
            
            else:
                log_abs = numpy.log(first.rl**2 + first.im**2) / 2
                arg = first.arg()
                
                r = numpy.exp(log_abs * second)
                f = arg * second
                return Numerical(
                    r * numpy.cos(f),
                    r * numpy.sin(f)
                )
        
        else:
            raise TypeError('Wrong arguments type: ' + str(type(first)) + ', ' + str(type(second)))
    
    def __abs__(self):
        return numpy.sqrt(self.rl ** 2 + self.im ** 2)
    
    def arg(self):
        if(self.rl > 0):
            return numpy.arctan( self.im / self.rl )
        
        elif(self.rl < 0):
            if(self.im >= 0):
                return numpy.pi - numpy.arctan( abs(self.im / self.rl) )
            
            else:
                return -numpy.pi + numpy.arctan( abs(self.im / self.rl) )
            
        else:
            if(self.im > 0):
                return numpy.pi / 2
            
            elif(self.im < 0):
                return -numpy.pi / 2
            
            else:
                raise ValueError('arg(0) is not defined')
            
class UnknownFunctionError(Exception):
    def __init__(self, f):
        self.function = f
        
    def __str__(self):
        return 'Can`t compute the following function: ' + str(cp(self.function))
    
    def get_function(self):
        return self.function
        
class Function:
    def __init__(self, params):  # class Function item describes a certain function
                                 # it has the information on function`s name
                                 # and what arguments it has
        params = self.check_params(params)
        
        self.name = params[0]
        self.args = copy_list(params[1:])
        
        self.check_number_of_args()
        self.check_values_of_args()
        
    def make_up_a_name(reduced_name, partials):
        out = ''
        
        for p in partials:
            out += str(p) + '\0'
        
        out += reduced_name
        return out
        
    def check_number_of_args(self):
        n_args = self.number_of_arguments_required()
        if(n_args == None):  # if the function name is unknown
            pass
        else:
            if (n_args != len(self.args)):  # if the number of arguments given is incorrect
                if(n_args == 0):
                    raise TypeError('<' + self.name + '> function must have zero arguments, but the number of arguments given is ' + str(len(self.args)))
                elif(n_args == 1):
                    raise TypeError('<' + self.name + '> function must have one argument, but the number of arguments given is ' + str(len(self.args)))
                else:
                    raise TypeError('<' + self.name + '> function must have ' + str(n_args) + ' arguments, but the number of arguments given given is ' + str(len(self.args)))        
                
    def check_name(name):
        reduced_name = name.split('\0')[-1]
        
        if not isinstance(reduced_name, str):  # name is supposed to tell the name of the function
                                       # so it must be a string
            raise TypeError('Function name must be a string')
        
        if (len(reduced_name) == 0):
            raise TypeError('Function name cannot be an empty string')
        
        forbidden_characters = '\0\r\n ()[]{}<>.,^*-=+/\'\"'
        for char in reduced_name:
            if(char in forbidden_characters):
                try:
                    name = unicodedata.name(char)
                except Exception:
                    name = 'chr(' + str(ord(char)) + ')'
                    
                raise TypeError('Function name cannot contain the following character: ' + name)
        
        numbers = '0123456789'
        if(reduced_name[0] in numbers):
            raise TypeError('Function name cannot start with a digit')
        
    def check_params(self, params):
        if isinstance(params, str):        # sometimes we want to make a function constant, like pi
                                           # in that case, since there is no arguments, the params
                                           # should look like (pi_constant_name,)
                                           # but it is convenient to add a shortcut, so that
                                           # params could be just pi_constant_name
            params = (params,)
        elif isinstance(params, tuple):
            pass
        
        else:
            raise TypeError('Wrong argument type: expected str or tuple, got ' + str(type(params)))
        
        
        for item in params[1:]:
            if not isinstance(item, Complex):  # all of the function`s arguments
                                               # must be Complex
                raise TypeError('Arguments of a function must be Complex')
            
        Function.check_name(params[0])
        
        return params
        
    def check_values_of_args(self):
        not_defined_for_zero_arg = (
            ln_function_name,
            arg_function_name,
            csc_function_name,
            cot_function_name,
            csch_function_name,
            coth_function_name,
            arcsec_function_name,
            arccsc_function_name,
            arcsech_function_name,
            arccsch_function_name
        )
        if(self.name in not_defined_for_zero_arg):
            if(self.args[0].is_zero()):
                raise ValueError(self.name + '(0) is not defined')
        
        not_defined_for_half_pi = (
            sec_function_name,
            tan_function_name
        )
        if(self.name in not_defined_for_half_pi):
            if(self.args[0] == half_pi_constant):
                raise ValueError(self.name + '(pi/2) is not defined')
        
        not_defined_for_one = (
            arctanh_function_name,
        )
        if(self.name in not_defined_for_one):
            if(self.args[0] == cp(1)):
                raise ValueError(self.name + '(1) is not defined')
        
    def number_of_arguments_required(self):
        name = self.name
        args = self.args
        
        if(name in (
            pi_constant_name,
            e_constant_name
            )):
            if(len(args) != 0):
                return 0
            
        elif(name in (
            sin_function_name, cos_function_name, tan_function_name,
            sec_function_name, csc_function_name, cot_function_name,
            tanh_function_name, sinh_function_name, cosh_function_name,
            coth_function_name, sech_function_name, csch_function_name,
            arcsin_function_name, arccos_function_name, arctan_function_name,
            arcsec_function_name, arccsc_function_name, arccot_function_name,
            arcsinh_function_name, arccosh_function_name, arctanh_function_name,
            arcsech_function_name, arccsch_function_name, arccoth_function_name,
            sqrt_function_name,
            ln_function_name,
            exp_function_name,
            rl_function_name, im_function_name, arg_function_name, abs_function_name,
            sign_function_name
            )):
            return 1
            
        elif(name in (
            pow_function_name,
            log_function_name
            )):
            if(len(args) != 2):
                return 2
        
        return None
    
    def str_partials(partials):
        out = ''
        
        if ColorfulPrint:
            out += colorama.Style.NORMAL + colorama.Fore.GREEN + 'D' + colorama.Style.RESET_ALL
        else:
            out += 'D'
        if(len(partials) == 1):
            if(partials[0] != 1):
                out += str(partials[0])
        
        else:
            if ColorfulPrint:
                out += colorama.Style.DIM + colorama.Fore.GREEN + '|' + colorama.Style.RESET_ALL
            else:
                out += '|'
                
            if ColorfulPrint:
                out += colorama.Fore.WHITE
            for i in range(len(partials)):
                out += str(partials[i])
                
                if(i != len(partials) - 1):
                    out += ' '
            if ColorfulPrint:
                out += colorama.Style.RESET_ALL
                
            if ColorfulPrint:
                out += colorama.Style.DIM + colorama.Fore.GREEN + '|' + colorama.Style.RESET_ALL
            else:
                out += '|'
        
        if ColorfulPrint:
            out += colorama.Fore.BLACK + '< ' + colorama.Style.RESET_ALL
        else:
            out += '< '
        
        return out
    
    def str_body_for_sqrt(args, elevation):
        out = ''
        tab = ' ' * elevation * TabulationStep
        
        if ColorfulPrint:
            out += colorama.Fore.MAGENTA + sqrt_function_name + colorama.Style.RESET_ALL
        
        else:
            out += sqrt_function_name
        
        if(args[0].is_constant()):
            if ColorfulPrint:
                out += colorama.Style.NORMAL + colorama.Fore.MAGENTA + '[' + colorama.Style.RESET_ALL
                
            else:
                out += '['   
                
            out += ' '
            out += args[0].show(0, one_line=True)
            out += ' '
            
            if ColorfulPrint:
                out += colorama.Style.NORMAL + colorama.Fore.MAGENTA + ']' + colorama.Style.RESET_ALL
                
            else:
                out += ']'                   
            
        else:
            if ColorfulPrint:
                out += colorama.Style.NORMAL + colorama.Fore.MAGENTA + '[' + colorama.Style.RESET_ALL
                
            else:
                out += '['   
                
            out += '\n'
            
            out += args[0].show(elevation + 2)
            out += '\n'
            
            out += tab + '    '
            if ColorfulPrint:
                out += colorama.Style.NORMAL + colorama.Fore.MAGENTA + ']' + colorama.Style.RESET_ALL
                
            else:
                out += ']'
        
        return out
    
    def str_body_for_inverse_sqrt(args, elevation):
        out = ''
        tab = ' ' * elevation * TabulationStep
        
        if ColorfulPrint:
            out += colorama.Fore.WHITE + '1' + colorama.Fore.BLACK + ' / ' +\
                colorama.Fore.MAGENTA + sqrt_function_name + colorama.Style.RESET_ALL
        
        else:
            out += '1 / ' + sqrt_function_name
        
        if(args[0].is_constant()):
            if ColorfulPrint:
                out += colorama.Style.NORMAL + colorama.Fore.MAGENTA + '[' + colorama.Style.RESET_ALL
                
            else:
                out += '['                   
            
            out += ' '
            out += args[0].show(0, one_line=True)
            out += ' '
            
            if ColorfulPrint:
                out += colorama.Style.NORMAL + colorama.Fore.MAGENTA + ']' + colorama.Style.RESET_ALL
                
            else:
                out += ']'                   
            
        else:
            if ColorfulPrint:
                out += colorama.Style.NORMAL + colorama.Fore.MAGENTA + '[' + colorama.Style.RESET_ALL
                
            else:
                out += '['
            out += '\n'
            
            out += args[0].show(elevation + 2)
            out += '\n'
            
            out += tab + '    '
            if ColorfulPrint:
                out += colorama.Style.NORMAL + colorama.Fore.MAGENTA + ']' + colorama.Style.RESET_ALL
                
            else:
                out += ']'
        
        return out
    
    def str_body_for_reciprocal(args, elevation):
        out = ''
        tab = ' ' * elevation * TabulationStep
        
        if ColorfulPrint:
            out += colorama.Fore.WHITE + '1' + colorama.Fore.BLACK + ' / ' + colorama.Style.RESET_ALL
        
        else:
            out += '1 / '
        
        if(args[0].is_constant()):
            if ColorfulPrint:
                out += colorama.Fore.BLACK + '( ' + colorama.Style.RESET_ALL
            else:
                out += '( '
                
            out += args[0].show(0, one_line=True)
            
            if ColorfulPrint:
                out += colorama.Fore.BLACK + ' )' + colorama.Style.RESET_ALL
            else:
                out += ' )'
            
        else:
            if ColorfulPrint:
                out += colorama.Fore.BLACK + '(\n' + colorama.Style.RESET_ALL
            else:
                out += '(\n'
            
            out += args[0].show(elevation + 2)
            out += '\n'
            
            out += tab + '    '
            if ColorfulPrint:
                out += colorama.Fore.BLACK + ')' + colorama.Style.RESET_ALL
            else:
                out += ')'
        
        return out
    
    def str_body_for_pos_int_pow(args, elevation, p):
        out = ''
        tab = ' ' * elevation * TabulationStep
        
        if(args[0].is_constant()):
            if ColorfulPrint:
                out += colorama.Fore.BLACK + '( ' + colorama.Style.RESET_ALL
            else:
                out += '( '
                
            out += args[0].show(0, one_line=True)
            
            if ColorfulPrint:
                out += colorama.Fore.BLACK + ' ) ^ ' + colorama.Style.RESET_ALL
            else:
                out += ' ) ^ '
            
            if ColorfulPrint:
                out += colorama.Fore.WHITE + str(p) + colorama.Style.RESET_ALL
            else:
                out += str(p)
            
        else:
            out += '(\n'
            
            out += args[0].show(elevation + 2)
            out += '\n'
            
            out += tab + '    ) ^ ' + str(p)
        
        return out
    
    def str_body_for_neg_int_pow(args, elevation, p):
        out = ''
        tab = ' ' * elevation * TabulationStep
        
        if(args[0].is_constant()):
            if ColorfulPrint:
                out += colorama.Fore.WHITE + '1' + colorama.Fore.BLACK + ' / ( ' + colorama.Style.RESET_ALL
            else:
                out += '1 / ( '
                
            out += args[0].show(0, one_line=True)
            
            if ColorfulPrint:
                out += colorama.Fore.BLACK + ' ) ^ ' + colorama.Style.RESET_ALL
            else:
                out += ' ) ^ '
                
            if ColorfulPrint:
                out += colorama.Fore.WHITE + str(-p) + colorama.Style.RESET_ALL
            else:
                out += str(-p)
            
        else:
            out += '1 / (\n'
            
            out += args[0].show(elevation + 2)
            out += '\n'
            
            out += tab + '    ) ^ ' + str(-p)
        
        return out
    
    def str_body_for_general_function(args, reduced_name, elevation):
        out = ''
        tab = ' ' * elevation * TabulationStep
        
        if ColorfulPrint:
            out += colorama.Fore.MAGENTA + str(reduced_name) + colorama.Style.RESET_ALL
        
        else:
            out += str(reduced_name)
        
        if(len(args) == 0):
            if ColorfulPrint:
                return colorama.Fore.CYAN + str(reduced_name) + colorama.Style.RESET_ALL
        
            else:
                return str(reduced_name)
        
        elif((len(args) == 1) and args[0].is_constant()):
            if ColorfulPrint:
                out += colorama.Style.NORMAL + colorama.Fore.MAGENTA + '[' + colorama.Style.RESET_ALL
                
            else:
                out += '['                
            
            out += ' '
            out += args[0].show(0, one_line=True)
            out += ' '
            
            if ColorfulPrint:
                out += colorama.Style.NORMAL + colorama.Fore.MAGENTA + ']' + colorama.Style.RESET_ALL
                
            else:
                out += ']'   
            
        else:
            if ColorfulPrint:
                out += colorama.Style.NORMAL + colorama.Fore.MAGENTA + '[' + colorama.Style.RESET_ALL + '\n'
                
            else:
                out += '[\n'
            
            for i in range(len(args)):
                out += args[i].show(elevation + 2)
                
                if(i != len(args) - 1):
                    if ColorfulPrint:
                        out += colorama.Fore.BLACK + ', \n' + colorama.Style.RESET_ALL
                    else:
                        out += ', \n'
                        
                out += '\n'
                
            out += tab + '    '
            if ColorfulPrint:
                out += colorama.Style.NORMAL + colorama.Fore.MAGENTA + ']' + colorama.Style.RESET_ALL
                
            else:
                out += ']'
        
        return out
    
    def decode_name(self):
        if('\0' in self.name):
            lst = self.name.split('\0')
            
            reduced_name = lst[-1]
            partials = list(map(int, lst[:-1]))
                
            return reduced_name, partials
        
        else:
            return self.name, [0] * len(self.args)    
        
    def serial_index(f):
        out = f.name
        
        for arg in f.args:
            out += str(arg)
        
        return out
        
    def __eq__(first, second):
        if not isinstance(second, Function):
            raise TypeError('Wrong Function.__eq__ argument type: expected Function, got ' + str(type(second)))
        
        if(first.name != second.name):
            return False
        if(len(first.args) != len(second.args)):
            return False
        
        for i in range(len(first.args)):
            if(first.args[i] != second.args[i]):
                return False
        
        return True
    def __ne__(self, other):
        return not (self == other)
        
    def show(self, elevation):
        out = ''
        
        reduced_name, partials = self.decode_name()
        print_d_op = (True in map(lambda c: c > 0, partials))
        
        if print_d_op:
            out += Function.str_partials(partials)
        
        if((reduced_name == pow_function_name) and (self.args[1] == cp((1, 2)))):
            out += Function.str_body_for_sqrt(self.args, elevation)
            
        elif((reduced_name == pow_function_name) and (self.args[1] == cp((-1, 2)))):
            out += Function.str_body_for_inverse_sqrt(self.args, elevation)
                
        elif((reduced_name == pow_function_name) and (self.args[1] == cp(-1))):
            out += Function.str_body_for_reciprocal(self.args, elevation)
                
        elif((reduced_name == pow_function_name) and self.args[1].is_real_int()):
            if(len(self.args[1].rl.terms) == 0):
                p = 0
            
            else:
                p = self.args[1].rl.terms[0].rat.num
            
            if(p > 0):
                out += Function.str_body_for_pos_int_pow(self.args, elevation,  p)
                
            else:
                out += Function.str_body_for_neg_int_pow(self.args, elevation, p)
        
        else:
            out += Function.str_body_for_general_function(self.args, reduced_name, elevation)
                
        if(print_d_op):
            if ColorfulPrint:
                out += colorama.Fore.BLACK + ' >' + colorama.Style.RESET_ALL
            else:
                out += ' >'
        
        return out
    
    def __str__(self):
        return self.show(0)
    
    def wolfram_lang(self):
        alg_to_wolf_dict = {
            e_constant_name: 'E',
            pi_constant_name: 'Pi',
            tau_constant_name: '2 * Pi',
            sin_function_name: 'Sin',
            cos_function_name: 'Cos',
            sec_function_name: 'Sec',
            csc_function_name: 'Csc',
            tan_function_name: 'Tan',
            cot_function_name: 'Cot',
            sinh_function_name: 'Sinh',
            cosh_function_name: 'Cosh',
            sech_function_name: 'Sech',
            csch_function_name: 'Csch',
            tanh_function_name: 'Tanh',
            coth_function_name: 'Coth',
            arcsin_function_name: 'ArcSin',
            arccos_function_name: 'ArcCos',
            arcsec_function_name: 'ArcSec',
            arccsc_function_name: 'ArcCsc',
            arctan_function_name: 'ArcTan',
            arccot_function_name: 'ArcCot',
            arcsinh_function_name: 'ArcSinh',
            arccosh_function_name: 'ArcCosh',
            arcsech_function_name: 'ArcSech',
            arccsch_function_name: 'ArcCsch',
            arctanh_function_name: 'ArcTanh',
            arccoth_function_name: 'ArcCoth',
            sqrt_function_name: 'Sqrt',
            ln_function_name: 'Log',
            exp_function_name: 'Exp',
            rl_function_name: 'Re',
            im_function_name: 'Im',
            arg_function_name: 'Arg',
            abs_function_name: 'Abs'
        }
        
        out = ''
        
        reduced_name, partials = self.decode_name()
            
        if(self.name == pow_function_name):
            out += '(' + self.args[0].wolfram_lang() + ') ^' +\
                   '(' + self.args[1].wolfram_lang() + ')'
        
        elif(self.name == log_function_name):
            out += 'Log[' + self.args[1].wolfram_lang() + ', ' +\
                            self.args[0].wolfram_lang() + ']'
        
        else:
            includes_partials = False
            for p in partials:
                if(p > 0):
                    includes_partials = True
                    break
            
            if includes_partials:
                out += 'Derivative['
                for i in range(len(partials)):
                    out += str(partials[i])
                    
                    if(i != len(partials) - 1):
                        out += ','
                out += ']'
                
                try:
                    out += '[' + alg_to_wolf_dict[reduced_name] + ']'
                
                except KeyError:
                    out += '[' + reduced_name + ']'
            
            else:
                try:
                    out += alg_to_wolf_dict[reduced_name]
                
                except KeyError:
                    out += reduced_name
            
            if(len(self.args) > 0):
                out += '['
                for i in range(len(self.args)):
                    out += self.args[i].wolfram_lang()
                    
                    if(i != len(self.args) - 1):
                        out += ', '
                out += ']'
        
        return out
    
    def compileable(self):
        out = '('
        
        out += "\'" + self.name.replace('\0', '\\r') + "\'"
        
        for arg in self.args:
            out += ', '
            
            out += arg.compileable()
        
        out += ')'
        return out
    
    def af_pow(self):
        base = self.args[0]
        power = self.args[1]
        
        if base.is_zero():
            return cp(0)
        if((base == cp(1)) or (power == cp(0))):
            return cp(1)
        if(power == cp(1)):
            return base
        
        if(base.is_real_int() and power.is_real_int()):
            n = base.rl.terms[0].rat.num
            p = power.rl.terms[0].rat.num
            
            if(p > 0):
                return cp(n ** p)
            else:
                return cp( (1, n**(-p)) )
                
        if((len(base.rl.terms) == 1) and (len(base.im.terms) == 0)):
            if((len(power.rl.terms) == 1) and (len(power.im.terms) == 0)):
                if(len(power.rl.terms[0].irt) == 0):
                    def siutable(base):
                        if(len(base.rl.terms[0].irt) == 0):
                            return True
                        
                        else:
                            sign = cp(1)
                            for func in base.rl.terms[0].irt:
                                s = func.get_sign()
                                
                                if isinstance(s, type(None)):
                                    return False
                                
                                else:
                                    sign *= s
                            
                            return ( (sign == cp(1)) or (base.rl.terms[0].rat > 0) )
                        
                    if siutable(base):
                        base_num_primes = prime_factors(base.rl.terms[0].rat.num)
                        base_den_primes = prime_factors(base.rl.terms[0].rat.den)
                        # pow( 2 * 3 * pi, power ) = pow(2, power) * pow(3, power) * pow(pi, power)
                        
                        prat = power.rl.terms[0].rat
                        
                        primes = []
                        for n in base_num_primes:
                            primes.append( (n[0], Rational( n[1], 1) * prat) )
                        for d in base_den_primes:
                            primes.append( (d[0], Rational(-d[1], 1) * prat) )
                    
                        update = False
                        out = cp(1)
                        for fun in base.rl.terms[0].irt:
                            out *= cp( (1, (pow_function_name, cp(fun), power)) )
                        for prime in primes:
                            int_part, frac_part = prime[1].proper_fraction()
                            if(int_part != 0):
                                update = True
                            
                            out *= cp( ((prime[0] ** int_part, 1), (pow_function_name, cp(prime[0]), cp(frac_part))) )
                        
                        if update:
                            return out
                    
        if((len(base.rl.terms) == 1) and (len(base.im.terms) == 0)):  # if base is single-term
            if(len(base.rl.terms[0].irt) == 1):  # if the term is single-function
                if(base.rl.terms[0].irt[0].name == pow_function_name):
                    # pow( pow(a, b), c ) = pow(a, b * c)
                    # (this is, of course, not always true)
                    # for example, ( (-1) ** 2 ) ** (1/2) != (-1) ** 1
                    
                    rat = base.rl.terms[0].rat
                    n = base.rl.terms[0].irt[0].args[0]
                    p1 = base.rl.terms[0].irt[0].args[1]
                    p2 = power                        
                    if( (n.is_real_and_positive() and (p1.is_always_real() and p2.is_always_real())) or
                        (n.is_always_real() and (p1.is_real_int() and p2.is_real_int())) ):
                        return cp( (1, (pow_function_name, n, p1 * p2)) ) * cp( (1, (pow_function_name, cp(rat), p2)) )
                
                elif(base.rl.terms[0].irt[0].name == e_constant_name):
                    return cp( (1, (exp_function_name, power * base.rl.terms[0].rat)) )
                
                elif(base.rl.terms[0].irt[0].name == exp_function_name):
                    if power.is_real_int():
                        # pow(exp(a), b) = exp(a * b)
                        # if b is real and integer.
                        # if b = 1/2, for example, that`s not true.
                        # pow(exp(2 pi i), 1/2) = pow(1, 1/2) = 1
                        # exp(2 pi i * 1/2) = exp(pi i) = -1
                        return cp( (1, (exp_function_name, power * base.rl.terms[0].irt[0].args[0])) )                    
                        
            elif(len(base.rl.terms[0].irt) > 1):  # if the term contains several functions,
                                                  # such as 2 * pi * ln(3)
                if( (base.rl.terms[0].rat > 0) and power.is_always_real()):
                    siutable = True
                    for func in base.rl.terms[0].irt:
                        if not func.is_real_and_positive():
                            siutable = False
                            break
                        
                    if siutable:
                        term = base.rl.terms[0]
                        
                        out = cp( (1, (pow_function_name, cp( (term.rat.num, term.rat.den) ), power)) )
                        for f in term.irt:
                            out *= cp( (1, (pow_function_name, cp(f), power)) )
                        out.simplify()
                        
                        return out
        
        if((not base.has_functions()) and (power.is_real_int())):
            p = power.rl.terms[0].rat.num
            
            if(p > 0):
                return base ** p
            
            else:
                '''inv = ~base
                
                # base ** p = 1 / (base ** -p) = (~base) ** (-p) / ( (base * ~base) ** (-p) ) =
                # in order to calculate base ** p, we can multiply both the top
                # and bottom of the fraction by (~base) ** (-p)
                return ( inv ** (-p) ) / ( (base * inv) ** (-p) )'''
                
                inv = ~base
                
                if(len(base.rl.terms) == 0):
                    abs_2 = base.im.terms[0].rat ** 2
                    
                elif(len(base.im.terms) == 0):
                    abs_2 = base.rl.terms[0].rat ** 2
                    
                else:
                    abs_2 = base.rl.terms[0].rat ** 2 + base.im.terms[0].rat ** 2 # abs(base) ** 2
                
                return ( inv ** (-p) ) / ( abs_2 ** (-p) )
        
        if(len(base.rl.terms) + len(base.im.terms) > 1):
            rats = []
            for term in base.rl.terms:
                rats.append(term.rat)
            for term in base.im.terms:
                rats.append(term.rat)
            
            com = common_divider(rats)
            if( (com != Rational(1, 1)) and (com > 0) ):
                return cp( (1, (pow_function_name, cp(com), power)) ) * cp( (1, (pow_function_name, base / com, power)) )
        
        return cp(self)
    
    def af_exp(self):
        P = self.args[0]
        
        if(P == cp(0)):
            return cp(1)
        if(P == cp(1)):
            return e_constant_name
        
        if(len(P.rl.terms) >= 1):
            # this functions transforms expressions of type
            # exp(power + ln(num)) to exp(power) * num
            num = cp(1)
            power = cp(0, 1) * cp(P.im)
            update = False
            
            for term in P.rl.terms:
                if(len(term.irt) == 1):
                    if(term.irt[0].name == ln_function_name):
                        # exp( ln(a) ) = a
                        
                        if(term.rat == 1):
                            num *= term.irt[0].args[0]
                        else:
                            num *= cp( (1, (pow_function_name, term.irt[0].args[0],
                                                  cp(term.rat))) )
                        update = True
                    
                    else:
                        power += cp(term)
                else:
                    power += cp(term)
                        
            if update:
                return num * cp( (1, (exp_function_name, power)) )
        
        if(P.rl.is_always_real() and P.im.is_always_real()):
                if ExpandComplexExponentials:
                    # exp(f i) = cos(f) + i sin(f)
                    
                    if P.im.might_be_nonzero():
                        return cp( (1, [(cos_function_name, cp(P.im)), (exp_function_name, cp(P.rl))]),
                                   (1, [(sin_function_name, cp(P.im)), (exp_function_name, cp(P.rl))]) )
        
        return cp(self)
    
    def af_ln(self):
        N = self.args[0]
        
        if(N.is_zero()):
            raise ValueError('ln(0) is not defined')
        
        if(N == cp(1)):
            return cp(0)
        
        elif(N == cp(-1)):
            if TauMode:
                return cp( 0, ((1, 2), tau_constant_name) )
            else:
                return cp( 0, (1, pi_constant_name) )
        
        elif(N == cp(0, 1)):
            if TauMode:
                return cp( 0, ((1, 4), tau_constant_name) )
            else:
                return cp( 0, ((1, 2), pi_constant_name) )                
            
        elif(N == e_constant):
            return cp(1)
            
        if((len(N.rl.terms) == 1) and (len(N.im.terms) == 0)):  # if N is single-term
            term = N.rl.terms[0]
            
            if(len(N.rl.terms[0].irt) == 0):
                if(term.rat == 1):
                    return cp(0)
                
                elif(term.rat > 0):
                    pnum = prime_factors(term.rat.num)
                    pden = prime_factors(term.rat.den)
                    # if N is rational, it is decomposed into primes
                    # so that ln(2 * 3) = ln(2) + ln(3)
                    
                    out = cp(0)
                    for pr in pnum:
                        out += cp( ((pr[1], 1), (ln_function_name, cp(pr[0]))) )
                    for pr in pden:
                        out += cp( ((-pr[1], 1), (ln_function_name, cp(pr[0]))) )
                    
                    return out
                
            elif(len(N.rl.terms[0].irt) == 1):
                if(term.irt[0].name == pow_function_name):
                    if(term.rat == 1):
                        if(N.is_real_and_positive()):
                            # ln(pi ** 3) = 3 * ln(pi)                        
                            return cp( (1, (ln_function_name, term.irt[0].args[0])) ) * term.irt[0].args[1]
                    
                    else:
                        if( (term.rat > 0) and N.is_real_and_positive() ):
                            # ln(2 * pi ** 3) = ln(2) + 3 * ln(pi)
                            return (cp( (1, (ln_function_name, cp(term.rat))) ) +
                                    cp( (1, (ln_function_name, term.irt[0].args[0])) ) * term.irt[0].args[1])
                
                elif(term.irt[0].name == exp_function_name):
                    # ln(exp(Z)) = Z
                    return cp( (1, (ln_function_name, cp(term.rat))) ) + term.irt[0].args[0]
                
                else:
                    if(term.rat != 1):
                        if(term.rat > 0):
                            # ln(3 * pi) = ln(3) + ln(pi)
                            return (cp( (1, (ln_function_name, cp(term.rat))) ) +
                                    cp( (1, (ln_function_name, cp(term.irt[0]))) ))                            
                
            else:
                if(term.rat > 0):
                    # ln( 4/7 * a * b * c * ... ) =
                    # = ln(4) - ln(7) + ln(a) + ln(b) + ln(c) ...
                    out = (cp( (( 1, 1 ), (ln_function_name, cp( (term.rat.num, 1) ))) ) +
                           cp( ((-1, 1 ), (ln_function_name, cp( (term.rat.den, 1) ))) ))
                    for t in term.irt:
                        out += cp( (1, (ln_function_name, cp(t))) )                 
                    
                    return out
            
        return cp(self)
        
    def af_rl(self):
        if(self.args[0].is_zero()):
            return cp(0)
        
        elif(self.args[0].rl.is_zero() and self.args[0].im.might_be_nonzero()):
            return cp( (-1, (im_function_name, cp(self.args[0].im))) )  # Rl(bi) = -Im(b)
        
        else:
            if self.args[0].is_always_real():
                return self.args[0]
            if self.args[0].is_always_imag():
                return cp(0)
            
            if((len(self.args[0].rl.terms) + len(self.args[0].im.terms)) > 1):
                out = cp(0)
                
                # Rl(a + bi) = Rl(a) - Im(b)
                for term in self.args[0].rl.terms:
                    out += cp( (1, (rl_function_name, cp(term))) )
                for term in self.args[0].im.terms:
                    out += cp( (-1, (im_function_name, cp(term))) )
                
                return out
        
        return cp(self)
    
    def af_im(self):
        if(self.args[0].is_zero()):
            return cp(0)
        
        elif(self.args[0].rl.is_zero() and self.args[0].im.might_be_nonzero()):
            return cp( (1, (rl_function_name, cp(self.args[0].im))) )  # Im(bi) = Rl(b)
        
        else:
            if self.args[0].is_always_real():
                return cp(0)
            if self.args[0].is_always_imag():
                return self.args[0] / cp(0, 1)
            
            if((len(self.args[0].rl.terms) + len(self.args[0].im.terms)) > 1):
                out = cp(0)
                
                # Im(a + bi) = Im(a) + Rl(b)
                for term in self.args[0].rl.terms:
                    out += cp( (1, (im_function_name, cp(term))) )
                for term in self.args[0].im.terms:
                    out += cp( (1, (rl_function_name, cp(term))) )
                
                return out
        
        return cp(self)
    
    def af_sign(self):
        N = self.args[0]
        
        if N.rl.is_zero():
            if N.im.is_zero():
                return cp(0)
                
            else:
                return cp( (1, (sign_function_name, cp(N.im))) ) * cp(0, 1)
            
        else:
            if N.im.is_zero():
                if(len(N.rl.terms) == 1):
                    sign = N.rl.terms[0].get_sign()
                    
                    if isinstance(sign, Complex):
                        return sign
            
            else:
                if not N.has_functions():
                    return N / abs(N)
        
        
        if((len(N.rl.terms) == 1) and (len(N.im.terms) == 0)):
            s = N.rl.terms[0].get_sign()
            
            if isinstance(s, Complex):
                return s
        
        if((len(N.rl.terms) == 0) and (len(N.im.terms) == 1)):
            # sign(i z) = i sign(z)
            return cp( 0, Function( (sign_function_name, cp(N.im)) ) )
        
        return cp(self)
    
    def af_abs(self):
        N = self.args[0]
        
        if(N.is_real_and_positive()):
            return N
        
        if(N.rl.is_always_real() and N.im.is_always_real()):
            # abs(a + b i) = sqrt( a**2 + b**2 )
            return cp( (1, (pow_function_name,
                                 cp( (1, (pow_function_name, cp(N.rl), cp(2))) )  + cp( (1, (pow_function_name, cp(N.im), cp(2))) ),
                                 cp((1, 2)))) )
        
        if(N.rl.is_zero() and N.im.might_be_nonzero()):
            # abs(N i) = abs(N)
            return cp( (1, (abs_function_name, cp(N.im))) )
        
        if(N.rl.might_be_nonzero() and N.im.is_zero()):
            if(len(N.rl.terms) == 1):
                term = N.rl.terms[0]
                
                if(len(term.irt) == 0):
                    return abs(term.rat)
                
                if(term.rat < 0):
                    # abs(-N) = abs(N)
                    return cp( (1, (abs_function_name, -N)) )
                
                if( ((len(term.irt) > 0) and (term.rat != 1)) or (len(term.irt) > 1) ):
                    out = abs(term.rat)
                    
                    for fun in term.irt:
                        out *= cp( (1, (abs_function_name, cp(fun))) )
                    
                    # abs( a * b * c... ) = abs(a) * abs(b) * abs(c)...
                    return out
        
        return cp(self)
    
    def af_arg(self):
        N = self.args[0]
        
        if(N.is_always_real()):
            s = N.get_sign()
            
            if(s == cp(1)):
                return cp(0)
            
            elif(s == cp(-1)):
                return cp(pi_constant_name)
        
        if(N.rl.is_zero() and N.im.is_always_real()):
            s = N.im.get_sign()
            
            if(s == cp(1)):
                if TauMode:
                    return cp( ((1, 4), tau_constant_name) )
                
                else:
                    return cp( ((1, 2), pi_constant_name) )
            
            elif(s == cp(-1)):
                if TauMode:
                    return cp( ((-1, 4), tau_constant_name) )
                
                else:
                    return cp( ((-1, 2), pi_constant_name) )
            
            else:
                if(N.im.is_zero()):
                    raise ValueError('arg(0) is not defined')
        
        if ExpressArg:
            if not N.has_functions():
                if(N.rl.is_zero()):
                    if(N.im.is_zero()):
                        return cp(0)
                    
                    else:
                        if(N.im.terms[0].rat.get_sign() < 0):
                            if TauMode:
                                return cp( ((-1, 4), tau_constant_name) )
                            else:
                                return cp( ((-1, 2), pi_constant_name) )
                        
                        else:
                            if TauMode:
                                return cp( ((1, 4), tau_constant_name) )
                            else:
                                return cp( ((1, 2), pi_constant_name) )
                
                else:
                    if(N.im.is_zero()):
                        if(N.rl.terms[0].rat.get_sign() < 0):
                            if TauMode:
                                return cp( ((-1, 2), tau_constant_name) )
                            else:
                                return cp( ((-1, 1), pi_constant_name) )
                        
                        else:
                            return cp(0)
                    
                    else:
                        real = N.rl.terms[0].rat
                        imag = N.im.terms[0].rat
                        
                        if(real > 0):
                            return cp( (1, (arctan_function_name, cp(imag / real))) )
                        
                        else:
                            if(imag > 0):
                                return cp([ pi_constant_name, (1, (arctan_function_name, cp(imag / real))) ])
                            
                            else:
                                return cp([ (-1, pi_constant_name), (1, (arctan_function_name, cp(imag / real))) ])
    
        return cp(self)
    
    def af_sin(self):
        if(self.args[0].is_zero()):
            return cp(0)
        
        else:
            if(len(self.args[0].rl.terms) > 0):
                if(self.args[0].rl.terms[0].rat < Rational(0)):  # sin( -5/2 x ) = -sin( 5/2 x )
                    return -Function( (sin_function_name, -self.args[0]) ).alternate_form()
            
            a = self.args[0]
            
            if( (len(a.rl.terms) == 1) and (len(a.im.terms) == 0) ):
                if(len(a.rl.terms[0].irt) == 1):
                    rat = None
                    
                    if TauMode:
                        if(a.rl.terms[0].irt[0].name == tau_constant_name):
                            rat = a.rl.terms[0].rat * 2                           
                    else:
                        if(a.rl.terms[0].irt[0].name == pi_constant_name):
                            rat = a.rl.terms[0].rat
                
                    if isinstance(rat, Rational):
                        rat = rat % 2
                        
                        if(rat > Rational(1, 2)):
                            # sin( a ) = sin( pi - a )
                            return cp( (1, (sin_function_name, pi_constant * (Rational(1) - rat))) )
                        
                        
                        if(rat == Rational(1, 2)):
                            return cp( (1, 1) )
                        
                        elif(rat == Rational(1, 3)):
                            return cp( ((1, 2), (sqrt_function_name, cp(3))) )
                        
                        elif(rat == Rational(1, 4)):
                            return cp( ((1, 2), (sqrt_function_name, cp(2))) )
                        
                        elif(rat == Rational(1, 6)):
                            return cp( (1, 2) )
                        
                        else:
                            if TauMode:
                                return cp( (1, (sin_function_name, cp(rat / 2) * tau_constant)) )
                            
                            else:
                                return cp( (1, (sin_function_name, cp(rat) * pi_constant)) )
        
        return cp(self)
    
    def af_cos(self):
        if(self.args[0].is_zero()):
            return cp(1)
        else:
            if(len(self.args[0].rl.terms) > 0):
                if(self.args[0].rl.terms[0].rat < Rational(0)):  # cos( -5/2 x ) = cos( 5/2 x )
                    return Function( (cos_function_name, -self.args[0]) ).alternate_form()
            
            a = self.args[0]
            
            if( (len(a.rl.terms) == 1) and (len(a.im.terms) == 0) ):
                if(len(a.rl.terms[0].irt) == 1):
                    rat = None
                    
                    if TauMode:
                        if(a.rl.terms[0].irt[0].name == tau_constant_name):
                            rat = a.rl.terms[0].rat * 2                           
                    else:
                        if(a.rl.terms[0].irt[0].name == pi_constant_name):
                            rat = a.rl.terms[0].rat
                
                    if isinstance(rat, Rational):
                        rat = rat % 2
                        
                        if(rat > Rational(1, 2)):
                            # cos( a ) = -cos( pi - a )
                            return cp( (-1, (cos_function_name, pi_constant * (Rational(1) - rat))) )
                        
                        
                        if(rat == Rational(1, 2)):
                            return cp(0)
                        
                        elif(rat == Rational(1, 3)):
                            return cp( (1, 2) )
                        
                        elif(rat == Rational(1, 4)):
                            return cp( ((1, 2), (sqrt_function_name, cp(2))) )
                        
                        elif(rat == Rational(1, 6)):
                            return cp( ((1, 2), (sqrt_function_name, cp(3))) )
                        
                        else:
                            if TauMode:
                                return cp( (1, (cos_function_name, cp(rat / 2) * tau_constant)) )
                            
                            else:
                                return cp( (1, (cos_function_name, cp(rat) * pi_constant)) )
                            
        return cp(self)
    
    def af_sec(self):
        # sec(x) = 1 / cos(x)
        
        if(self.args[0].is_zero()):
            return cp(1)
        
        else:
            if(len(self.args[0].rl.terms) > 0):
                if(self.args[0].rl.terms[0].rat < Rational(0)):  # 1 / cos( -5/2 x ) = 1 / cos( 5/2 x )
                    return Function( (sec_function_name, -self.args[0]) ).alternate_form()
                
            a = self.args[0]
            
            if( (len(a.rl.terms) == 1) and (len(a.im.terms) == 0) ):
                if(len(a.rl.terms[0].irt) == 1):
                    rat = None
                    
                    if TauMode:
                        if(a.rl.terms[0].irt[0].name == tau_constant_name):
                            rat = a.rl.terms[0].rat * 2                         
                    else:
                        if(a.rl.terms[0].irt[0].name == pi_constant_name):
                            rat = a.rl.terms[0].rat
                
                    if isinstance(rat, Rational):
                        rat = rat % 2
                        
                        if(rat > Rational(1, 2)):
                            # 1 / cos( a ) = -1 / cos( pi - a )
                            return cp( (-1, (sec_function_name, pi_constant * (Rational(1) - rat))) )
                        
                        
                        if(rat == Rational(1, 2)):
                            raise ValueError('sec(pi/2) is not defined')
                        
                        elif(rat == Rational(1, 3)):
                            return cp( (2, 1) )
                        
                        elif(rat == Rational(1, 4)):
                            return cp( (1, (sqrt_function_name, cp(2))) )
                        
                        elif(rat == Rational(1, 6)):
                            return cp( ((2, 3), (pow_function_name, cp(3), cp((1, 2))) ))
                        
                        else:
                            if TauMode:
                                return cp( (1, (sec_function_name, cp(rat / 2) * tau_constant)) )
                            
                            else:
                                return cp( (1, (sec_function_name, cp(rat) * pi_constant)) )
    
        return cp(self)
                                
    def af_csc(self):
        if(self.args[0].is_zero()):
            raise ValueError('csc(0) is not defined')
        
        else:
            if(len(self.args[0].rl.terms) > 0):
                if(self.args[0].rl.terms[0].rat < Rational(0)):  # 1 / sin( -5/2 x ) = -1 / sin( 5/2 x )
                    return -Function( (csc_function_name, -self.args[0]) ).alternate_form()
            
            a = self.args[0]
            
            if( (len(a.rl.terms) == 1) and (len(a.im.terms) == 0) ):
                if(len(a.rl.terms[0].irt) == 1):
                    rat = None
                    
                    if TauMode:
                        if(a.rl.terms[0].irt[0].name == tau_constant_name):
                            rat = a.rl.terms[0].rat * 2                         
                    else:
                        if(a.rl.terms[0].irt[0].name == pi_constant_name):
                            rat = a.rl.terms[0].rat
                
                    if isinstance(rat, Rational):
                        rat = rat % 2
                        
                        if(rat > Rational(1, 2)):
                            # 1 / sin( a ) = 1 / sin( pi - a )
                            return cp( (1, (csc_function_name, pi_constant * (Rational(1) - rat))) )
                        
                        
                        if(rat == Rational(1, 2)):
                            return cp( (1, 1) )
                        
                        elif(rat == Rational(1, 3)):
                            return cp( ((2, 3), (sqrt_function_name, cp(3))) )
                        
                        elif(rat == Rational(1, 4)):
                            return cp( (1, (sqrt_function_name, cp(2))) )
                        
                        elif(rat == Rational(1, 6)):
                            return cp( (2, 1) )
                        
                        else:
                            if TauMode:
                                return cp( (1, (csc_function_name, cp(rat / 2) * tau_constant)) )
                            
                            else:
                                return cp( (1, (csc_function_name, cp(rat) * pi_constant)) )
        
        return cp(self)    
    
    def af_tan(self):
        if(self.args[0].is_zero()):
            return cp(0)
        
        else:
            if(len(self.args[0].rl.terms) > 0):
                if(self.args[0].rl.terms[0].rat < Rational(0)):  # tan( -5/2 x ) = -tan( 5/2 x )
                    return -Function( (tan_function_name, -self.args[0]) ).alternate_form()
            
            a = self.args[0]
            
            if( (len(a.rl.terms) == 1) and (len(a.im.terms) == 0) ):
                if(len(a.rl.terms[0].irt) == 1):
                    rat = None
                    
                    if TauMode:
                        if(a.rl.terms[0].irt[0].name == tau_constant_name):
                            rat = a.rl.terms[0].rat * 2                                
                    else:
                        if(a.rl.terms[0].irt[0].name == pi_constant_name):
                            rat = a.rl.terms[0].rat
                
                    if isinstance(rat, Rational):
                        rat = rat % 2
                        
                        if(rat > Rational(1, 2)):
                            # tan( a ) = -tan( pi - a )
                            return cp( (-1, (tan_function_name, pi_constant * (Rational(1) - rat))) )
                        
                        if(rat == Rational(1, 2)):
                            raise ValueError('tan(pi/2) is undefined')
                        
                        elif(rat == Rational(1, 3)):
                            return cp( (1, (pow_function_name, cp(3), cp((1, 2))) ))
                        
                        elif(rat == Rational(1, 4)):
                            return cp(1)
                        
                        elif(rat == Rational(1, 6)):
                            return cp( ((1, 3), (sqrt_function_name, cp(3))) )
                        
                        else:
                            if TauMode:
                                return cp( (1, (tan_function_name, cp(rat / 2) * tau_constant)) )
                            
                            else:
                                return cp( (1, (tan_function_name, cp(rat) * pi_constant)) )
        
        return cp(self)
    
    def af_cot(self):
        if(self.args[0].is_zero()):
            raise ValueError('cot(0) is undefined')
        
        else:
            if(len(self.args[0].rl.terms) > 0):
                if(self.args[0].rl.terms[0].rat < Rational(0)):  # cot( -5/2 x ) = -cot( 5/2 x )
                    return -Function( (cot_function_name, -self.args[0]) ).alternate_form()
            
            a = self.args[0]
            
            if( (len(a.rl.terms) == 1) and (len(a.im.terms) == 0) ):
                if(len(a.rl.terms[0].irt) == 1):
                    rat = None
                    
                    if TauMode:
                        if(a.rl.terms[0].irt[0].name == tau_constant_name):
                            rat = a.rl.terms[0].rat * 2                                
                    else:
                        if(a.rl.terms[0].irt[0].name == pi_constant_name):
                            rat = a.rl.terms[0].rat
                
                    if isinstance(rat, Rational):
                        rat = rat % 2
                        
                        if(rat > Rational(1, 2)):
                            # cot( a ) = -cot( pi - a )
                            return cp( (-1, (cot_function_name, pi_constant * (Rational(1) - rat))) )
                        
                        if(rat == Rational(1, 2)):
                            return cp(0)
                        
                        elif(rat == Rational(1, 3)):
                            return cp( ((1, 3), (sqrt_function_name, cp(3))) )
                        
                        elif(rat == Rational(1, 4)):
                            return cp(1)
                        
                        elif(rat == Rational(1, 6)):
                            return cp( (1, (pow_function_name, cp(3), cp((1, 2))) ))
                        
                        else:
                            if TauMode:
                                return cp( (1, (cot_function_name, cp(rat / 2) * tau_constant)) )
                            
                            else:
                                return cp( (1, (cot_function_name, cp(rat) * pi_constant)) )
        
        return cp(self)
    
    def af_arcsin(self):
        f = self.args[0]
        
        if(f.is_zero()):
            return cp(0)            
        
        if((len(f.rl.terms) == 1) and (len(f.im.terms) == 0)):
            if(f.rl.terms[0].rat < 0):
                return cp( (-1, (arcsin_function_name, -f)) )
            
            else:
                if TauMode:
                    if(f == cp((1, 2))):
                        return cp( ((1, 12), tau_constant_name) )
                    
                    if(f == cp(((1, 2), (pow_function_name, cp(2), cp((1, 2)))))):
                        return cp( ((1, 8), tau_constant_name) )
                    
                    if(f == cp(((1, 2), (pow_function_name, cp(3), cp((1, 2)))))):
                        return cp( ((1, 6), tau_constant_name) )
                    
                    if(f == cp(1)):
                        return cp( ((1, 4), tau_constant_name) )
                    
                else:
                    if(f == cp((1, 2))):
                        return cp( ((1, 6), pi_constant_name) )
                    
                    if(f == cp(((1, 2), (pow_function_name, cp(2), cp((1, 2)))))):
                        return cp( ((1, 4), pi_constant_name) )
                    
                    if(f == cp(((1, 2), (pow_function_name, cp(3), cp((1, 2)))))):
                        return cp( ((1, 3), pi_constant_name) )
                    
                    if(f == cp(1)):
                        return cp( ((1, 2), pi_constant_name) )
        
        return cp(self)
    
    def af_arccos(self):
        f = self.args[0]
        
        if(f.is_zero()):
            if TauMode:
                return cp( ((1, 4), tau_constant_name) )
            
            else:
                return cp( ((1, 2), pi_constant_name) )
                    
        if((len(f.rl.terms) == 1) and (len(f.im.terms) == 0)):
            if(f.rl.terms[0].rat < 0):
                return pi_constant - cp( (1, (arccos_function_name, -f)) )       
            
            else:
                if TauMode:
                    if(f == cp((1, 2))):
                        return cp( ((1, 6), tau_constant_name) )
                    
                    if(f == cp(((1, 2), (pow_function_name, cp(2), cp((1, 2)))))):
                        return cp( ((1, 8), tau_constant_name) )
                    
                    if(f == cp(((1, 2), (pow_function_name, cp(3), cp((1, 2)))))):
                        return cp( ((1, 12), tau_constant_name) )
                    
                    if(f == cp(1)):
                        return cp(0)
                
                else:
                    if(f == cp((1, 2))):
                        return cp( ((1, 3), pi_constant_name) )
                    
                    if(f == cp(((1, 2), (pow_function_name, cp(2), cp((1, 2)))))):
                        return cp( ((1, 4), pi_constant_name) )
                    
                    if(f == cp(((1, 2), (pow_function_name, cp(3), cp((1, 2)))))):
                        return cp( ((1, 6), pi_constant_name) )
                    
                    if(f == cp(1)):
                        return cp(0)
        
        return cp(self)
    
    def af_arcsec(self):
        f = self.args[0]
        
        if(f.is_zero()):
            raise ValueError('arcsec(0) is not defined')
                    
        if((len(f.rl.terms) == 1) and (len(f.im.terms) == 0)):
            if(f.rl.terms[0].rat < 0):
                # arcsec(-f) = arccos(-1/f) = pi - arccos(1/f) = pi - arcsec(f)
                return pi_constant - cp( (1, (arcsec_function_name, -f)) )       
            
            else:
                if TauMode:
                    if(f == cp((2, 1))):
                        return cp( ((1, 6), tau_constant_name) )
                    
                    if(f == cp( (1, (pow_function_name, cp(2), cp((1, 2))) ))):
                        return cp( ((1, 8), tau_constant_name) )
                    
                    if(f == cp( ((2, 3), (pow_function_name, cp(3), cp((1, 2))) ))):
                        return cp( ((1, 12), tau_constant_name) )
                    
                    if(f == cp(1)):
                        return cp(0)
                
                else:
                    if(f == cp((2, 1))):
                        return cp( ((1, 3), pi_constant_name) )
                    
                    if(f == cp( (1, (pow_function_name, cp(2), cp((1, 2))) ))):
                        return cp( ((1, 4), pi_constant_name) )
                    
                    if(f == cp( ((2, 3), (pow_function_name, cp(3), cp((1, 2))) ))):
                        return cp( ((1, 6), pi_constant_name) )
                    
                    if(f == cp(1)):
                        return cp(0)
        
        return cp(self)
    
    def af_arccsc(self):
        f = self.args[0]
        
        if(f.is_zero()):
            raise ValueError('arccsc(0) is not defined')
        
        if((len(f.rl.terms) == 1) and (len(f.im.terms) == 0)):
            if(f.rl.terms[0].rat < 0):
                # arccsc(-f) = arcsin(-1/f) = -arcsin(1/f) = -arccsc(f)
                return cp( (-1, (arccsc_function_name, -f)) )
            
            else:
                if TauMode:
                    if(f == cp((2, 1))):
                        return cp( ((1, 12), tau_constant_name) )
                    
                    if(f == cp( (1, (pow_function_name, cp(2), cp((1, 2))) ))):
                        return cp( ((1, 8), tau_constant_name) )
                    
                    if(f == cp( ((2, 3), (pow_function_name, cp(3), cp((1, 2))) ))):
                        return cp( ((1, 6), tau_constant_name) )
                    
                    if(f == cp(1)):
                        return cp( ((1, 4), tau_constant_name) )
                    
                else:
                    if(f == cp((2, 1))):
                        return cp( ((1, 6), pi_constant_name) )
                    
                    if(f == cp( (1, (pow_function_name, cp(2), cp((1, 2))) ))):
                        return cp( ((1, 4), pi_constant_name) )
                    
                    if(f == cp( ((2, 3), (pow_function_name, cp(3), cp((1, 2))) ))):
                        return cp( ((1, 3), pi_constant_name) )
                    
                    if(f == cp(1)):
                        return cp( ((1, 2), pi_constant_name) )
        
        return cp(self)    
    
    def af_arctan(self):
        f = self.args[0]
        
        if(f.is_zero()):
            return cp(0)
        
        if((len(f.rl.terms) == 1) and (len(f.im.terms) == 0)):
            if(f.rl.terms[0].rat < 0):
                return cp( (-1, (arctan_function_name, -f)) )
            
            else:
                if TauMode:
                    if(f == cp( ((1, 3), (pow_function_name, cp(3), cp((1, 2))))) ):
                        return cp( ((1, 12), tau_constant_name))
                    
                    if(f == cp(1)):
                        return cp( ((1, 8), tau_constant_name) )
                    
                    if(f == cp( (1, (pow_function_name, cp(3), cp((1, 2)))))):
                        return cp( ((1, 6), tau_constant_name) )
                
                else:
                    if(f == cp( ((1, 3), (pow_function_name, cp(3), cp((1, 2))))) ):
                        return cp( ((1, 6), pi_constant_name))
                    
                    if(f == cp(1)):
                        return cp( ((1, 4), pi_constant_name) )
                    
                    if(f == cp( (1, (pow_function_name, cp(3), cp((1, 2)))))):
                        return cp( ((1, 3), pi_constant_name) )
            
        return cp(self)
    
    def af_arccot(self):
        f = self.args[0]
        
        if(f.is_zero()):
            if TauMode:
                return cp( ((1, 4), tau_constant_name) )
            
            else:
                return cp( ((1, 2), pi_constant_name) )
        
        if((len(f.rl.terms) == 1) and (len(f.im.terms) == 0)):
            if(f.rl.terms[0].rat < 0):
                return cp( (-1, (arccot_function_name, -f)) )
            
            else:
                if TauMode:
                    if(f == cp( ((1, 3), (pow_function_name, cp(3), cp((1, 2))))) ):
                        return cp( ((1, 6), tau_constant_name))
                    
                    if(f == cp(1)):
                        return cp( ((1, 8), tau_constant_name) )
                    
                    if(f == cp( (1, (pow_function_name, cp(3), cp((1, 2)))))):
                        return cp( ((1, 12), tau_constant_name) )
                    
                else:
                    if(f == cp( ((1, 3), (pow_function_name, cp(3), cp((1, 2))))) ):
                        return cp( ((1, 3), pi_constant_name))
                    
                    if(f == cp(1)):
                        return cp( ((1, 4), pi_constant_name) )
                    
                    if(f == cp( (1, (pow_function_name, cp(3), cp((1, 2)))))):
                        return cp( ((1, 6), pi_constant_name) )
        
        return cp(self)
    
    def af_sinh(self):
        f = self.args[0]
        
        if(f.is_zero()):
            return cp(0)
        if((len(f.rl.terms) == 1) and (len(f.im.terms) == 0)):
            if(f.rl.terms[0].rat < 0):
                return cp( (-1, (sinh_function_name, -f)) )
            
        return cp(self)
            
    def af_cosh(self):
        f = self.args[0]
        
        if(f.is_zero()):
            return cp(1)
        if((len(f.rl.terms) == 1) and (len(f.im.terms) == 0)):
            if(f.rl.terms[0].rat < 0):
                return cp( (1, (cosh_function_name, -f)) )      
            
        return cp(self)
            
    def af_sech(self):
        f = self.args[0]
        
        if(f.is_zero()):
            return cp(1)
        if((len(f.rl.terms) == 1) and (len(f.im.terms) == 0)):
            if(f.rl.terms[0].rat < 0):
                return cp( (1, (sech_function_name, -f)) )
        
        return cp(self)
    
    def af_csch(self):
        f = self.args[0]
        
        if(f.is_zero()):
            raise ValueError('csch(0) is not defined')
        if((len(f.rl.terms) == 1) and (len(f.im.terms) == 0)):
            if(f.rl.terms[0].rat < 0):
                return cp( (-1, (csch_function_name, -f)) )
        
        return cp(self)
    
    def af_tanh(self):
        f = self.args[0]
                
        if(f.is_zero()):
            return cp(0)
        
        if((len(f.rl.terms) == 1) and (len(f.im.terms) == 0)):
            if(f.rl.terms[0].rat < 0):
                return cp( (-1, (tanh_function_name, -f)) )
        
        return cp(self)
    
    def af_coth(self):
        f = self.args[0]
        
        if(f.is_zero()):
            raise ValueError('coth(0) is not defined')
        
        if((len(f.rl.terms) == 1) and (len(f.im.terms) == 0)):
            if(f.rl.terms[0].rat < 0):
                return cp( (-1, (coth_function_name, -f)) )        
        
        return cp(self)
    
    def af_arcsinh(self):
        f = self.args[0]
        
        if(f.is_zero()):
            return cp(0)            
        
        if((len(f.rl.terms) == 1) and (len(f.im.terms) == 0)):
            if(f.rl.terms[0].rat < 0):
                return cp( (-1, (arcsinh_function_name, -f)) )
        
        return cp(self)
    
    def af_arccosh(self):
        f = self.args[0]
        
        if(f == cp(0)):
            return cp(0, ((1, 2), pi_constant_name))
        if(f == cp(1)):
            return cp(0)
        
        return cp(self)
    
    def af_arcsech(self):
        f = self.args[0]
        
        if(f == cp(0)):
            raise ValueError('arcsech(0) is not defined')
        
        elif(f == cp(1)):
            return cp(0)
        
        return cp(self)
    
    def af_arccsch(self):
        f = self.args[0]
        
        if(f == cp(0)):
            raise ValueError('arccsch(0) is not defined')
        
        if( (len(f.rl.terms) == 1) and (len(f.im.terms) == 0) ):
            if(f.rl.terms[0].rat < 0):
                return cp( (-1, (arccsch_function_name, -f)) )
            
        return cp(self)
    
    def af_arctanh(self):
        f = self.args[0]
        
        if(f.is_zero()):
            return cp(0)
        if(f == cp(1)):
            raise ValueError('arctanh(1) is not defined')
        
        if((len(f.rl.terms) == 1) and (len(f.im.terms) == 0)):
            if(f.rl.terms[0].rat < 0):
                return cp( (-1, (arctanh_function_name, -f)) )   
        
        return cp(self)
    
    def af_arccoth(self):
        f = self.args[0]
        
        if(f.is_zero()):
            return cp( 0, ((1, 2), pi_constant_name) )
        
        if((len(f.rl.terms) == 1) and (len(f.im.terms) == 0)):
            if(f.rl.terms[0].rat < 0):
                return cp( (-1, (arccoth_function_name, -f)) )        
        
        return cp(self)
            
    def alternate_form(self):
        for i in range(len(self.args)):
            self.args[i].simplify()
        
        if(self.name == pow_function_name):  # N ** P
            return self.af_pow()
            
        elif(self.name == sqrt_function_name):
            # sqrt(x) is changed to pow(x, 1/2)
            return cp( (1, (pow_function_name, self.args[0], cp((1, 2)))) )
        
        elif(self.name == log_function_name):
            # log(a, b) is changed to ln(a) / ln(b)
            return cp( (1, [
                (ln_function_name, self.args[0]),
                (pow_function_name, cp(
                    (1, (ln_function_name, self.args[1]))
                    ), cp(-1))
            ]) )
        
        elif(self.name == pi_constant_name):
            if TauMode:
                return cp( ((1, 2), tau_constant_name) )
            
        elif(self.name == tau_constant_name):
            if not TauMode:
                return cp( ((2, 1), pi_constant_name) )
                
        elif(self.name == exp_function_name):
            return self.af_exp()
            
        elif(self.name == ln_function_name):
            return self.af_ln()
                        
        elif(self.name == rl_function_name):
            return self.af_rl()
                
        elif(self.name == im_function_name):
            return self.af_im()
                        
        elif(self.name == sign_function_name):
            return self.af_sign()
                
        elif(self.name == abs_function_name):
            return self.af_abs()
                    
        elif(self.name == arg_function_name):
            return self.af_arg()
                
        elif(self.name == sin_function_name):
            return self.af_sin()
        
        elif(self.name == cos_function_name):
            return self.af_cos()
        
        elif(self.name == sec_function_name):
            return self.af_sec()
        
        elif(self.name == csc_function_name):
            return self.af_csc()
        
        elif(self.name == tan_function_name):
            return self.af_tan()
        
        elif(self.name == cot_function_name):
            return self.af_cot()
                    
        elif(self.name == arcsin_function_name):
            return self.af_arcsin()
                
        elif(self.name == arccos_function_name):
            return self.af_arccos()
        
        elif(self.name == arcsec_function_name):
            return self.af_arcsec()
        
        elif(self.name == arccsc_function_name):
            return self.af_arccsc()
        
        elif(self.name == arctan_function_name):
            return self.af_arctan()
        
        elif(self.name == arccot_function_name):
            return self.af_arccot()
                    
        elif(self.name == sinh_function_name):
            return self.af_sinh()
                
        elif(self.name == cosh_function_name):
            return self.af_cosh()
                
        elif(self.name == sech_function_name):
            return self.af_sech()
        
        elif(self.name == csch_function_name):
            return self.af_csch()
        
        elif(self.name == tanh_function_name):
            return self.af_tanh()
        
        elif(self.name == coth_function_name):
            return self.af_coth()
                
        elif(self.name == arcsinh_function_name):
            return self.af_arcsinh()
    
        elif(self.name == arccosh_function_name):
            return self.af_arccosh()
            
        elif(self.name == arcsech_function_name):
            return self.af_arcsech()
            
        elif(self.name == arccsch_function_name):
            return self.af_arccsch()
                
        elif(self.name == arctanh_function_name):
            return self.af_arctanh()
        
        elif(self.name == arccoth_function_name):
            return self.af_arccoth()
        
        return cp(self)
    
    def get_sign(self):
        always_positive = [pi_constant_name, e_constant_name, abs_function_name]
        positive_for_real = [cosh_function_name, exp_function_name]
        positive_for_positive = [sinh_function_name, arctan_function_name, arctan_function_name]
        
        if(self.name in always_positive):
            return cp(1)
        if(self.name in positive_for_real):
            if(self.args[0].is_always_real()):
                return cp(1)
        if(self.name in positive_for_positive):
            if(cp( (1, (sign_function_name, self.args[0])) ).simplified() == cp(1)):
                return cp(1)
        
        
        if(self.name == pow_function_name):
            if(self.args[0].is_real_and_positive() and self.args[1].is_always_real()):
                return cp(1)
        
        if(self.name == ln_function_name):
            if(self.args[0].is_zero()):
                raise ValueError('ln(0) is not defined')
            
            if(self.args[0].is_real_rational()):
                r = self.args[0].rl.terms[0].rat
                
                if(r == 1):
                    return None
                elif(r > 1):
                    return cp(1)
                elif(0 < r < 1):
                    return cp(-1)
        
        return None
    
    def is_always_real(self):
        if(self.name == pi_constant_name):
            return True
        if(self.name == e_constant_name):
            return True
        if((self.name == rl_function_name) or
           (self.name == im_function_name) or
           (self.name == arg_function_name) or
           (self.name == abs_function_name)):
            return True
        
        if(self.name == pow_function_name):
            if(self.args[0].is_real_and_positive() and
               self.args[1].is_always_real()):
                return True
        
        if(self.name == ln_function_name):
            if(self.args[0].is_real_and_positive()):
                return True
        
        # if the argument of the following functions is real, the output is real, too
        real_for_real = [exp_function_name, sin_function_name, cos_function_name, tan_function_name, cot_function_name, sinh_function_name, cosh_function_name, tanh_function_name, coth_function_name, arctan_function_name]
        if(self.name in real_for_real):
            if self.args[0].is_always_real():
                return True
        
        return False
    
    def is_real_and_positive(self):
        if not self.is_always_real():
            return False
        
        s = self.get_sign()
        if not isinstance(s, int):
            return False
        if(s != 1):
            return False
        
        return True
    
    def partial(self, i):
        if( (self.name == pow_function_name) and (i == 0) ):
            return self.args[1] * cp( (1, (pow_function_name, self.args[0], self.args[1] - cp(1))) )
        
        elif( (self.name == pow_function_name) and (i == 1) ):
            return cp( (1, (ln_function_name, self.args[0])) ) * cp(self)
        
        elif( (self.name == exp_function_name) and (i == 0) ):
            return cp(self)
        
        elif( (self.name == ln_function_name) and (i == 0) ):
            # d/dx[ln(x)] = 1 / x
            return cp( (1, (pow_function_name, self.args[0], cp(-1))) )
        
        elif( (self.name == log_function_name) and (i == 0) ):
            # d/dx[log2(x)] = 1 / ln(2) x
            return cp( (1, (pow_function_name, self.args[0] * cp( (1, (ln_function_name, self.args[1])) ), cp(-1))) )
        
        elif( (self.name == log_function_name) and (i == 1) ):
            # d/dx[logx(b)] = -ln(b) / (x ln(x)^2)
            return cp( (-1, [(ln_function_name, self.args[0]),
                             (pow_function_name, self.args[1], cp(-1)),
                             (pow_function_name, cp( (1, (ln_function_name, self.args[1])) ), cp(-2))]) )
        
        elif( (self.name == sin_function_name) and (i == 0) ):
            # d/dx[sin x] = cos x
            return cp(Function( (cos_function_name, self.args[0]) ))
        
        elif( (self.name == cos_function_name) and (i == 0) ):
            # d/dx[cos x] = -sin x
            return -cp(Function( (sin_function_name, self.args[0]) ))
        
        elif( (self.name == sec_function_name) and (i == 0) ):
            # d/dx[sec x] = sec(x) tan(x)
            return cp( (1, [(sec_function_name, self.args[0]),
                            (tan_function_name, self.args[0])]) )
        
        elif( (self.name == csc_function_name) and (i == 0) ):
            # d/dx[csc x] = -csc(x) cot(x)
            return cp( (-1, [(csc_function_name, self.args[0]),
                             (cot_function_name, self.args[0])]) )
        
        elif( (self.name == tan_function_name) and (i == 0) ):
            # d/dx[tan x] = 1 / cos^2 x
            return cp( (1, (pow_function_name, cp(Function( (sec_function_name, self.args[0]) )), cp(2))) )
        
        elif( (self.name == cot_function_name) and (i == 0) ):
            # d/dx[cot x] = -1 / sin^2 x
            return cp( (-1, (pow_function_name, cp(Function( (csc_function_name, self.args[0]) )), cp(2))) )
        
        elif( (self.name == arcsin_function_name) and (i == 0) ):
            # d/dx[arcsin x] = 1 / sqrt(1 - x**2)
            return cp( (1, (pow_function_name, cp(1) - cp( (1, (pow_function_name, self.args[0], cp(2))) ), cp((-1, 2)))) )
        
        elif( (self.name == arccos_function_name) and (i == 0) ):
            # d/dx[arccos x] = -1 / sqrt(1 - x**2)
            return cp( (-1, (pow_function_name, cp(1) - cp( (1, (pow_function_name, self.args[0], cp(2))) ), cp((-1, 2)))) )
        
        elif( (self.name == arcsec_function_name) and (i == 0) ):
            # d/dx[arcsec x] = 1 / (x^2 * sqrt(1 - 1 / x^2))
            
            y = cp(1) - cp( (1, (pow_function_name, self.args[0], cp(-2))) )
            
            return cp( (1, [( pow_function_name, self.args[0], cp(-2) ),
                            ( pow_function_name, y, cp((-1, 2)) )]) )
        
        elif( (self.name == arccsc_function_name) and (i == 0) ):
            # d/dx[arccsc x] = -1 / (x^2 * sqrt(1 - 1 / x^2))
            
            y = cp(1) - cp( (1, (pow_function_name, self.args[0], cp(-2))) )
            
            return cp( (-1, [( pow_function_name, self.args[0], cp(-2) ),
                             ( pow_function_name, y, cp((-1, 2)) )]) )
        
        elif( (self.name == arctan_function_name) and (i == 0) ):
            # d/dx[arctan x] = 1 / (1 + x**2)
            return cp( (1, (pow_function_name, cp(1) + self.args[0] ** 2, cp(-1))) )
        
        elif( (self.name == arccot_function_name) and (i == 0) ):
            # d/dx[arccot x] = -1 / (1 + x**2)
            return cp( (-1, (pow_function_name, cp(1) + self.args[0] ** 2, cp(-1))) )
        
        elif( (self.name == sinh_function_name) and (i == 0) ):
            # d/dx[sinh x] = cosh x
            return cp(Function( (cosh_function_name, self.args[0]) ))
        
        elif( (self.name == cosh_function_name) and (i == 0) ):
            # d/dx[cosh x] = -sinh x
            return cp(Function( (sinh_function_name, self.args[0]) ))
        
        elif( (self.name == sech_function_name) and (i == 0) ):
            # d/dx[sech x] = -sech(x) tanh(x)
            return cp( (-1, [(sech_function_name, self.args[0]),
                             (tanh_function_name, self.args[0])]) )
        
        elif( (self.name == csch_function_name) and (i == 0) ):
            # d/dx[csch x] = -csch(x) coth(x)
            return cp( (-1, [(csch_function_name, self.args[0]),
                             (coth_function_name, self.args[0])]) )
        
        elif( (self.name == tanh_function_name) and (i == 0) ):
            # d/dx[tanh x] = 1 / cosh^2 x
            return cp( (1, (pow_function_name, cp(Function( (cosh_function_name, self.args[0]) )), cp(-2))) )
        elif( (self.name == coth_function_name) and (i == 0) ):
            # d/dx[coth x] = -1 / sinh^2 x
            return cp( (-1, (pow_function_name, cp(Function( (sinh_function_name, self.args[0]) )), cp(-2))) )
        
        elif( (self.name == arcsinh_function_name) and (i == 0) ):
            # d/dx[arcsin x] = 1 / sqrt(1 + x**2)
            return cp( (1, (pow_function_name, cp( (1, (pow_function_name, self.args[0], cp(2))) ) + cp(1), cp((-1, 2)))) )
        
        elif( (self.name == arccosh_function_name) and (i == 0) ):
            # d/dx[arccos x] = 1 / (sqrt(x - 1) * sqrt(x + 1) )
            return cp( (1, [
                (pow_function_name, self.args[0] - cp(1), cp( ((-1, 2)) )),
                (pow_function_name, self.args[0] + cp(1), cp( ((-1, 2)) ))
            ]) )
        
        elif( (self.name == arcsech_function_name) and (i == 0) ):
            # d/dx[arcsech x] = -1 / (x^2 * sqrt(1/x - 1) * sqrt(1/x + 1))
            
            one_over_x = cp( (1, (pow_function_name, self.args[0], cp(-1))) )
            
            return cp( (-1, [
                (pow_function_name, self.args[0], cp(-2)),
                (pow_function_name, one_over_x - cp(1), cp((-1, 2))),
                (pow_function_name, one_over_x + cp(1), cp((-1, 2)))
                             ]) )
        
        elif( (self.name == arccsch_function_name) and (i == 0) ):
            # d/dx[arccsch x] = -1 / (x^2 * sqrt(1/x^2 + 1))
            
            return cp( (-1, [
                (pow_function_name, self.args[0], cp(-2)),
                (pow_function_name, cp(1) + cp( (1, (
                    pow_function_name, self.args[0], cp(-2))) ), cp((-1, 2)))
            ]) )
        
        elif( (self.name == arctanh_function_name) and (i == 0) ):
            # d/dx[arctan x] = 1 / (1 - x**2)
            return cp( (1, (pow_function_name, cp(1) - self.args[0] ** 2, cp(-1))) )
        
        elif( (self.name == arccoth_function_name) and (i == 0) ):
            # d/dx[arccot x] = 1 / (1 - x**2)
            return cp( (1, (pow_function_name, cp(1) - self.args[0] ** 2, cp(-1))) )   
        
        else:
            reduced_name, partials = self.decode_name()
            
            partials[i] += 1
            
            return cp(Function( (Function.make_up_a_name(reduced_name, partials), *self.args) ))
    
    def differentiated(self, var):
        if(self.name == var):
            return cp(1)
        
        else:
            out = cp(0)
            
            # df(a, b, c)/dx = pf/pa * da/dx + pf/pb * db/dx + pf/pc * dc/dx
            # where d/dz denotes full derivative and p/px means partial derivative.
            for i in range(len(self.args)):
                out += self.partial(i) * self.args[i].differentiated(var)
            
            return out
        
    def numerical(self, substitution_list):
        for subst in substitution_list:  # trying to find itself in substitution_list
            if(subst[0] == self):
                return subst[1]
            
            
        args = []
        for arg in self.args:
            args.append(arg.numerical(substitution_list))
        
        if(self.name == pow_function_name):
            if( (args[0].rl > 0) and (args[0].im == 0) and (args[1].im == 0) ):
                return numpy.power(args[0], args[1])
            
            else:
                return args[0] ** args[1]
        
        elif(self.name == sqrt_function_name):
            if( (args[0].rl > 0) and (args[0].im == 0) ):
                return numpy.sqrt(args[0].rl)
            
            else:
                return args[0] ** 0.5
        
        elif(self.name == exp_function_name):
            return cexp(args[0])
        
        elif(self.name == ln_function_name):
            if( (args[0].im == 0) and (args[0].rl > 0) ):
                return Numerical(numpy.log(args[0].rl), 0.0)
            
            else:
                return cln(args[0])
        
        elif(self.name == log_function_name):
            if( (args[0].im == 0) and (args[0].rl > 0) and
                (args[1].im == 0) and (args[1].rl > 0) ):
                return Numerical(numpy.log(args[0].rl) / numpy.log(args[1].rl), 0.0)
            
            else:
                return cln(args[0]) / cln(args[1])
        
        elif(self.name == rl_function_name):
            return Numerical(args[0].rl, 0)
        
        elif(self.name == im_function_name):
            return Numerical(args[0].im, 0)
        
        elif(self.name == abs_function_name):
            return Numerical(abs(args[0]), 0.0)
        
        elif(self.name == arg_function_name):
            return Numerical((args[0]).arg(), 0.0)
        
        elif(self.name == sign_function_name):
            return csign(args[0])
            
        elif(self.name == sin_function_name):
            if(args[0].im == 0):
                return Numerical( numpy.sin(args[0].rl), 0.0 )
            
            else:
                return csin(args[0])
        
        elif(self.name == cos_function_name):
            if(args[0].im == 0):
                return Numerical( numpy.cos(args[0].rl), 0.0 )
            
            else:
                return ccos(args[0])
            
        elif(self.name == sec_function_name):
            if(args[0].im == 0):
                return Numerical( 1 / numpy.cos(args[0].rl), 0.0 )
            
            else:
                return csec(args[0])
            
        elif(self.name == csc_function_name):
            if(args[0].im == 0):
                return Numerical( 1 / numpy.sin(args[0].rl), 0.0 )
            
            else:
                return ccsc(args[0])            
        
        elif(self.name == tan_function_name):
            if(args[0].im == 0):
                return Numerical( numpy.tan(args[0].rl), 0.0 )
            
            else:
                return ctan(args[0])
            
        elif(self.name == cot_function_name):
            if(args[0].im == 0):
                return Numerical( 1 / numpy.tan(args[0].rl), 0.0 )    
            
            else:
                return ccot(args[0])
        
        elif(self.name == sinh_function_name):
            if(args[0].im == 0):
                return Numerical( numpy.sinh(args[0].rl), 0.0 )
            
            else:
                return csinh(args[0])
            
        elif(self.name == cosh_function_name):
            if(args[0].im == 0):
                return Numerical( numpy.cosh(args[0].rl), 0.0 )
            
            else:
                return ccosh(args[0])
        
        elif(self.name == sech_function_name):
            if(args[0].im == 0):
                return Numerical( 1 / numpy.cosh(args[0].rl), 0.0 )
            
            else:
                return csech(args[0])
            
        elif(self.name == csch_function_name):
            if(args[0].im == 0):
                return Numerical( 1 / numpy.sinh(args[0].rl), 0.0 )
            
            else:
                return ccsch(args[0])
            
        elif(self.name == tanh_function_name):
            if(args[0].im == 0):
                return Numerical( numpy.tanh(args[0].rl), 0.0 )
            
            else:
                return ctanh(args[0])
        
        elif(self.name == coth_function_name):
            if(args[0].im == 0):
                return Numerical( 1 / numpy.tanh(args[0].rl), 0.0 )
            
            else:
                return ccoth(args[0])
            
        elif(self.name == arcsin_function_name):
            if( (args[0].im == 0) and (-1 <= args[0].rl <= 1) ):
                return Numerical( numpy.arcsin(args[0].rl), 0.0 )
            
            else:
                return carcsin(args[0])
            
        elif(self.name == arccos_function_name):
            if( (args[0].im == 0) and (-1 <= args[0].rl <= 1) ):
                return Numerical( numpy.arccos(args[0].rl), 0.0 )
            
            else:
                return carccos(args[0])
            
        elif(self.name == arcsec_function_name):
            if( (args[0].im == 0) and (abs(args[0].rl) >= 1) ):
                return Numerical( numpy.arccos(1 / args[0].rl), 0.0 )
            
            else:
                return carcsec(args[0])
            
        elif(self.name == arccsc_function_name):
            if( (args[0].im == 0) and (abs(args[0].rl) >= 1) ):
                
                return Numerical( numpy.arcsin(1 / args[0].rl), 0.0 )
            
            else:
                return carccsc(args[0])
        
        elif(self.name == arctan_function_name):
            if(args[0].im == 0):
                return Numerical( numpy.arctan(args[0].rl), 0.0 )
            
            else:
                return carctan(args[0])
            
        elif(self.name == arccot_function_name):
            if(args[0].im == 0):
                # arccot(x) = arctan(1 / x)
                return Numerical( numpy.arctan(1 / args[0].rl), 0.0 )
            
            else:
                return carccot(args[0])
            
        elif(self.name == arcsinh_function_name):
            if(args[0].im == 0):
                return Numerical( numpy.arcsinh(args[0].rl), 0.0 )
            
            else:
                return carcsinh(args[0])
            
        elif(self.name == arccosh_function_name):
            if( (args[0].im == 0) and (args[0].rl >= 1) ):
                return Numerical( numpy.arccosh(args[0].rl), 0.0 )
            
            else:
                return carccosh(args[0])
            
        elif(self.name == arcsech_function_name):
            if( (args[0].im == 0) and (0 < args[0].rl <= 1) ):
                return Numerical( numpy.arccosh(1 / args[0].rl), 0.0 )
            
            else:
                return carcsech(args[0])
            
        elif(self.name == arccsch_function_name):
            if( (args[0].im == 0) and (args[0].rl != 0) ):
                return Numerical( numpy.arcsinh(1 / args[0].rl), 0.0 )
            
            else:
                return carccsch(args[0])
        
        elif(self.name == arctanh_function_name):
            if( (args[0].im == 0) and (abs(args[0].rl) < 1) ):
                return Numerical( numpy.arctanh(args[0].rl), 0.0 )
            
            else:
                return carctanh(args[0])
            
        elif(self.name == arccoth_function_name):
            if( (args[0].im == 0) and (abs(args[0].rl) > 1) ):
                # arccoth(x) = arctanh(1 / x)
                return Numerical( numpy.arctanh(1 / args[0].rl), 0.0 )
            
            else:
                return carccoth(args[0])
            
        elif(self.name == pi_constant_name):
            return Numerical(numpy.pi, 0.0)
        
        elif(self.name == tau_constant_name):
            return Numerical(2 * numpy.pi, 0.0)
    
        elif(self.name == e_constant_name):
            return Numerical(numpy.e, 0.0)
        
        else:
            raise UnknownFunctionError(self)
    
    
class Term:
    def __init__(self, newrat, newirt):  # object of class Term is a product of
                                         # a rational coefficient (which in turn is a class Rational item)
                                         # and a number of irrational factors
                                         #    (they are class Function items)
        if not(isinstance(newrat, Rational) and
               isinstance(newirt, list)):
            raise TypeError('Wrong arguments` type: expected <class \'Algebraica.Rational\'>, <class \'list\'>, got ' + str(type(newrat)) + ', ' + str(type(newirt)))
        
        self.rat = newrat.copy()
        self.irt = copy_list(newirt)
        
        self.sort()
        self.group()
        
    def serial_index(t):
        out = ''
        
        for f in t.irt:
            out += f.serial_index()
        
        return out
    
    def __mul__(first, second):
        if(isinstance(second, Term)):
            
            return Term(
                first.rat * second.rat,
                first.irt + second.irt
            )
        
        elif (isinstance(second, int) or isinstance(second, float) or isinstance(second, Rational)):
            return Term(
                first.rat * second,
                first.irt
            )
        
        elif isinstance(second, str):
            return Term(
                first.rat,
                first.irt + [Function((second))]
            )
        
        else:
            raise TypeError('Wrong argument type: ' + str(type(second)))
    def __imul__(first, second):
        return first * second
    
    def __truediv__(first, second):
        if(isinstance(second, Term)):
            
            return Term(
                first.rat * second.rat,
                first.irt + second.irt
            )
        
        elif ((isinstance(second, int) or isinstance(second, float)) or isinstance(second, Rational)):
            return Term(
                first.rat / second,
                first.irt
            )
        
        else:
            raise TypeError('Wrong argument type: ' + str(type(second)))
        
    def __eq__(first, second):
        if isinstance(second, Term):
            return ((first.rat == second.rat) and compare_lists(first.irt, second.irt))
        
        else:
            raise TypeError('Wrong argument type: ' + str(type(second)))
    def __ne__(first, second):
        return not (first == second)
        
    def show(self, elevation, one_line=False):
        out = ''
        
        tab = ' ' * elevation * TabulationStep
        
        show_rat = (self.rat != 1) or (len(self.irt) == 0)
        if show_rat:
            out += tab + str(self.rat)
            
        if(len(self.irt) > 0):
            if show_rat:
                if ColorfulPrint:
                    out += colorama.Fore.BLACK + ' * ' + colorama.Style.RESET_ALL
                
                else:
                    out += ' * '
                
            for j in range(len(self.irt)):
                if((j != 0) or not show_rat):
                    out += tab
                    
                out += self.irt[j].show(elevation)
                
                if(j != len(self.irt) - 1):
                    if ColorfulPrint:
                        if(one_line):
                            out += colorama.Fore.BLACK + ' * ' + colorama.Style.RESET_ALL
                        else:
                            out += colorama.Fore.BLACK + ' *\n' + colorama.Style.RESET_ALL

                        
                    else:
                        if(one_line):
                            out += ' * '
                        else:
                            out += ' *\n'
        
        return out
    def __str__(self):
        return self.show(0)
    
    def wolfram_lang(self):
        out = str(self.rat)
        
        for i in range(len(self.irt)):
            out += ' * ' + self.irt[i].wolfram_lang()
        
        return out
    
    def compileable(self):
        out = '('
        
        out += self.rat.compileable()
        
        if(len(self.irt) == 0):
            pass
        
        elif(len(self.irt) == 1):
            out += ', ' + self.irt[0].compileable()
        
        else:
            out += ', ['
            for i in range(len(self.irt)):
                out += self.irt[i].compileable()
                
                if(i != len(self.irt) - 1):
                    out += ', '
            out += ']'
        
        out += ')'
        return out
    
    def sort(self):
        self.irt.sort(key=Function.serial_index)
    
    def group(self):
        if(len(self.irt) == 0):
            return
        
        new_irt = []
        
        i = 0
        while(i < len(self.irt) - 1):
            f = self.irt[i]
            
            p = 0
            while((i < len(self.irt)) and (self.irt[i] == f)):
                p += 1
                i += 1
            
            if(p == 1):
                new_irt.append(f)
            else:
                new_irt.append(Function( (pow_function_name, cp(f), cp(p)) ))
        if(i == len(self.irt) - 1):
            new_irt.append(self.irt[i])
        
        self.irt = copy_list(new_irt)
    
    def collect_exponentials(self):
        P = cp(0)
        new_irt = []
        
        for fun in self.irt:
            if(fun.name == exp_function_name):
                P += fun.args[0]
            
            elif(fun.name == e_constant_name):
                P += cp(1)
            
            else:
                new_irt.append(fun)
                
        if(P.might_be_nonzero()):
            new_irt.append(Function( (exp_function_name, P) ))
        
        self.irt = copy_list(new_irt)
        
    '''def collect_powers(self):   # this is the old verion of the function
                                   # now it`s obsolete
        new_irt = []
        
        irt = self.irt
        while(len(irt) > 0):
            fun = irt.pop()
            
            if(fun.name == pow_function_name):
                base = fun.args[0]
                power = fun.args[1]
                
                i = 0
                while(i < len(irt)):
                    if(irt[i].name == pow_function_name):
                        if(irt[i].args[0] == base):
                            power += irt.pop(i).args[1]
                        
                        else:
                            i += 1
                    
                    else:
                        i += 1
                
                new_irt.append(Function( (pow_function_name, base, power) ))
            
            else:
                new_irt.append(fun)
        
        self.irt = copy_list(new_irt)'''
    
    def collect_powers(self):
        out = cp(self.rat)
        
        while(len(self.irt) > 0):
            func = self.irt[0]
            if(func.name == pow_function_name):
                base = func.args[0]
            else:
                base = cp(func)
            power = cp(0)
            
            i = 0
            while(i < len(self.irt)):
                if(cp(self.irt[i]) == base):
                    power += cp(1)
                    self.irt.pop(i)
                    
                elif( (self.irt[i].name == pow_function_name) and
                      (self.irt[i].args[0] == base) ):
                    power += self.irt[i].args[1]
                    self.irt.pop(i)
                
                else:
                    i += 1
            
            if(power == cp(1)):
                out *= base
            else:
                out *= cp(Function( (pow_function_name,
                                     base,
                                     power) ))
        
        return out
        
    def collect_trig(self):
        def increment_pows(pow_sin, pow_cos, name, increment=1):
            if(name == sin_function_name):
                pow_sin += increment
            
            elif(name == cos_function_name):
                pow_cos += increment
                
            elif(name == sec_function_name):
                pow_cos -= increment
            
            elif(name == csc_function_name):
                pow_sin -= increment
            
            elif(name == tan_function_name):
                pow_sin += increment
                pow_cos -= increment
            
            elif(name == cot_function_name):
                pow_sin -= increment
                pow_cos += increment
            
            return pow_sin, pow_cos
        
        def sin_to_pow(p, arg):
            if(p == 1):
                return Function( (sin_function_name, arg) )
            
            else:
                return Function( (pow_function_name,
                                  cp( (1, (sin_function_name, arg)) ),
                                  cp(p)) )
        
        def cos_to_pow(p, arg):
            if(p == 1):
                return Function( (cos_function_name, arg) )
            
            else:
                return Function( (pow_function_name,
                                  cp( (1, (cos_function_name, arg)) ),
                                  cp(p)) )

        def sec_to_pow(p, arg):
            if(p == 1):
                return Function( (sec_function_name, arg) )
            
            else:
                return Function( (pow_function_name,
                                  cp( (1, (sec_function_name, arg)) ),
                                  cp(p)) )
        
        def csc_to_pow(p, arg):
            if(p == 1):
                return Function( (csc_function_name, arg) )
            
            else:
                return Function( (pow_function_name,
                                  cp( (1, (csc_function_name, arg)) ),
                                  cp(p)) )

        def tan_to_pow(p, arg):
            if(p == 1):
                return Function( (tan_function_name, arg) )
            
            else:
                return Function( (pow_function_name,
                                  cp( (1, (tan_function_name, arg)) ),
                                  cp(p)) )
            
        def cot_to_pow(p, arg):
            if(p == 1):
                return Function( (cot_function_name, arg) )
            
            else:
                return Function( (pow_function_name,
                                  cp( (1, (cot_function_name, arg)) ),
                                  cp(p)) )
        
        def get_psin_and_pcos(irt, arg, name, pow_sin, pow_cos):
            i = 0
            while(i < len(irt)):
                if((irt[i].name in trig_names) and
                   (irt[i].args[0] == arg)):
                    f = irt.pop(i)
                    
                    pow_sin, pow_cos = increment_pows(pow_sin, pow_cos, f.name)
                
                elif( (irt[i].name == pow_function_name) and
                      (len(irt[i].args[0].rl.terms) == 1) and
                      (len(irt[i].args[0].im.terms) == 0) and
                      (len(irt[i].args[0].rl.terms[0].irt) == 1) and
                      (irt[i].args[0].rl.terms[0].irt[0].name in trig_names) and
                      (irt[i].args[0].rl.terms[0].irt[0].args[0] == arg) and
                      irt[i].args[1].is_real_int() ):
                    
                    f = irt.pop(i)
                    
                    if(len(f.args[1].rl.terms) == 0):
                        p = 0
                    else:
                        p = f.args[1].rl.terms[0].rat.num
                    name = f.args[0].rl.terms[0].irt[0].name
                    
                    pow_sin, pow_cos = increment_pows(pow_sin, pow_cos, name, p)
                
                else:
                    i += 1
            
            return pow_sin, pow_cos        
        
        def collect_trig_with_specific_arg(irt, arg, name, psin0, pcos0):
            pow_sin, pow_cos = get_psin_and_pcos(irt, arg, name, psin0, pcos0)
            
            if(pow_sin > 0):
                if(pow_cos > 0):
                    return [sin_to_pow(pow_sin, arg),
                            cos_to_pow(pow_cos, arg)]
                
                elif(pow_cos < 0):
                    if(abs(pow_cos) > abs(pow_sin)):
                        return [sec_to_pow(abs(pow_cos) - pow_sin, arg),
                                tan_to_pow(pow_sin, arg)]
                    
                    elif(abs(pow_cos) < abs(pow_sin)):
                        return [sin_to_pow(pow_sin - abs(pow_cos), arg),
                                tan_to_pow(abs(pow_cos), arg)]
                    
                    else:  # pow_sin = -pow_cos
                        return [tan_to_pow(pow_sin, arg)]
                
                else:  # pow_cos = 0
                    return [sin_to_pow(pow_sin, arg)]
            
            elif(pow_sin < 0):
                if(pow_cos > 0):
                    if(abs(pow_cos) > abs(pow_sin)):
                        return [cos_to_pow(pow_cos - abs(pow_sin), arg),
                                cot_to_pow(abs(pow_sin), arg)]
                    
                    elif(abs(pow_cos) < abs(pow_sin)):
                        return [csc_to_pow(abs(pow_sin) - pow_cos, arg),
                                cot_to_pow(pow_cos, arg)]
                    
                    else:  # pow_cos = -pow_sin
                        return [cot_to_pow(abs(pow_sin), arg)]
                
                elif(pow_cos < 0):
                    return [sec_to_pow(abs(pow_cos), arg),
                            csc_to_pow(abs(pow_sin), arg)]
                
                else:  # pow_cos = 0
                    return [csc_to_pow(abs(pow_sin), arg)]
            
            else:  # pow_sin = 0
                if(pow_cos > 0):
                    return [cos_to_pow(pow_cos, arg)]
                
                elif(pow_cos < 0):
                    return [sec_to_pow(abs(pow_cos), arg)]
                
                else:  # pow_sin = pow_cos = 0
                    return []
            
        
        
        new_irt = []
        
        trig_names = (sin_function_name,
                      cos_function_name,
                      sec_function_name, 
                      csc_function_name,
                      tan_function_name,
                      cot_function_name)
        
        while(len(self.irt) > 0):
            func = self.irt.pop(0)
            
            if(func.name in trig_names):
                psin, pcos = increment_pows(0, 0, func.name)
                
                new_irt += collect_trig_with_specific_arg(self.irt, func.args[0], func.name, psin, pcos)
            
            
            elif( (func.name == pow_function_name) and
                  (len(func.args[0].rl.terms) == 1) and
                  (len(func.args[0].im.terms) == 0) and
                  (len(func.args[0].rl.terms[0].irt) == 1) and
                  (func.args[0].rl.terms[0].irt[0].name in trig_names) and
                  func.args[1].is_real_int() ):
                
                if(len(func.args[1].rl.terms) == 0):
                    p = 0
                else:
                    p = func.args[1].rl.terms[0].rat.num
                name = func.args[0].rl.terms[0].irt[0].name
                    
                psin, pcos = increment_pows(0, 0, name, increment=p)
                
                new_irt += collect_trig_with_specific_arg(self.irt, func.args[0].rl.terms[0].irt[0].args[0], func.name, psin, pcos)
                self.rat *= func.args[0].rl.terms[0].rat ** p
            
            else:
                new_irt.append(func)
        
        self.irt = copy_list(new_irt)
    
    def collect_htrig(self):
        def increment_pows(pow_sinh, pow_cosh, name, increment=1):
            if(name == sinh_function_name):
                pow_sinh += increment
            
            elif(name == cosh_function_name):
                pow_cosh += increment
                
            elif(name == sech_function_name):
                pow_cosh -= increment
            
            elif(name == csch_function_name):
                pow_sinh -= increment
            
            elif(name == tanh_function_name):
                pow_sinh += increment
                pow_cosh -= increment
            
            elif(name == coth_function_name):
                pow_sinh -= increment
                pow_cosh += increment
            
            return pow_sinh, pow_cosh
        
        def sinh_to_pow(p, arg):
            if(p == 1):
                return Function( (sinh_function_name, arg) )
            
            else:
                return Function( (pow_function_name,
                                  cp( (1, (sinh_function_name, arg)) ),
                                  cp(p)) )
        
        def cosh_to_pow(p, arg):
            if(p == 1):
                return Function( (cosh_function_name, arg) )
            
            else:
                return Function( (pow_function_name,
                                  cp( (1, (cosh_function_name, arg)) ),
                                  cp(p)) )

        def sech_to_pow(p, arg):
            if(p == 1):
                return Function( (sech_function_name, arg) )
            
            else:
                return Function( (pow_function_name,
                                  cp( (1, (sech_function_name, arg)) ),
                                  cp(p)) )
        
        def csch_to_pow(p, arg):
            if(p == 1):
                return Function( (csch_function_name, arg) )
            
            else:
                return Function( (pow_function_name,
                                  cp( (1, (csch_function_name, arg)) ),
                                  cp(p)) )

        def tanh_to_pow(p, arg):
            if(p == 1):
                return Function( (tanh_function_name, arg) )
            
            else:
                return Function( (pow_function_name,
                                  cp( (1, (tanh_function_name, arg)) ),
                                  cp(p)) )
            
        def coth_to_pow(p, arg):
            if(p == 1):
                return Function( (coth_function_name, arg) )
            
            else:
                return Function( (pow_function_name,
                                  cp( (1, (coth_function_name, arg)) ),
                                  cp(p)) )
        
        def get_psinh_and_pcosh(irt, arg, name, pow_sinh, pow_cosh):
            i = 0
            while(i < len(irt)):
                if((irt[i].name in htrig_names) and
                   (irt[i].args[0] == arg)):
                    f = irt.pop(i)
                    
                    pow_sinh, pow_cosh = increment_pows(pow_sinh, pow_cosh, f.name)
                
                elif( (irt[i].name == pow_function_name) and
                      (len(irt[i].args[0].rl.terms) == 1) and
                      (len(irt[i].args[0].im.terms) == 0) and
                      (len(irt[i].args[0].rl.terms[0].irt) == 1) and
                      (irt[i].args[0].rl.terms[0].irt[0].name in htrig_names) and
                      (irt[i].args[0].rl.terms[0].irt[0].args[0] == arg) and
                      irt[i].args[1].is_real_int() ):
                    
                    f = irt.pop(i)
                    
                    if(len(f.args[1].rl.terms) == 0):
                        p = 0
                    else:
                        p = f.args[1].rl.terms[0].rat.num
                    name = f.args[0].rl.terms[0].irt[0].name
                    
                    pow_sinh, pow_cosh = increment_pows(pow_sinh, pow_cosh, name, p)
                
                else:
                    i += 1
            
            return pow_sinh, pow_cosh        
        
        def collect_htrig_with_specific_arg(irt, arg, name, psinh0, pcosh0):
            pow_sinh, pow_cosh = get_psinh_and_pcosh(irt, arg, name, psinh0, pcosh0)
            
            if(pow_sinh > 0):
                if(pow_cosh > 0):
                    return [sinh_to_pow(pow_sinh, arg),
                            cosh_to_pow(pow_cosh, arg)]
                
                elif(pow_cosh < 0):
                    if(abs(pow_cosh) > abs(pow_sinh)):
                        return [sech_to_pow(abs(pow_cosh) - pow_sinh, arg),
                                tanh_to_pow(pow_sinh, arg)]
                    
                    elif(abs(pow_cosh) < abs(pow_sinh)):
                        return [sinh_to_pow(pow_sinh - abs(pow_cosh), arg),
                                tanh_to_pow(abs(pow_cosh), arg)]
                    
                    else:  # pow_sinh = -pow_cosh
                        return [tanh_to_pow(pow_sinh, arg)]
                
                else:  # pow_cosh = 0
                    return [sinh_to_pow(pow_sinh, arg)]
            
            elif(pow_sinh < 0):
                if(pow_cosh > 0):
                    if(abs(pow_cosh) > abs(pow_sinh)):
                        return [cosh_to_pow(pow_cosh - abs(pow_sinh), arg),
                                coth_to_pow(abs(pow_sinh), arg)]
                    
                    elif(abs(pow_cosh) < abs(pow_sinh)):
                        return [csch_to_pow(abs(pow_sinh) - pow_cosh, arg),
                                coth_to_pow(pow_cosh, arg)]
                    
                    else:  # pow_cosh = -pow_sinh
                        return [coth_to_pow(abs(pow_sinh), arg)]
                
                elif(pow_cosh < 0):
                    return [sech_to_pow(abs(pow_cosh), arg),
                            csch_to_pow(abs(pow_sinh), arg)]
                
                else:  # pow_cosh = 0
                    return [csch_to_pow(abs(pow_sinh), arg)]
            
            else:  # pow_sinh = 0
                if(pow_cosh > 0):
                    return [cosh_to_pow(pow_cosh, arg)]
                
                elif(pow_cosh < 0):
                    return [sech_to_pow(abs(pow_cosh), arg)]
                
                else:  # pow_sinh = pow_cosh = 0
                    return []
            
        
        
        new_irt = []
        
        htrig_names = (sinh_function_name,
                      cosh_function_name,
                      sech_function_name, 
                      csch_function_name,
                      tanh_function_name,
                      coth_function_name)
        
        while(len(self.irt) > 0):
            func = self.irt.pop(0)
            
            if(func.name in htrig_names):
                psinh, pcosh = increment_pows(0, 0, func.name)
                
                new_irt += collect_htrig_with_specific_arg(self.irt, func.args[0], func.name, psinh, pcosh)
            
            
            elif( (func.name == pow_function_name) and
                  (len(func.args[0].rl.terms) == 1) and
                  (len(func.args[0].im.terms) == 0) and
                  (len(func.args[0].rl.terms[0].irt) == 1) and
                  (func.args[0].rl.terms[0].irt[0].name in htrig_names) and
                  func.args[1].is_real_int() ):
                
                if(len(func.args[1].rl.terms) == 0):
                    p = 0
                else:
                    p = func.args[1].rl.terms[0].rat.num
                name = func.args[0].rl.terms[0].irt[0].name
                    
                psinh, pcosh = increment_pows(0, 0, name, increment=p)
                
                new_irt += collect_htrig_with_specific_arg(self.irt, func.args[0].rl.terms[0].irt[0].args[0], func.name, psinh, pcosh)
                self.rat *= func.args[0].rl.terms[0].rat ** p
            
            else:
                new_irt.append(func)
        
        self.irt = copy_list(new_irt)        
    
    def change_sign(self):
        self.rat.change_sign()    
    
    def copy(self):
        return Term(self.rat, self.irt)
    
    def eq_irt(self, other):
        return compare_lists(self.irt, other.irt)
    
    def might_be_nonzero(self):
        return self.rat.might_be_nonzero()
    
    def has_functions(self):
        if not self.might_be_nonzero():
            return False
        
        return (len(self.irt) > 0)
    
    def is_constant(self):
        for fun in self.irt:
            if(len(fun.args) > 0):
                return False
        
        return True
    
    def is_always_real(self):
        for fun in self.irt:
            if not fun.is_always_real():
                return False
        
        return True
    
    def get_sign(self):
        sign = cp(self.rat.get_sign())
        
        for func in self.irt:
            try:
                sign *= func.get_sign()
            
            except TypeError:  # if the sign of the function is not known
                return None
        
        return sign
    
    def is_real_and_positive(self):
        return self.get_sign() == cp(1)

    def expanded(self):
        out = cp( (self.rat.num, self.rat.den) )
        
        for fun in self.irt:
            out *= fun.alternate_form()
            
            '''if(cp(fun).numerical() != fun.alternate_form().numerical()):
                print('Incorrectly simplifying ', cp(fun).compileable())
                input()'''
        
        return out
    
    def transformed(self, tr):
        out = cp( (self.rat.num, self.rat.den) )
        
        for fun in self.irt:
            for i in range(len(fun.args)):
                fun.args[i].transform(tr)
            
            out *= tr(fun)
        
        return out    
    
    def differentiated(self, var):
        out = cp(0)
        
        for i in range(len(self.irt)):
            functions = []
            for j in range(len(self.irt)):
                if(i != j):
                    functions.append(self.irt[j])
            
            C = Irrational(Term(self.rat, functions))
            coef = Complex(C, make_irrational(0))
            
            out += self.irt[i].differentiated(var) * coef
        
        return out
    
    def numerical(self, substitution_list):
        out = Numerical(self.rat.num / self.rat.den, 0.0)
        
        for func in self.irt:
            out *= func.numerical(substitution_list)
        
        return out
    
    
class Irrational:
    def __init__(self, *new_terms):  # object of class Irrational is simply a sum
                                     # of a number of class Term objects
        for term in new_terms:
            if not isinstance(term, Term):
                raise TypeError('Wrong type of Irrational.__init__ arguments: expected <class \'Algebraica.Term\'>, got ' + str(type(term)))
            
        self.terms = []
        for i in range(len(new_terms)):
            self.terms.append(new_terms[i].copy())  # to avoid awkward situations with
                                                    # pointers, we need to utilize .copy()
        
        self.sort()
        self.reduce()                                 
        
    def __str__(self):
        return '{\n' + self.show(1) + '\n}'
    
    def wolfram_lang(self):
        out = ''
        
        for i in range(len(self.terms)):
            out += self.terms[i].wolfram_lang()
            
            if(i != len(self.terms) - 1):
                out += ' + '
        
        return out
    
    def compileable(self):
        if(len(self.terms) == 0):
            return '0'
        
        elif(len(self.terms) == 1):
            return self.terms[0].compileable()
        
        else:
            out = '['
            
            for i in range(len(self.terms)):
                out += self.terms[i].compileable()
                
                if(i != len(self.terms) - 1):
                    out += ', '
            
            out += ']'
            
        return out
    
    def __add__(first, second):
        if isinstance(second, Irrational):
            return Irrational(*((first.terms + second.terms).copy()))
        
        elif(isinstance(second, int) or
             isinstance(second, float) or
             isinstance(second, Rational) or
             isinstance(second, Term) or
             isinstance(second, str)):
            irrational_second = make_irrational(second)
            
            return first + irrational_second
        
        else:
            raise TypeError('Wrong Irrational.__add__ argument type: ' + str(type(second)))
    def __iadd__(first, second):
        return first + second
    
    def __sub__(first, second):
        if isinstance(second, Irrational):
            return Irrational(*((first.terms + (-second).terms).copy()))
        
        elif(isinstance(second, int) or
             isinstance(second, float) or
             isinstance(second, Rational) or
             isinstance(second, Term) or
             isinstance(second, str)):
            irrational_second = make_irrational(second)
            
            return first - irrational_second        
        
        else:
            raise TypeError('Wrong Irrational.__sub__ argument type: ' + str(type(second)))
    def __isub__(first, second):
        return first - second
    
    def __mul__(first, second):
        if(isinstance(second, Irrational)):
            
            out = Irrational()
            
            for t1 in first.terms:
                for t2 in second.terms:
                    t = t1 * t2
                    out.terms.append(t)
                
            out.reduce()
            return out
        
        elif(isinstance(second, int) or
             isinstance(second, float) or
             isinstance(second, Rational) or
             isinstance(second, Term) or
             isinstance(second, str)):
            res_terms = copy_list(first.terms)
            
            for i in range(len(first.terms)):
                res_terms[i] *= second
            
            return Irrational(*res_terms)
        
        else:
            raise TypeError('Wrong Irrational.__mul__ argument type: ' + str(type(second)))
    def __imul__(first, second):
        return first * second
    
    def __eq__(first, second):
        if isinstance(second, Irrational):
            first.reduce()
            second.reduce()
            
            if(len(first.terms) == len(second.terms)):
                for i in range(len(first.terms)):
                    if not(first.terms[i] == second.terms[i]):
                        return False
                
                return True
            
            else:
                return False
        
        else:
            raise TypeError('Wrong Irrational.__eq__ argument type: expected Irrational, got ' + str(type(second)))
    
    def __neg__(self):
        res = self.copy()
        
        for i in range(len(res.terms)):
            res.terms[i].change_sign()
        
        res.reduce()
        return res
    
    def __truediv__(first, second):  # note that Irrational cannot be divided by another Irrational
        if (isinstance(second, int) or isinstance(second, float) or isinstance(second, Rational)):
            res = first.copy()
            
            for i in range(len(res.terms)):
                res.terms[i] /= second
            
            res.reduce()
            return res
        
        else:
            raise TypeError('Wrong Irrational.__truediv__ argument type: ' + str(type(second)))
    def __idiv__(first, second):
        return first / second

    def copy(self):
        return Irrational(*self.terms)
    
    def show(self, elevation, ol=False):
        if self.might_be_nonzero():
            out = ''
            
            if(ol):
                out = ' ' * TabulationStep * elevation
        
            for i in range(len(self.terms)):
                if(ol):
                    out += self.terms[i].show(0, one_line=True)
                else:
                    out += self.terms[i].show(elevation)
                    
                if(i != len(self.terms) - 1):
                    if ColorfulPrint:
                        if(ol):
                            out += colorama.Style.BRIGHT + colorama.Fore.BLACK + ' + ' + colorama.Style.RESET_ALL
                            
                        else:
                            out += colorama.Style.BRIGHT + colorama.Fore.BLACK + ' +\n' + colorama.Style.RESET_ALL
                    
                    else:
                        if(ol):
                            out += ' + '
                            
                        else:
                            out += ' +\n'
            
            return out
        
        else:
            return ' ' * elevation * TabulationStep + '0'
    
    def might_be_nonzero(self):
        for term in self.terms:
            if(term.rat.num != 0):
                return True
        
        return False
    def is_zero(self):
        return not self.might_be_nonzero()
    
    def has_functions(self):
        if self.is_zero():
            return False
        
        for term in self.terms:
            if term.has_functions():
                return True
            
        return False
    
    def is_constant(self):
        for term in self.terms:
            if not term.is_constant():
                return False
        
        return True
    
    def is_always_real(self):
        for term in self.terms:
            if not term.is_always_real():
                return False
        
        return True
    
    def is_real_and_positive(self):
        if self.is_zero():
            return False
        
        for term in self.terms:
            if not term.is_real_and_positive():
                return False
        
        return True
    def is_real_and_negative(self):
        return (-self).is_real_and_positive()
    
    def is_rational(self):
        if(len(self.terms) == 0):
            return True
        
        elif(len(self.terms) == 1):
            if(len(self.terms[0].irt) == 0):
                return True
            
            else:
                return False
            
        else:
            return False
        
    def get_sign(self):
        if(len(self.terms) == 1):
            return self.terms[0].get_sign()
        
        else:
            return None
    
    def sort(self):
        self.terms.sort(key=Term.serial_index)
    
    def reduce(self):
        if(len(self.terms) > 1):
            newterms = []
                     
            while(len(self.terms) > 1):
                term = self.terms.pop(0)
                if not term.might_be_nonzero():
                    continue

                new_rat = term.rat
                while((len(self.terms) > 0) and term.eq_irt(self.terms[0])):
                    new_rat += self.terms.pop(0).rat
                
                new_rat.reduce()
                newterms.append(Term(new_rat, term.irt))
                
            if(len(self.terms) > 0):
                if(self.terms[0].might_be_nonzero()):
                    newterms.append(self.terms[0])
            
            self.terms = copy_list(newterms)
        
        elif(len(self.terms) == 1):
            if(self.terms[0].might_be_nonzero()):
                pass
            else:
                self.terms = []
            
    def reduced(self):
        res = self.copy()
        res.reduce()
        
        return res
    
    def differentiated(self, var):
        out = cp(0)
        
        for term in self.terms:
            out += term.differentiated(var)
        
        return out
    
    def numerical(self, substitution_list):
        out = Numerical(0.0, 0.0)
        
        for term in self.terms:
            out += term.numerical(substitution_list)
        
        return out
    

class Complex:
    def __init__(self, new_real, new_imag, new_is_simplified=False):  # object of class Complex is just M + N i,
                                             # where M and N are class Irrational objects and
                                             # i is the square root of negative one
        if not(isinstance(new_real, Irrational) and
               isinstance(new_imag, Irrational)):
            
            raise TypeError('Wrong Complex.__init__ arguments` type: ' + str(type(new_real)) + ', ' + str(type(new_imag)))
        
        self.rl = new_real.copy()  # name "rl" is misleading, since rl is often a complex number.
        self.im = new_imag.copy()  # since "im" is also complex, self.im is not always the imaginary part of self.
        self.is_simplified = new_is_simplified
        
    def __str__(self):
        if AutoSimplify:
            self.simplify()
            
        return '<<<\n' + self.show(1) + '\n>>>'
    
    def wolfram_lang(self):
        if AutoSimplify:
            self.simplify()        
        
        if(self.rl.is_zero()):
            if(self.im.is_zero()):
                return '0'
            
            else:
                return '(' + self.im.wolfram_lang() + ') I'
        
        else:
            if(self.im.is_zero()):
                return self.rl.wolfram_lang()
            
            else:
                return self.rl.wolfram_lang() + ' + ' +\
                       '(' + self.im.wolfram_lang() + ') I'
    
    def compileable(self):
        return 'cp( ' + self.rl.compileable() + ', ' + self.im.compileable() + ' )'
    
    def __add__(first, second):
        if isinstance(second, Complex):
            
            return Complex(
                first.rl + second.rl,
                first.im + second.im
            )
        
        elif(isinstance(second, int) or
             isinstance(second, float)):
            
            return Complex(
                first.rl + second,
                first.im
            )
        
        else:
            raise TypeError('Wrong Complex.__add__ argument type: ' + str(type(second)))
    def __iadd__(first, second):
        return Complex.__add__(first, second)
    
    def __mul__(first, second):
        if(isinstance(second, Complex)):
            return Complex(first.rl * second.rl - first.im * second.im,
                           first.rl * second.im + first.im * second.rl)
        
        elif(isinstance(second, int) or isinstance(second, float) or isinstance(second, Rational)):
            return Complex(first.rl * second,
                          first.im * second)
        
        else:
            # if the second input was not recognized,
            # it makes sense to try to make it Complex
            
            if(isinstance(second, tuple)):
                return first * cp(*second)
            
            else:
                return first * cp(second)
    def __imul__(self, other):
        return Complex.__mul__(self, other)
    
    def __sub__(first, second):
        if isinstance(second, Complex):
            
            return Complex(
                first.rl - second.rl,
                first.im - second.im
            )
        
        elif(isinstance(second, int) or
             isinstance(second, float)):
            
            return Complex(
                first.rl - second,
                first.im
            )
        
        else:
            raise TypeError('Wring Complex.__sub__ argument type: ' + str(type(second)))
    def __isub__(first, second):
        return Complex.__sub__(first, second)
    
    def __invert__(self):  # turns (a + b i) to (a - b i)
        if(self.rl.is_always_real() and self.im.is_always_real()):
            return Complex(
                self.rl,
                -self.im,
                self.is_simplified
            )
        
        else:
            return Complex(
                make_irrational( ('Rl', self) ),
                -make_irrational( ('Im', self) ),
                self.is_simplified
            )
    
    def __eq__(first, second):
        if isinstance(second, Complex):
            return ((first.rl == second.rl) and
                    (first.im == second.im))
        
        elif isinstance(second, type(None)):
            return False
        
        else:
            raise TypeError('Wrong Complex.__eq__ argument type: ' + str(type(second)))
    def __ne__(first, second):
        return not Complex.__eq__(first, second)
    
    def __truediv__(first, second):
        if isinstance(second, Complex):
            C = Complex(
                irrational((1, 1), [(pow_function_name, second * (~second), cp(-1))]),
                irrational((0, 1), [])
            )
            
            return (first * (~second)) * C
        
        elif(isinstance(second, int) or isinstance(second, float) or isinstance(second, Rational)):
            return Complex(
                first.rl / second,
                first.im / second,
                first.is_simplified
            )
        
        else:
            raise TypeError('Wrong Complex.__truediv__ argument type: ' + str(type(second)))
    
    def __neg__(self):
        return Complex(-self.rl, -self.im, self.is_simplified)
    
    def __pow__(first, second):
        if(isinstance(second, int)):
            if(second == 0):
                return cp(1)
            
            elif(second > 0):
                out = cp(1)
                
                for i in range(second):
                    out *= first
                
                return out
            
            elif(second < 0):
                return cp( (-1, (pow_function_name, first, cp(-second))) )
            
        elif(isinstance(second, Complex)):
            abs_val = abs(first)
            cmp_arg = first.arg()
            
            ln_first = cp(
                (1, [(ln_function_name, abs_val)])
            ) + cp(0, 1) * cmp_arg
            
            return cp((1, [(exp_function_name, ln_first * second)]), 0)
        
        else:
            raise TypeError('Wrong Complex.__pow__ argument type: ' + str(type(second)))
        
    def __abs__(self):
        return cp( (1, (abs_function_name, self)) )
    
    def has_functions(self):
        return self.rl.has_functions() or self.im.has_functions()
    
    def is_constant(self):
        return self.rl.is_constant() and self.im.is_constant()
    
    def is_always_real(self):
        return (self.rl.is_always_real() and
                self.im.is_zero())    
    
    def is_real_and_positive(self):
        return (self.rl.is_real_and_positive() and
                self.im.is_zero())
    def is_real_and_negative(self):
        return (self.rl.is_real_and_negative() and
                self.im.is_zero())
    
    def is_even(self):
        return (self.is_real_int() and (self.rl.terms[0].rat.num % 2 == 0))
    
    def is_real_rational(self):
        if not self.is_always_real():
            return False
        
        if(len(self.rl.terms) == 0):
            return True
        elif(len(self.rl.terms) == 1):
            if(len(self.rl.terms[0].irt) == 0):
                return True
            
            else:
                return False
            
        else:
            return False
    
    def is_real_int(self):
        if not self.is_always_real():
            return False
        
        if(len(self.rl.terms) == 0):
            return True
        elif(len(self.rl.terms) == 1):
            if(len(self.rl.terms[0].irt) == 0):
                if(self.rl.terms[0].rat.den == 1):
                    return True
                
                else:
                    return False
            
            return False
            
        else:
            return False
    
    def is_always_imag(self):
        return self.rl.is_zero() and self.im.is_always_real()
    
    def might_be_nonzero(self):
        return self.rl.might_be_nonzero() or self.im.might_be_nonzero()
    def is_zero(self):
        return not self.might_be_nonzero()
    
    def get_sign(self):
        if(self.rl.is_zero()):
            if(self.im.is_zero()):
                return cp(0)
            
            else:
                # sign(i z) = i sign(z)
                return cp( 0, Function( (sign_function_name, cp(self.im)) ) )
        
        else:
            if(self.im.is_zero()):
                if(len(self.rl.terms) == 1):
                    return self.rl.terms[0].get_sign()
            
            else:
                return None

    def arg(self):
        if self.rl.is_zero():
            if self.im.is_zero():
                return cp(0)
            
            else:
                return cp( ((1, 2), [pi_constant_name, (sign_function_name, cp( (1, (im_function_name, self)) ) )]) )
        
        elif self.im.is_zero():
            return cp(0)
        
        else:
            return cp( (1, (arg_function_name, self)) )
    
    def show(self, elevation, one_line=False):
        out = ''
        tab = ' ' * elevation * TabulationStep
        
        rl_nz = self.rl.might_be_nonzero()
        im_nz = self.im.might_be_nonzero()
        
        if(rl_nz):
            if(self.rl.is_constant()):
                if(im_nz):
                    out += tab
                    
                    if ColorfulPrint:
                        out += colorama.Fore.BLACK + '{ ' + colorama.Style.RESET_ALL
                    else:
                        out += '{ '
                        
                    out += self.rl.show(0, ol=True)
                    
                    if ColorfulPrint:
                        out += colorama.Fore.BLACK + ' }' + colorama.Style.RESET_ALL
                    else:
                        out += ' }'
                    
                else:
                    out += self.rl.show(elevation, ol=True)
            
            else:
                if(im_nz):
                    if not one_line:
                        out += tab
                        
                        if ColorfulPrint:
                            out += colorama.Fore.BLACK + '{\n' + colorama.Style.RESET_ALL
                        else:
                            out += '{\n'
                    
                    if(one_line):
                        out += self.rl.show(0, ol=one_line)
                    else:
                        out += self.rl.show(elevation + 1, ol=one_line)
                else:
                    out += self.rl.show(elevation, ol=one_line)
                
                if(im_nz):
                    if not one_line:
                        out += '\n' + tab 
                        
                        if ColorfulPrint:
                            out += colorama.Fore.BLACK + '}' + colorama.Style.RESET_ALL
                        else:
                            out += '}'
            
        if (rl_nz and im_nz):
            if one_line:
                if ColorfulPrint:
                    out += colorama.Fore.BLACK + ' + ' + colorama.Style.RESET_ALL
                else:
                    out += ' + '
            else:
                out += '\n'
                out += tab
                if ColorfulPrint:
                    out += colorama.Fore.BLACK + ' + ' + colorama.Style.RESET_ALL
                else:
                    out += ' + '                
                out += '\n'        
        
        if(im_nz):
            if(self.im.is_constant()):
                out += tab
                if ColorfulPrint:
                    out += colorama.Fore.BLACK + '{ ' + colorama.Style.RESET_ALL
                else:
                    out += '{ '
                
                out += self.im.show(0, ol=True)
                
                if ColorfulPrint:
                    out += colorama.Fore.BLACK + ' } ' + colorama.Style.RESET_ALL
                else:
                    out += ' } '
                
                if(ColorfulPrint):
                    out += colorama.Fore.CYAN + 'i' + colorama.Style.RESET_ALL
                
                else:
                    out += 'i'
                
            else:
                if one_line:
                    if ColorfulPrint:
                        out += colorama.Fore.BLACK + '{ ' + colorama.Style.RESET_ALL
                    else:
                        out += '{ '
                    out += self.im.show(elevation + 1, ol=True)
                    if ColorfulPrint:
                        out += colorama.Fore.BLACK + ' } ' + colorama.Style.RESET_ALL
                    else:
                        out += ' } '
                    
                    if(ColorfulPrint):
                        out += colorama.Fore.CYAN + 'i' + colorama.Style.RESET_ALL
                    
                    else:
                        out += 'i'                    
                else:
                    out += tab
                    if ColorfulPrint:
                        out += colorama.Fore.BLACK + '{\n' + colorama.Style.RESET_ALL
                    else:
                        out += '{\n'
                    out += self.im.show(elevation + 1, ol=False)
                    out += '\n' + tab
                    if ColorfulPrint:
                        out += colorama.Fore.BLACK + '} ' + colorama.Style.RESET_ALL
                    else:
                        out += '} '
                    
                    if(ColorfulPrint):
                        out += colorama.Fore.CYAN + 'i' + colorama.Style.RESET_ALL
                    
                    else:
                        out += 'i'                    
            
        if not(rl_nz or im_nz):
            out += self.rl.show(elevation, ol=one_line)
        
        return out
    
    def transformed_once(self, tr):
        new_rl = cp(0)
        new_im = cp(0)
        
        for term in self.rl.terms:
            new_rl += term.transformed(tr)
            
        for term in self.im.terms:
            new_im += term.transformed(tr)
        
        return new_rl + new_im * cp(0, 1)        
    def transformed(self, tr):
        old = self
        res = self.transformed_once(tr)

        while(res != old):
            old = res.copy()
            res = res.transformed_once(tr)
        
        return res
    def transform(self, tr):
        res = self.transformed(tr)
        
        self.rl = res.rl.copy()
        self.im = res.im.copy()
        self.is_simplified = res.is_simplified
    
    def copy(self):
        return Complex(self.rl,
                       self.im,
                       new_is_simplified = self.is_simplified)
    
    def set_equal_to(self, other):
        self.rl = other.rl
        self.im = other.im
        self.is_simplified = other.is_simplified
    
    def collect_powers(self):
        new_rl = cp(0)
        new_im = cp(0)
        
        for term in self.rl.terms:
            new_rl += term.collect_powers()
            
        for term in self.im.terms:
            new_im += term.collect_powers()
        
        self.set_equal_to(new_rl + new_im * cp(0, 1))
    
    def simplified_once(self):
        self.collect_powers()
        
        new_rl = cp(0)
        new_im = cp(0)
        
        for term in self.rl.terms:
            term.collect_exponentials()
            term.collect_trig()
            term.collect_htrig()
            new_rl += term.expanded()
            
        for term in self.im.terms:
            term.collect_exponentials()
            term.collect_trig()
            term.collect_htrig()
            new_im += term.expanded()
        
        return new_rl + new_im * cp(0, 1)
    def simplified(self):
        if(self.is_simplified):
            return self
        
        else:
            old = self
            res = self.simplified_once()

            while(res != old):
                old = res.copy()
                res = res.simplified_once()
            
            res.is_simplified = True
            return res
    def simplify(self):
        self.set_equal_to(self.simplified())
        
    def differentiated(self, var):
        if AutoSimplify:
            self.simplify()
        
        return self.rl.differentiated(var) + self.im.differentiated(var) * cp(0, 1)
    
    def differentiate(self, var):
        D = self.differentiated(var)
        
        self.rl = D.rl.copy()
        self.im = D.im.copy()
        self.is_simplified = False
        
    def numerical(self, *args):
        if AutoSimplify:
            self.simplify()
            
        def check_sb_list(substitution_list):
            if not isinstance(substitution_list, list):
                raise TypeError('Substitution list must be a list')
                
            for item in substitution_list:
                if not isinstance(item, tuple):
                    raise TypeError('Elements of substitution list must be tuples')
                
                if(len(item) != 2):
                    raise TypeError('Elements of substitution list must have length 2')
                
                if not(isinstance(item[0], Function) and
                       isinstance(item[1], Numerical)):
                    raise TypeError('Elements of substitution list must be of type (Function, Numerical)')
        
        if(len(args) == 0):
            substitution_list = []
        
        elif(len(args) == 1):
            substitution_list = args[0]
        
        elif(len(args) == 2):
            if(isinstance(args[0], Function) and
               isinstance(args[1], Numerical)):
                substitution_list = [args]
            
            else:
                raise TypeError('Complex.numerical can take two arguments that are of type Function and Numerical, not ' + str(type(args[0])) + ' and ' + str(type(args[1])))
        
        else:
            raise TypeError('Complex.numerical takes only up to two additional arguments, while the number of arguments given is ' + str(len(args)))
        
        check_sb_list(substitution_list)
        return self.rl.numerical(substitution_list) + self.im.numerical(substitution_list) * Numerical(0.0, 1.0)
    


def copy_list(lst):  # this function is made to copy a multi-dimensional list
                     # it is simply a recursive version of .copy()
    out = []
    
    for item in lst:
        if(type(item) == type([])):
            out.append(copy_list(item))
            
        else:
            try:
                out.append(item.copy())
            except Exception:
                out.append(item)
    
    return out

def expand_plist(primes, search_range):
    min_p = max(primes) + 1
    max_p = min_p + search_range
    
    print('Finding primes between ' + str(min_p) + ' and ' + str(max_p))
    primes += list(range(min_p, max_p))  # min_p >= 2
    
    i = 0
    while(i < len(primes)):
        change = False
        
        prime = primes[i]
        n = 2 * prime
        j = 0
        while True:
            if(j >= len(primes)):
                break
            
            while(primes[j] > n):
                n += prime
            
            if(primes[j] == n):
                del primes[j]
                change = True
            else:
                j += 1
        
        if not change:
            i += 1
    
global prime_list
prime_list = [2]
print('Loading primes')
expand_plist(prime_list, InitialMaxPrime)

def prime_factors(n):  # note that this function may also return non-primes
                       # if n is too large
    
    global prime_list
    factors = []
    
    abs_n = abs(n)
    i = 0
    while True:
        if(i > len(prime_list) - 1):
            if(prime_list[-1] > HighestAllowedPrime):
                break
            else:
                expand_plist(prime_list, PrimeNumberSearchStep)
        prime = prime_list[i]
            
        if(prime > abs_n):
            break
        
        power = 0
        while(abs_n % prime == 0):
            abs_n = abs_n // prime
            power += 1
        
        if(power > 0):
            factors.append( (prime, power) )
        
        i += 1
            
    if(n < 0):
        factors.append( (-1, 1) )
    if(abs_n != 1):
        factors.append( (abs_n, 1) )
    
    return factors
    
def irrational(frac, raw_irrational_stuff):  # Irrational.__init__ is pretty inconvenient,
                                             # so that`s why this shortcut was created
    if(isinstance(frac[0], int) and
       isinstance(frac[1], int)):
        irrational_stuff = []
        for item in raw_irrational_stuff:
            irrational_stuff.append(Function(item))
        
        return Irrational(Term(Rational(frac[0], frac[1]), irrational_stuff)).reduced()

    else:
        raise TypeError
    
    
def make_irrational(*inp):  # you, user, want things to be easy to use and convenient
                            # constructing Irrational by __init__ would be extremely inconvenient,
                            # so that`s why make_irrational function was created
                            # it simply considers many possible ways you might want to
                            # express an irrational number
    
    if(len(inp) == 1):
        a = inp[0]
        
        if(isinstance(a, int)):  # if OUT is a real integer
                                 # such as a = 7
                                 # OUT = 7/1 + 0/1 i
            return irrational((a, 1), [])
        
        elif(isinstance(a, float)):  # if a is real
                                     # for example, a = 1.6
                                     # OUT = 8/5 + 0/1 i
            return irrational(float_to_frac(a), [])
        
        elif(isinstance(a, str)):   # if OUT is a constant
                                    # like a = pi_constant_name
                                    # OUT = 1/1 * pi + 0/1 i
            return irrational((1, 1), [a])
        
        elif(isinstance(a, Rational)):
            return irrational((a.num, a.den), [])        
        
        elif(isinstance(a, Function)):
            return irrational((1, 1), [tuple([a.name] + list(a.args))])
        
        elif(isinstance(a, Term)):
            return Irrational(a)
        
        elif(isinstance(a, Irrational)):
            return a
        
        elif(isinstance(a, tuple)):
            if(isinstance(a[0], str)):  # if OUT is a tuple denoting function,
                                        # such as a = (pow_function_name, cp(3), cp(1/2))
                                        # OUT = (3/1 + 0/1 i) ** (1/2 + 0/1 i)
                                        # then a[0] is the name of the function
                                        #  and a[1:] are its` arguments
                for i in range(1, len(a)):
                    if not isinstance(a[i], Complex):
                        raise TypeError('Tuple denoting a Function has entries of a wrong type: expected Complex, got ' + str(type(a[i])))
                
                return irrational((1, 1), [a])
            
            else:
                if(len(a) == 2):
                    if(isinstance(a[0], int) and isinstance(a[1], int)):  # if a is a tuple denoting fraction
                                                                          # like a = (2, 3)
                                                                          # OUT = 2/3 + 0/1 i
                        return irrational((a[0], a[1]), [])
                    
                    elif(isinstance(a[0], int) and isinstance(a[1], str)):    # if a is a tuple representing a constant with an integer coefficient
                                                                              # such as a = (6, 'e')
                                                                              # OUT = 6 * e + 0 i
                        return irrational((a[0], 1), [(a[1],)])                    
                    
                    elif(isinstance(a[0], tuple) and isinstance(a[1], str)):  # if a is a tuple representing a constant with a rational coefficient
                                                                              # such as a = ((2, 3), 'pi')
                                                                              # OUT = 2/3 * pi + 0 i
                        return irrational((a[0][0], a[0][1]), [(a[1],)])
                    
                    elif(isinstance(a[0], int) and isinstance(a[1], tuple)):
                                             # if a is a tuple describing a function with an integer coefficient
                                             # for example, a = (2, ('sqrt', cp(5)))
                                             # OUT = 2 * sqrt[5] + 0 i
                        if(isinstance(a[1][0], str)):
                            for i in range(1, len(a[1])):
                                if not isinstance(a[1][i], Complex):
                                    raise TypeError('Tuple denoting arguments of a Function has wrong type of entries: expected Complex, got ' + str(type(a[1][i])))
                        else:
                            raise TypeError('Tuple denoting a Function does not contain the name of the function')
                        
                        return irrational((a[0], 1), [a[1]])                    
                    
                    elif(isinstance(a[0], tuple) and isinstance(a[1], tuple)):
                                             # if a is a tuple describing a function with a rational coefficient
                                             # for example, a = ((4, 9), ('ln', cp(7)))
                                             # OUT = 4/9 * ln[7] + 0 i
                        if(isinstance(a[1][0], str)):
                            for i in range(1, len(a[1])):
                                if not isinstance(a[1][i], Complex):
                                    raise TypeError('Tuple denoting arguments of a Function has wrong type of entries: expected Complex, got ' + str(type(a[1][i])))
                        else:
                            raise TypeError('Tuple denoting a Function does not contain the name of the function')
                        
                        return irrational((a[0][0], a[0][1]), [a[1]])
                    
                    elif (isinstance(a[0], int) and isinstance(a[1], list)):
                                             # if a is a tuple representing a product of functions with an integer coefficient
                                             # for instance, a = (4, [pi_constant_name, (sqrt_function_name, cp(2))])
                                             # OUT = 4 * pi * sqrt[2] + 0/1 i
                        return irrational((a[0], 1), a[1])           
                    
                    elif(isinstance(a[0], tuple) and isinstance(a[1], list)):
                                             # if a is a tuple representing a product of functions with a rational coefficient
                                             # for instance, a = ((2, 3), [pi_constant_name, (sqrt_function_name, cp(2))])
                                             # OUT = 2/3 * pi * sqrt[2] + 0/1 i
                        return irrational((a[0][0], a[0][1]), a[1])
                    
                    else:
                        raise ValueError('make_irrational cannot interpret the tuple given: ' + str(a))
                
                else:
                    raise ValueError('Wrong make_irrational tuple argument length: ' + str(len(a)))
        
        elif(isinstance(a, list)):  # if a is a list, then each element of the
                                    # list represents a term of a sum
                                    # each of these terms is to be read by conventions
            out = irrational((0, 1), [])
            
            for term in a:
                out += make_irrational(term)
                
            return out
    
        else:
            raise TypeError('Wrong make_irrational argument type: ' + str(type(a)))
        
    else:
        raise ValueError('Wrong number of make_irrational arguments: ' + str(len(inp)))
    
    
def float_to_frac(f):
    N = f
    D = 1
    while((N != int(N)) and (D < 10**6)):
        N *= 10
        D *= 10
        
    return int(N), int(D)

def compare_lists(A, B):
    if(len(A) != len(B)):
        return False
    
    for i in range(len(A)):
        if(type(A[i]) != type(B[i])):
            return False
        
        if(isinstance(A[i], list)):
            if not(compare_lists(A[i], B[i])):
                return False
        else:
            if (A[i] != B[i]):
                return False
    # this function compares two multi-dimensional lists
    # the two lists are considered equal, if all their elements match
    # for example, [1, [2, 3]] == [1, [2, 3]]
    #              [1, [2, 4]] != [1, [2, 4]]
    #              [1, [2, 3]] != [1]
    
    return True

def common_divider(rats):
    pset = set()
    for rat in rats:
        pf = rat.prime_factors()
        for item in pf:
            pset.add(item[0])
    
    num = 1
    den = 1
    for prime in pset:
        lst = []
        
        power = 10**10
        for rat in rats:
            pf = rat.prime_factors()
            
            pw = 0
            for p in pf:
                if(p[0] == prime):
                    pw = p[1]
            
            if(abs(pw) < abs(power)):
                power = pw
        
        if(power > 0):
            num *= prime**power
        if(power < 0):
            den *= prime**(-power)
    
    return Rational(num, den)

def sign(n): # n is a float or an integer
    if(isinstance(n, float) or isinstance(n, int)):
        if(n == 0):
            return 0
        
        else:
            return n / abs(n)
    
    else:
        raise TypeError('Wrong data type: ' + str(type(n)))
    
def Str(n, length): # n is a float or an int,
               # length is an int that determines the length of the output
    if isinstance(n, int):
        n = float(n)
    
    elif isinstance(n, float):
        pass
    
    else:
        raise TypeError('Wrong input type: ' + str(type(n)))
    
    if(n == 0):
        lg = 0
    else:
        lg = int(numpy.log(abs(n)) // numpy.log(10))
    coef = abs(n) / 10**lg
    
    
    if(lg == 0):
        str_coef_len = length
        if(n < 0):
            str_coef_len -= 1
        if(str_coef_len <= 3):
            str_coef_len = 6
    
    else:
        str_coef_len = length - 1
        if(n < 0):
            str_coef_len -= 1
        str_coef_len -= len(str(lg))
        if(str_coef_len <= 3):
            str_coef_len = 6
    
    str_coef = str(coef)
    if(len(str_coef) == str_coef_len):
        pass
    if(len(str_coef) > str_coef_len):
        str_coef = str_coef[:str_coef_len]
    if(len(str_coef) < str_coef_len):
        str_coef += '0' * (str_coef_len - len(str_coef))
    
    out = ''
    if(n < 0):
        if ColorfulPrint:
            out += colorama.Fore.CYAN + '-' + colorama.Style.RESET_ALL
        else:
            out += '-'
        
    if ColorfulPrint:
        out += colorama.Fore.CYAN + str_coef + colorama.Style.RESET_ALL
    else:
        out += str_coef
    
    if(lg != 0):
        if ColorfulPrint:
            out += colorama.Style.NORMAL + colorama.Fore.CYAN + 'e' + colorama.Style.BRIGHT + colorama.Fore.CYAN + str(lg) + colorama.Style.RESET_ALL
        else:
            out += 'e' + str(lg)
    
    return out


def cpower(b, p): # b is a class Numerical object
    if isinstance(b, Numerical):
        return b ** p
    
    else:
        raise TypeError('Wrong data type: ' + str(type(b)))
    
def csqrt(n): # n is a class Numerical object
    if isinstance(n, Numerical):
        return n ** 0.5
    
    else:
        raise TypeError('Wrong data type: ' + str(type(n)))    
    
def cexp(n): # n is a class Numerical object
    if isinstance(n, Numerical):
        amplitude = numpy.exp(n.rl)
        
        return Numerical(
            amplitude * numpy.cos(n.im),
            amplitude * numpy.sin(n.im)
        )
    else:
        raise TypeError('Wrong data type: ' + str(type(n)))

def cln(n): # n is a class Numerical object
    if isinstance(n, Numerical):
        return Numerical(numpy.log(abs(n)), n.arg())
    
    else:
        raise TypeError('Wrong data type: ' + str(type(n)))
    
def clog(a, b): # a and b are class Numerical objects
                # returns log of a with base b, Ln(a) / Ln(b)
    if(isinstance(a, Numerical) and isinstance(b, Numerical)):
        return cln(a) / cln(b)
    
    else:
        raise TypeError('Wrong arguments type: ' + str(type(a)) + ', ' + str(type(b)))
    
def csign(n): # n is a class Numerical object
    if isinstance(n, Numerical):
        return n / abs(n)
    
    else:
        raise TypeError('Wrong data type: ' + str(type(n)))
    
def csin(n): # n is a class Numerical object
    if isinstance(n, Numerical):
        return Numerical(
                    numpy.sin(n.rl) * numpy.cosh(n.im),
                    numpy.cos(n.rl) * numpy.sinh(n.im)
                )
    
    else:
        raise TypeError('Wrong data type: ' + str(type(n)))
    
def ccos(n): # n is a class Numerical object
    if isinstance(n, Numerical):
        return Numerical(
                            numpy.cos(n.rl) * numpy.cosh(n.im),
                            -numpy.sin(n.rl) * numpy.sinh(n.im)
                        )
    
    else:
        raise TypeError('Wrong data type: ' + str(type(n)))
    
def ccsc(n): # n is a class Numerical object
    if isinstance(n, Numerical):
        return Numerical(
                    numpy.sin(n.rl) * numpy.cosh(n.im),
                    numpy.cos(n.rl) * numpy.sinh(n.im)
                ).reciprocal()
    
    else:
        raise TypeError('Wrong data type: ' + str(type(n)))
    
def csec(n): # n is a class Numerical object
    if isinstance(n, Numerical):
        return Numerical(
                            numpy.cos(n.rl) * numpy.cosh(n.im),
                            -numpy.sin(n.rl) * numpy.sinh(n.im)
                        ).reciprocal()
    
    else:
        raise TypeError('Wrong data type: ' + str(type(n)))    
    
def ctan(n): # n is a class Numerical object
    if isinstance(n, Numerical):
        sa = numpy.sin(n.rl)
        ca = numpy.cos(n.rl)
        sb = numpy.sinh(n.im)
        cb = numpy.cosh(n.im)
        
        coef = (ca ** 2 * cb ** 2 +
                sa ** 2 * sb ** 2)
        
        return Numerical(
            sa * ca / coef,
            sb * cb / coef
        )
    
    else:
        raise TypeError('Wrong data type: ' + str(type(n)))
    
def ccot(n): # n is a class Numerical object
    if isinstance(n, Numerical):
        sa = numpy.sin(n.rl)
        ca = numpy.cos(n.rl)
        sb = numpy.sinh(n.im)
        cb = numpy.cosh(n.im)
        
        coef = (sa ** 2 * cb ** 2 +
                ca ** 2 * sb ** 2)
        
        return Numerical(
            sa * ca / coef,
            -sb * cb / coef
        )
    
    
    else:
        raise TypeError('Wrong data type: ' + str(type(n)))    

def carcsin(n): # n is a class Numerical object
    if isinstance(n, Numerical):
        i = Numerical( 0.0, 1.0 )
        
        return -i * cln(i * n + csqrt(-n**2 + 1))
    
    else:
        raise TypeError('Wrong data type: ' + str(type(n)))    

def carccos(n): # n is a class Numerical object
    if isinstance(n, Numerical):
        return Numerical(numpy.pi/2, 0) - carcsin(n)
    
    else:
        raise TypeError('Wrong data type: ' + str(type(n)))
    
def carcsec(n): # n is a class Numerical object
    if isinstance(n, Numerical):
        return carccos(n.reciprocal())
    
    else:
        raise TypeError('Wrong data type: ' + str(type(n)))    

def carccsc(n): # n is a class Numerical object
    if isinstance(n, Numerical):
        return carcsin(n.reciprocal())
        
    else:
        raise TypeError('Wrong data type: ' + str(type(n)))

def carctan(n): # n is a class Numerical object
    if isinstance(n, Numerical):
        return Numerical(
            Numerical(
                1 - n.rl**2 - n.im**2,
                2 * n.rl
            ).arg() / 2,
            numpy.log(n.rl**2 + (n.im + 1)**2) / 2 -
            numpy.log(4 * n.rl**2 + (n.rl**2 + n.im**2 - 1)**2) / 4
        )
    
    else:
        raise TypeError('Wrong data type: ' + str(type(n)))    

def carccot(n): # n is a class Numerical object
    if isinstance(n, Numerical):
        return Numerical(
            Numerical(
                n.rl**2 + n.im**2 - 1,
                2 * n.rl
            ).arg() / 2,
            numpy.log(n.rl**2 + (n.im - 1)**2) / 2 -
            numpy.log(4 * n.rl**2 + (n.rl**2 + n.im**2 - 1)**2) / 4
        )
    
    else:
        raise TypeError('Wrong data type: ' + str(type(n)))

def csinh(n): # n is a class Numerical object
    if isinstance(n, Numerical):
        return Numerical(
            numpy.sinh(n.rl) * numpy.cos(n.im),
            numpy.cosh(n.rl) * numpy.sin(n.im)
            )
    
    else:
        raise TypeError('Wrong data type: ' + str(type(n)))
    
def ccosh(n): # n is a class Numerical object
    if isinstance(n, Numerical):
        return Numerical(
            numpy.cosh(n.rl) * numpy.cos(n.im),
            numpy.sinh(n.rl) * numpy.sin(n.im)
        )
    
    else:
        raise TypeError('Wrong data type: ' + str(type(n)))
    
def csech(n): # n is a class Numerical object
    if isinstance(n, Numerical):
        return ccosh(n).reciprocal()
    
    else:
        raise TypeError('Wrong data type: ' + str(type(n)))
    
def ccsch(n): # n is a class Numerical object
    if isinstance(n, Numerical):
        return csinh(n).reciprocal()
    
    else:
        raise TypeError('Wrong data type: ' + str(type(n)))    
    
def ctanh(n): # n is a class Numerical object
    if isinstance(n, Numerical):
        sa = numpy.sinh(n.rl)
        ca = numpy.cosh(n.rl)
        sb = numpy.sin(n.im)
        cb = numpy.cos(n.im)
        
        coef = (ca ** 2 * cb ** 2 +
                sa ** 2 * sb ** 2)
        
        return Numerical(
            sa * ca / coef,
            sb * cb / coef
        )        
    
    else:
        raise TypeError('Wrong data type: ' + str(type(n)))
    
def ccoth(n): # n is a class Numerical object
    if isinstance(n, Numerical):
        sa = numpy.sinh(n.rl)
        ca = numpy.cosh(n.rl)
        sb = numpy.sin(n.im)
        cb = numpy.cos(n.im)
        
        coef = (sa ** 2 * cb ** 2 +
                ca ** 2 * sb ** 2)
        
        return Numerical(
            sa * ca / coef,
            -sb * cb / coef
        )
    
    else:
        raise TypeError('Wrong data type: ' + str(type(n)))
    
def carcsinh(n): # n is a class Numerical object
    if isinstance(n, Numerical):
        return cln(n + csqrt(n**2 + 1))
    
    else:
        raise TypeError('Wrong data type: ' + str(type(n)))
    
def carccosh(n): # n is a class Numerical object
    if isinstance(n, Numerical):
        return cln(n + csqrt(n + 1) * csqrt(n - 1))
    
    else:
        raise TypeError('Wrong data type: ' + str(type(n)))
    
def carcsech(n): # n is a class Numerical object
    if isinstance(n, Numerical):
        return carccosh(n.reciprocal())
    
    else:
        raise TypeError('Wrong data type: ' + str(type(n)))
    
def carccsch(n): # n is a class Numerical object
    if isinstance(n, Numerical):
        return carcsinh(n.reciprocal())
    
    else:
        raise TypeError('Wrong data type: ' + str(type(n)))    
    
def carctanh(n): # n is a class Numerical object
    if isinstance(n, Numerical):
        return ( cln(n + 1) - cln(-n + 1) ) / 2
    
    else:
        raise TypeError('Wrong data type: ' + str(type(n)))
    
def carccoth(n): # n is a class Numerical object
    if isinstance(n, Numerical):
        return (cln(n.reciprocal() + 1) - cln(-n.reciprocal() + 1)) / 2
    
    else:
        raise TypeError('Wrong data type: ' + str(type(n)))    


def set_names_to_secure():
    global e_constant_name
    global pi_constant_name
    global tau_constant_name
    global sin_function_name
    global cos_function_name
    global sec_function_name
    global csc_function_name
    global tan_function_name
    global cot_function_name
    global sinh_function_name
    global cosh_function_name
    global sech_function_name
    global csch_function_name
    global tanh_function_name
    global coth_function_name
    global arcsin_function_name
    global arccos_function_name
    global arcsec_function_name
    global arccsc_function_name
    global arctan_function_name
    global arccot_function_name
    global arcsinh_function_name
    global arccosh_function_name
    global arcsech_function_name
    global arccsch_function_name
    global arctanh_function_name
    global arccoth_function_name
    global ln_function_name
    global log_function_name
    global pow_function_name
    global sqrt_function_name
    global exp_function_name
    global rl_function_name
    global im_function_name
    global arg_function_name
    global abs_function_name
    global sign_function_name
    
    e_constant_name =         '__e__'
    pi_constant_name =        '__pi__'
    tau_constant_name =        '__tau__'
    sin_function_name =       '__sin__'
    cos_function_name =       '__cos__'
    sec_function_name =       '__sec__'
    csc_function_name =       '__csc__'
    tan_function_name =       '__tan__'
    cot_function_name =       '__cot__'
    sinh_function_name =      '__sinh__'
    cosh_function_name =      '__cosh__'
    sech_function_name =      '__sech__'
    csch_function_name =      '__csch__'
    tanh_function_name =      '__tanh__'
    coth_function_name =      '__coth__'
    arcsin_function_name =    '__arcsin__'
    arccos_function_name =    '__arccos__'
    arcsec_function_name =    '__arcsec__'
    arccsc_function_name =    '__arccsc__'
    arctan_function_name =    '__arctan__'
    arccot_function_name =    '__arccot__'
    arcsinh_function_name =   '__arcsinh__'
    arccosh_function_name =   '__arccosh__'
    arcsech_function_name =   '__arcsech__'
    arccsch_function_name =   '__arcscsh__'
    arctanh_function_name =   '__arctanh__'
    arccoth_function_name =   '__arccoth__'
    ln_function_name =        '__ln__'
    log_function_name =       '__log__'
    pow_function_name =       '__pow__'
    sqrt_function_name =      '__sqrt__'
    exp_function_name =       '__exp__'
    rl_function_name =        '__Rl__'
    im_function_name =        '__Im__'
    arg_function_name =       '__arg__'
    abs_function_name =       '__abs__'
    sign_function_name =      '__sign__'
    
def set_names_to_default():
    global e_constant_name
    global pi_constant_name
    global tau_constant_name
    global sin_function_name
    global cos_function_name
    global sec_function_name
    global csc_function_name
    global tan_function_name
    global cot_function_name
    global sinh_function_name
    global cosh_function_name
    global sech_function_name
    global csch_function_name
    global tanh_function_name
    global coth_function_name
    global arcsin_function_name
    global arccos_function_name
    global arcsec_function_name
    global arccsc_function_name
    global arctan_function_name
    global arccot_function_name
    global arcsinh_function_name
    global arccosh_function_name
    global arcsech_function_name
    global arccsch_function_name
    global arctanh_function_name
    global arccoth_function_name
    global ln_function_name
    global log_function_name
    global pow_function_name
    global sqrt_function_name
    global exp_function_name
    global rl_function_name
    global im_function_name
    global arg_function_name
    global abs_function_name
    global sign_function_name
    
    e_constant_name =         'e'
    pi_constant_name =        'pi'
    tau_constant_name =       'tau'
    sin_function_name =       'sin'
    cos_function_name =       'cos'
    csc_function_name =       'csc'
    sec_function_name =       'sec'
    tan_function_name =       'tan'
    cot_function_name =       'cot'
    sinh_function_name =      'sinh'
    cosh_function_name =      'cosh'
    sech_function_name =      'sech'
    csch_function_name =      'csch'
    tanh_function_name =      'tanh'
    coth_function_name =      'coth'
    arcsin_function_name =    'arcsin'
    arccos_function_name =    'arccos'
    arcsec_function_name =    'arcsec'
    arccsc_function_name =    'arccsc'
    arctan_function_name =    'arctan'
    arccot_function_name =    'arccot'
    arcsinh_function_name =   'arcsinh'
    arccosh_function_name =   'arccosh'
    arcsech_function_name =   'arcsech'
    arccsch_function_name =   'arccsch'
    arctanh_function_name =   'arctanh'
    arccoth_function_name =   'arccoth'
    ln_function_name =        'ln'
    log_function_name =       'log'
    pow_function_name =       'pow'
    sqrt_function_name =      'sqrt'
    exp_function_name =       'exp'
    rl_function_name =        'Rl'
    im_function_name =        'Im'
    arg_function_name =       'arg'
    abs_function_name =       'abs'
    sign_function_name =      'sign'
    
def return_to_default():
    set_names_to_default()
    
    global TabulationStep
    TabulationStep = 4
    
    global TauMode
    TauMode = False
    
    global AutoSimplify
    AutoSimplify = True
    
    global ExpandComplexExponentials
    ExpandComplexExponentials = False
    
    global ExpressArg
    ExpressArg = True
    
    global NumericalNumbersPrintLength
    NumericalNumbersPrintLength = 12
    
    global ColorfulPrint
    ColorfulPrint = True
    
    
def show_names():
    print('Names of constants and functions:')
    print('    Euler`s number                           ' + e_constant_name)
    print('    pi constant                              ' + pi_constant_name)
    print('    tau constant                             ' + tau_constant_name)
    print()
    
    print('    exponential function                     ' + exp_function_name)
    print('    natural logarithm function               ' + ln_function_name)
    print('    arbitrary base logarithm function        ' + log_function_name)
    print('    power function                           ' + pow_function_name)
    print('    square root function                     ' + sqrt_function_name)
    print()
    
    print('    real part function                       ' + rl_function_name)
    print('    imaginary part function                  ' + im_function_name)
    print('    complex argument function                ' + arg_function_name)
    print('    absolute value function                  ' + abs_function_name)
    print('    complex signum function                  ' + sign_function_name)
    print()
    
    print('    sine function                            ' + sin_function_name)
    print('    cosine function                          ' + cos_function_name)
    print('    secant function                          ' + sec_function_name)
    print('    cosecant function                        ' + csc_function_name)
    print('    tangent function                         ' + tan_function_name)
    print('    cotangent function                       ' + cot_function_name)
    print()
    
    print('    hyperbolic sine function                 ' + sinh_function_name)
    print('    hyperbolic cosine function               ' + cosh_function_name)
    print('    hyperbolic secant function               ' + sech_function_name)
    print('    hyperbolic cosecant function             ' + csch_function_name)
    print('    hyperbolic tangent function              ' + tanh_function_name)
    print('    hyperbolic cotangent function            ' + coth_function_name)
    print()
    
    print('    inverse sine function                    ' + arcsin_function_name)
    print('    inverse cosine function                  ' + arccos_function_name)
    print('    inverse secant function                  ' + arcsec_function_name)
    print('    inverse cosecant function                ' + arccsc_function_name)
    print('    inverse tangent function                 ' + arctan_function_name)
    print('    inverse cotangent function               ' + arccot_function_name)
    print()
    
    print('    inverse hyperbolic sine function         ' + arcsinh_function_name)
    print('    inverse hyperbolic cosine function       ' + arccosh_function_name)
    print('    inverse hyperbolic secant function       ' + arcsech_function_name)
    print('    inverse hyperbolic cosecant function     ' + arccsch_function_name)
    print('    inverse hyperbolic tangent function      ' + arctanh_function_name)
    print('    inverse hyperbolic cotangent function    ' + arccoth_function_name)
    print()
    
def function(*args, partials=None):
    def check_partials(partials, n_args):
        if not isinstance(partials, list):
            raise TypeError('Partial derivative list must be a list, not ' + str(type(partials)))
            
        if(len(partials) != n_args):
            raise TypeError('Wrong length of partial derivative list for a function with ' +\
                            str(n_args) + 'arguments: ' + str(len(partials)))
        
        for item in partials:
            if not isinstance(item, int):
                raise TypeError('Orders of partial derivatives must be represented by int, not by ' + str(type(item)))
            
            if(item < 0):
                raise TypeError('Order of partial derivative cannot be negative: ' + str(item))
        
    if(len(args) == 0):
        raise TypeError('function(*args) must take at least one argument')
    
    if isinstance(partials, type(None)):
        return Function(args)
    
    else:
        if(len(args) > 1):
            check_partials(partials, len(args) - 1)
            
            return Function( (Function.make_up_a_name(args[0], partials), *args[1:]) )
        
        else:
            raise TypeError('Cannot differentiate a function without arguments')

def numerical(*args):
    if(len(args) == 0):
        raise TypeError('numerical(*args) must take at least one argument')
    
    elif(len(args) == 1):
        if(isinstance(args[0], float) or
           isinstance(args[0], int)):
            return Numerical(args[0], 0)
        
        else:
            raise TypeError('numerical(arg) takes int or float as an argument, not ' + str(type(args[0])))
    
    elif(len(args) == 2):
        if( (isinstance(args[0], float) or isinstance(args[0], int)) and 
            (isinstance(args[1], float) or isinstance(args[1], int))):
            return Numerical(args[0], args[1])
        
        else:
            raise TypeError('numerical(rl, im) takes ints or floats as arguments, not ' +\
                            str(type(args[0])) + ' and ' + str(type(args[1])))    
    
    else:
        raise TypeError('Too many arguments for numerical(args*): ' + str(len(args)))
    
def cp(*args):  # this little function is here to make your life easier
                # it reduces the amount of work you have to do in order
                # to create a Complex
    
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
        raise ValueError('Expected 1 or 2 arguments, got ' + str(len(args)))
    
    return c



pi_constant = cp(pi_constant_name)
tau_constant = cp(tau_constant_name)
half_pi_constant = cp( ((1, 2), pi_constant_name) )
e_constant = cp(e_constant_name)
