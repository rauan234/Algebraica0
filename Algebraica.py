import math



global TabulationStep
TabulationStep = 4
# this refers to Complex.__str__ function

global TauMode
TauMode = False
# if true, tau = 2 pi will be used instead of pi

global AutoSimplify
AutoSimplify = False
# if true, all class Complex are simplified before each print

global e_constant_name
global pi_constant_name
global tau_constant_name
global sin_function_name
global cos_function_name
global tan_function_name
global ctg_function_name
global sinh_function_name
global cosh_function_name
global tanh_function_name
global ctgh_function_name
global arcsin_function_name
global arccos_function_name
global arctan_function_name
global arcctg_function_name
global arcsinh_function_name
global arccosh_function_name
global arctanh_function_name
global arcctgh_function_name
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
tan_function_name =       'tan'
ctg_function_name =       'ctg'
sinh_function_name =      'sinh'
cosh_function_name =      'cosh'
tanh_function_name =      'tanh'
ctgh_function_name =      'ctgh'
arcsin_function_name =    'arcsin'
arccos_function_name =    'arccos'
arctan_function_name =    'arctan'
arcctg_function_name =    'arcctg'
arcsinh_function_name =   'arcsinh'
arccosh_function_name =   'arccosh'
arctanh_function_name =   'arctanh'
arcctgh_function_name =   'arcctgh'
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
        if(self.den == 1):
            return str(self.num)
        
        elif(self.num == 0):
            return '0'
        
        else:
            return str(self.num) + '/' + str(self.den)  
    
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
        return not first == second
    
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
        if((self.den == 1) or (self.num == 1)):  # in this case the fraction
                                                 # cannot be reduced
            return
        
        if(self.den < 0):
            self.den *= -1
            self.num *= -1
        
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
    
class Function:
    def __init__(self, params):  # class Function item describes a certain function
                                 # it has the information on function`s name
                                 # and what arguments it has
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
        
        
        if not isinstance(params[0], str):  # params[0] is supposed to tell the name of the function
                                            # so that`s why it must be a string
            raise TypeError('Function name must be a string')
        if (len(params[0]) == 0):
            raise ValueError('Function name cannot be an empty string')
        
        for item in params[1:]:
            if not isinstance(item, Complex):  # all of the function`s arguments
                                               # must be Complex
                raise TypeError('Arguments of a function must be Complex')
        
        self.name = params[0]
        self.args = copy_list(params[1:])
        
        n_args = self.number_of_arguments_required()
        if(n_args == None):  # if the function name is unknown
            pass
        else:
            if (n_args != len(self.args)):  # if the number of arguments given is incorrect
                if(n_args == 0):
                    raise RuntimeError('<' + self.name + '> function must have zero arguments, but the number of arguments given is ' + str(len(self.args)))
                elif(n_args == 1):
                    raise RuntimeError('<' + self.name + '> function must have one argument, but the number of arguments given is ' + str(len(self.args)))
                else:
                    raise RuntimeError('<' + self.name + '> function must have ' + str(n_args) + ' arguments, but the number of arguments given given is ' + str(len(self.args)))
        
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
            sin_function_name, cos_function_name, tan_function_name, ctg_function_name,
            sinh_function_name, cosh_function_name, tanh_function_name, ctgh_function_name,
            arcsin_function_name, arccos_function_name, arctan_function_name, arcctg_function_name,
            arcsinh_function_name, arccosh_function_name, arctanh_function_name, arcctgh_function_name,
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
        tab = ' ' * elevation * TabulationStep
        
        out += str(self.name)
        
        if(len(self.args) == 0):
            return self.name
        
        elif((len(self.args) == 1) and self.args[0].is_constant()):
            out += '[ '
            out += self.args[0].show(0, one_line=True)
            out += ' ]'
            
        else:
            out += '[\n'
            
            for i in range(len(self.args)):
                out += self.args[i].show(elevation + 2)
                
                if(i != len(self.args) - 1):
                    out += ', \n'
                out += '\n'
                
            out += tab + '    ]'

        return out
    
    def __str__(self):
        return self.show(0)
    
    def alternate_form(self):
        for i in range(len(self.args)):
            self.args[i].simplify()
        
        if(self.name == pow_function_name):  # N ** P
            N = self.args[0]
            P = self.args[1]
            
            if((N == cp(1)) or (P == cp(0))):
                return cp(1)
            if(P == cp(1)):
                return N
            
            if(N.im.might_be_nonzero() and N.rl.is_zero()):
                # pow( a * i, b) = pow(a, b) * pow(i, b)
                return cp(0, ((1, 1), (pow_function_name, N / cp(0, 1), P) )) * (cp(0, 1) ** P)
            
            if(N.is_real_int() and P.is_real_int()):
                n = N.rl.terms[0].rat.num
                p = P.rl.terms[0].rat.num
                
                if(p > 0):
                    return cp(n ** p)
                else:
                    return cp( (1, n**(-p)) )            
                    
            if((len(N.rl.terms) == 1) and (len(N.im.terms) == 0)):
                if((len(P.rl.terms) == 1) and (len(P.im.terms) == 0)):
                    if(len(P.rl.terms[0].irt) == 0):
                        N_num_primes = prime_factors(N.rl.terms[0].rat.num)
                        N_den_primes = prime_factors(N.rl.terms[0].rat.den)
                        # pow( 2 * 3 * pi, P ) = pow(2, P) * pow(3, P) * pow(pi, P)
                        
                        prat = P.rl.terms[0].rat
                        
                        primes = []
                        for n in N_num_primes:
                            primes.append( (n[0], Rational( n[1], 1) * prat) )
                        for d in N_den_primes:
                            primes.append( (d[0], Rational(-d[1], 1) * prat) )
                    
                        out = cp(1)
                        for fun in N.rl.terms[0].irt:
                            out *= cp( ((1, 1), (pow_function_name, cp(fun), P)) )
                        for prime in primes:
                            int_part, frac_part = prime[1].proper_fraction()
                            out *= cp( ((prime[0] ** int_part, 1), (pow_function_name, cp(prime[0]), cp(frac_part))) )
                        
                        if(cp(self) != out):
                            return out
            
            if((len(P.rl.terms) + len(P.im.terms)) > 1):  # if P is multi-term
                out = cp(1)
                
                # pow( a, b + c ) = pow(a, b) * pow(a, c)
                for term in P.rl.terms:
                    out *= cp( ((1, 1), (pow_function_name, N, cp(term))) )
                for term in P.im.terms:
                    out *= cp( ((1, 1), (pow_function_name, N, cp(term) * cp(0, 1))) )                    
                
                return out
            
            if((len(N.rl.terms) == 1) and (len(N.im.terms) == 0)):  # if N is single-term
                if(len(N.rl.terms[0].irt) == 1):  # if the term is single-function
                    if(N.rl.terms[0].irt[0].name == pow_function_name):
                        # pow( pow(a, b), c ) = pow(a, b * c)
                        # (this is, of course, not always true)
                        
                        if(N.rl.terms[0].irt[0].args[0].is_real_and_positive()):
                            n = N.rl.terms[0].irt[0].args[0]
                            p1 = N.rl.terms[0].irt[0].args[1]
                            p2 = P
                            
                            return cp( ((1, 1), (pow_function_name, n, p1 * p2)) )
                            
                elif(len(N.rl.terms[0].irt) > 1):  # if the term contains several functions,
                                                   # such as 2 * pi * ln(3)
                    term = N.rl.terms[0]
                    
                    out = cp( ((1, 1), (pow_function_name, cp( (term.rat.num, term.rat.den) ), P)) )
                    for f in term.irt:
                        out *= cp( ((1, 1), (pow_function_name, cp(f), P)) )
                    out.simplify()
                    
                    if(cp(self) != out):
                        return out
                    
            if((len(N.rl.terms) == 1) and (len(N.im.terms) == 0)):
                if(len(N.rl.terms[0].irt) == 1):
                    if(N.rl.terms[0].irt[0].name == e_constant_name):
                        return cp( ((1, 1), (exp_function_name, P * N.rl.terms[0].rat)) )
            
        elif(self.name == sqrt_function_name):
            # sqrt(x) is changed to pow(x, 1/2)
            return cp( ((1, 1), (pow_function_name, self.args[0], cp((1, 2)))) )
        
        elif(self.name == log_function_name):
            # log(a, b) is changed to ln(a) / ln(b)
            return cp( ((1, 1), [
                (ln_function_name, self.args[0]),
                (pow_function_name, cp(
                    ((1, 1), (ln_function_name, self.args[1]))
                    ), cp(-1))
            ]) )
        
        elif(self.name == pi_constant_name):
            if TauMode:
                return cp( ((1, 2), tau_constant_name) )
            
        elif(self.name == tau_constant_name):
            if not TauMode:
                return cp( ((2, 1), pi_constant_name) )
                
        elif(self.name == exp_function_name):
            P = self.args[0]
            
            if(P == cp(0)):
                return cp(1)
            if(P == cp(1)):
                return e_constant_name
            
            if((len(P.rl.terms) + len(P.im.terms)) > 1):  # if the argument is multi-term
                out = cp(1)
                
                # exp(a + b) is changed to exp(a) * exp(b)
                power = P
                for term in power.rl.terms:
                    out *= cp( ((1, 1), (exp_function_name, cp(term))) )
                for term in power.im.terms:
                    out *= cp( ((1, 1), (exp_function_name, cp(term) * cp(0, 1))) )                    
                
                return out
            
            if((len(P.rl.terms) == 1) and (len(P.im.terms) == 0)):
                if(len(P.rl.terms[0].irt) == 1):
                    if(P.rl.terms[0].irt[0].name == ln_function_name):
                        # exp( ln(a) ) = a
                        
                        if(P.rl.terms[0].rat == 1):
                            return P.rl.terms[0].irt[0].args[0]
                        else:
                            return cp( ((1, 1), (pow_function_name, P.rl.terms[0].irt[0].args[0],
                                                  cp(P.rl.terms[0].rat))) )    
            
        elif(self.name == ln_function_name):
            N = self.args[0]
            
            if(N == cp(1)):
                    return cp(0)
            elif(N == e_constant):
                    return cp(1)
                
            if((len(N.rl.terms) == 1) and (len(N.im.terms) == 0)):  # if N is single-term
                term = N.rl.terms[0]
                
                if(len(N.rl.terms[0].irt) == 0):
                    if(term.rat == 1):
                        return cp(0)
                    
                    else:
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
                            # ln(pi ** 3) = 3 * ln(pi)
                            return cp( ((1, 1), (ln_function_name, term.irt[0].args[0])) ) * term.irt[0].args[1]
                        
                        else:
                            # ln(2 * pi ** 3) = ln(2) + 3 * ln(pi)
                            return (cp( ((1, 1), (ln_function_name, cp(term.rat))) ) +
                                    cp( ((1, 1), (ln_function_name, term.irt[0].args[0])) ) * term.irt[0].args[1])
                    
                    else:
                        if(term.rat != 1):
                            # ln(3 * pi) = ln(3) + ln(pi)
                            return (cp( ((1, 1), (ln_function_name, cp(term.rat))) ) +
                                    cp( ((1, 1), (ln_function_name, cp(term.irt[0]))) ))                            
                    
                else:
                    # ln( 4/7 * a * b * c * ... ) =
                    # = ln(4) - ln(7) + ln(a) + ln(b) + ln(c) ...
                    out = (cp( (( 1, 1 ), (ln_function_name, cp( (term.rat.num, 1) ))) ) +
                           cp( ((-1, 1 ), (ln_function_name, cp( (term.rat.den, 1) ))) ))
                    for t in term.irt:
                        out += cp( ((1, 1), (ln_function_name, cp(t))) )                 
                    
                    return out
                        
        elif(self.name == rl_function_name):
            if(self.args[0].is_zero()):
                return cp(0)
            
            elif(self.args[0].rl.is_zero() and self.args[0].im.might_be_nonzero()):
                return cp( ((-1, 1), (im_function_name, cp(self.args[0].im))) )  # Rl(bi) = -Im(b)
            
            else:
                if self.args[0].is_always_real():
                    return self.args[0]
                if self.args[0].is_always_imag():
                    return cp(0)
                
                if((len(self.args[0].rl.terms) + len(self.args[0].im.terms)) > 1):
                    out = cp(0)
                    
                    # Rl(a + bi) = Rl(a) - Im(b)
                    for term in self.args[0].rl.terms:
                        out += cp( ((1, 1), (rl_function_name, cp(term))) )
                    for term in self.args[0].im.terms:
                        out += cp( ((-1, 1), (im_function_name, cp(term))) )
                    
                    return out
                
        elif(self.name == im_function_name):
            if(self.args[0].is_zero()):
                return cp(0)
            
            elif(self.args[0].rl.is_zero() and self.args[0].im.might_be_nonzero()):
                return cp( ((1, 1), (rl_function_name, cp(self.args[0].im))) )  # Im(bi) = Rl(b)
            
            else:
                if self.args[0].is_always_real():
                    return cp(0)
                if self.args[0].is_always_imag():
                    return self.args[0] / cp(0, 1)
                
                if((len(self.args[0].rl.terms) + len(self.args[0].im.terms)) > 1):
                    out = cp(0)
                    
                    # Im(a + bi) = Im(a) + Rl(b)
                    for term in self.args[0].rl.terms:
                        out += cp( ((1, 1), (im_function_name, cp(term))) )
                    for term in self.args[0].im.terms:
                        out += cp( ((1, 1), (rl_function_name, cp(term))) )
                    
                    return out                
                        
        elif(self.name == sign_function_name):
            N = self.args[0]
            if(N.is_real_rational()):
                return cp(N.rl.terms[0].rat.get_sign())
            
            if(N.is_always_real()):
                if((len(N.rl.terms) == 1) and (len(N.im.terms) == 0)):
                    if(len(N.rl.terms[0].irt) == 1):
                        s = N.rl.terms[0].irt[0].sign()
                        
                        if(s != None):
                            return cp(s)
                        
                    elif(len(N.rl.terms[0].irt) > 1):
                        out = cp(N.rl.terms[0].rat.get_sign())
                        
                        # sign( -3/2 * a * b * c ... ) =
                        # = sign(-3/2) * sign(a) * sign(b) * sign(c) ...
                        for fun in N.rl.terms[0].irt:
                            s = fun.sign()
                            
                            if(s == None):
                                out *= cp( ((1, 1), (sign_function_name, cp(fun))) )
                            else:
                                out *= s

                        return out
                    
        elif(self.name == arg_function_name):
            N = self.args[0]
            
            if(N.rl.is_rational() and N.im.is_rational()):
                if(N.im.is_zero()):
                    return cp(0)
                elif(N.rl.is_zero()):
                    return cp(N.im.terms[0].rat.get_sign()) * half_pi
                
                else:
                    return cp( ((1, 1), (arctan_function_name, cp(N.im.terms[0].rat / N.rl.terms[0].rat) )) )
                        
        elif(self.name == sin_function_name):
            if(len(self.args[0].rl.terms) == 0):
                return cp(0)
            else:
                if(self.args[0].rl.terms[0].rat < Rational(0)):  # sin( -5/2 x ) = -sin( 5/2 x )
                    return -Function( (sin_function_name, -self.args[0]) ).alternate_form()
                
                else:
                    if(self.args[0] == cp(0)):
                        return cp(0)
                    
                    elif(self.args[0] == pi_constant ):
                        return cp(0)            
                    
                    elif(self.args[0] == cp( ((1, 2), pi_constant_name)) ):
                        return cp( (1, 1) )
                    
                    elif(self.args[0] == cp( ((1, 3), pi_constant_name)) ):
                        return cp( ((1, 2), (sqrt_function_name, cp(3))) )  
                    
                    elif(self.args[0] == cp( ((2, 3), pi_constant_name)) ):
                        return cp( ((1, 2), (sqrt_function_name, cp(3))) )              
                    
                    elif(self.args[0] == cp( ((1, 4), pi_constant_name)) ):
                        return cp( ((1, 2), (sqrt_function_name, cp(2))) )
                    
                    elif(self.args[0] == cp( ((1, 6), pi_constant_name)) ):
                        return cp( (1, 2) )       
        
        elif(self.name == cos_function_name):
            return Function( (sin_function_name, half_pi - self.args[0]) )
        
        elif(self.name == tan_function_name):
            if(len(self.args[0].rl.terms) == 0):
                return cp(0)
            else:
                if(self.args[0].rl.terms[0].rat < Rational(0)):
                    return -Function( (tan_function_name, -self.args[0]) ).alternate_form()
                
                else:            
                    if(self.args[0] == cp(0)):
                        return cp(0)
                    
                    elif(self.args[0] == pi_constant):
                        return cp(0)            
                    
                    elif(self.args[0] == cp( ((1, 3), pi_constant_name)) ):
                        return cp( ((1, 1), (sqrt_function_name, cp(3))) )  
                    
                    elif(self.args[0] == cp( ((2, 3), pi_constant_name)) ):
                        return cp( ((-1, 1), (sqrt_function_name, cp(3))) )              
                    
                    elif(self.args[0] == cp( ((1, 4), pi_constant_name)) ):
                        return cp(1)
                    
                    elif(self.args[0] == cp( ((1, 6), pi_constant_name)) ):
                        return cp( ((1, 3), (sqrt_function_name, cp(3))) )
                    
        elif(self.name == arcsin_function_name):
            f = self.args[0]
            
            if(f.is_zero()):
                return cp(0)            
            
            if((len(f.rl.terms) == 1) and (len(f.im.terms) == 0)):
                if(f.rl.terms[0].rat < 0):
                    return cp( ((-1, 1), (arcsin_function_name, -f)) )
                
                else:
                    if(f == cp((1, 2))):
                        return cp( ((1, 6), pi_constant_name) )
                    if(f == cp(((1, 2), (pow_function_name, cp(2), cp((1, 2)))))):
                        return cp( ((1, 4), pi_constant_name) )
                    if(f == cp(((1, 2), (pow_function_name, cp(3), cp((1, 2)))))):
                        return cp( ((1, 3), pi_constant_name) )
                    if(f == cp(1)):
                        return half_pi
                
        elif(self.name == arccos_function_name):
            f = self.args[0]
            
            if(f.is_zero()):
                return half_pi
                        
            if((len(f.rl.terms) == 1) and (len(f.im.terms) == 0)):
                if(f.rl.terms[0].rat < 0):
                    return pi_constant - cp( ((1, 1), (arccos_function_name, -f)) )       
                
                else:
                    if(f == cp((1, 2))):
                        return cp( ((1, 3), pi_constant_name) )
                    if(f == cp(((1, 2), (pow_function_name, cp(2), cp((1, 2)))))):
                        return cp( ((1, 4), pi_constant_name) )
                    if(f == cp(((1, 2), (pow_function_name, cp(3), cp((1, 2)))))):
                        return cp( ((1, 6), pi_constant_name) )
                    if(f == cp(1)):
                        return cp(0)              
        
        elif(self.name == arctan_function_name):
            f = self.args[0]
            
            if(f.is_zero()):
                return cp(0)
            
            if((len(f.rl.terms) == 1) and (len(f.im.terms) == 0)):
                if(f.rl.terms[0].rat < 0):
                    return cp( ((-1, 1), (arctan_function_name, -f)) )
                
                else:
                    if(f == cp( ((1, 3), (pow_function_name, cp(3), cp((1, 2))))) ):
                        return cp( ((1, 6), pi_constant_name))
                    if(f == cp(1)):
                        return cp( ((1, 4), pi_constant_name) )
                    if(f == cp( ((1, 1), (pow_function_name, cp(3), cp((1, 2)))))):
                        return cp( ((1, 3), pi_constant_name) )
                    
        elif(self.name == sinh_function_name):
            f = self.args[0]
            
            if(f.is_zero()):
                return cp(0)
            if((len(f.rl.terms) == 1) and (len(f.im.terms) == 0)):
                if(f.rl.terms[0].rat < 0):
                    return cp( ((-1, 1), (sinh_function_name, -f)) )
                
        elif(self.name == cosh_function_name):
            f = self.args[0]
            
            if(f.is_zero()):
                return cp(1)
            if((len(f.rl.terms) == 1) and (len(f.im.terms) == 0)):
                if(f.rl.terms[0].rat < 0):
                    return cp( ((1, 1), (cosh_function_name, -f)) )      
        
        elif(self.name == tanh_function_name):
            f = self.args[0]
                    
            if(f.is_zero()):
                return cp(0)
            if((len(f.rl.terms) == 1) and (len(f.im.terms) == 0)):
                if(f.rl.terms[0].rat < 0):
                    return cp( ((-1, 1), (tanh_function_name, -f)) )              
                
        return cp(self)
    
    def sign(self):
        always_positive = [pi_constant_name, e_constant_name]
        
        if(self.name in always_positive):
            return 1
        
        return None
    
    def is_always_real(self):
        if(self.name == pi_constant_name):
            return True
        if(self.name == e_constant_name):
            return True
        if((self.name == rl_function_name) or
           (self.name == im_function_name)):
            return True
        
        if(self.name == pow_function_name):
            if(self.args[0].is_always_real() and
               self.args[1].is_always_real()):
                return True
            
        # if the argument of the following functions is real, the output is real too
        lst = [exp_function_name, ln_function_name, sin_function_name, cos_function_name, tan_function_name, 'cotan', sinh_function_name, cosh_function_name, tanh_function_name, 'cotanh']
        if(self.name in lst):
            if self.args[0].is_always_real():
                return True
        
        return False
    
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
                out += ' * '
                
            for j in range(len(self.irt)):
                if((j != 0) or not show_rat):
                    out += tab
                    
                out += self.irt[j].show(elevation)
                
                if(j != len(self.irt) - 1):
                    if(one_line):
                        out += ' * '
                    else:
                        out += ' *\n'
        
        return out
    def __str__(self):
        return self.show(0)
    
    def sort(self):
        self.irt.sort(key=function_serial_index)
    
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
    
    def is_real_and_positive(self):
        sign = self.rat.get_sign()
        
        for fun in self.irt:
            s = fun.sign()
            
            if(s == None):
                return False
            else:
                sign *= s
        
        return (sign == 1)

    def expanded(self):
        out = cp( (self.rat.num, self.rat.den) )
        
        for fun in self.irt:
            out *= fun.alternate_form()
        
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
        for term in self.terms:
            if not term.is_real_and_positive():
                return False
        
        return True    
    
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
    
    def sort(self):
        self.terms.sort(key=term_serial_index)
    
    def reduce(self):
        if(len(self.terms) > 1):
            newterms = []
            
            while(len(self.terms) > 1):
                term = self.terms.pop()
                if not term.might_be_nonzero():
                    continue

                new_rat = term.rat
                while(len(self.terms) and term.eq_irt(self.terms[0])):
                    new_rat += self.terms.pop().rat
                
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
    
    

class Complex:
    def __init__(self, new_real, new_imag, new_is_simplified=False):  # object of class Complex is just M + N i,
                                             # where M and N are class Irrational objects and
                                             # i is the square root of negative one
        if not(isinstance(new_real, Irrational) and
               isinstance(new_imag, Irrational)):
            
            raise TypeError('Wrong Complex.__init__ arguments` type: ' + str(type(new_real)) + ', ' + str(type(new_imag)))
        
        self.rl = new_real.copy()
        self.im = new_imag.copy()
        self.is_simplified = new_is_simplified
        
    def __str__(self):
        if AutoSimplify:
            self.simplify()
            
        return '<<<\n' + self.show(1) + '\n>>>'
    
    def __add__(first, second):
        if isinstance(second, Complex):
            return Complex(
                first.rl + second.rl,
                first.im + second.im
            )
        else:
            raise TypeError('Wrong Complex.__add__ argument type: expected Complex, got ' + str(type(second)))
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
        else:
            raise TypeError('Wring Complex.__sub__ argument type: ' + str(type(second)))
    def __isub__(first, second):
        return Complex.__sub__(first, second)
    
    def __invert__(self):  # turns (a + b i) to (a - b i)
        return Complex(
            make_irrational( ('Rl', self) ),
            make_irrational( ('Im', self) ),
            self.is_simplified
        )
    
    def __eq__(first, second):
        if isinstance(second, Complex):
            return ((first.rl == second.rl) and
                    (first.im == second.im))
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
                inv = cp(1)
                
                for i in range(-second):
                    inv *= first
                
                return cp(((1, 1), [(pow_function_name, inv, cp(-1))]))
            
        elif(isinstance(second, Complex)):
            abs_val = abs(first)
            cmp_arg = first.arg()
            
            ln_first = cp(
                ((1, 1), [(ln_function_name, abs_val)])
            ) + cp(0, 1) * cmp_arg
            
            return cp(((1, 1), [(exp_function_name, ln_first * second)]), 0)
        
        else:
            raise TypeError('Wrong Complex.__pow__ argument type: ' + str(type(second)))
        
    def __abs__(self):
        return cp( ((1, 1), (abs_function_name, self)) )
    
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

    def arg(self):
        if self.rl.is_zero():
            if self.im.is_zero():
                return cp(0)
            
            else:
                return cp( ((1, 2), [pi_constant_name, (sign_function_name, cp( ((1, 1), (im_function_name, self)) ) )]) )
        
        elif self.im.is_zero():
            return cp(0)
        
        else:
            return cp( ((1, 1), (arg_function_name, self)) )
    
    def show(self, elevation, one_line=False):
        out = ''
        tab = ' ' * elevation * TabulationStep
        
        rl_nz = self.rl.might_be_nonzero()
        im_nz = self.im.might_be_nonzero()
        
        if(rl_nz):
            if(self.rl.is_constant()):
                if(im_nz):
                    out += tab
                    out += '{ '
                    out += self.rl.show(0, ol=True)
                    out += ' }'
                    
                else:
                    out += self.rl.show(elevation, ol=True)
            
            else:
                if(im_nz):
                    if not one_line:
                        out += tab
                        out += '{\n'
                    
                    if(one_line):
                        out += self.rl.show(0, ol=one_line)
                    else:
                        out += self.rl.show(elevation + 1, ol=one_line)
                else:
                    out += self.rl.show(elevation, ol=one_line)
                
                if(im_nz):
                    if not one_line:
                        out += '\n' + tab + '}'
            
        if (rl_nz and im_nz):
            if one_line:
                out += ' + '
            else:
                out += '\n'
                out += tab + '+\n'        
        
        if(im_nz):
            if(self.im.is_constant()):
                out += tab + '{ '
                out += self.im.show(0, ol=True)
                out += ' } i'
                
            else:
                if one_line:
                    out += '{ '
                    out += self.im.show(elevation + 1, ol=True)
                    out += ' } i'
                else:
                    out += tab + '{\n'
                    out += self.im.show(elevation + 1, ol=False)
                    out += '\n' + tab + '} i'
            
        if not(rl_nz or im_nz):
            out += self.rl.show(elevation, ol=one_line)
        
        return out
    
    def copy(self):
        return Complex(self.rl,
                       self.im,
                       new_is_simplified = self.is_simplified)
    
    def simplified_once(self):
        new_rl = cp(0)
        new_im = cp(0)
        
        for term in self.rl.terms:
            new_rl += term.expanded()
        for term in self.im.terms:
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
        res = self.simplified()
        
        self.rl = res.rl
        self.im = res.im
        self.is_simplified = res.is_simplified
    


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

def prime_factors(n):
    factors = []
    
    imax = n
    for i in range(2, imax + 1):
        p = 0
        while(n % i == 0):
            n = n // i
            p += 1
        
        if(p > 0):
            factors.append( (i, p) )
    
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
                    
                    elif(isinstance(a[0], tuple) and isinstance(a[1], str)):  # if a is a tuple representing a constant with a coefficient
                                                                              # such as a = ((2, 3), pi_constant_name)
                                                                              # OUT = 2/3 * pi + 0 i
                        return irrational((a[0][0], a[0][1]), [(a[1],)])
                    
                    elif(isinstance(a[0], tuple) and isinstance(a[1], tuple)):
                                             # if a is a tuple describing a function with a coefficient
                                             # for example, a = ((4, 9), (ln_function_name, cp(7)))
                                             # OUT = 4/9 * ln[7] + 0/1 i
                        if(isinstance(a[1][0], str)):
                            for i in range(1, len(a[1])):
                                if not isinstance(a[1][i], Complex):
                                    raise TypeError('Tuple denoting arguments of a Function has wrong type of entries: expected Complex, got ' + str(type(a[1][i])))
                        else:
                            raise TypeError('Tuple denoting a Function does not contain the name of the function')
                        
                        return irrational((a[0][0], a[0][1]), [a[1]])
                    
                    elif (isinstance(a[0], tuple) and isinstance(a[1], list)):
                                             # if a is a tuple representing a product of functions with a coefficient
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

def function_serial_index(f):
    out = f.name
    
    for arg in f.args:
        out += str(arg)
    
    return out

def term_serial_index(t):
    
    out = ''
    
    for f in t.irt:
        out += function_serial_index(f)
    
    return out

def sign(inp):
    return inp / abs(inp)

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

def set_names_to_secure():
    global e_constant_name
    global pi_constant_name
    global tau_constant_name
    global sin_function_name
    global cos_function_name
    global tan_function_name
    global ctg_function_name
    global sinh_function_name
    global cosh_function_name
    global tanh_function_name
    global ctgh_function_name
    global arcsin_function_name
    global arccos_function_name
    global arctan_function_name
    global arcctg_function_name
    global arcsinh_function_name
    global arccosh_function_name
    global arctanh_function_name
    global arcctgh_function_name
    global ln_function_name
    global log_function_name
    global pow_function_name
    global sqrt_function_name
    global exp_function_name
    global rl_function_name
    global im_function_name
    global arg_function_name
    global sign_function_name
    
    e_constant_name =         '__e__'
    pi_constant_name =        '__pi__'
    tau_constant_name =        '__tau__'
    sin_function_name =       '__sin__'
    cos_function_name =       '__cos__'
    tan_function_name =       '__tan__'
    ctg_function_name =       '__ctg__'
    sinh_function_name =      '__sinh__'
    cosh_function_name =      '__cosh__'
    tanh_function_name =      '__tanh__'
    ctgh_function_name =      '__ctgh__'
    arcsin_function_name =    '__arcsin__'
    arccos_function_name =    '__arccos__'
    arctan_function_name =    '__arctan__'
    arcctg_function_name =    '__arcctg__'
    arcsinh_function_name =   '__arcsinh__'
    arccosh_function_name =   '__arccosh__'
    arctanh_function_name =   '__arctanh__'
    arcctgh_function_name =   '__arcctgh__'
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
    global tan_function_name
    global ctg_function_name
    global sinh_function_name
    global cosh_function_name
    global tanh_function_name
    global ctgh_function_name
    global arcsin_function_name
    global arccos_function_name
    global arctan_function_name
    global arcctg_function_name
    global arcsinh_function_name
    global arccosh_function_name
    global arctanh_function_name
    global arcctgh_function_name
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
    tan_function_name =       'tan'
    ctg_function_name =       'ctg'
    sinh_function_name =      'sinh'
    cosh_function_name =      'cosh'
    tanh_function_name =      'tanh'
    ctgh_function_name =      'ctgh'
    arcsin_function_name =    'arcsin'
    arccos_function_name =    'arccos'
    arctan_function_name =    'arctan'
    arcctg_function_name =    'arcctg'
    arcsinh_function_name =   'arcsinh'
    arccosh_function_name =   'arccosh'
    arctanh_function_name =   'arctanh'
    arcctgh_function_name =   'arcctgh'
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
    AutoSimplify = False
    
def show_names():
    print('Names of constants and varibles:')
    print('    e_constant_name =         ' + e_constant_name)
    print('    pi_constant_name =        ' + pi_constant_name)
    print('    tau_constant_name =       ' + tau_constant_name)
    print('    sin_function_name =       ' + sin_function_name)
    print('    cos_function_name =       ' + cos_function_name)
    print('    tan_function_name =       ' + tan_function_name)
    print('    ctg_function_name =       ' + ctg_function_name)
    print('    sinh_function_name =      ' + sinh_function_name)
    print('    cosh_function_name =      ' + cosh_function_name)
    print('    tanh_function_name =      ' + tanh_function_name)
    print('    ctgh_function_name =      ' + ctgh_function_name)
    print('    arcsin_function_name =    ' + arcsinh_function_name)
    print('    arccos_function_name =    ' + arccosh_function_name)
    print('    arctan_function_name =    ' + arctanh_function_name)
    print('    arcctg_function_name =    ' + arcctgh_function_name)
    print('    arcsinh_function_name =   ' + arcsinh_function_name)
    print('    arccosh_function_name =   ' + arccosh_function_name)
    print('    arctanh_function_name =   ' + arctanh_function_name)
    print('    arcctgh_function_name =   ' + arcctgh_function_name)
    print('    ln_function_name =        ' + ln_function_name)
    print('    log_function_name =       ' + log_function_name)
    print('    pow_function_name =       ' + pow_function_name)
    print('    sqrt_function_name =      ' + sqrt_function_name)
    print('    exp_function_name =       ' + exp_function_name)
    print('    rl_function_name =        ' + rl_function_name)
    print('    im_function_name =        ' + im_function_name)
    print('    arg_function_name =       ' + arg_function_name)
    print('    abs_function_name =       ' + abs_function_name)
    print('    sign_function_name =      ' + sign_function_name)
    
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
half_pi = cp( ((1, 2), pi_constant_name) )
e_constant = cp(e_constant_name)
