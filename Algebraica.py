import math



TabulationStep = 4

MultiplyPowers = False
# if true, (x ** 2) ** 1/2 = x



class Rational:
    def __init__(self, *args):  # object of class Rational represents a fraction
                                # it has an integer numerator and an integer denominator
        if(len(args) == 1):
            f = args[0]
            
            if isinstance(f, int):
                self.num = f
                self.den =  1
            
            elif isinstance(f, float):
                self.num, self.den = float_to_frac(f)
            
            else:
                raise RuntimeError
        
        elif(len(args) == 2):
            numerator = args[0]
            denominator = args[1]
            
            if not(isinstance(numerator, int) and isinstance(denominator, int)):
                raise RuntimeError
            
            self.num = numerator
            self.den = denominator
        
        else:
            raise RuntimeError
        
        self.reduce()
        
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
            raise RuntimeError
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
            raise RuntimeError
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
            raise RuntimeError
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
            raise RuntimeError
    def __idiv__(first, second):
        return first / second
    
    def __eq__(first, second):
        if isinstance(second, Rational):
            return ((first.num == second.num) and (first.den == second.den))
        elif isinstance(second, int):
            return ((first.num == first.den * second))
        else:
            raise RuntimeError
    def __ne__(first, second):
        return not first == second
    
    def __gt__(first, second):
        return (first.num * second.den > second.num * first.den)
    def __lt__(first, second):
        return (first.num * second.den < second.num * first.den)
    def __ge__(first, second):
        return (first.num * second.den >= second.num * first.den)
    def __le__(first, second):
        return (first.num * second.den <= second.num * first.den)    
    
    def __neg__(self):
        return Rational(-self.num, self.den)
    
    def __invert__(self):
        return Rational(self.den, self.num)
    
    def reduce(self):
        if(self.den < 0):
            self.den *= -1
            self.num *= -1
        
        gcd = math.gcd(abs(self.num), abs(self.den))
        
        self.num = self.num // gcd
        self.den = self.den // gcd
    
    def change_sign(self):
        self.num *= -1
        
    def get_sign(self):
        return sign(self.num) * sign(self.den)
        
    def copy(self):
        return Rational(self.num, self.den)
    
    def nonzero(self):
        return self.num != 0
    
    def proper_fraction(self):
        c = self.num // self.den
        f = Rational(
            self.num % self.den,
            self.den
        )
        if(c < 0):
            f += c
            c = 0
        
        return c, f
    
class Function:
    def __init__(self, params):  # class Function item describes a certain function
                                 # it has the information on function`s name
                                 # and what arguments it has
        if not isinstance(params, tuple):  # sometimes we want to make a function constant, like pi
                                           # in that case, since there is no arguments, the params
                                           # should look like ('pi',)
                                           # but it is convenient to add a shortcut, so that
                                           # params could be just 'pi'
            params = (params,)
        if not isinstance(params[0], str):  # params[0] is supposed to tell the name of the function
                                            # so that`s why it must be a string
            raise RuntimeError
        
        for item in params[1:]:
            if not isinstance(item, Complex):  # all of the function`s arguments
                                               # must be Complex
                raise RuntimeError
        
        self.name = params[0]
        self.args = params[1:]
        
    def __eq__(first, second):
        if not isinstance(second, Function):
            raise RuntimeError
        
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
        
        if(self.name == 'pow'):
            if((self.args[0] == cp(1)) or (self.args[1] == cp(0))):
                return cp(1)
            
            if(self.args[0].rl.nonzero() and self.args[0].im.is_zero()):
                if(len(self.args[0].rl.terms) == 1):
                    if(len(self.args[0].rl.terms[0].irt) == 1):
                        if(self.args[0].rl.terms[0].irt[0].name == 'pow'):
                            if(MultiplyPowers):
                                n = self.args[0].rl.terms[0].irt[0].args[0]
                                p1 = self.args[0].rl.terms[0].irt[0].args[1]
                                p2 = self.args[1]
                                
                                return cp( ((1, 1), ('pow', n, p1 * p2)) )
                            
                    elif(len(self.args[0].rl.terms[0].irt) > 1):
                        term = self.args[0].rl.terms[0]
                        p = self.args[1]
                        
                        out = cp( ((1, 1), ('pow', cp( (term.rat.num, term.rat.den) ), p)) )
                        for f in term.irt:
                            out *= cp( ((1, 1), ('pow', cp(f), p)) )
                        out.simplify()
                        
                        if(cp(self) != out):
                            return out
            
            if(self.args[0].im.nonzero() and self.args[0].rl.is_zero()):
                return cp(0, ((1, 1), ('pow', self.args[0] / cp(0, 1), self.args[1]) )) * (cp(0, 1) ** self.args[1])
            
            if((len(self.args[0].rl.terms) == 1) and (len(self.args[0].im.terms) == 1)):
                if((len(self.args[0].rl.terms[0].irt) == 0) and (len(self.args[0].im.terms[0].irt) == 0)):
                    return self.args[0] ** self.args[1]
            
            if(self.args[0].is_rational() and self.args[1].is_rational()):
                num_primes = prime_factors(self.args[0].rl.terms[0].rat.num)
                den_primes = prime_factors(self.args[0].rl.terms[0].rat.den)
                
                P = self.args[1].rl.terms[0].rat.copy()
                
                primes = []
                for n in num_primes:
                    primes.append( (n[0], Rational( n[1], 1) * P) )
                for d in den_primes:
                    primes.append( (d[0], Rational(-d[1], 1) * P) )
                    
                out = cp(1)
                for prime in primes:
                    new_pow, leftover = prime[1].proper_fraction()
                    out *= cp( ((prime[0] ** new_pow, 1), ('pow', cp(prime[0]), cp(leftover))) )
                
                if(cp(self) != out):
                    return out
            
            if((len(self.args[1].rl.terms) + len(self.args[1].im.terms)) > 1):
                out = cp(1)
                
                power = self.args[1]
                for term in power.rl.terms:
                    out *= cp( ((1, 1), ('pow', self.args[0], cp(term))) )
                for term in power.im.terms:
                    out *= cp( ((1, 1), ('pow', self.args[0], cp(term) * cp(0, 1))) )                    
                
                return out
            
            if(self.args[0].is_rational() and self.args[1].is_rational()):
                if(len(self.args[0].rl.terms) == 1):
                    if(self.args[1].rl.terms[0].rat.den == 1):
                        n = self.args[0].rl.terms[0].rat.num
                        p = self.args[1].rl.terms[0].rat.num
                        
                        if(p > 0):
                            return cp(n ** p)
                        else:
                            return cp( (1, n**(-p)) )
            
        elif(self.name == 'sqrt'):
            return cp( ((1, 1), ('pow', self.args[0], cp((1, 2)))) )
        
        elif(self.name == 'log'):
            return cp( ((1, 1), [
                ('ln', self.args[0]),
                ('pow', cp(
                    ((1, 1), ('ln', self.args[1]))
                    ), cp(-1))
            ]) )
                
        elif(self.name == 'exp'):
            if(self.args[0] == cp(0)):
                return cp(1)
            if(self.args[0] == cp(1)):
                return EulersNumber
            
            if((len(self.args[0].rl.terms) + len(self.args[0].im.terms)) > 1):
                out = cp(1)
                
                power = self.args[0]
                for term in power.rl.terms:
                    out *= cp( ((1, 1), ('exp', cp(term))) )
                for term in power.im.terms:
                    out *= cp( ((1, 1), ('exp', cp(term) * cp(0, 1))) )                    
                
                return out
            
            if(len(self.args[0].rl.terms) == 1):
                if(len(self.args[0].rl.terms[0].irt) == 1):
                    if(self.args[0].rl.terms[0].irt[0].name == 'ln'):
                        if(self.args[0].rl.terms[0].rat == 1):
                            return self.args[0].rl.terms[0].irt[0].args[0]
                        else:
                            return cp( ((1, 1), ('pow', self.args[0].rl.terms[0].irt[0].args[0],
                                                  cp(self.args[0].rl.terms[0].rat))) )    
            
        elif(self.name == 'ln'):
            if(self.args[0] == cp(1)):
                    return cp(0)
                
            elif(self.args[0] == EulersNumber):
                    return cp(1)
                
            elif(len(self.args[0].rl.terms) == 1):
                if(len(self.args[0].im.terms) == 0):
                    if(len(self.args[0].rl.terms[0].irt) == 0):
                        p = self.args[0].rl.terms[0]
                        
                        if(p.rat == 1):
                            return cp(0)
                        
                        else:
                            if(p.rat.den == 1):
                                return cp( (( 1, 1 ), ('ln', cp( (p.rat.num, 1) ))) )     
                            
                            elif(p.rat.num == 1):
                                return cp( ((-1, 1 ), ('ln', cp( (p.rat.den, 1) ))) )                           
                                
                            else:
                                return (cp( (( 1, 1 ), ('ln', cp( (p.rat.num, 1) ))) ) +
                                        cp( ((-1, 1 ), ('ln', cp( (p.rat.den, 1) ))) ))
                        
                    elif(len(self.args[0].rl.terms[0].irt) == 1):
                        p = self.args[0].rl.terms[0]
                        
                        if(p.irt[0].name == 'pow'):
                            return (cp( ((1, 1), ('ln', cp(p.rat))) ) +
                                    cp( ((1, 1), ('ln', p.irt[0].args[0])) ) * p.irt[0].args[1])
                        
                        else:
                            if(p.rat == 1):
                                return cp( ((1, 1), ('ln', cp(p.irt[0]))) )                                 
                            
                            else:
                                return (cp( ((1, 1), ('ln', cp(p.rat))) ) +
                                        cp( ((1, 1), ('ln', cp(p.irt[0]))) ))                            
                        
                    else:
                        term = self.args[0].rl.terms[0]
                        
                        out = (cp( (( 1, 1 ), ('ln', cp( (term.rat.num, 1) ))) ) +
                               cp( ((-1, 1 ), ('ln', cp( (term.rat.den, 1) ))) ))
                        for t in term.irt:
                            out += cp( ((1, 1), ('ln', cp(t))) )                 
                        
                        return out
                    
                    if(len(self.args[0].rl.terms[0].irt) == 0):
                        term = self.args[0].rl.terms[0]
                        
                        if(term.rat.den != 1):
                            out = (cp( (( 1, 1 ), ('ln', cp( (term.rat.num, 1) ))) ) +
                                   cp( ((-1, 1 ), ('ln', cp( (term.rat.den, 1) ))) ))
                            return out
                        
        elif(self.name == 'Rl'):
            if(self.args[0].is_zero()):
                return cp(0)
            
            else:
                if self.args[0].is_always_real():
                    return self.args[0]
                if self.args[0].is_always_imag():
                    return cp(0)
                
                if((len(self.args[0].rl.terms) + len(self.args[0].im.terms)) > 1):
                    out = cp(0)
                    
                    for term in self.args[0].rl.terms:
                        out += cp( ((1, 1), ('Rl', cp(term))) )
                    for term in self.args[0].im.terms:
                        # out += cp( ((1, 1), ('Rl', cp(term) * cp(0, 1) )) )
                        out += cp( ((-1, 1), ('Im', cp(term))) )
                        # these two are equivalent
                    
                    return out
                
        elif(self.name == 'Im'):
            if(self.args[0].is_zero()):
                return cp(0)
            
            else:
                if self.args[0].is_always_real():
                    return cp(0)
                if self.args[0].is_always_imag():
                    return self.args[0] / cp(0, 1)
                
                if((len(self.args[0].rl.terms) + len(self.args[0].im.terms)) > 1):
                    out = cp(0)
                    
                    for term in self.args[0].rl.terms:
                        out += cp( ((1, 1), ('Im', cp(term))) )
                    for term in self.args[0].im.terms:
                        # out += cp( ((1, 1), ('Rl', cp(term) * cp(0, 1) )) )
                        out += cp( ((1, 1), ('Rl', cp(term))) )
                        # these two are equivalent
                    
                    return out                
                        
        elif(self.name == 'sign'):
            if(self.args[0].is_rational()):
                return cp(self.args[0].rl.terms[0].rat.get_sign())
            
            if(self.args[0].is_always_real()):
                if(len(self.args[0].rl.terms) == 1):
                    if(len(self.args[0].rl.terms[0].irt) == 1):
                        s = self.args[0].rl.terms[0].irt[0].sign()
                        
                        if(s != None):
                            return cp(s)
                        
                    elif(len(self.args[0].rl.terms[0].irt) > 1):
                        out = cp(self.args[0].rl.terms[0].rat.get_sign())
                        
                        for fun in self.args[0].rl.terms[0].irt:
                            out *= cp( ((1, 1), ('sign', cp(fun))) )

                        return out
                        
        elif(self.name == 'sin'):
            if(len(self.args[0].rl.terms) == 0):
                return cp(0)
            else:
                if(self.args[0].rl.terms[0].rat < Rational(0)):
                    return -Function( ('sin', -self.args[0]) ).alternate_form()
                
                else: 
                    if(self.args[0] == cp(0)):
                        return cp(0)
                    
                    elif(self.args[0] == Pi ):
                        return cp(0)            
                    
                    elif(self.args[0] == cp((1, 2)) * Pi ):
                        return cp( (1, 1) )
                    
                    elif(self.args[0] == cp((1, 3)) * Pi ):
                        return cp( ((1, 2), ('sqrt', cp(3))) )  
                    
                    elif(self.args[0] == cp((2, 3)) * Pi ):
                        return cp( ((1, 2), ('sqrt', cp(3))) )              
                    
                    elif(self.args[0] == cp((1, 4)) * Pi ):
                        return cp( ((1, 2), ('sqrt', cp(2))) )
                    
                    elif(self.args[0] == cp((1, 6)) * Pi ):
                        return cp( (1, 2) )       
        
        elif(self.name == 'cos'):
            return Function( ('sin', Half_pi - self.args[0]) )
        
        elif(self.name == 'tan'):
            if(len(self.args[0].rl.terms) == 0):
                return cp(0)
            else:
                if(self.args[0].rl.terms[0].rat < Rational(0)):
                    return -Function( ('tan', -self.args[0]) ).alternate_form()
                
                else:            
                    if(self.args[0] == cp(0)):
                        return cp(0)
                    
                    elif(self.args[0] == Pi ):
                        return cp(0)            
                    
                    elif(self.args[0] == cp((1, 3)) * Pi ):
                        return cp( ((1, 1), ('sqrt', cp(3))) )  
                    
                    elif(self.args[0] == cp((2, 3)) * Pi ):
                        return cp( ((-1, 1), ('sqrt', cp(3))) )              
                    
                    elif(self.args[0] == cp((1, 4)) * Pi ):
                        return cp(1)
                    
                    elif(self.args[0] == cp((1, 6)) * Pi ):
                        return cp( ((1, 3), ('sqrt', cp(3))) )
        
        return cp(self)
    
    def sign(self):
        always_positive = ['pi', 'e']
        
        if(self.name in always_positive):
            return 1
        
        return None
    
    def is_always_real(self):
        if(self.name == 'pi'):
            return True
        if(self.name == 'e'):
            return True
        if((self.name == 'Rl') or
           (self.name == 'Im')):
            return True
        
        if(self.name == 'pow'):
            if(self.args[0].is_always_real() and
               self.args[1].is_always_real()):
                return True
            
        # if the argument of the following functions is real, the output is real too
        lst = ['exp', 'ln', 'sin', 'cos', 'tan', 'cotan', 'sinh', 'cosh', 'tanh', 'cotanh']
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
            raise RuntimeError
        
        self.rat = newrat
        self.irt = newirt
        
        self.sort()
        self.group()
    
    def __mul__(first, second):
        if(isinstance(first, Term) and
           isinstance(second, Term)):
            
            return Term(
                first.rat * second.rat,
                first.irt + second.irt
            )
        
        elif (isinstance(first, Term) and
              (isinstance(second, int) or isinstance(second, float))):
            return Term(
                first.rat * second,
                first.irt
            )
        
        else:
            raise RuntimeError
        
    def __eq__(first, second):
        return ((first.rat == second.rat) and compare_lists(first.irt, second.irt))
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
                new_irt.append(Function( ('pow', cp(f), cp(p)) ))
        if(i == len(self.irt) - 1):
            new_irt.append(self.irt[i])
        
        self.irt = copy_list(new_irt)
    
    def change_sign(self):
        self.rat.change_sign()    
    
    def copy(self):
        return Term(self.rat.copy(), copy_list(self.irt))
    
    def eq_irt(self, other):
        return compare_lists(self.irt, other.irt)
    
    def nonzero(self):
        return self.rat.nonzero()
    
    def has_functions(self):
        if not self.nonzero():
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
                raise RuntimeError
            
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
        
        else:
            S = make_irrational(second)
            return Irrational(*((first.terms + S.terms).copy()))
        
    def __iadd__(first, second):
        return first + second
    
    def __sub__(first, second):
        if isinstance(second, Irrational):
            return Irrational(*((first.terms + (-second).terms).copy()))
        
        else:
            S = make_irrational(second)
            return Irrational(*((first.terms + (-S).terms).copy()))
    def __isub__(first, second):
        return first - second
    
    def __mul__(first, second):
        if(isinstance(first, Irrational) and
           isinstance(second, Irrational)):
            
            out = Irrational()
            
            for t1 in first.terms:
                for t2 in second.terms:
                    t = t1 * t2
                    out.terms.append(t)
                
            out.reduce()
            return out  
        
        else:
            return first * make_irrational(second)
    def __imul__(first, second):
        return first * second
    
    def __eq__(first, second):
        if not isinstance(second, Irrational):
            raise RuntimeError
        first.reduce()
        second.reduce()
        
        if(len(first.terms) != len(second.terms)):
            return False
        
        for i in range(len(first.terms)):
            if not(first.terms[i] == second.terms[i]):
                return False
        
        return True
    
    def __neg__(self):
        res = self.copy()
        
        for i in range(len(res.terms)):
            res.terms[i].change_sign()
        
        res.reduce()
        return res
    
    def __truediv__(first, second):
        if (isinstance(second, int) or isinstance(second, float)):
            res = first.copy()
            
            for i in range(len(res.terms)):
                res.terms[i] *= second
            
            res.reduce()
            return res
        
        else:
            raise RuntimeError
    def __idiv__(first, second):
        return first / second

    def copy(self):
        return Irrational(*self.terms)
    
    def show(self, elevation, ol=False):
        if self.nonzero():
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
    
    def nonzero(self):
        for term in self.terms:
            if(term.rat.num != 0):
                return True
        
        return False
    def is_zero(self):
        return not self.nonzero()
    
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
    
    def sort(self):
        self.terms.sort(key=term_serial_index)
    
    def reduce(self):
        if(len(self.terms) > 1):
            newterms = []
            
            while(len(self.terms) > 1):
                term = self.terms.pop()
                if not term.nonzero():
                    continue

                new_rat = term.rat
                while(len(self.terms) and term.eq_irt(self.terms[0])):
                    new_rat += self.terms.pop().rat
                
                new_rat.reduce()
                newterms.append(Term(new_rat, term.irt))
                
            if(len(self.terms) > 0):
                if(self.terms[0].nonzero()):
                    newterms.append(self.terms[0])
            
            self.terms = copy_list(newterms)
        
        elif(len(self.terms) == 1):
            if(self.terms[0].nonzero()):
                pass
            else:
                self.terms = []
            
    def reduced(self):
        res = self.copy()
        res.reduce()
        
        return res
    
    

class Complex:
    def __init__(self, new_real, new_imag, new_is_simpfd=False):  # object of class Complex is just M + N i,
                                             # where M and N are class Irrational objects and
                                             # i is the square root of negative one
        if not(isinstance(new_real, Irrational) and
               isinstance(new_imag, Irrational)):
            
            raise RuntimeError
        
        self.rl = new_real.copy()
        self.im = new_imag.copy()
        self.is_simplified = new_is_simpfd
        
    def __str__(self):
        return '<<<\n' + self.show(1) + '\n>>>'
    
    def __add__(first, second):
        return Complex(
            first.rl + second.rl,
            first.im + second.im
        )
    
    def __iadd__(first, second):
        return Complex.__add__(first, second)
    
    def __mul__(first, second):
        if(isinstance(second, Complex)):
            return Complex(first.rl * second.rl - first.im * second.im,
                           first.rl * second.im + first.im * second.rl)
        
        elif(isinstance(first, Complex) and 
             (isinstance(second, int) or
             isinstance(second, float))):
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
        return Complex(
            first.rl - second.rl,
            first.im - second.im
        )
    
    def __isub__(first, second):
        return Complex.__sub__(first, second)
    
    def __invert__(self):
        res = self.copy()
        
        res.im *= -1
        
        return res
    
    def __eq__(first, second):
        return ((first.rl == second.rl) and
                (first.im == second.im))
    
    def __ne__(first, second):
        return not Complex.__eq__(first, second)
    
    def __truediv__(first, second):
        C = Complex(
            irrational((1, 1), [('pow', second * (~second), cp(-1))]),
            irrational((0, 1), [])
        )
        
        return (first * (~second)) * C
    
    def __neg__(self):
        return Complex(-self.rl, -self.im)
    
    def __pow__(first, second):
        if(isinstance(first, Complex) and isinstance(second, int)):
            if(second == 0):
                return cm(1)
            
            elif(second > 0):
                out = first.copy()
                
                for i in range(1, second):
                    out *= first
                
                return out
            
            elif(second < 0):
                inv = first.copy()
                
                for i in range(1, -second):
                    inv *= first
                
                return cp(((1, 1), [('pow', inv, cp(-1))]))
            
        elif(isinstance(first, Complex) and isinstance(second, Complex)):
            abs_val = abs(first)
            cmp_arg = first.arg()
            
            ln_first = cp(
                ((1, 1), [('ln', abs_val)])
            ) + cp(0, 1) * cmp_arg
            
            return cp(((1, 1), [('exp', ln_first * second)]), 0)
        
        else:
            raise RuntimeError
        
    def __abs__(self):
        return cp(((1, 1), [('pow', self * (~self), cp((1,2)))]))
    
    
    def is_real(self):
        return self.rl.nonzero() and self.im.is_zero()
    
    def is_imag(self):
        return self.rl.is_zero() and self.im.nonzero()
    
    def has_functions(self):
        return self.rl.has_functions() or self.im.has_functions()
    
    def is_constant(self):
        return self.rl.is_constant() and self.im.is_constant()
    
    def is_rational(self):
        if self.im.nonzero():
            return False
        
        if(len(self.rl.terms) != 1):
            return False
        if(len(self.rl.terms[0].irt) > 0):
            return False
        
        return True
    
    def is_always_real(self):
        return (self.rl.is_always_real() and
                self.im.is_zero())
    
    def is_always_imag(self):
        return self.rl.is_zero() and self.im.is_always_real()

    def arg(self):
        if self.rl.is_zero():
            if self.im.is_zero():
                return cp(0)
            
            return cp( ((1, 2), ['pi', ('sign', cp( ((1, 1), ('Im', self)) ) )]) )
        
        return cp(((1, 1), [
            ('arctan', self.imag() * cp(((1, 1), [
                ('pow', self.real(), cp(-1))
            ])))
        ]))
    
    def show(self, elevation, one_line=False):
        out = ''
        tab = ' ' * elevation * TabulationStep
        
        rl_nz = self.rl.nonzero()
        im_nz = self.im.nonzero()
        
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
                    out += tab
                    out += '{\n'
                    
                    if(one_line):
                        out += self.rl.show(0, ol=one_line)
                    else:
                        out += self.rl.show(elevation + 1, ol=one_line)
                else:
                    out += self.rl.show(elevation, ol=one_line)
                
                if(im_nz):
                    out += '\n' + tab + '}'
            
        if (rl_nz and im_nz):
            out += '\n'
            out += tab + '+\n'        
        
        if(im_nz):
            if(self.im.is_constant()):
                out += tab + '{ '
                out += self.im.show(0, ol=True)
                out += ' } i'
                
            else:
                out += tab + '{\n'
                out += self.im.show(elevation + 1, ol=one_line)
                out += '\n' + tab + '} i'
            
        if not(rl_nz or im_nz):
            out += self.rl.show(elevation, ol=one_line)
        
        return out
    
    def nonzero(self):  # if this function returns True, it does not always mean that the number is nonzero
                        # what it actually means is that the number MIGHT be nonzero
        return self.rl.nonzero() or self.im.nonzero()
    def is_zero(self):
        return not self.nonzero()
    
    def copy(self):
        return Complex(self.rl.copy(),
                       self.im.copy(),
                       new_is_simpfd = self.is_simplified)
    
    def real(self):  # note that the 'real' part can have complex functions inside it
                     # and thus actually have an imaginary component
        return cp(self.rl)
    
    def imag(self):
        return cp(self.im)
    
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
        raise RuntimeError
    
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
                                    # like a = 'pi'
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
                                        # such as a = ('pow', cp(3), cp(1/2))
                                        # OUT = (3/1 + 0/1 i) ** (1/2 + 0/1 i)
                                        # then a[0] is the name of the function
                                        #  and a[1:] are its` arguments
                for i in range(1, len(a)):
                    if not isinstance(a[i], Complex):
                        raise RuntimeError
                
                return irrational((1, 1), [a])
            
            else:
                if(len(a) == 2):
                    if(isinstance(a[0], int) and isinstance(a[1], int)):  # if a is a tuple denoting fraction
                                                                          # like a = (2, 3)
                                                                          # OUT = 2/3 + 0/1 i
                        return irrational((a[0], a[1]), [])
                    
                    elif(isinstance(a[0], tuple) and isinstance(a[1], str)):  # if a is a tuple representing a constant with a coefficient
                                                                              # such as a = ((2, 3), 'pi')
                                                                              # OUT = 2/3 * pi + 0 i
                        return irrational((a[0][0], a[0][1]), [(a[1],)])
                    
                    elif(isinstance(a[0], tuple) and isinstance(a[1], tuple)):
                                             # if a is a tuple describing a function with a coefficient
                                             # for example, a = ((4, 9), ('ln', cp(7)))
                                             # OUT = 4/9 * ln[7] + 0/1 i
                        if(isinstance(a[1][0], str)):
                            for i in range(1, len(a[1]) - 1):
                                if not isinstance(a[1][i], Complex):
                                    raise RuntimeError
                        else:
                            raise RuntimeError
                        
                        return irrational((a[0][0], a[0][1]), [a[1]])
                    
                    elif (isinstance(a[0], tuple) and isinstance(a[1], list)):
                                             # if a is a tuple representing a product of functions with a coefficient
                                             # for instance, a = ((2, 3), ['pi', ('sqrt', cp(2))])
                                             # OUT = 2/3 * pi * sqrt[2] + 0/1 i
                        return irrational((a[0][0], a[0][1]), a[1])
                    
                    else:
                        raise RuntimeError
                
                else:
                    raise RuntimeError
        
        elif(isinstance(a, list)):  # if a is a list, then each element of the
                                    # list represents a term of a sum
                                    # each of these terms is to be read by conventions
            out = irrational((0, 1), [])
            
            for term in a:
                out += make_irrational(term)
                
            return out
    
        else:
            raise RuntimeError
        
    else:
        raise RuntimeError
    
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
    
    return True
    
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
        raise RuntimeError
    
    return c



EulersNumber = cp('e')
Pi = cp('pi')
Half_pi = cp( ((1, 2), 'pi') )
