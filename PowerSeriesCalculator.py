import os 
import sys
import math
#import matplotlib.pyplot as plt
#import numpy as np
#check if sympy is installed
try:
    from sympy import *
    from sympy.solvers import solve
    from sympy.abc import x, y, n, C, c, h, e, E, k
    from sympy.plotting import plot
    print('sympy loading succeed.')
    print('----------------------\n')
    
except ModuleNoteFoundError:
    print('Sympy missing, installing.')
    os.system('pip install sympy')
    print('sympy loading succeed.')
    print('----------------------\n')


'''
Right now I have all the recursive approaches finished, what left are non-recursive versions
with one annoying bug in s1-s2=0 case
'''
class PowerSeriesCalculator(object):
    def __init__(self, P_x, Q_x):
        #Cn, Hn, En are calculated via their recursive formulas at a given n
        #we use lists to store these values
        self.cList = []
        self.hList = []
        self.eList = []

        #y2cList serves for the second solution y2(x)
        #since y2(x) don't share the same recursive formla as y1(x)
        self.y2cList = []

        self.aCoeff = []
        self.bCoeff = []

        self.cList.append(1)   #c0 = 1 by default
        self.y2cList.append(1) #c0 = 1 by default
        self.hList.append(0)   #h0 = 0 by default
        self.eList.append(1)   #e0 = 1 by default
        
        #maximum loop allowed to compute series coefficients for input P(x) and Q(x)
        #because in some cases P(x) and Q(x) could have infinite long coefficients. i.e. e^x, sin(x) etc.
        self.MAXDEPTH = 20 
        self.loopDepth = 4 #loopDepth controls how many output terms the calculator will compute

        #indicator for complex-root-case
        self.complexCase = False

        #essential variables for calculation
        self.P_x = P_x
        self.Q_x = Q_x

        self.aTerm = 0 #could be 0 or the last value in aCoefficient list
        self.bTerm = 0 #could be 0 or the last value in bCoefficient list

        self.largestRoot = 0
        self.s2 = 0

        self.alpha = 0
        self.omega = 0
        self.sigma = 0

        #Index values for recursive formula calculation
        self.cIndex = 0
        self.hIndex = 0
        self.k_value_without_n = 0

        self.cPart = None


    '''
    This function is for the one-time calculation
    for example, after powerSeries = PowerSeriesCalculator(Px, Qx)
    we can use result = powerSeries.run(4) to return the 4-terms calculation result
    '''
    def run(self, loopDepth):
        self.loopDepth = int(loopDepth)
        self.MAXDEPTH = self.loopDepth * 4

        #calculate a and b
        a = x * self.P_x
        a = simplify(a)
        b = x ** 2 * self.Q_x
        b = simplify(b)
        print("a is: ", a)
        print("b is: ", b)
        
        #check if a and b are both constants
        #if they are, then use CauchyEuler method
        if a.is_constant() and b.is_constant():
            print('\n******************')
            print("*CauchyEuler Case*")
            print('******************\n')

            self.CauchyEuler(a, b)
         
        #regular singular case   
        else:
            print('\n***********************')
            print("*Regular Singular Case*")
            print('***********************\n')
            self.RegSing(a, b)
        


    def CauchyEuler(self, a, b):
        formula = x ** 2 + (a - 1) * x + b
        roots = solve(formula, x)
        print("Roots:", roots)
        #in the case that the roots are two same number
        if len(roots) == 1:
            y1 = x ** roots[0]
            y2 = x ** roots[0] * ln(x)
        #when the roots are either two different real number or two complex numbers    
        else:    
            #roots are two different real numbers
            if roots[0].is_real and roots[1].is_real:
                y1 = x ** roots[0]
                y2 = x ** roots[1]

            #roots are two complex numbers
            else:
                #the equation is r = sigma +/- i * w
                #in order to extract sigma and iw
                #we plus root1 with root2
                #the result will remove iw and left 2 * sigma
                sigma = (roots[0] + roots[1])/2 #genius
                insideContenty1 = abs((roots[0] - sigma) * 1/I) * ln(x)
                insideContenty2 = abs((roots[1] - sigma) * 1/I) * ln(x)
                y1 = x**sigma * sin(insideContenty1)
                y2 = x**sigma * cos(insideContenty2)
        print("y1(x) is: ", y1)
        print("y2(x) is: ", y2)
        print("y(x) =","c1 *",y1,"+","c2 *",y2)
        return y1, y2



    def RegSing(self, a, b):
        #we use termcalculator method to calculate terms of a and b
        #here aCoeff and bCoeff are the terms with zeros included
        self.aCoeff = self.termCalculator(a)
        a0 = self.aCoeff[0]

        self.bCoeff = self.termCalculator(b)
        b0 = self.bCoeff[0]

        #calculate the largest root s1
        function = self.FsCalculator(n, a0, b0)
        print('function:', function)
        roots = self.FsRootSolver(function)
        print('roots:', roots)

        #extract the roots of the Fs function
        #
        if roots[0].is_real:
            print("real case")
            #s1 should be the largest root
            #thus s2 should be the minimum root respectively
            self.largestRoot = max(roots)
            self.s2 = min(roots)
        else:
            self.complexCase = True
            self.omega = abs(roots[0] - roots[1])/2 #this is w(omega) without I
            self.largestRoot = self.omega * I
            self.sigma = (roots[0] + roots[1])/2 #this is sigma
        print('largestRoot:',self.largestRoot)

        #non-recursive approach can deal with all the cases
        #but not as fast as recursive approach does
        #however, recursive approach only works under some restricted conditions
        #that is, only when a and b are polynomial
        
        if a.is_polynomial() and b.is_polynomial(): 
            print('\n>>>a and b are polynomial, use the recursive approach.<<<\n')
            self.recursiveApproach()
        else:
            print('\n>>>a and b are not polynomial, use the non-recursive approach.<<<\n')
            self.nonRecursiveApproach()
        
        #self.recursiveApproach()
        #self.nonRecursiveApproach()

    def nonRecursiveApproach(self):
        #self.listOptimizer(self.aCoeff)
        #self.listOptimizer(self.bCoeff)
        print('Coefficients of a:', self.aCoeff)
        print('Coefficients of b:', self.bCoeff)

        if self.complexCase:
            print('Non-recursive with complex roots')

            results = self.calculateComplexCPartNonRec() #minus one because the 0-term(c0) is already given
            print("y1(x) =", x**self.sigma, "* [",results[0],"]")
            self.functionGrapher(x**self.sigma * results[0])

            print("y2(x) =", x**self.sigma, "* [",results[1],"]")
            self.functionGrapher(x**self.sigma * results[1])

        else:
            print('Non-recursive with real number roots')
            self.calculateCPartNonRec()
            print("y(x) =", x**self.largestRoot, "* [",self.cPart,"]")
            print('s2 is: ', self.s2)
            self.functionGrapher(x**self.largestRoot * self.cPart)

            '''
            This part is so buggy
            and it has to be checked by professor
            The question is: do we actually have a second solution if p(x) or q(x) is not polynomial?
            '''
            #this k value would be n-1 or n-2 or anything
            #depending on the a0 and b0
            #it should be the maximum length of the coefficient list minus one
            self.listOptimizer(self.aCoeff)
            self.listOptimizer(self.bCoeff)
    
            maxLengthOfList = max(len(self.aCoeff), len(self.bCoeff))
            self.k_value_without_n = maxLengthOfList - 1

            #there are only three cases to consider, either aCoeff or bCoeff has shorter length
            #or they have the same length
            if len(self.aCoeff) < maxLengthOfList:
                self.aTerm = 0
                self.bTerm = self.bCoeff[-1]
            elif len(self.bCoeff) < maxLengthOfList:
                self.bTerm = 0
                self.aTerm = self.aCoeff[-1]
            else:
                self.aTerm = self.aCoeff[-1]
                self.bTerm = self.bCoeff[-1]
            self.calculateSecondSolution()
        return True



    def recursiveApproach(self):
        #this k value would be n-1 or n-2 or anything
        #depending on the a0 and b0
        #it should be the maximum length of the coefficient list minus one
        self.listOptimizer(self.aCoeff)
        self.listOptimizer(self.bCoeff)
    
        maxLengthOfList = max(len(self.aCoeff), len(self.bCoeff))
        self.k_value_without_n = maxLengthOfList - 1

        #there are only three cases to consider, either aCoeff or bCoeff has shorter length
        #or they have the same length
        if len(self.aCoeff) < maxLengthOfList:
            self.aTerm = 0
            self.bTerm = self.bCoeff[-1]
        elif len(self.bCoeff) < maxLengthOfList:
            self.bTerm = 0
            self.aTerm = self.aCoeff[-1]
        else:
            self.aTerm = self.aCoeff[-1]
            self.bTerm = self.bCoeff[-1]
    
        #these part below are executing the formula
        #calculating both y1 and y2
        insideBracketPart = (n - self.k_value_without_n + self.largestRoot) * self.aTerm + self.bTerm
        indexOfC = n - self.k_value_without_n #this is the index of Ck

        '''in complex roots case, we have to add one more step
        that is to expand the f(s) function again to eliminate I^2
        since I^2 is actually equals to c given that s^2 + c = 0
        we could extract I^2 by simply: s^2 - f(s)
        '''

        if self.complexCase:
            Fs = self.FsCalculator(self.largestRoot + n, self.aCoeff[0], self.bCoeff[0])
            Fs_expanded = expand(Fs)
            Fs = Fs_expanded.subs(I**2, -1)
            print('f(s):',Fs)
            outsideBracketPart = -1/Fs
        else:
            outsideBracketPart = -1/self.FsCalculator(self.largestRoot + n, self.aCoeff[0], self.bCoeff[0])

        combinePart = insideBracketPart * outsideBracketPart #combine part is the part that without Ck
        combinePart = simplify(combinePart)
        print("Recursive formula for y1(x): Cn =", combinePart, "* C" + str(indexOfC).replace(" ", ""))

        if self.complexCase:
            print('Complex')
            self.cPart = self.calculateComplexCPartRec(combinePart) 
            print("y1(x) =", x**self.sigma, "* [",self.cPart[0],"]")
            self.functionGrapher(x**self.sigma * self.cPart[0])

            print("y2(x) =", x**self.sigma, "* [",self.cPart[1],"]")
            self.functionGrapher(x**self.sigma * self.cPart[1])

        else:
            self.cPart = self.calculateCPartRec(combinePart, self.k_value_without_n, self.cList)
            print("y1(x) =", x**self.largestRoot, "* [",self.cPart,"]")
            self.functionGrapher(x**self.largestRoot * self.cPart)
            self.calculateSecondSolution()
            #try:
            #    self.calculateSecondSolution()
            #except:
            #    print("Failed to determine y2(x)")

    def calculateCPartRec(self, function, cIndex, coefficientList):
        print('Calculate normal C recursive')
        output = coefficientList[0]
        #start the loop
        term = 1
        count = 1
        while count <= self.loopDepth-1:
            #print("term: ", term)
            #print(term - cIndex)

            if (term - cIndex) < 0:
                cComponentNum = 0
            else:
                cComponentNum = coefficientList[term - cIndex]
            #print("cComponentNum: ", cComponentNum)
            result = (function * cComponentNum).subs(n, term)

            coefficientList.append(result)
            #we can treat output as a sympy object
            #which means we normally multiply the x's
            #then convert the output to str
            cTerm = result * x ** term
            #print(cTerm)
            output += cTerm
            #print(output)
            if result != 0:
                count += 1
            term+=1
        return output

    def calculateCPartNonRec(self):
        #Cn = 1/f(s1+n) * sum_of{[(k+1)a(n-k)+b(n-k)]*Ck}(k=0-->n-1)
        n = 1
        #insidePart = ((k + self.largestRoot) * self.aCoeff[n - k] + self.bCoeff[n - k]) * self.cList[k]
        for Cindex in range(1, self.loopDepth): #calculate Cns
            #print("n =", n)
            Fs = self.FsCalculator(self.largestRoot + n, self.aCoeff[0], self.bCoeff[0])
            #print('Fs is: ', Fs)
            #initialize sumPart
            sumPart = 0

                
            #for each Cn, execute the sum function
            #be careful that the original implementation may have index Error
            #theoretically, the index of aCoeffList and bCoeffList and cList should not exceeding user's desired loop times
            #which means that we only need to calculate a, b coefficients in that range
            for kValue in range(n):  
                sumPart += ((kValue + self.largestRoot) * self.aCoeff[n - kValue] + self.bCoeff[n - kValue]) * self.cList[kValue] 
            #determine cn and store it in cList
            cn = -1/Fs * sumPart
            self.cList.append(cn)
            n += 1

        print('CList is: ', self.cList)
        
        #now we calculate y(x)
        #y(x) = x**largestRoot * sum_of(Cn * x**n)
        self.cPart = self.cList[0] #the first term is c0
        for index in range(1, self.loopDepth):
            self.cPart += self.cList[index]* x**index

    def calculateComplexCPartRec(self, function):
        print('Calculate Complex recursive')
        print(function)
        outputy1 = self.cList[0]
        outputy2 = self.cList[0]
        #start the loop
        for term in range(1, self.loopDepth):
            print('term: ', term)
            print('Index: ', term - self.k_value_without_n)
            cComponentNum = self.cList[term - self.k_value_without_n]
            result = (function * cComponentNum).subs(n, term)
            result_expanded = expand(result)
            result = result_expanded.subs(I**2, -1) #now we get cn(s1)

            #calculate the real part and imagine part
            #realPart = (result1 + result2)/2
            #imaginePart = (result1 - result2)/(2* I)
            #I'll use the method of sympy
            #to extract the real and imaginary part
            realPart = re(result)
            imaginePart = im(result)
            self.cList.append(result)
            print('-----------------', result)

            #multiply the x's
            #then convert the output to str
            #note that we have two solutions
            #that is, y1(x) and y2(x)
            y1_cTerm = (realPart * cos(self.omega * ln(x)) - imaginePart * sin(self.omega * ln(x))) * x ** term
            y2_cTerm = (realPart * sin(self.omega * ln(x)) + imaginePart * cos(self.omega * ln(x))) * x ** term
            outputy1 += y1_cTerm
            outputy2 += y2_cTerm
            #print(output)
            term += 1

        return [outputy1, outputy2]


    def calculateComplexCPartNonRec(self):
        print('Calculate Complex Non-recursive')
        outputy1 = self.cList[0]
        outputy2 = self.cList[0]
        n = sympify('n')
        #n = 1
        Fs = self.FsCalculator(self.largestRoot + n, self.aCoeff[0], self.bCoeff[0])
        Fs_expanded = expand(Fs)
        Fs = Fs_expanded.subs(I**2, -1)
        print('f(s):',Fs)
        #outsideBracketPart = -1/Fs

        for nValue in range(1, self.loopDepth): #calculate Cns
            #initialize sumPart
            sumPart = 0
            for kValue in range(nValue):  #k value is from 0 to n-1
                #print('nValue - kValue: ', nValue, kValue)
                sumPart += ((kValue + self.largestRoot) * self.aCoeff[nValue - kValue] + self.bCoeff[nValue - kValue]) * self.cList[kValue] 
                #print('each part: ', sumPart)
            #print('Sum part is: ', sumPart)
            result = -1/Fs * sumPart
            result = result.subs(n, nValue)
            result_expanded = expand(result)
            result = result_expanded.subs(I**2, -1)
            #print('---------------', result)

            realPart = re(result)
            imaginePart = im(result)
            self.cList.append(result)
            print('real Part:', realPart)
            print('imaginary Part:', imaginePart)

            #multiply the x's
            #then convert the output to str
            #note that we have two solutions
            #that is, y1(x) and y2(x)
            y1_cTerm = (realPart * cos(self.omega * ln(x)) - imaginePart * sin(self.omega * ln(x))) * x ** nValue
            y2_cTerm = (realPart * sin(self.omega * ln(x)) + imaginePart * cos(self.omega * ln(x))) * x ** nValue
            outputy1 += y1_cTerm
            outputy2 += y2_cTerm

        return [outputy1, outputy2]




    def calculateSecondSolution(self):
        '''
        calculate recursive formula for y2
        there are three cases: s1 = s2, s1 - s2 is an interger, s1 - s2 is a non-integer
        '''
        s1s2Difference = self.largestRoot - self.s2
        print("s2: ", self.s2)

        if s1s2Difference != 0:
            print("difference is: ", s1s2Difference)
            if s1s2Difference.is_integer:
                print("s1 - s2 is an integer")
                alpha = self.enRecursiveCalculator(s1s2Difference)
               
                #extract user specified number of en terms
                enPart = self.eList[0]
                for index in range(1, self.loopDepth):
                    enPart += self.eList[index] * x**index

                y2 = alpha * x**self.largestRoot * self.cPart * ln(x) + x**self.s2 * enPart
                print("y2(x) =", y2)
                self.functionGrapher(y2)

            else:
                print("s1 - s2 is a non-integer")
                #use the formula y2 = (x - x0)**s2 * sum_of(Cn(s2) * (x - x0)**n)
                y2insideBracketPart = (n - self.k_value_without_n + self.s2) * self.aTerm + self.bTerm
                y2indexOfC = n - self.k_value_without_n #this is the index of Ck
                y2outsideBracketPart = -1/self.FsCalculator(self.s2 + n, self.aCoeff[0], self.bCoeff[0])
                y2combinePart = y2insideBracketPart * y2outsideBracketPart
                y2combinePart = simplify(y2combinePart)
                print("Recursive formula for y2(x): Cn =", y2combinePart, "* C" + str(y2indexOfC).replace(" ", ""))
                y2cPart = self.calculateCPartRec(y2combinePart, self.k_value_without_n, self.y2cList) #minus one because the 0-term(c0) is already given
                print("y2(x) =", x**self.s2, "* [",y2cPart,"]")
                self.functionGrapher(x**self.s2 * y2cPart)
                    
                
        else:
            print("s1 = s2")
            #in this case, the second solution y2(x) is
            #y2(x) = y1(x) * ln(x-x0) + (x - x0)**s1 * sum_of_(hn(x - x0)**n)
            #where the recursive formula of hn = -1/n**2 * sum_of_{[(k + s1) * a{n-k} + b{n-k}] * hk + a{n-k} * ck} - 2/n * cn
            hnSeries = self.hnRecursiveCalculator()
            y2 = x**self.largestRoot * self.cPart * ln(x) + x**self.largestRoot * hnSeries #is this implementation correct?
            print("y2(x) =", y2)
            self.functionGrapher(y2)


    def calculateHPart(self, function):
        output = self.hList[0]
        #start the loop
        term = 1
        count = 1
        while count <= self.loopDepth:
            #print("term: ", term)
            #print(term - hIndex)

            if (term - self.k_value_without_n) < 0: #if the index of hList is negative, then it will be zero 
                tempFunction = function.subs(c, 0)
                tempFunction = tempFunction.subs(h, 0)
            else:
                tempFunction = function.subs(c, self.cList[term - self.k_value_without_n])
                print(self.cList[term - self.k_value_without_n])
                print("hList is: ", self.hList)
                #print(tempFunction)
                #print(term - self.k_value_without_n)
                tempFunction = tempFunction.subs(h, self.hList[term - self.k_value_without_n])

            #print(self.loopDepth)
            #print('What we want:', term, 'What we have:', len(self.cList))
            tempFunction = tempFunction.subs(C, self.cList[term-1])
            result = tempFunction.subs(n, term)

            self.hList.append(result)
            #we can treat output as a sympy object
            #which means we normally multiply the x's
            #then convert the output to str
            hTerm = result * x ** term
            #print(hTerm)
            output += hTerm
            #print(output)
            if result != 0:
                count += 1
            term+=1
        return output


    '''
    if we multiply our function with list[n] where n is an unknown integer
    then we will have an IndexError
    so the solution is that we can exclude Ck and Hk in the formula when printing
    then use a sympy variable to replace it
    once we reach the stage of substituting values
    we can then assign the sympy variable with list[n] where n is now known
    '''
    def hnRecursiveCalculator(self):
        #C refers to Cn, and c refers to ck, h refers to hk
        formula = -1/n**2 * (((n - self.k_value_without_n + self.largestRoot) * self.aTerm + self.bTerm) * h + self.aTerm * c) - 2/n * C
        #print("aaaaaaaaaaaaa", self.aTerm)
        print("k: ", self.k_value_without_n)
        #do some tricks to print out the k value of Ck and Hk
        formulaForPrint = str(formula).replace('C', 'C' + str(n).replace(' ', ''))
        formulaForPrint = formulaForPrint.replace('c', 'C' + str(n-self.k_value_without_n).replace(' ', ''))
        formulaForPrint = formulaForPrint.replace('h', 'H' + str(n-self.k_value_without_n).replace(' ', ''))
        print("Hn =", formulaForPrint)
        print(formula)

        return self.calculateHPart(formula)


    def calculateEnPart(self, startIndex, endIndex):
        count = len(self.eList)
        n = startIndex
        while count <= self.loopDepth and n <= endIndex:
            en = -1/self.FsCalculator(self.s2 + n, self.aCoeff[0], self.bCoeff[0])
            sum = 0
            for k in range(n): #k will be range from 0 to n - 1
                if (n - k) >= len(self.aCoeff):
                    aTermTemp = 0
                else:
                    aTermTemp = self.aCoeff[n - k]

                if (n - k) >= len(self.bCoeff):
                    bTermTemp = 0
                else:
                    bTermTemp = self.bCoeff[n - k]

                sum += ((k + self.s2) * aTermTemp + bTermTemp) * self.eList[k]
            en = en * sum #combine these two parts
            en = simplify(en)

            if en != 0:
                count += 1
            self.eList.append(en)
            n += 1


    '''
    In this case when s1-s2 is an integer, the second solution is:
    y2(x) = alpha * y1(x) * log(x) + x**s2 * sum_of{en * x**n| n range from zero to infinity}
    where alpha = -1/m * sum_of{[(k+s2) * am-k + bm-k] * ek}
    '''
    def enRecursiveCalculator(self, m):
        #determine the recursive formula for en first
        #en = -1/F(s2+n) * sum_of{[(k+s2) * an-k + bn-k] * ek}
        enFormula = -1/self.FsCalculator(self.s2 + n, self.aCoeff[0], self.bCoeff[0]) * ((n - m - 1 + self.s2) * self.aTerm + self.bTerm) * e
        enFormula = simplify(enFormula)
        #enFormulaForPrint = str(enFormula).replace('e', 'e' + str(n - m - 1).replace(' ', ''))
        #why did I use n-m-1 instead of k_value?
        #why n-m-1 doesn't work?
        enFormulaForPrint = str(enFormula).replace('e', 'e' + str(n - self.k_value_without_n).replace(' ', ''))
        print('WWWWWWWWWWWWWWWWWW', n - m - 1, self.k_value_without_n)
        print(enFormula)
        print("en = ", enFormulaForPrint)

        #calculate en coefficients, starting from second term
        #since the first term e0 = 1 is given
        self.calculateEnPart(1, m - 1)

        #now we calculate alpha 
        #where alpha = -1/m * sum_of{[(k+s2) * am-k + bm-k] * ek}
        alpha = -1/m * ((n - self.k_value_without_n + self.s2) * self.aTerm + self.bTerm) * e
        alpha = simplify(alpha)
        alphaForPrint = str(alpha).replace('e', 'e' + str(self.k_value_without_n).replace(' ', ''))

        #note that here we didn't consider negative index case, which might cause some bug
        #but we do so because k_value_without_n is usually positive
        alpha = alpha.subs(e, self.eList[self.k_value_without_n])
        print("alpha = ", alphaForPrint, "=", alpha)

        #em = -1/m * sum_of{[am-k * ek} - 1/m * alpha
        emFormula = -1/m * (self.aTerm * e) - 1/m * alpha
        print(emFormula)
        emFormula = simplify(emFormula)
        emFormulaForPrint = str(emFormula).replace('e', 'e' + str(n - self.k_value_without_n).replace(' ', ''))
        print("em = ", emFormulaForPrint)
        emFormula = emFormula.subs(e, self.eList[self.k_value_without_n])
        print("e" + str(m) + "=", emFormula)
        #add this value to our eList
        #note that e2 cannot be determined by en formula
        #since it will be dividing by zero
        #but em formula gives us the rest of the e values
        #self.eList[m] = emFormula
        self.eList.append(emFormula)
        self.listOptimizer(self.eList)

        #calculate the rest of the terms of en
        #starting from m + 1
        self.calculateEnPart(m + 1, self.MAXDEPTH)
        print(self.eList)
        return alpha


    '''
    This function extracts coefficients of a given function
    '''
    def termCalculator(self, inputFunction):
        #the formula for calculating terms is
        #Cn = f(nth derivative)(x0)/n!
        output = []
        for derivative in range(self.MAXDEPTH):
            tempTerm = diff(inputFunction, x, derivative)

            #plug in x0
            cn = tempTerm.subs(x, 0)/factorial(derivative)
            output.append(cn)

        return output


    '''
    This function aims to remove the zero terms of the series
    for the benefit of our recursive approach.
    Input: List object
    Output: Shorter (probably) list
    Algorithm:
    start from the end of the list and clear all zeros
    if we encounter a nonzero number, terminate
    '''
    def listOptimizer(self, inputList):
        #print('loop start')
        for coefficient in range(len(inputList) - 1, -1, -1):
            #print(coefficient)
            #print(inputList)
            if inputList[coefficient] == 0:
                inputList.pop()
            else:
                break    
        return True


    '''
    This function 
    '''
    def FsCalculator(self, s, a0, b0):
        formula = s ** 2 + (a0 - 1) * s + b0
        result = simplify(formula)
        return result

    

    def FsRootSolver(self, function):
        try:
            roots = list(solve(function, n))       
        except:
            print("This f(s) equation does not have solution.\n")
        return roots


    '''
    This grapher will ask for domain of x, and graph the function in this range
    User has three tries for any unexpected errors, for instance, wrong input type, enter an unintended number by mistype,
    or the graphing of the function is not possible.
    After three tries, the function will terminate.
    User can enter 'n' to change the domain and regraph the function,
    Otherwise if they are satisfied they can enter 'y' to finish graphing.
    '''
    def functionGrapher(self, function):
        notTerminate = True
        counter = 3
        while notTerminate:
            try:
                print("Please Enter the domain")
                lowerBond = input("Lower Bond: ")
                upperBond = input("Upper Bond: ")
                
                plot(function, (x, lowerBond, upperBond))
                while True:
                    choice = input("Is this the graph you want or do you want to change the domain? Enter 'y' to finish graphing, or enter 'n' to change the domain and regraph the function: ")
                    if choice.lower() == 'y':
                        notTerminate = False
                        break
                    elif choice.lower() == 'n':
                        break
                    else:
                        print('Unknown input, please enter again.')
            except:
                if counter <= 1:
                    notTerminate = False
                else:
                    counter -= 1
                    print('Something unexpected happend when tring to graph the function, please enter again. This function will automatically exit after '+str(counter)+' times of such error.')
                

    def clean(self):
        print('Cleanning...')
        self.cList = []
        self.hList = []
        self.eList = []
        self.y2cList = []
        self.aCoeff = []
        self.bCoeff = []
        self.cList.append(1) #c0 = 1 by default
        self.y2cList.append(1)
        self.hList.append(0) #h0 = 0 by default
        self.eList.append(1) #e0 = 1 by default
        self.MAXDEPTH = 0 
        self.loopDepth = 0
        self.complexCase = False
        self.P_x = 0
        self.Q_x = 0
        self.aTerm = 0 
        self.bTerm = 0 
        self.largestRoot = 0
        self.s2 = 0
        self.alpha = 0
        self.omega = 0
        self.sigma = 0
        self.cIndex = 0
        self.hIndex = 0
        self.k_value_without_n = 0
        self.cPart = None
        print('Cleanning Finished.')



def main():
    onSwitch = True
    print("Welcome!\n")
    print("IMPORTANT Note:")
    print('-------------------------------------------------------------------------------------------------------')
    print('Make sure to use "*" symbol for any mutiplication in your input (such as 2 * x instead of 2x)')
    print('Also, please use "**" symbol for power, such as x**2 which means x^2')
    print('-------------------------------------------------------------------------------------------------------\n')

    while onSwitch:
        P_x = input("Please Enter P(x): ")
        if P_x.lower() == "quit":
            print("See you next time.")
            break
        P_x = sympify(P_x)
        Q_x = input("Please Enter Q(x): ")
        #repetition of the if statement is essential
        #since we have to make sure we can quit whenever we enter "quit"
        if Q_x.lower() == "quit":
            print("See you next time.")
            break        
        Q_x = sympify(Q_x)

        loopDepth = int(input("Please Enter how many terms (including zero terms) you want to have: "))
        calculator = PowerSeriesCalculator(P_x, Q_x)
        calculator.run(loopDepth)
        
        #-----------------------------
        #after the calculation
        #ask user if loop again or not
        choice = input("Do you want to solve more questions?(y/n): ")
        if choice.lower() == 'n' or choice.lower() == 'quit':
            print("See you next time.")
            onSwitch = False
        elif choice.lower() == 'y':
            pass
        else:
            print("Ooops, I don't get what you mean, so I'll just continue.")
            print('Hint: You can type "quit" to quit')


def test1():
    '''function1 = x**2 + 2 * x + 1
    function2 = sin(x)
    function3 = cos(x)
    #termCalculator(function1)
    #termCalculator(function2)
    termCalculator(function3)'''
    test1 = False
    test2 = False
    test3 = False
    test4 = False
    test5 = False
    test6 = False
    test7 = True
    test8 = False
    test9 = False

    #s1 - s2 = 1/2 (non-integer case)
    #y1 = x * (1 + x/3 + x**2/15 + x**3/105)
    #y2 = x**1/2 * (1 + x/2 + x**2/8 + x**3/48)
    '''
    PASSED
    '''
    if test1:
        Px = -(1+x)/(2*x)
        Qx = 1/(2*x**2)
        test = PowerSeriesCalculator(Px,Qx)
        test.run(4)
    
    #s1 = s2
    #y1 = x * (1 - x + x**2/2 - x**3/6) = x * e^-x
    '''
    ERROR
    List index of h out of range
    '''
    if test2:
        Px = -(1-x)/x
        Qx = 1/x**2
        test = PowerSeriesCalculator(Px,Qx)
        test.run(4)

    #the complex non-recursive case given on the last meeting
    '''
    RUN WITHOUT ERROR
    but the correctness is unknown
    '''
    if test3:
        Px = 1/(x - x**2)
        Qx = (1 + x)/x**2
        test = PowerSeriesCalculator(Px,Qx)
        test.run(4)

    #Cauchy Euler case
    '''
    PASSED
    '''
    if test4:
        Px = 3/x
        Qx = 2/x**2
        test = PowerSeriesCalculator(Px,Qx)
        test.run(4)

    #the example complex non-recursive case 
    '''
    PASSED
    '''
    if test5:
        Px = (1+x)/x
        Qx = 1/x**2
        test = PowerSeriesCalculator(Px,Qx)
        test.run(4)
    
    #s1 - s2 = int m case
    '''
    RUN WITHOUT ERROR
    '''
    if test6:
        Px = 1/x
        Qx = -(1+x)/x**2
        test = PowerSeriesCalculator(Px,Qx)
        test.run(4)

    #non-recursive real case
    '''
    ERROR
    hn List out of range
    '''
    if test7:
        Px = 3/(x - x**2)
        Qx = 1/x**2
        test = PowerSeriesCalculator(Px,Qx)
        test.run(4)

        '''
        OUTPUT:
        aTerms: [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
        bTerms: [1]
        '''

    '''
    RUN WITHOUT ERROR
    but why would test.run(4) return such a long solution?
    '''
    if test8:
        Px = sin(x)/x
        Qx = 1/x**2
        test = PowerSeriesCalculator(Px,Qx)
        test.run(4)

    if test9:
        Px = 2
        Qx = 1
        test = PowerSeriesCalculator(Px,Qx)
        test.run(4)

        
test1()
#main()
