import os 
import sys
import math
#import matplotlib.pyplot as plt
#import numpy as np
#check if sympy is installed
try:
    from sympy import *
    from sympy.solvers import solve
    from sympy.abc import x, y, n, C, c, h, e, E
    from sympy.plotting import plot
    print('sympy loading succeed.')
    
except ModuleNoteFoundError:
    print('Sympy missing, installing.')
    os.system('pip install sympy')



#things to optimize:
#global aCoeff and bCoeff
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
        
        #calculate a and b
        a = x * P_x
        a = simplify(a)
        b = x ** 2 * Q_x
        b = simplify(b)
        print("a is: ", a)
        print("b is: ", b)
        
        #check if a and b are both constants
        #if they are, then use CauchyEuler method
        if a.is_constant() and b.is_constant():
            print("CauchyEuler Case")
            CauchyEuler(a, b)
         
        #regular singular case   
        else:
            print("RegSing Case")
            RegSing(a, b)
        
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


            
def CauchyEuler(a, b):
    formula = x ** 2 + (a - 1) * x + b
    roots = solve(formula, x)
    print("Roots:", roots)
    #in the case that the roots are two same number
    if len(roots) == 1:
        y1 = x ** roots[0]
        y2 = x ** roots[0] * ln(x)
    #when the roots are either two different real number or two complex numbers    
    else:    
        #in the case that the roots are two different real numbers
        if roots[0].is_real and roots[1].is_real:
            y1 = x ** roots[0]
            y2 = x ** roots[1]

        #in the case of roots are two complex numbers
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
    return 1



def RegSing(a, b):
    #Cn, Hn, En are calculated via their recursive formulas at a given n
    #create global lists that store all the values of Cn, Hn, En
    #initialize cList, hList, eList in every run
    global cList
    global hList
    global eList
    global y2cList
    cList = []
    hList = []
    eList = []
    #y2cList serves for the second solution y2(x)
    #since y2(x) don't share the same recursive formla as y1(x)
    y2cList = []

    c0 = 1 #by default
    cList.append(c0)
    y2cList.append(c0)

    h0 = 0 #by default
    hList.append(h0)

    e0 = 1
    eList.append(e0)

    complexCase = False
    global loopDepth
    loopDepth = int(input("Please Enter how many terms (including zero terms) you want to have: "))
    global MAXDEPTH
    MAXDEPTH = loopDepth * 4
    
    #we use termcalculator method to calculate terms of a and b
    #here aCoeff and bCoeff are the terms with zeros included
    global aCoeff
    global bCoeff
    aCoeff = termCalculator(a)
    a0 = aCoeff[0]

    bCoeff = termCalculator(b)
    b0 = bCoeff[0]

    #calculate the largest root s1
    function = FsCalculator(n, a0, b0)
    print('function:', function)
    roots = FsRootSolver(function)
    print('roots:', roots)

    if roots[0].is_real:
        print("real case")
        #s1 should be the largest root
        #thus s2 should be the minimum root respectively
        largestRoot = max(roots)
        s2 = min(roots)
    else:
        complexCase = True
        omega = abs(roots[0] - roots[1])/2 #this is w(omega) without I
        print('omega: ', omega)
        largestRoot = omega * I
        sigma = (roots[0] + roots[1])/2 #this is sigma
        print('sigma: ', sigma)
        s2 = 0
    print('largestRoot:',largestRoot)
    #nonRecursiveApproach(largestRoot, sigma)
    if a.atoms(sin, cos, sinh, cosh, tan, cot) or b.atoms(sin, cos, sinh, cosh, tan, cot): #check if P(x) or Q(x) have infinite many terms, if they do, use non-recursive approach
        nonRecursiveApproach(largestRoot, sigma)
        #print(cList)
    else:
        recursiveApproach(largestRoot, complexCase, s2, omega, sigma)



def nonRecursiveApproach(largestRoot, sigma):
    print('Non-recursive approach')
    #Cn = 1/f(s1+n) * sum_of{[(k+1)a(n-k)+b(n-k)]*Ck}(k=0-->n-1)
    n = 1
    for Cindex in range(1, loopDepth): #calculate Cns
        #initialize sumPart
        sumPart = 0
        for kValue in range(n): #for each Cn, execute the sum function
            #be careful that the original implementation may have index Error
            #theoretically, the index of aCoeffList and bCoeffList and cList should not exceeding user's input number
            #that is we only need to calculate a, b coefficients in that range
            sumPart += ((kValue + largestRoot) * aCoeff[n - kValue] + bCoeff[n - kValue]) * cList[kValue] 
        #determine cn and store it in cList
        cn = -1/FsCalculator(largestRoot + n, aCoeff[0], bCoeff[0]) * sumPart
        cList.append(cn)
        n += 1
        
    #now we calculate y(x)
    #y(x) = x**largestRoot * sum_of(Cn * x**n)
    cPart = cList[0] #the first term is c0
    for index in range(1, loopDepth):
        cPart += cList[index]* x**index
    #print('F(s) =', FsCalculator(largestRoot + n, aCoeff[0], bCoeff[0]))
    print('cList: ', cList)
    print("y(x) =", x**sigma, "* [",cPart,"]")
    functionGrapher(x**sigma * cPart)
    return True



def recursiveApproach(largestRoot, complexCase, s2, omega, sigma):
    #this k value would be n-1 or n-2 or anything
    #depending on the a0 and b0
    #it should be the maximum length of the coefficient list minus one
    listOptimizer(aCoeff)
    listOptimizer(bCoeff)
    
    maxLengthOfList = max(len(aCoeff), len(bCoeff))
    k_value_without_n = maxLengthOfList - 1

    #there are only three cases to consider, either aCoeff or bCoeff has shorter length
    #or they have the same length
    if len(aCoeff) < maxLengthOfList:
        aTerm = 0
        bTerm = bCoeff[-1]
    elif len(bCoeff) < maxLengthOfList:
        bTerm = 0
        aTerm = aCoeff[-1]
    else:
        aTerm = aCoeff[-1]
        bTerm = bCoeff[-1]
    
    #these part below are executing the formula
    #calculating both y1 and y2
    
    insideBracketPart = (n - k_value_without_n + largestRoot) * aTerm + bTerm
    indexOfC = n - k_value_without_n #this is the index of Ck

    #in complex roots case, we have to add one more step
    #that is to expand the f(s) function again to eliminate I^2
    #since I^2 is actually equals to c given that s^2 + c = 0
    #we could extract I^2 by simply: s^2 - f(s)
    if complexCase:
        Fs = FsCalculator(largestRoot + n, aCoeff[0], bCoeff[0])
        Fs_expanded = expand(Fs)
        Fs = Fs_expanded.subs(I**2, -1)
        print('f(s):',Fs)
        outsideBracketPart = -1/Fs
    else:
        outsideBracketPart = -1/FsCalculator(largestRoot + n, aCoeff[0], bCoeff[0])


           
    combinePart = insideBracketPart * outsideBracketPart #combine part is the part that without Ck
    combinePart = simplify(combinePart)
    print("Recursive formula for y1(x): Cn =", combinePart, "* C" + str(indexOfC).replace(" ", ""))

    if complexCase:
        print('Complex Case 1')
        cPart = calculateComplexCPart(combinePart, k_value_without_n, omega, sigma, loopDepth-1) #minus one because the 0-term(c0) is already given
        print(cList)
        print("y1(x) =", x**sigma, "* [",cPart[0],"]")
        functionGrapher(x**sigma * cPart[0])

        print("y2(x) =", x**sigma, "* [",cPart[1],"]")
        functionGrapher(x**sigma * cPart[1])

    else:
        cPart = calculateCPart(combinePart, k_value_without_n, cList) #minus one because the 0-term(c0) is already given
        print("y1(x) =", x**largestRoot, "* [",cPart,"]")
        functionGrapher(x**largestRoot * cPart)

        #calculate recursive formula for y2
        #determine if (s1 - s2) is a integer or a non-integer
        #then use the formula y2 = (x - x0)**s2 * sum_of(Cn(s2) * (x - x0)**n)
        s1s2Difference = largestRoot - s2
        print("s2: ", s2)

        if s1s2Difference != 0:
            print("difference is: ", s1s2Difference)
            if s1s2Difference.is_integer:
                print("s1 - s2 is an integer")
                alpha = enRecursiveCalculator(aTerm, bTerm, k_value_without_n, s2, s1s2Difference)
               
                #extract user specified number of en terms
                enPart = eList[0]
                for index in range(1, loopDepth):
                    enPart += eList[index] * x**index

                y2 = alpha * x**largestRoot * cPart * ln(x) + x**s2 * enPart
                print("y2(x) =", y2)
                functionGrapher(y2)

            else:
                print("s1 - s2 is a non-integer")
                y2insideBracketPart = (n - k_value_without_n + s2) * aTerm + bTerm
                y2indexOfC = n - k_value_without_n #this is the index of Ck
                y2outsideBracketPart = -1/FsCalculator(s2 + n, aCoeff[0], bCoeff[0])
                y2combinePart = y2insideBracketPart * y2outsideBracketPart
                y2combinePart = simplify(y2combinePart)
                print("Recursive formula for y2(x): Cn =", y2combinePart, "* C" + str(y2indexOfC).replace(" ", ""))
                y2cPart = calculateCPart(y2combinePart, k_value_without_n, loopDepth-1, y2cList) #minus one because the 0-term(c0) is already given
                print("y2(x) =", x**s2, "* [",y2cPart,"]")
                functionGrapher(x**s2 * y2cPart)
                    
                
        else:
            print("s1 = s2")
            #in this case, the second solution y2(x) is
            #y2(x) = y1(x) * ln(x-x0) + (x - x0)**s1 * sum_of_(hn(x - x0)**n)
            #where the recursive formula of hn = -1/n**2 * sum_of_{[(k + s1) * a{n-k} + b{n-k}] * hk + a{n-k} * ck} - 2/n * cn
            hnSeries = hnRecursiveCalculator(aCoeff, bCoeff, k_value_without_n, loopDepth-1, largestRoot)
            y2 = x**largestRoot * cPart * ln(x) + x**largestRoot * hnSeries
            print("y2(x) =", y2)
            functionGrapher(y2)
                

        #try:
        #    #calculate recursive formula for y2
        #    #determine if (s1 - s2) is a integer or a non-integer
        #    #then use the formula y2 = (x - x0)**s2 * sum_of(Cn(s2) * (x - x0)**n)
        #    s1s2Difference = largestRoot - s2
        #    print("s2: ", s2)

        #    if s1s2Difference != 0:
        #        print("difference is: ", s1s2Difference)
        #        if s1s2Difference.is_integer:
        #            print("s1 - s2 is an integer")
        #            enRecursiveCalculator(aTerm, bTerm, k_value_without_n, s2, s1s2Difference)
        #        else:
        #            print("s1 - s2 is a non-integer")
        #            y2insideBracketPart = (n - k_value_without_n + s2) * aTerm + bTerm
        #            y2indexOfC = n - k_value_without_n #this is the index of Ck
        #            y2outsideBracketPart = -1/FsCalculator(s2 + n, aCoeff[0], bCoeff[0])
        #            y2combinePart = y2insideBracketPart * y2outsideBracketPart
        #            y2combinePart = simplify(y2combinePart)
        #            print("Recursive formula for y2(x): Cn =", y2combinePart, "* C" + str(y2indexOfC).replace(" ", ""))
        #            y2cPart = calculateCPart(y2combinePart, k_value_without_n, y2cList) #minus one because the 0-term(c0) is already given
        #            print("y2(x) =", x**s2, "* [",y2cPart,"]")
        #            functionGrapher(x**s2 * y2cPart)
                
        #    else:
        #        print("s1 = s2")
        #        #in this case, the second solution y2(x) is
        #        #y2(x) = y1(x) * ln(x-x0) + (x - x0)**s1 * sum_of_(hn(x - x0)**n)
        #        #where the recursive formula of hn = -1/n**2 * sum_of_{[(k + s1) * a{n-k} + b{n-k}] * hk + a{n-k} * ck} - 2/n * cn
        #        hnSeries = hnRecursiveCalculator(aTerm, bTerm, k_value_without_n, loopDepth-1, largestRoot)
        #        y2 = x**largestRoot * cPart * ln(x) + x**largestRoot * hnSeries
        #        print("y2(x) =", y2)
        #        functionGrapher(y2)
              
        #except:
        #    print("Failed to determine y2(x)")



def termCalculator(inputFunction):
    #the formula for calculating terms is
    #Cn = f(nth derivative)(x0)/n!
    output = []
    for derivative in range(MAXDEPTH):
        tempTerm = diff(inputFunction, x, derivative)
        #plug in x0
        cn = tempTerm.subs(x, 0)/factorial(derivative)
        output.append(cn)
    #print('loop end')
    return output


    
def listOptimizer(inputList):
    #start from the end of the list and clear all zeros
    #if we encounter a nonzero number, stop
    #print('loop start')
    for coefficient in range(len(inputList) - 1, -1, -1):
        #print(coefficient)
        #print(inputList)
        if inputList[coefficient] == 0:
            inputList.pop()
        else:
            break    
    return True


    
def FsCalculator(s, a0, b0):
    formula = s ** 2 + (a0 - 1) * s + b0
    result = simplify(formula)
    return result



def FsRootSolver(function):
    try:
        roots = list(solve(function, n))       
    except:
        print("This f(s) equation does not have solution.\n")

    return roots



def calculateCPart(function, cIndex, coefficientList):
    output = coefficientList[0]
    #start the loop
    term = 1
    count = 1
    while count <= loopDepth-1:
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



def calculateComplexCPart(function, cIndex, omega, sigma, loopDepth):
    outputy1 = cList[0]
    outputy2 = cList[0]
    #start the loop
    for term in range(loopDepth):
        term += 1
        cComponentNum = cList[term - cIndex]
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
        cList.append(result)

        #multiply the x's
        #then convert the output to str
        #note that we have two solutions
        #that is, y1(x) and y2(x)
        y1_cTerm = (realPart * cos(omega * ln(x)) - imaginePart * sin(omega * ln(x))) * x ** term
        y2_cTerm = (realPart * sin(omega * ln(x)) + imaginePart * cos(omega * ln(x))) * x ** term
        outputy1 += y1_cTerm
        outputy2 += y2_cTerm
        #print(output)
    return [outputy1, outputy2]



def calculateHPart(function, hIndex, loopDepth):
    output = hList[0]
    #start the loop
    term = 1
    count = 1
    while count <= loopDepth:
        #print("term: ", term)
        #print(term - hIndex)

        if (term - hIndex) < 0: #if the index of hList is negative, then it will be zero 
            tempFunction = function.subs(c, 0)
            tempFunction = tempFunction.subs(h, 0)
        else:
            tempFunction = function.subs(c, cList[term - hIndex])
            tempFunction = tempFunction.subs(h, hList[term - hIndex])

        tempFunction = tempFunction.subs(C, cList[term])
        result = tempFunction.subs(n, term)

        hList.append(result)
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



#def calculateEnPart(startTerm, function, eIndex):
#    term = startTerm
#    count = startTerm
#    while count < loopDepth:
#        #print("term: ", term)
#        #print(term - eIndex)

#        if (term - eIndex) < 0: #if the index of eList is negative, then it will be zero 
#            tempFunction = function.subs(e, 0)
#        else:
#            tempFunction = function.subs(e, eList[term - eIndex])

#        result = tempFunction.subs(n, term)
#        #if a number is divided by zero
#        #we set the number to zero
#        if result == zoo:
#            print("zoo")
#            result = 0

#        eList.append(result)

#        count += 1
#        term+=1
#    print("en coefficients are: ", eList)
#    return 1



def calculateEnPart(startIndex, endIndex, s2):
    count = len(eList)
    n = startIndex
    while count <= loopDepth and n <= endIndex:
        en = -1/FsCalculator(s2 + n, aCoeff[0], bCoeff[0])
        sum = 0
        for k in range(n): #k will be range from 0 to n - 1
            if (n - k) >= len(aCoeff):
                aTerm = 0
            else:
                aTerm = aCoeff[n - k]

            if (n - k) >= len(bCoeff):
                bTerm = 0
            else:
                bTerm = bCoeff[n - k]

            sum += ((k + s2) * aTerm + bTerm) * eList[k]
        en = en * sum #combine these two parts
        en = simplify(en)

        if en != 0:
            count += 1
        eList.append(en)
        n += 1



def hnRecursiveCalculator(aTerm, bTerm, k_value_without_n, loopDepth, largestRoot):
    '''if we multiply our function with list[n] where n is an unknown integer
    then we will have an IndexError
    so the solution is that we can exclude Ck and Hk in the formula when printing
    then use a sympy variable to replace it
    once we reach the stage of substituting values
    we can then assign the sympy variable with list[n] where n is now known'''
    
    #C refers to Cn, and c refers to ck, h refers to hk
    formula = -1/n**2 * (((n - k_value_without_n + largestRoot) * aTerm + bTerm) * h + aTerm * c) - 2/n * C
    #do some tricks to print out the k value of Ck and Hk
    formulaForPrint = str(formula).replace('C', 'C' + str(n).replace(' ', ''))
    formulaForPrint = formulaForPrint.replace('c', 'C' + str(n-k_value_without_n).replace(' ', ''))
    formulaForPrint = formulaForPrint.replace('h', 'H' + str(n-k_value_without_n).replace(' ', ''))
    print("Hn =", formulaForPrint)

    return calculateHPart(formula, k_value_without_n, loopDepth)



def enRecursiveCalculator(aTerm, bTerm, k_value_without_n, s2, m):
    '''
    In this case when s1-s2 is an integer, the second solution is:
    y2(x) = alpha * y1(x) * log(x) + x**s2 * sum_of{en * x**n| n range from zero to infinity}
    '''

    ##determine the recursive formula for en first
    ##en = -1/F(s2+n) * sum_of{[(k+s2) * an-k + bn-k] * ek}
    #enFormula = -1/FsCalculator(s2 + n, aCoeff[0], bCoeff[0]) * ((n - m - 1 + s2) * aTerm + bTerm) * e
    #enFormula = simplify(enFormula)
    #enFormulaForPrint = str(enFormula).replace('e', 'e' + str(n - m - 1).replace(' ', ''))
    #print("en = ", enFormulaForPrint)

    #calculate en coefficients, starting from second term
    #since the first term e0 = 1 is given
    calculateEnPart(1, m - 1, s2)

    #now we calculate alpha 
    #where alpha = -1/m * sum_of{[(k+s2) * am-k + bm-k] * ek}
    alpha = -1/m * ((n - k_value_without_n + s2) * aTerm + bTerm) * e
    alpha = simplify(alpha)
    alphaForPrint = str(alpha).replace('e', 'e' + str(k_value_without_n).replace(' ', ''))

    #note that here we didn't consider negative index case, which might cause some bug
    #but we do so because k_value_without_n is usually positive
    alpha = alpha.subs(e, eList[k_value_without_n])
    print("alpha = ", alphaForPrint, "=", alpha)

    #em = -1/m * sum_of{[am-k * ek} - 1/m * alpha
    emFormula = -1/m * (aTerm * e) - 1/m * alpha
    emFormula = simplify(emFormula)
    emFormulaForPrint = str(emFormula).replace('e', 'e' + str(n - k_value_without_n).replace(' ', ''))
    print("em = ", emFormulaForPrint)
    emFormula = emFormula.subs(e, eList[k_value_without_n])
    print("e" + str(m) + "=", emFormula)
    #add this value to our eList
    #note that e2 cannot be determined by en formula
    #since it will be dividing by zero
    #but em formula gives us the rest of the e values
    eList.append(emFormula)
    listOptimizer(eList)

    #calculate the rest of the terms of en
    #starting from m + 1
    calculateEnPart(m + 1, MAXDEPTH, s2)
    print('En coefficients are: ', eList)
    return alpha


def modifyMAXDEPTH():
    userInputDepth = input("Please Enter the new MAXDEPTH you want to set: ")
    MAXDEPTH = userInputDepth



def functionGrapher(function):
    print("Please Enter the domain")
    lowerBond = input("Lower Bond: ")
    upperBond = input("Upper Bond: ")
    
    plot(function, (x, lowerBond, upperBond))



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
    test6 = True

    if test1:
        #test the case that P(x) = -(1+x)/(2*x), Q(x) = 1/(2*x**2)
        #in this case, we input a and b
        a = x* (-(1+x)/(2*x))
        b = x**2 * (1/(2*x**2))
        RegSing(a, b)
    
    if test2:
        #test the case that P(x) = -(1-x)/x, Q(x) = 1/x**2
        a = x * (-(1-x)/x)
        b = x**2 * (1/x**2)
        RegSing(a,b)

    if test3:
        function = x * sin(x) + 2
        print(function.atoms(sin, cos))
        if function.atoms(sin, cos):
            print('Yes')

    if test4:
        #plot(x**2, (x, -10, 10))
        plot(x*E**(-x), (x, -2, 10))

    if test5:
        a = -(1+x)/2
        b = 1/2
        RegSing(a, b)

    if test6:
        a = x* (1/x)
        b = x**2 * (-(1+x)/x**2)
        RegSing(a, b)


        
#test1()
main()