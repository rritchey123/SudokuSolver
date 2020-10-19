import SudokuBoard
import Variable
import Domain
import Trail
import Constraint
import ConstraintNetwork
import time
import random
from collections import defaultdict 
import copy

class BTSolver:

    # ==================================================================
    # Constructors
    # ==================================================================

    def __init__ ( self, gb, trail, val_sh, var_sh, cc ):
        self.network = ConstraintNetwork.ConstraintNetwork(gb)
        self.hassolution = False
        self.gameboard = gb
        self.trail = trail

        self.varHeuristics = var_sh
        self.valHeuristics = val_sh
        self.cChecks = cc

    # ==================================================================
    # Consistency Checks
    # ==================================================================

    # Basic consistency check, no propagation done
    def assignmentsCheck ( self ):
        for c in self.network.getConstraints():
            if not c.isConsistent():
                return False
        return True

    """
        Part 1 TODO: Implement the Forward Checking Heuristic -------------------------------------------------------------

        This function will do both Constraint Propagation and check
        the consistency of the network

        (1) If a variable is assigned then eliminate that value from
            the square's neighbors.

        Note: remember to trail.push variables before you assign them
        Return: a tuple of a dictionary and a bool. The dictionary contains all MODIFIED variables, mapped to their MODIFIED domain.
                The bool is true if assignment is consistent, false otherwise.
    """
    def forwardChecking ( self ):
        #get all assigned variables
        assignedVars=[]
        returnDict=defaultdict()
        for c in self.network.constraints:
            for v in c.vars:
                if v.isAssigned() and v not in assignedVars:
                    assignedVars.append(v)
                    
        
        for v in assignedVars:
            curVal=v.getAssignment()
            for neighbor in self.network.getNeighborsOfVariable(v): #forsure got all the neighbors
                
                #only want to pay attention to neighbors that are changeable and have our val in their domain
                if neighbor.getDomain().contains(curVal) and neighbor.isChangeable:
                    
                    self.trail.push(neighbor) #keep track of the change you made
                    
                    neighbor.removeValueFromDomain(curVal) #correctly removes the value from the domain
                    returnDict[neighbor]=neighbor.getDomain()

                    if neighbor.size()==0: #also think this is good cuz if we take out only value, we need to backtraack
                        return (returnDict,False) #do we return the returnDict with it?
                    #look at what the return dict is used for
                    
        return (returnDict,self.assignmentsCheck())

    # =================================================================
	# Arc Consistency
	# =================================================================
    def arcConsistency( self ):
        assignedVars = []
        for c in self.network.constraints:
            for v in c.vars:
                if v.isAssigned():
                    assignedVars.append(v)
        while len(assignedVars) != 0:
            av = assignedVars.pop(0)
            for neighbor in self.network.getNeighborsOfVariable(av):
                if neighbor.isChangeable and not neighbor.isAssigned() and neighbor.getDomain().contains(av.getAssignment()):
                    neighbor.removeValueFromDomain(av.getAssignment())
                    if neighbor.domain.size() == 1:
                        neighbor.assignValue(neighbor.domain.values[0])
                        assignedVars.append(neighbor)

    
    """
        Part 2 TODO: Implement both of Norvig's Heuristics

        This function will do both Constraint Propagation and check
        the consistency of the network

        (1) If a variable is assigned then eliminate that value from
            the square's neighbors.
            
            ----simple forward checking-----
            
        (2) If a constraint has only one possible place for a value
            then put the value there.

        Note: remember to trail.push variables before you assign them
        Return: a pair of a dictionary and a bool. The dictionary contains all variables 
		        that were ASSIGNED during the whole NorvigCheck propagation, and mapped to the values that they were assigned.
                The bool is true if assignment is consistent, false otherwise.
    """
    def norvigCheck ( self ):
        tup = self.forwardChecking()
        #tup[0]=The dictionary contains all MODIFIED variables, mapped to their MODIFIED domain.
        #tup[1]=true if consistent, else false
        if tup[1]==False:
            return ({},False)
        
        returnDict=defaultdict()
        # make sure we are only looking at unassigned variables
        for c in self.network.constraints:
            counter=defaultdict(int)
            varDict=defaultdict()
            for var in c.vars:
                for val in var.domain.values:
                    counter[val]+=1
                    varDict[val]=var
            for val, tempCount in counter.items():
                if tempCount==1 and not varDict[val].isAssigned():
                    var=varDict[val]
                    self.trail.push(var)
                    var.assignValue(val)
                    var.setDomain(Domain.Domain(val)) #find way to alter domain
                    returnDict[var]= val
        return (returnDict,self.assignmentsCheck())

    """
         Optional TODO: Implement your own advanced Constraint Propagation

         Completing the three tourn heuristic will automatically enter
         your program into a tournament.
     """
    def getTournCC ( self ):
        return False

    # ==================================================================
    # Variable Selectors
    # ==================================================================

    # Basic variable selector, returns first unassigned variable
    def getfirstUnassignedVariable ( self ):
        for v in self.network.variables:
            if not v.isAssigned():
                return v

        # Everything is assigned
        return None

    """
        Part 1 TODO: Implement the Minimum Remaining Value Heuristic ------------------------------------------------------

        Return: The unassigned variable with the smallest domain
    """
    def getMRV ( self ):
        
        #seems to always return the bottom right of the board
        mrv=None
        mrvCount=10000000
        for v in self.network.variables: #for every variable
            if not v.isAssigned() and v.isChangeable() and v.size() < mrvCount: #if unassigned and changeable
                mrvCount=v.size()
                mrv=v
        return mrv

    """
        Part 2 TODO: Implement the Minimum Remaining Value Heuristic
                       with Degree Heuristic as a Tie Breaker

        Return: The unassigned variable with the smallest domain and affecting the  most unassigned neighbors.
                If there are multiple variables that have the same smallest domain with the same number of unassigned neighbors, add them to the list of Variables.
                If there is only one variable, return the list of size 1 containing that variable.
    """
    def MRVwithTieBreaker ( self ):
        ''' cnt is always 9
            -------variables domain is not changing------
        '''
        arr=[]
        cnt=1000000
        for v in self.network.variables: #for all variables
            if not v.isAssigned() and v.isChangeable() and v.size() == cnt:
                #if unassigned and has smallest domain alone
                arr.append(v)
            elif not v.isAssigned() and v.isChangeable() and v.size() < cnt:
                #if unassigned and has smallest domain as well, then append
                arr=[v]
                cnt=v.size()
        if len(arr)==0:
            return [None]
       
        maxDegree=-1
        retArr=[]
        for v in arr: #for all variables with the smallest domain
            curDegree=0
            for neighbor in self.network.getNeighborsOfVariable(v): #loop through neighbors
                if not neighbor.isAssigned() and v.isChangeable(): #if theneighbor is unassigned
                    curDegree+=1 #increment curDegree
            
            if curDegree == maxDegree: # if the degree of the curent var (v) is greater than max Degree
                retArr.append(v)
            elif curDegree > maxDegree: #if the degree of the current var (v) is the same, append to return array
                maxDegree=curDegree 
                retArr=[v]
        return retArr

    """
         Optional TODO: Implement your own advanced Variable Heuristic

         Completing the three tourn heuristic will automatically enter
         your program into a tournament.
     """
    def getTournVar ( self ):
        return None

    # ==================================================================
    # Value Selectors
    # ==================================================================

    # Default Value Ordering
    def getValuesInOrder ( self, v ):
        values = v.domain.values
        return sorted( values )

    """
        Part 1 TODO: Implement the Least Constraining Value Heuristic ----------------------------------------------------

        The Least constraining value is the one that will knock the least
        values out of it's neighbors domain.

        Return: A list of v's domain sorted by the LCV heuristic
                The LCV is first and the MCV is last
    """
    def getValuesLCVOrder ( self, v ):
        '''
        -------------------------------------------------need to implement the tie breaker-------------------------------
        if values have the same number, order those values in numerical order
        '''
        
        leDict=defaultdict()
        
        for val in v.getValues(): #for all values
            tempCount=0
            for neighbor in self.network.getNeighborsOfVariable(v): #for all neighbors
                if neighbor.getDomain().contains(val) and neighbor.isChangeable() and not neighbor.isAssigned(): #if val in domain of neighbot
                    tempCount+=1 #increment counter
            leDict[val]=tempCount #after all neighbors, store val as key and count as value
        
        #now we need to sort the dictionary by values then return the corresponding keys
        
        returnList=[i for i,j in sorted(leDict.items(), key = lambda kv:(kv[1], kv[0]))]
        

        return returnList

    """
         Optional TODO: Implement your own advanced Value Heuristic

         Completing the three tourn heuristic will automatically enter
         your program into a tournament.
     """
    def getTournVal ( self, v ):
        return None

    # ==================================================================
    # Engine Functions
    # ==================================================================

    def solve ( self, time_left=600):
        if time_left <= 60:
            return -1

        start_time = time.time()
        if self.hassolution:
            return 0

        # Variable Selection
        v = self.selectNextVariable()

        # check if the assigment is complete
        if ( v == None ):
            # Success
            self.hassolution = True
            return 0

        # Attempt to assign a value
        for i in self.getNextValues( v ):

            # Store place in trail and push variable's state on trail
            self.trail.placeTrailMarker()
            self.trail.push( v )

            # Assign the value
            v.assignValue( i )

            # Propagate constraints, check consistency, recur
            if self.checkConsistency():
                elapsed_time = time.time() - start_time 
                new_start_time = time_left - elapsed_time
                if self.solve(time_left=new_start_time) == -1:
                    return -1
                
            # If this assignment succeeded, return
            if self.hassolution:
                return 0

            # Otherwise backtrack
            self.trail.undo()
        
        return 0

    def checkConsistency ( self ):
        if self.cChecks == "forwardChecking":
            return self.forwardChecking()[1]

        if self.cChecks == "norvigCheck":
            return self.norvigCheck()[1]

        if self.cChecks == "tournCC":
            return self.getTournCC()

        else:
            return self.assignmentsCheck()

    def selectNextVariable ( self ):
        if self.varHeuristics == "MinimumRemainingValue":
            return self.getMRV()

        if self.varHeuristics == "MRVwithTieBreaker":
            return self.MRVwithTieBreaker()[0]

        if self.varHeuristics == "tournVar":
            return self.getTournVar()

        else:
            return self.getfirstUnassignedVariable()

    def getNextValues ( self, v ):
        if self.valHeuristics == "LeastConstrainingValue":
            return self.getValuesLCVOrder( v )

        if self.valHeuristics == "tournVal":
            return self.getTournVal( v )

        else:
            return self.getValuesInOrder( v )

    def getSolution ( self ):
        return self.network.toSudokuBoard(self.gameboard.p, self.gameboard.q)
