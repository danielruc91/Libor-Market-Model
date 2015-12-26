
# coding: utf-8

# In[188]:

import bisect

class RateStructure:
    
    def __init__( self, time, RateAtTime ):
        self.t = time
        self.rates = RateAtTime.copy()
        
    def __getitem__( self, i ):
        return self.rates[ i ]
    
    def getRates( self ):
        return self.rates.copy()
    
    def __iter__( self ):
        return iter( self.rates )
    
    def __str__( self ):
        return "At time " + str( self.t ) + ", the rate structure is " + str( self.rates )
    
    # For comparision methods, we have the following assumption:
    #    for <, <=, >, >= operators, we use time as the key
    #                                that is, the result is equivalent to directly compare the two time members
    #                                for the purpose of sorting
    #    for ==, != operators, use both rate and time member as the key
    
    def __lt__( self, other ):
        return self.t < other.t
    
    def __le__( self, other ):
        return self.t <= other.t
    
    def __eq__( self, other ):
        return self.t == other.t and self.rates == other.rates
    
    def __nq__( self, other ):
        return not( self == other )
    
    def __gt__( self, other ):
        return not( self <= other )
    
    def __ge__( self, other ):
        return not( self < other )
    
    def __len__( self ):
        return len( self.rates )


# In[189]:

class SimResult:
    
    def __init__( self ):
        self.__result = []
        
    def __len__( self ):
        return len( self.__result )
        
    def add( self, aRateStruct ):
        """Add one simulation of RateStructure to data members.
        Data member would be maintained as an ordered list. It's recommended that RateStructure is 
        added in ascending order.
        """
        if( self.__result == [] or self.__result[ -1 ] < aRateStruct ):
            self.__result.append( aRateStruct )
        else:
            bisect.insort_right( self.__result, aRateStruct )
        
    def __getitem__( self, i ):
        return self.__result[ i ]
    
    def get( self, time, i ):
        """Access the ith forward rate at time. If time is not in data, linearly interpolate a value """
        assert( self.__sorted( ) )
        # __result is assumed to be an ordered list, since the only method to add elements is add and it would add
        # elements in order. If it turns out that __result is not sorted, then some procedure must have modified the
        # the data structure illegaly
        j = bisect.bisect_right( self.__result, RateStructure( time, [] ) )
        return self.__interpolate( time, j, i )
            
    def __interpolate( self, time, timeInd, rateInd ):
        """Linearly interpolate two forwards rates to get a forward rate at the desired time point.
        accepts 3 argument, time( double ), timeInd and rateInd( int )
        return the interpolated result
        """
        if( timeInd == 0 ):
            if( time == self.__result[ 0 ].t ):
                return self.__result[ 0 ][ rateInd ]
            else:
                raise Exception( "Trying to access " + str( time ) + " while the earlist time point is " +                                 str( self.__result[ 0 ].t ) )
                
        elif( timeInd == len( self.__result ) ):
            if( time == self.__result[ -1 ].t ):
                return self.__result[ -1 ][ rateInd ]
            else:
                raise Exception( "Trying to access " + str( time ) + " while the oldest time point is " +                                 str( self.__result[ 0 ].t ) )
                
        else:
            ratio = ( time - self.__result[ timeInd - 1 ].t ) / ( self.__result[ timeInd ].t - self.__result[ timeInd - 1 ].t )
            diff = self.__result[ timeInd ][ rateInd ] - self.__result[ timeInd - 1 ][ rateInd ]
            return self.__result[ timeInd - 1 ][ rateInd ]+ ratio * diff
    
    def __sorted( self ):
        return all( self.__result[ i ] <= self.__result[ i + 1 ] for i in range( len( self.__result ) - 1 ) )
    
    def __iter__( self ):
        return iter( self.__result )


# In[190]:

from abc import ABCMeta
from abc import abstractmethod
import bisect
import math
from scipy.integrate import quad
from scipy.optimize import minimize

class volatility( metaclass = ABCMeta ):
    """Volatility generating class"""
    
    @abstractmethod
    def get( self, curTime, criticalTimePoint ):
        ...
        
    @abstractmethod
    def getI( self, curTime, criticalTimePoint, i ):
        ...
    
    @abstractmethod
    def calibrate( self, tarVol, criticalTimePoint, tau ):
        ...
        
class vol2( volatility ):
    """Volatility generator representing function 2"""
    def __init__( self, volList ):
        self.__vol = volList.copy()
    
    def setVol( self, volList ):
        self.__vol = volList.copy()
    
    def get( self, curTime, criticalTimePoint ):
        assert( len( self.__vol ) == len( criticalTimePoint ) )
        k = bisect.bisect_left( criticalTimePoint, curTime )
        return [ self.__vol[ i ] for i in range( 0, len( self.__vol ) - k ) ]
    
    def getI( self, curTime, criticalTimePoint, i ):
        return self.__vol[ i ]
    
    def calibrate( self, tarVol, criticalTimePoint, tau ):
        assert( len( tarVol ) == len( self.__vol ) == len( criticalTimePoint ) == len( tau ) )
        for i in range( 0, len( tarVol ) ):
            temp = sum( [ self.__vol[ j ] ** 2 * tau[ i - 1 - j ] for j in range( 0, i ) ] )
            temp = criticalTimePoint[ i ] * tarVol[ i ] ** 2 - temp
            self.__vol[ i ] = math.sqrt( temp / tau[ i ] )
    
    def __str__( self ):
        return str( self.__vol )
    

class vol6( volatility ):
    """Volatility generator representing function 6"""
    def __init__( self, a, b, c, d ):
        self.a, self.b, self.c, self.d = a, b, c, d
    
    def get( self, curTime, criticalTimePoint ):
        k = bisect.bisect_left( criticalTimePoint, curTime )
        return [ self.__volFun( ele - curTime ) for ele in criticalTimePoint[ k : ] ]
    
    def getI( self, curTime, criticalTimePoint, i ):
        return self.__volFun( criticalTimePoint[ i ] - curTime )
    
    def __str__( self ):
        return "sigma(t) = [ a( T - t ) + d ] * exp( -b * ( T - t ) ) + c, a = " + str( self.a ) +                  " b = " + str( self.b ) + " c = " + str( self.c ) + " d = " + str( self.d )
    
    def calibrate( self, tarVol, criticalTimePoint, tau ):
        assert( len( tarVol ) == len( criticalTimePoint ) == len( tau ) )
        return minimize( self.__errToTarVol, [ self.a, self.b, self.c, self.d ], ( tarVol, criticalTimePoint ),                         bounds = [ (0.0, None), (0.0,None), (None, None), (None,None) ], method = 'L-BFGS-B' )
        
    def __volFun( self, t ):
        return ( self.a * t + self.d ) * math.exp( -self.b * t) + self.c
    
    def __varFun( self, t ):
        return self.__volFun( t ) ** 2
    
    def __integralFun( self, lower, upper ):
        return quad( self.__varFun, lower, upper )[ 0 ]
    
    def __toBSVar( self, ctp ):
        """ctp represents criticalTimePoint"""
        temp = [ self.__integralFun( 0, ctp[ 0 ] ) ]
        for i in range( len( ctp ) - 1 ):
            temp.append( self.__integralFun( ctp[ i ], ctp[ i + 1 ] ) + temp[ -1 ] )
        
        return [ ele[ 0 ] / ele[ 1 ] for ele in zip( temp, ctp ) ]
    
    def __errToTarVol( self, x, tarVol, ctp ):
        """ctp represents criticalTimePoint
        It is worth noticing that when ever __errToTarVol is invoked, a,b,c,d are set as the x. This might not be a favourable
        feature in large applications as it relies on the user of this client. However, we adopted this design for the 
        following reason:
            1. __errToTarVol is a private function and (currently is invoked only in calibrating, this means we really want to
            set the 4 parameters
            2. Our appliaction is a bit small and this design fits our need in a much straightforward way.
        """
        self.a, self.b, self.c, self.d = x
        temp = self.__toBSVar( ctp )
        return sum( [ ( ele[ 0 ] ** 2 - ele[ 1 ] ) ** 2 for ele in zip( tarVol, temp ) ] )
        

class vol7( volatility ):
    def __init__( self, a, b, c, d, phiList ):
        self.a, self.b, self.c, self.d = a, b, c, d
        self.__phi = phiList.copy()
        
    def setPhi( self, phiList ):
        self.__phi = phiList.copy()
    
    def get( self, curTime, criticalTimePoint ):
        assert( len( criticalTimePoint ) == len( self.__phi ) )
        k = bisect.bisect_left( criticalTimePoint, curTime )
        return [ self.__phi[ i ] * self.__volFun( criticalTimePoint[ i ] - curTime ) for i in range( k, len( self.__phi) ) ]
    
    def getI( self,  curTime, criticalTimePoint, i ):
        return self.__phi[ i ] * self.__volFun( criticalTimePoint[ i ] - curTime )
    
    def calibrate( self, tarVol, criticalTimePoint, tau ):
        """To calibrate a vol7 object, first calibrate a,b,c,d 4 parameters in the vol6 way,
        then calibrate __phi values to make an exact match"""
        assert( len( tarVol ) == len( criticalTimePoint ) == len( tau ) == len( self.__phi ) )
        
        minimize( self.__errToTarVol, [ self.a, self.b, self.c, self.d ], ( tarVol, criticalTimePoint ),                  bounds = [ (0.0, None), (0.0,None), (None, None), (None,None) ], method = 'L-BFGS-B' )
        temp = self.__toBSVar( criticalTimePoint )
        
        assert( len( temp ) == len( criticalTimePoint ) )
        
        for i in range( len( self.__phi ) ):
            self.__phi[ i ] = tarVol[ i ] * math.sqrt( criticalTimePoint[ i ] / temp[ i ] )
    
    def __str__( self ):
        return "sigma(t) = phi(i) * [ a( T(i) - t ) + d ] * exp( -b * ( T(i) - t ) ) + c, a = " + str( self.a ) +                  ", b = " + str( self.b ) + ", c = " + str( self.c ) + ", d = " + str( self.d ) + "\nphi = " +                 str( self.__phi )
    
    def __volFun( self, t ):
        return ( self.a * t + self.d ) * math.exp( -self.b * t) + self.c
    
    def __varFun( self, t ):
        return self.__volFun( t ) ** 2
    
    def __integralFun( self, lower, upper ):
        return quad( self.__varFun, lower, upper )[ 0 ]
    
    def __toBSVar( self, ctp ):
        """ctp represents criticalTimePoint"""
        temp = [ self.__integralFun( 0, ctp[ 0 ] ) ]
        for i in range( len( ctp ) - 1 ):
            temp.append( self.__integralFun( ctp[ i ], ctp[ i + 1 ] ) + temp[ -1 ] )
        
        return [ ele[ 0 ] / ele[ 1 ] for ele in zip( temp, ctp ) ]
    
    def __errToTarVol( self, x, tarVol, ctp ):
        """ctp represents criticalTimePoint"""
        self.a, self.b, self.c, self.d = x
        temp = self.__toBSVar( ctp )
        return sum( [ ( ele[ 0 ] ** 2 - ele[ 1 ] ) ** 2 for ele in zip( tarVol, temp ) ] )


# In[191]:

import numpy as np

class correlator( metaclass = ABCMeta ):
    
    @abstractmethod
    def get( self, curTime, criticalTimePoint ):
        ...
        
    @abstractmethod
    def calibrate( self, tarCorr, criticalTimePoint, tau ):
        ...


class parametricCorrelator( correlator ):
    
    def get( self, curTime, criticalTimePoint ):
        k = bisect.bisect_left( criticalTimePoint, curTime )
        n = len( criticalTimePoint ) - k
        corrMatrix = np.identity( n )
        
        for i in range( n ):
            for j in range( i + 1, n ):
                temp = self.corrFun( criticalTimePoint[ i ], criticalTimePoint[ j ] )
                corrMatrix[ i ][ j ] = temp
                corrMatrix[ j ][ i ] = temp
        
        return corrMatrix
    
    def calibrate( self, tarCorr, criticalTimePoint, tau ):
        assert( isinstance( tarCorr, numpy.matrix ) and tarCorr.shape[ 0 ] ==                tarCorr.shape[ 1 ] == len( criticalTimePoint ) )
        minimize( self.__errToTarCorr, self.getPara(), ( tarCorr, criticalTimePoint ),                  method = 'L-BFGS-B' )
        
    def __errToTarCorr( self, x, tarCorr, ctp ):
        self.setPara( x )
        return numpy.sum( numpy.square( tarCorr - self.get( 0, ctp ) ) )
        
    @abstractmethod
    def getPara( self ):
        ...
    
    @abstractmethod
    def setPara( self, x ):
        ...
        
    @abstractmethod
    def getBound( self ):
        ...
    
    @abstractmethod
    def corrFun( self, Ti, Tj ):
        ...

class Corr1( parametricCorrelator ):
    """Correlation calculator representing function 1
        rho(i,j) = exp( -beta | T(i) - T(j) |)
    """
    
    def __init__( self, beta ):
        assert( beta >= 0 )
        self.beta = beta
        
    def getPara( self ):
        return [ self.beta ]
    
    def setPara( self, x ):
        beta = x
    
    def getBound( self ):
        return [ ( 0.0, None ) ]
        
    def corrFun( self, Ti, Tj ):
        return math.exp( -self.beta * abs( Ti - Tj ) )
    
class Corr2( parametricCorrelator ):
    """Correlation calculator representing function 2
        rho(i,j) = rho(inf) + ( 1 - rho( inf ) ) * exp( -beta | T(i) - T(j) |)
        rho( inf ) represents the lowest possible correlation
    """
    def __init__( self, beta, rhoInf ):
        assert( beta >= 0 and rhoInf > 0.0 and rhoInf < 1.0 )
        self.beta, self.rhoInf = beta, rhoInf
    
    def getPara( self ):
        return [ self.beta, self.rhoInf ]
    
    def setPara( self, x ):
        self.beta, self.rhoInf = x
    
    def getBound( self ):
        return [ ( 0.0, None ), ( 0.0, 1.0 ) ]
        
    def corrFun( self, Ti, Tj ):
        return self.rhoInf + ( 1 - self.rhoInf ) * math.exp( -self.beta * abs( Ti - Tj ) )
    


# In[192]:

class corrToBSIV:
    def __init__( self ):
        self.__theList = []
    
    def append( self, i, j, coef ):
        self.__theList.append( ( i, j, coef ) )
    
    def __call__( self, corrMatrix ):
        return sum( [ ele[ 2 ] * corrMatrix[ ele[ 0 ], ele[ 1 ] ] for ele in self.__theList ] )

def prodIfNotNone( a, b ):
    return a * b if ( a is not None and b is not None ) else None


# In[193]:

class LiborMarketModel( metaclass = ABCMeta ):
    """LiborMarketModel is a base class for all model simulators. It's an abstract base class and all derived class
    shall implement the 'simulate' function.
    A libor market model shall have 4 elements, a volatility calculator, a correlation calculator, a RateStructure
    indicating the initial forwards, and a list of price sensitive time points. These 4 elements would define the 
    specification of a market model.
    However, it's not sufficient to simulate with only these 4 elements. Therefore, it's implemented as an abstract 
    base class and derived class shall implement the 'simulate' function.
    Besides, a model should be calibrated to market data before put into use. Therefore, it's required that derived 
    classes implement a 'calibrate' function, and invoke this function at initilization. 
    """
    
    def __init__( self, iniState, volCalc, corrCalc, criticalTimePoint, tau, numeraireIndex ):
        """
        iniState             :    initial state of term structures, should belong to RateStructure class
        volCalc              :    volatility calculator, should be derived class of volatility
        corrCalc             :    correlation calculator, should be derived class of correlator
        criticalTimePoint    :    list of reset dates
        tau                  :    list of forward rate period
        numeraireIndex       :    the index of numeraire forward rate
        capletBSVol          :    list of BS implied volatility for caplets, should have the same length as iniState
        swaptionBSVol        :    list of swaption BS implied volatilites, in the form ( T0, Tn, iv )
                                  largest Tn shall not exceed the largest crititcal time point
        """
        # I miss c++ when I have to do type checks
        assert( isinstance( iniState, RateStructure ) )
        assert( isinstance( volCalc, volatility ) )
        assert( isinstance( corrCalc, correlator) ) 
        # Numeraire forward rate should be one of the forward rates
        assert( numeraireIndex >= 0 and numeraireIndex < len( iniState ) )
        # length of Critical Time Point shall equal to len to iniState
        assert( len( criticalTimePoint ) == len( iniState ) == len( tau ) )
        
        # Initialize: initial state, critical time point, numeraire index
        self.__rate0, self.__ctp, self.__tau, self.__numInd = iniState, criticalTimePoint, tau, numeraireIndex
        
        # Initialze: volatility calculator, correlation calculator
        self.__vol, self.__corr = volCalc, corrCalc
    
    def getDrift( self, rateStruct ):
        """Drift calculator. rateStruct is a RateStructure, containing time and forward rates information.
        Return a list of float, which is the drift for different forward rates. For dead forward rates, use 'None'
        to hold the position.
        """
        assert( isinstance( rateStruct, RateStructure ) )
        assert rateStruct.t <= self.__ctp[ self.__numInd ], "Numeraire Forward Rate is dead"
        
        """Cache the result of volatility and correlation matrix. Drift is usually called in simulations and volatility 
        and correlation matrix would be used in generating random variables
        """
        self.tempVol = self.__vol.get( rateStruct.t, self.__ctp )
        self.tempCor = self.__corr.get( rateStruct.t, self.__ctp )
        
        assert( len( self.tempVol ) == self.tempCor.shape[ 0 ] == self.tempCor.shape[ 1 ] )
        
        return self.driftGivenVolCor( rateStruct, self.tempVol, self.tempCor )
    
    def driftGivenVolCor( self, rateStruct, volTerm, corMatrix ):
        
        assert rateStruct.t <= self.__ctp[ self.__numInd ], "Numeraire Forward Rate is dead"
        
        k = bisect.bisect_left( self.__ctp, rateStruct.t )
        
        #Use 'None' to represent dead forward rates. A forward rate is not considered 'dead' on the reset date.
        drift = [ None for i in range( k ) ]
        
        for i in range( k, self.__numInd ):
            temp = 0.0
            for j in range( i + 1, self.__numInd + 1 ):
                temp1 = rateStruct[ j ] * self.__tau[ j ]
                temp -= temp1 * corMatrix[ i ][ j ] * volTerm[ j - k ] / ( 1 + temp1 )
            
            drift.append( volTerm[ i - k ] * temp )
        
        for i in range( max( k, self.__numInd ), len( rateStruct ) ):
            temp = 0.0
            for j in range( self.__numInd + 1, i + 1 ):
                temp1 = rateStruct[ j ] * self.__tau[ j ]
                temp += temp1 * corMatrix[ i - k, j - k ] * volTerm[ j - k ] / ( 1 + temp1 )
            
            drift.append( volTerm[ i - k ] * temp )
            
        assert len( rateStruct ) == len( drift ) == len( self.__ctp ), str(rateStruct)
        
        return drift
        
    
    @abstractmethod
    def simulate( self, finishTime ):
        """Simulate from the time 0 to finishTime.
        finishTime is a float, representing the furtherest time to go in the future.
        Return a SimResult
        """
        ...
    
    def calibrate( self, capletBSVol, swaptionBSVol ):
        self.__vol.calibrate( capletBSVol, self.__ctp, self.__tau )
        
        zcpPrc = [ 1.0 / ( 1.0 + self.__tau[ i ] * self.__rate0[ i ] )  for i in range( len( self.__tau ) ) ]
        for i in range( 1, len( zcpPrc ) ):
            zcpPrc[ i ] *= zcpPrc[ i - 1 ]
        
        ivCalc   = [ self.__getIVCalc( ele, zcpPrc ) for ele in swaptionBSVol ]
        
        tarBSVol = [ ele[ 0 ] * ele[ 2 ] for ele in swaptionBSVol ]
        
        minimize( self.__errToBSSswaptionVol, self.__corr.getPara(), ( tarBSVol, ivCalc ),                  method = 'L-BFGS-B', bounds = self.__corr.getBound( ) )
    
    def __errToBSSswaptionVol( self, corrPara, tarBSVol, ivCalc ):
        self.__corr.setPara( corrPara )
        corrMatrix = self.__corr.get( 0.0, self.__ctp )
        return sum( [ ( ele[ 0 ] - ele[ 1 ]( corrMatrix ) ) ** 2 for ele in zip( tarBSVol, ivCalc ) ] )
    
    def __getIVCalc( self, swaptionIV, zcpPrc ):
        T0, Tn, bsIV = swaptionIV
        start = bisect.bisect_left( self.__ctp, T0 )         # start is included
        end   = bisect.bisect_right( self.__ctp, Tn )        # end is not included
        wht   = [ 0.0 for ele in zcpPrc ]
        temp  = sum( [ self.__tau[ i ] * zcpPrc[ i ] for i in range( start, end ) ] )
        for i in range( start, end ):
            wht[ i ] = self.__tau[ i ] * zcpPrc[ i ] / temp
        swaprate = sum( [ wht[ i ] * self.__rate0[ i ] for i in range( start, end ) ] )
        
        res = corrToBSIV( )
        for i in range( start, end ):
            for j in range( start, end ):
                c = lambda t: self.__vol.getI( t, self.__ctp, i ) * self.__vol.getI( t, self.__ctp, j ) 
                coef = quad( c, 0, T0 )[ 0 ]
                res.append( i, j, coef * wht[ i ] * wht[ j ] * self.__rate0[ i ] * self.__rate0[ j ] / ( swaprate ** 2 ) )
        
        return res
    
    def getIniState( self ):
        return self.__rate0
    
    def getInd( self ):
        return self.__numInd
    
    def getTau( self ):
        return self.__tau
    
        
class ShortJumper( LiborMarketModel ):
    """ShortJumper is a type of LiborMarketModel. It utilize Euler Scheme to simulate the evolvement of term structures.
    To receive accurate result in Euler Scheme, each time the model only move forward a small step in time. Indeed, it is 
    agreed that the time step should be between 0 and 1.0/12 (1 month):
        if time step is too large ( >1.0/12 ), it is automatically reset to 1.0 / 12;
        if time step is too small ( <=0 ), it is reset to 0.01;
        if time step is not provided during initilization, the default is 0.01
    """
    
    maxTimeStep = 1.0 / 12
    
    def __init__( self, iniState, volCalc, corrCalc, criticalTimePoint, tau, numeraireIndex, timeStep = 0.01 ):
        """Initialization of ShortJumper.
        To initialize a ShortJumper, there must be certain components:
            
            timeStep             :    time step of simulations
        """
        super(ShortJumper, self).__init__( iniState, volCalc, corrCalc, criticalTimePoint, tau, numeraireIndex )
        
        self.__dt     = ShortJumper.maxTimeStep if ( timeStep > ShortJumper.maxTimeStep ) else         ( timeStep if timeStep > 0 else 0.01 )
       
    def simulate( self, finishTime ):
        
        curRate, res = self.getIniState(), SimResult()
        res.add( curRate )
        
        while( curRate.t <= finishTime ):
            
            drift  = self.getDrift( curRate )
            
            mu     = self.__dt * ( np.array( [ ele for ele in drift if ele is not None ] ) - 0.5 * np.square( self.tempVol ) )
            cov    = self.__dt * np.diag( self.tempVol ) * self.tempCor * np.diag( self.tempVol )
            x      = np.exp( np.random.multivariate_normal( mu, cov ) )
            
            temp   = [ None for ele in drift if ele is None ] + x.tolist()
            
            curRate= RateStructure( curRate.t + self.__dt, [ prodIfNotNone( ele[0], ele[1] ) for ele in zip( curRate, temp ) ] )
            res.add( curRate )
            
        return res


class PredCorrect( LiborMarketModel ):
    """PredCorrect is a type of LiborMarketModel. It utilize Predictor-Corrector scheme to simulate the evolvement of 
    term structures.
    """
    
    def __init__( self, iniState, volCalc, corrCalc, criticalTimePoint, tau, numeraireIndex, timeStep = 0.25 ):
        """Initialization of PredCorrect.
        To initialize a PredCorrect, there must be certain components:
            
            timeStep             :    time step of simulations
        """
        super(PredCorrect, self).__init__( iniState, volCalc, corrCalc, criticalTimePoint, tau, numeraireIndex )
        
        self.__dt     = timeStep
       
    def simulate( self, finishTime ):
        
        curRate, res = self.getIniState(), SimResult()
        res.add( curRate )
        
        while( curRate.t <= finishTime ):
            
            drift   = self.getDrift( curRate )
            
            zero    = np.array( [ 0.0 for ele in drift if ele is not None ] )
            cov     = self.__dt * np.diag( self.tempVol ) * self.tempCor * np.diag( self.tempVol )
            x       = np.random.multivariate_normal( zero, cov )
            
            mu0     = self.__dt * ( np.array( [ ele for ele in drift if ele is not None ] ) - 0.5 * np.square( self.tempVol ) )
            temp    = [ ele for ele in drift if ele is None ] + np.exp( mu0 + x ).tolist()
            
            tmpRate = RateStructure( curRate.t, [ prodIfNotNone( ele[ 0 ], ele[ 1 ] ) for ele in zip( curRate, temp ) ] )
            
            drift   = self.driftGivenVolCor( tmpRate, self.tempVol, self.tempCor )
            mu1     = self.__dt * ( np.array( [ ele for ele in drift if ele is not None ] ) - 0.5 * np.square( self.tempVol ) )
            
            mu      = ( mu0 + mu1 ) / 2
            temp    = [ ele for ele in drift if ele is None ] + np.exp( mu + x ).tolist()
            
            curRate = RateStructure( curRate.t + self.__dt, [ prodIfNotNone( ele[0], ele[1] ) for ele in zip( curRate, temp ) ] )
            
            res.add( curRate )
            
        return res
        
class IterPredCorrect( LiborMarketModel ):
    """IterPredCorrect is a type of LiborMarketModel. It utilize Iterative Predictor-Corrector scheme to simulate 
    the evolvement of term structures.
    """
    
    def __init__( self, iniState, volCalc, corrCalc, criticalTimePoint, tau, numeraireIndex, timeStep = 0.25 ):
        """Initialization of IterPredCorrect.
        To initialize a IterPredCorrect, there must be certain components:
            
            timeStep             :    time step of simulations
        """
        super( IterPredCorrect, self).__init__( iniState, volCalc, corrCalc, criticalTimePoint, tau, numeraireIndex )
        
        self.__dt     = timeStep
       
    def simulate( self, finishTime ):
        
        tau, numInd, curRate, res = self.getTau(), self.getInd(), self.getIniState(), SimResult()
        res.add( curRate )
        
        while( curRate.t <= finishTime ):
            
            drift   = self.getDrift( curRate )
            
            zero    = np.array( [ 0.0 for ele in drift if ele is not None ] )
            cov     = self.__dt * np.diag( self.tempVol ) * self.tempCor * np.diag( self.tempVol )
            x       = np.random.multivariate_normal( zero, cov )
            
            k       = len( drift ) - len( zero )
            temp    = [ None for ele in drift ]                 # Create a list having the same length as drift
            
            for i in range( numInd, k - 1, -1 ):
                mu = 0
                for j in range( i + 1, numInd + 1 ):
                    temp1 = temp[ j ] * tau[ j ]
                    mu   -= temp1 * self.tempCor[ i - k ][ j - k ] * self.tempVol[ j - k ] / ( 1 + temp1 )
                mu = ( mu + drift[ i ] ) / 2
                temp[ i ] = curRate[ i ] * np.exp( self.__dt * ( mu - 0.5 * self.tempVol[ i - k ] ** 2 ) + x[ i - k ] )
            
            for i in range( numInd + 1, len( drift ) ):
                mu = 0
                for j in range( numInd, i ):
                    temp1 = temp[ j ] * tau[ j ]
                    mu   += temp1 * self.tempCor[ i - k ][ j - k ] * self.tempVol[ j - k ] / ( 1 + temp1 )
                mu = ( mu + drift[ i ] ) / 2
                temp[ i ] = curRate[ i ] * np.exp( self.__dt * ( mu - 0.5 * self.tempVol[ i - k ] ** 2 ) + x[ i - k ] )
            
            curRate = RateStructure( curRate.t + self.__dt, temp )
            
            res.add( curRate )
            
        return res  


# In[194]:

ctp = [ 0.25 * ( i + 1 ) for i in range( 10 )]
r0  = RateStructure( 0.0, [ 0.005, 0.01, 0.016, 0.02, 0.024, 0.027, 0.03, 0.033, 0.035, 0.037 ] )
dt  = [ 0.25 for i in range( 10 ) ]
ind = 5

tarVol = [0.38,0.385,0.396666667,0.405,0.412,0.411666667,0.408571429,0.405,0.398888889,0.392 ]

tarCor = [ (0.25, 1.0, 0.3 ), ( 0.25, 1.25, 0.32 ), (0.25, 1.5, 0.33), (0.5, 1.0, 0.31), (0.5, 1.5, 0.34 ) ]

vola = vol2( [ 0.2 for i in range( 10 ) ] )
volb = vol6( 0.2, 0.3, 0.4, 0.5 )
volc = vol7( 0.2, 0.3, 0.4, 0.5, [ 0.2 for i in range( 11 ) ] )

cora = Corr1( 0.9 )
corb = Corr2( 0.9, 0.9 )


# In[195]:

import random
"""
a = vol7( 0.2, 0.3, 0.4, 0.5, [ random.random() for i in range( 11 ) ] )
print( a )
print( a.get( 0.1, [0.25*i for i in range( 1, 12 )] ) )
a.calibrate( [0.34, 0.39, 0.42, 0.38, 0.35, 0.33, 0.32, 0.30, 0.28, 0.25, 0.24 ], [0.25*i for i in range( 1, 12 )], \
            [0.25 for i in range( 11 ) ] )
print( a )

print( a.get( 0.1, [0.25*i for i in range( 1, 12 )]))
print( a.get( 0.3, [0.25*i for i in range( 1, 12 )]))
"""

a = vol2( [ random.random() for i in range( 10 ) ] )
a.calibrate([0.38,0.385,0.396666667,0.405,0.412,0.411666667,0.408571429,0.405,0.398888889,0.392 ], [0.25*i for i in range( 1, 11 )], \
            [0.25 for i in range( 10 ) ] )


# In[196]:

testLMM = IterPredCorrect( r0, vola, cora, ctp, dt, ind, 0.25 )

testLMM.calibrate( tarVol, tarCor )


# In[197]:

print( testLMM.getDrift( r0 ) )


# In[198]:

aRes = testLMM.simulate( 1.5 )


# In[199]:

import os

os.chdir( r'C:\Users\Danie_000\Desktop')

outFile = open( 'temp.txt', 'w' )

for ele in aRes:
    print( ele, file = outFile )
    
outFile.close()


# In[200]:

bisect.bisect_left( ctp, 0.25 + 0.25 )


# In[201]:

print( aRes[ 2 ])


# In[202]:



testLMM.getDrift( aRes[ 2 ] )

