from __future__ import division

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 15:55:12 2018

@author: yik
"""
"""
This module is used to set up and visulize a model with the two scalar potential having the tree level form:
V(h,s) = lh(h^2-vh^2)^2+ls(s^2-vs^2)^2+lsh*s^2*h^2

The model is to be defined as a class of model(ms, tanb, sint), which needs three inputs as physical parameters

The higgs mass and vev are set default as 125GeV and 246Gev

Bare parameters entering the tree level potential can be called by the method model().info()
"""

import numpy as np
import matplotlib.pyplot as plt
#from cosmoTransitions import generic_potential
import generic_potential_cwd as generic_potential

"""        
        self.Y1 = 2*80.4/246
"""

"""
Adding two usage functions under the class 'model':
    
    prettyPrintTcTrans(self): 
        print critical temperatures of the phase transitions,
        and the correspoinding high and low T phases
    
    plotPhases2D(self, **plotArgs):
        plot the 2D phase diagram (v_s, v_h) at different temperatures.
        such a plot shows the phase trajectary as T changes.
"""

class model(generic_potential.generic_potential):

    def init(self, vh2, vs2, lh, ls, lsh, yd, v2re):
        
        self.vh2 = vh2
        self.vs2 = vs2
        self.lh = lh
        self.ls = ls
        self.lsh = lsh
        self.yd = yd
        #self.ch = ch
        #self.cs = cs
                
        
        self.Ndim = 2
        
        self.renormScaleSq = v2re

        self.Y1 = 2*80.4/246
        self.Y2 = 2*(91.2**2 - 80.4**2)**0.5/246
        self.Yt = 2**0.5*172.4/246 #just curious about 1/6

        self.nwt = 4
        self.nzt = 2
        
        self.nwl = 2
        self.nzl = 1
        
        self.nx = 3
        
        self.nt = 12
        self.nd = 1


    def test(self):
        print "hi:", self.lsh

            
    def V0(self, X):
        X = np.asanyarray(X)
        phi1,phi2 = X[...,0], X[...,1]
        #phi1 -> h, phi2 -> s
        r = 0.25 * (self.lh*(phi1**2-self.vh2)**2 + self.ls*(phi2**2-self.vs2)**2 + self.lsh*phi1**2*phi2**2)
        #r = self.lh*phi1**4 - 2*self.lh*self.vh2*phi1**2 + self.ls*phi2**4 - 2*self.ls*self.vs2*phi2**2 + self.lsh*phi1**2*phi2**2
        return r
  
    
    def boson_massSq(self, X, T):
        X = np.asanyarray(X)
        phi1,phi2 = X[...,0], X[...,1]

        ls = self.ls
        lh = self.lh
        lsh = self.lsh
        muh2 = lh*self.vh2
        mus2 = ls*self.vs2
        
        # Tong: Mass corrections. See paper 2.2. Does A have mass correction?
        
        ringh = (3.*self.Y1**2./16. + self.Y2**2./16. + self.lh/2 + self.Yt**2./4. + self.lsh/12.)*T**2.*phi1**0.
        rings = (self.ls/3. + self.lsh/6. + self.yd**2/12.)*T**2.*phi1**0.
        ringwl = 11.*self.Y1**2.*T**2.*phi1**0./6.
        ringbl = 11.*self.Y2**2.*T**2.*phi1**0./6.
        ringchi = (3.*self.Y1**2./16. + self.Y2**2./16. + self.lh/2. + self.Yt**2./4. + self.lsh/12.)*T**2.*phi1**0.
        

        # keep only the h and S dependent terms (take A=0)
        a = -self.lh*self.vh2 + 3.*self.lh*phi1**2. + 0.5*self.lsh*phi2**2 + ringh #mh^2
        b = -self.ls*self.vs2 + 3.*self.ls*phi2**2. + 0.5*self.lsh*phi1**2. + rings #ms^2
        c = -self.ls*self.vs2 + self.ls*phi2**2. + 0.5*self.lsh*phi1**2. #mA^2
        mab = self.lsh*phi2*phi1
        mbc = 0.
        mac = 0.
        #Eigenvalues
        A = c  #for calculation of eigenvalues of ((a,c),(c,b))
        B = 0.5*(a+b)
        C = 0.5*np.sqrt(a**2-2*a*b+b**2+4*mab**2)
        
        mwl = 0.25*self.Y1**2.*phi1**2. + ringwl
        mwt = 0.25*self.Y1**2.*phi1**2.
        
        mzgla = 0.25*self.Y1**2.*phi1**2. + ringwl
        mzglb = 0.25*self.Y2**2.*phi1**2. + ringbl
        mzgc = - 0.25*self.Y1*self.Y2*phi1**2.
        mzglA = .5*(mzgla + mzglb)
        mzglB = np.sqrt(.25*(mzgla-mzglb)**2. + mzgc**2.)
        
        mzgta = 0.25*self.Y1**2.*phi1**2.
        mzgtb = 0.25*self.Y2**2.*phi1**2.
        mzgtA = .5*(mzgta + mzgtb)
        mzgtB = np.sqrt(.25*(mzgta-mzgtb)**2. + mzgc**2.)
        
        mzl = mzglA + mzglB
        mzt = mzgtA + mzgtB
        mgl = mzglA - mzglB
        mgt = mzgtA - mzgtB
        
        mx = -muh2 + lh*phi1**2. + 0.5*lsh*phi2**2. + ringchi
 
        M = np.array([A, B+C, B-C, mwl, mwt, mzl, mzt, mgl, mgt, mx])
        

#        M = np.array([100, 100, 100, 100, 100, 100, 100, 100, 100])

        M = np.rollaxis(M, 0, len(M.shape))

        dof = np.array([1, 1, 1, self.nwl, self.nwt, self.nzl, self.nzt, self.nzl, self.nzt, self.nx])

        c = np.array([1.5, 1.5, 1.5, 5./6., 5./6., 5./6., 5./6., 5./6., 5./6., 1.5])#check Goldstones
               
        return M, dof, c
        
 
    def fermion_massSq(self, X):
        X = np.asanyarray(X)
        phi1, phi2 = X[...,0], X[...,1]

        mt = 0.5*self.Yt**2.*phi1**2.
        md = 2*self.yd**2.*phi2**2.
        M = np.array([mt, md])

        M = np.rollaxis(M, 0, len(M.shape))

        dof = np.array([self.nt, self.nd])

        return M, dof


    def approxZeroTMin(self):
        # There are generically two minima at zero temperature in this model,
        # and we want to include both of them.
        return [np.array([246., 0.1])]


    def forbidPhaseCrit(self, X):
        """
        forbid negative phases for both h and s
        """
        return any([np.array([X])[...,0] < -5.0, np.array([X])[...,1] < -5.0])
 
    
    def V0s0(self, phi):
        #r =  -0.5*self.m12*phi**2 + 0.25*self.l1*phi**4
        r =  -0.5*self.vh2*phi**2 + 0.25*self.lh*phi**4
        return r

       
    def info(self):
        print 'Bare parameters:'
        print 'vh^2=',self.vh2,',','vs^2=',self.vs2
        print 'lambh=',self.lh,',','lambs=',self.ls,',','lamb_sh=',self.lsh
        print 'physical parameters:'
        #print 'ms=',self.ms,',', 'tanb=',self.tanb,',','sint=',self.sint

        
    def prettyPrintTcTrans(self):
        if self.TcTrans is None:
            raise RuntimeError("self.TcTrans has not been set. "
                "Try running self.calcTcTrans() first.")
        if len(self.TcTrans) == 0:
            print("No transitions for this potential.\n")
        for trans in self.TcTrans:
            trantype = trans['trantype']
            if trantype == 1:
                trantype = 'First'
            elif trantype == 2:
                trantype = 'Second'
            print("%s-order transition at Tc = %0.4g" %
                  (trantype, trans['Tcrit']))
            print("High-T phase:\n  key = %s; vev = %s" %
                  (trans['high_phase'], trans['high_vev']))
            print("Low-T phase:\n  key = %s; vev = %s" %
                  (trans['low_phase'], trans['low_vev']))
            print("Energy difference = %0.4g = (%0.4g)^4" %
                  (trans['Delta_rho'], trans['Delta_rho']**.25))
            print("")

 
    def plotPhases2D(self, **plotArgs):
        import matplotlib.pyplot as plt
        if self.phases is None:
            self.getPhases()
        for key, p in self.phases.items():
            plt.plot(p.X[...,1], p.X[...,0], **plotArgs)
        plt.xlabel(R"$v_s(T)$")
        plt.ylabel(R"$v_h(T)$")
    
"""
Here defines some useful functions to visualize the tree level potential

v0h(m, tanb) used to show V-h plot at T=0 when s = s_vev. 
            m is the input model. tanb = s_vev/h_vev

vsh(m, box, T) is used to show the 1-loop effective V-hs plot at temperature T. 
            m: input model. 
            box = (xmin, xmax, ymin, ymax) to set the range of the plot
            T: temperature of the potential
"""

#V0 as function of h at s = s_vev
def v0h(m,tanb):
    plt.plot(np.arange(-300, 300, 2),m.V0([np.arange(-300, 300, 2),246*tanb]))

def vh(m, box, s, T, n=50):
    xmin, xmax = box
    X = np.linspace(xmin, xmax, n)
    Y = np.linspace(s, s, n)
    XY = np.zeros((n, 2))
    XY[...,0], XY[...,1] = X, Y    
    plt.plot(X, m.Vtot(XY, T))
    plt.show()
 
def vs(m, box, T, n=50):
    xmin, xmax = box
    X = np.linspace(xmin, xmax, n)
    Y = np.linspace(0, 0, n)
    XY = np.zeros((n, 2))
    XY[...,0], XY[...,1] = Y, X    
    plt.plot(X, m.Vtot(XY, T))
    plt.show()
    
def vh1T(m, box, s, T, n=50):
    xmin, xmax = box
    X = np.linspace(xmin, xmax, n)
    Y = np.linspace(s, s, n)
    XY = np.zeros((n, 2))
    XY[...,0], XY[...,1] = X, Y    
    plt.plot(X, m.V1T_from_X(XY, T))
    plt.show()

def vs1T(m, box, T, n=50):
    xmin, xmax = box
    X = np.linspace(xmin, xmax, n)
    Y = np.linspace(0, 0, n)
    XY = np.zeros((n, 2))
    XY[...,1], XY[...,0] = X, Y    
    plt.plot(X, m.V1T_from_X(XY, T))
    plt.show()
             
"""
make use of the Vtot method in generic_potential
"""

def vsh(m, box, T, n=50, clevs=200, cfrac=1., **contourParams):
    xmin,xmax,ymin,ymax = box
    X = np.linspace(xmin, xmax, n).reshape(n,1)*np.ones((1,n))
    Y = np.linspace(ymin, ymax, n).reshape(1,n)*np.ones((n,1))
    XY = np.zeros((n, n, 2))
    XY[...,0], XY[...,1] = X, Y
    Z = m.Vtot(XY, T)
    # Take the log to amplify variation around minimum
    #Z = np.log(Z)
    minZ, maxZ = min(Z.ravel()), max(Z.ravel())
    if minZ > 0:
        Z = np.log(Z)
    minZ, maxZ = min(Z.ravel()), max(Z.ravel())
    N1 = np.linspace(minZ, minZ+(maxZ-minZ)*cfrac, clevs)
    #plt.figure()
    cset1 = plt.contourf(X,Y,Z,N1, **contourParams)
    N2 = np.linspace(minZ, minZ+(maxZ-minZ)*cfrac, clevs/4)
    cset2 = plt.contour(X,Y,Z,N2, colors='white')
    for c in cset2.collections:
        c.set_linestyle('dotted')
        c.set_linewidth(0.8)
    plt.axis(box)
    plt.xlabel('h')
    plt.ylabel('S')
    plt.colorbar(cset1)
    #plt.show()
