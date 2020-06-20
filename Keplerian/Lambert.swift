//
//  Lambert.swift
//  KSPTelemetry
//
//  Created by Mike Muszynski on 5/24/20.
//  Copyright © 2020 Mike Muszynski. All rights reserved.
//

import Foundation

typealias LambertSolution = (Vector3D, Vector3D, Double)

class LambertSolver {
    
    enum ConfigurationError: Error {
        case invalidOrbits
    }
    
    var mu: Double
    var position1: Vector3D
    var velocity1: Vector3D?
    var position2: Vector3D
    var velocity2: Vector3D?
    var dt: Double
    var maxRevs: Int = 0
    var prograde = true
    var transferType = 1
    
    var debug = true
    
    init(position1: Vector3D, position2: Vector3D, dt: Double, mu: Double = CelestialBody.kerbol.gravitationalParameter) {
        self.position1 = position1
        self.position2 = position2
        self.dt = dt
        self.mu = mu
    }
    
    init(orbit1: Orbit, orbit2: Orbit, departureTime: KSPDate, travelTime: Double) throws {
        guard orbit1.centralBody == orbit2.centralBody else { throw ConfigurationError.invalidOrbits }
        let arrivalTime = departureTime.timeIntervalSinceReferenceDate + travelTime
    
        self.mu = orbit1.centralBody.gravitationalParameter
        (self.position1, self.velocity1) = orbit1.cartesian(atTime: Double(departureTime.timeIntervalSinceReferenceDate))
        (self.position2, self.velocity2) = orbit2.cartesian(atTime: Double(arrivalTime))
        
        self.dt = Double(travelTime)
        
        
        if debug {
            print("******INITIALIZING LAMBERT SOLVER*********")
            print("travel time:  \(dt)")
            print("mu:           \(mu)")
            print("departure:    [\(position1.x), \(position1.y), \(position1.z), \(velocity1?.x ?? 0), \(velocity1?.y ?? 0), \(velocity1?.z ?? 0)]")
            print("arrival:      [\(position2.x), \(position2.y), \(position2.z), \(velocity2?.x ?? 0), \(velocity2?.y ?? 0), \(velocity2?.z ?? 0)]")
            print("lambert([\(position1.x), \(position1.y), \(position1.z), \(velocity1?.x ?? 0), \(velocity1?.y ?? 0), \(velocity1?.z ?? 0)], [\(position2.x), \(position2.y), \(position2.z), \(velocity2?.x ?? 0), \(velocity2?.y ?? 0), \(velocity2?.z ?? 0)], \(dt), mu=\(mu))")
            print("******************************************")
        }
    }
    
    func solve() ->  [LambertSolution] {
        return mattsLambertSovlerFromCU()
    }
    
    private func mattsLambertSovlerFromCU(debug: Bool = false) -> [LambertSolution] {
        if debug {
            print("lambert([\(position1.x),\(position1.y),\(position1.z), \(velocity1!.x),\(velocity1!.y),\(velocity1!.z)], [\(position2.x),\(position2.y),\(position2.z),\(velocity2!.x),\(velocity2!.y),\(velocity2!.z)], \(dt), mu=\(mu))")
        }
        
        /*
         #------------------------------------------------------------------
         # lambert()
         #
         #Reference:
         #    http://ccar.colorado.edu/imd/2015/documents/LambertHandout.pdf
         #
         #Inputs:
         #    X_0: initial state vector; velocities required to calculate C3
         #        and vInf
         #    X_f: final state vector; velocities required to calculate C3
         #        and vInf
         #    dt_0: transfer time in seconds
         #keywords:
         #    DM: direction of motion (if you want to force it. It is
         #        nominaly calculated by the script.)
         #Outputs:
         #
         #-----------------------------------------------------------------
             r_0 = X_0[0:3]
             r_f = X_f[0:3]
             n0 = X_0[0:3]/norm(X_0[0:3])
             nf = X_f[0:3]/norm(X_f[0:3])
             r0 = norm(X_0[0:3])
             rf = norm(X_f[0:3])
             sqrtMu = sqrt(mu)
         
         */
        
        let r_0 = self.position1
        let r_f = self.position2
        let n0 = self.position1.normalized
        let nf = self.position2.normalized
        
        let r0 = self.position1.magnitude
        let rf = self.position2.magnitude
        
        let sqrtMu = sqrt(mu)
        let dt0 = self.dt
        
        /*
         
         ####################################################################
         #    Approximate delta nu.
         #
         #    This section assumes that the transfer orbit is in the
         #    ecliptic. Since we only use it to find DM, it's a good enough
         #    approximation for reasonable transfers.
         #
         ####################################################################
         */
        
        //nu_0 = arctan2(n0[1],n0[0])
        //nu_f = arctan2(nf[1],nf[0])
        //delta_nu = nu_f - nu_0

        var nu_0 = atan2(n0.y, n0.x)
        var nu_f = atan2(nf.y, nf.x)
        var delta_nu = nu_f - nu_0
        
        /*if delta_nu < 0:
            delta_nu += 2*pi
         if delta_nu > 2*pi:
            delta_nu -= 2*pi
         */
        
        if delta_nu < 0 {
            delta_nu += 2 * .pi
        } else if delta_nu > 2 * .pi {
            delta_nu -= 2 * .pi
        }
        
        /*
             #this delta nu is NOT an approximiation
             #since it uses the dot product
        if abs(delta_nu) < pi:
            DM = 1
        else:
            DM = -1*/
        let dm = Double(abs(delta_nu) < .pi ? 1 : -1)
        
        //cos_delta_nu = n0.dot(nf)
        let cos_delta_nu = n0 * nf
        
        //A = DM * sqrt(r0*rf*(1+cos_delta_nu))
        let A = dm * sqrt(r0 * rf * (1 + cos_delta_nu))
        
        /*
        if delta_nu == 0 or A == 0:
            print('Trajectory Cannot be Computed') */
        if delta_nu == 0 || A == 0 { return [] }
        
        //c2 = 1./2.
        //c3 = 1./6.
        var c2: Double = 1/2
        var c3: Double = 1/6
        
        /*

             if revs == 0:
                 psi = 0
                 psi_up = 4*pi**2
                 psi_low = -4*pi**2
             else:
                 psi_up = (2*(revs + 1)*pi)**2 - 0.01
                 psi_low = (2*revs*pi)**2 + 0.01
                 psi = (psi_up + psi_low)/2
         
         */
        
        var psi, psi_up, psi_low: Double
        if maxRevs == 0 {
            psi = 0
            psi_up = 4 * .pi * .pi
            psi_low = -4 * .pi * .pi
        } else {
            psi_up = pow(2 * Double(maxRevs + 1) * .pi, 2) - 0.01
            psi_low = pow(2 * Double(maxRevs) * .pi, 1) + 0.01
            psi = (psi_up + psi_low) / 2
        
            func tof(_ psi: Double) -> Double {
                //sqrtPsi = sqrt(psi)
                let sqrtPsi = sqrt(psi)

                //c2 = (1 - cos(sqrtPsi))/psi
                let c2 = (1 - cos(sqrtPsi)) / psi

                //sqrtC2 = sqrt(c2)
                let sqrtC2 = sqrt(c2)
                
                //c3 = (sqrtPsi - sin(sqrtPsi))/sqrtPsi**3
                let c3 = (sqrtPsi - sin(sqrtPsi)) / pow(sqrtPsi, 3)
                
                //y = r0 + rf + A*(psi*c3-1)/sqrtC2
                let y = r0 + rf + A * (psi * c3 - 1) / sqrtC2
                
                //xi = sqrt(y/c2)
                let xi = sqrt(y/c2)
                
                //sqrtY = sqrt(y)
                let sqrtY = sqrt(y)
                
                //TOF = (xi**3*c3+A*sqrtY)/sqrtMu
                let TOF = (pow(xi, 3) * c3 + A * sqrtY) / sqrtMu

                //return TOF
                return TOF
            }
            
            var lowBound = psi_low + 0.01
            var highBound = psi_up - 0.01
            var testPsi = stride(from: lowBound, to: highBound, by: (-lowBound + highBound) / 15)
            
            var testTOF = testPsi.map { (tof: tof($0), psi: $0) }
            var sortTOF = testTOF.sorted { $0.tof < $1.tof }
            
            var minimizingPsi: Double = sortTOF[0].psi
            
            while sortTOF[1].tof - sortTOF[0].tof > 1e-6 {
                //Here, matt used argmin and subtracted/added one to the index
                //Doesn't that potentially give a negative index?
                lowBound = sortTOF[0].psi
                highBound = sortTOF[1].psi
                testPsi = stride(from: lowBound, to: highBound, by: (lowBound + highBound / 15))
                
                testTOF = testPsi.map { (tof: tof($0), psi: $0) }
                sortTOF = testTOF.sorted { $0.tof < $1.tof }
                
                minimizingPsi = sortTOF[0].psi
            }
            
            let minimumTOF = tof(minimizingPsi)
            
            if minimumTOF > dt {
                fatalError("time of flight is required to be larger than minimum supplied. this should throw?")
            }
            
            
            if transferType == 3 {
                psi_up = minimizingPsi
            } else {
                psi_low = minimizingPsi
            }

            psi = (psi_low + psi_up)/2
        }
        
        /*
                 #define a function for calculating TOF from psi quickly
                 def TOF(psi):
                     sqrtPsi = sqrt(psi)
                     c2 = (1 - cos(sqrtPsi))/psi
                     sqrtC2 = sqrt(c2)
                     c3 = (sqrtPsi - sin(sqrtPsi))/sqrtPsi**3
                     y = r0 + rf + A*(psi*c3-1)/sqrtC2
                     xi = sqrt(y/c2)
                     sqrtY = sqrt(y)
                     TOF = (xi**3*c3+A*sqrtY)/sqrtMu
                     return TOF

                 #find psi associated with minimum TOF so we know where the
                 #cutoff between type III and type IV is
                 testPsi = linspace(psi_low + 0.01,psi_up - 0.01,15)
                 testTOF = TOF(testPsi)
                 sortTOF = sorted(testTOF)
                 while sortTOF[1] - sortTOF[0] > 1e-6:
                     testPsi = linspace(
                         testPsi[argmin(testTOF)-1],
                         testPsi[argmin(testTOF)+1],
                         15)
                     testTOF = TOF(testPsi)
                     sortTOF = sorted(testTOF)
                     minimizingPsi = testPsi[argmin(testTOF)]


                 minimumTOF = TOF(minimizingPsi)
                 
                 #if there are no solutions for this TOF, exit the function
                 #and return values so the thing calling this doesn't break
                 if minimumTOF > dt_0:
                     return \
                         {
                             'v_0' : nan,
                             'v_f': nan,
                             'DM': nan,
                             'vInfDepart': array([nan,nan,nan]),
                             'vInfArrive': array([nan,nan,nan]),
                             'magVInfArrive': nan,
                             'C3': nan,
                             'delta_nu': nan,
                             'i': nan,
                             'psi': nan,
                             'minTOF': minimumTOF,
                             'minimizingPsi': minimizingPsi
                         }


                 if transferType == 3:
                     psi_up = minimizingPsi
                 else:
                     psi_low = minimizingPsi

                 psi = (psi_low + psi_up)/2
         
         #initialize dt to make sure we always walk through the loop
            */
        
        //dt = dt_0+0.0001
        var dt = dt0 + 0.0001
        var iterations = 0
        
        /*
        if transferType == 3:
            def dtCheck(dt,dt_0):
                return dt >= dt_0
        else:
            def dtCheck(dt,dt_0):
                return dt <= dt_0
         */
        var dtCheck: (Double, Double) -> Bool
        if transferType == 3 {
            dtCheck = { (dt, dt_0) in
                return dt >= dt_0
            }
        } else {
            dtCheck = { (dt, dt_0) in
                return dt <= dt_0
            }
        }
        
            /*
         while abs(dt - dt_0) > 1e-6 and i<1000:
         i += 1

         if A > 0 and y < 0:
             while y < 0:
                 psi = psi + 0.1
                 y = r0 + rf + A*(psi*c3 - 1)/sqrtC2

         sqrtY = sqrt(y)
         xi = sqrtY/sqrtC2
         dt = (xi**3*c3+A*sqrtY)/sqrtMu


         if dtCheck(dt,dt_0):
             psi_low = psi
         else:
             psi_up = psi

         psi = (psi_up + psi_low)/2.
         */
        
        var y: Double = 0
        
        while abs(dt - dt0) > 1e-6 && iterations < 1000 {
            iterations += 1
            
            //sqrtC2 = sqrt(c2)
            let sqrtC2 = sqrt(c2)
            
            //y = r0 + rf + A*(psi*c3 - 1)/sqrtC2
            y = r0 + rf + A * (psi * c3 - 1) / sqrtC2
            
            if A > 0 && y < 0 {
                while y < 0 {
                    psi = psi + 0.1
                    y = r0 + rf + A * (psi * c3 - 1) / sqrtC2
                }
            }
            
            let sqrtY = sqrt(y)
            let xi = sqrtY / sqrtC2
            dt = (pow(xi, 3) * c3 + A * sqrtY) / sqrtMu
                    
            if dtCheck(dt, dt0) {
                psi_low = psi
            } else {
                psi_up = psi
            }
            
            psi = (psi_up + psi_low) / 2
            
             /*
             if psi > 1e-6:
                 sqrtPsi = sqrt(psi)
                 c2 = (1-cos(sqrtPsi))/psi
                 c3 = (sqrtPsi-sin(sqrtPsi))/sqrtPsi**3

             elif psi < -1e-6:
                 sqrtNegPsi = sqrt(-psi)
                 c2 = (1-cosh(sqrtNegPsi))/psi
                 c3 = (sinh(sqrtNegPsi)-sqrtNegPsi)/sqrtNegPsi**3
             
             else:
                 c2 = 1./2.
                 c3 = 1./6.

             */
            if psi > 1e-6 {
                let sqrtPsi = sqrt(psi)
                c2 = (1 - cos(sqrtPsi)) / psi
                c3 = (sqrtPsi - sin(sqrtPsi)) / pow(sqrtPsi, 3)
            } else if psi < -1e-6 {
                let sqrtPsi = sqrt(-psi)
                c2 = (1 - cosh(sqrtPsi))/psi
                c3 = (sinh(sqrtPsi) - sqrtPsi) / pow(sqrtPsi, 3)
            } else {
                c2 = 1/2
                c3 = 1/6
            }
            
        }
         /*

             if i == 1000: print("max iterations reached")

             f = 1 - y/r0
             g_dot = 1 - y/rf
             g = A*sqrtY/sqrtMu

             v_0 = (rf*nf - f*r0*n0)/g
             v_f = (g_dot*rf*nf - r0*n0)/g

             try:
                 minTOF = TOF(minimizingPsi)
             except:
                 minTOF = nan
                 minimizingPsi = nan
         
         */
        
        if iterations == 1000 { print("Warning: Maximum iterations reached") }
        let f = 1 - y / r0
        let g_dot = 1 - y / rf
        let g = A * sqrt(y) / sqrtMu
        
        let v_0 = (rf * nf - f * r0 * n0) / g
        let v_f = (g_dot * rf * nf - r0 * n0) / g
        
        return [(v_0, v_f, 0)]
        
        /*

             return \
                 {
                     'v_0' : v_0,
                     'v_f': v_f,
                     'DM': DM,
                     'vInfDepart': v_0 - X_0[3:6],
                     'vInfArrive': v_f - X_f[3:6],
                     'magVInfArrive': sqrt((v_f - X_f[3:6]).dot(v_f - X_f[3:6])),
                     'C3': (v_0 - X_0[3:6]).dot(v_0 - X_0[3:6]),
                     'delta_nu': delta_nu,
                     'i': i,
                     'psi': psi,
                     'minTOF': minTOF,
                     'minimizingPsi': minimizingPsi
                 }
         */
    }
    
}

func convergingSolver(_ range: ClosedRange<Double>, tolerance: Double, _ function: (Double)->Double) -> Double {
    if range.upperBound - range.lowerBound < tolerance { return range.upperBound + range.lowerBound / 2 }
    
    let distance = range.upperBound - range.lowerBound
    let lower = range.lowerBound...range.lowerBound + distance / 2
    let upper = lower.upperBound...range.upperBound
    
    if function(lower.upperBound) < 0 {
        return convergingSolver(upper, tolerance: tolerance, function)
    } else {
        return convergingSolver(lower, tolerance: tolerance, function)
    }
    
}
