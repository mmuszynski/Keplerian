//
//  Orbit.swift
//  KSPCockpitPanel
//
//  Created by Mike Muszynski on 7/22/17.
//  Copyright © 2017 Mike Muszynski. All rights reserved.
//

import Foundation
import simd

//http://www.bogan.ca/orbits/kepler/orbteqtn.html

fileprivate extension Double {
    static let twoPi = 2 * Double.pi
}

public class Orbit: Codable {
    
    //-MARK: Classcal Orbital Elements
    //The Keplerian Orbital elements
    public var semiMajorAxis: Double
    public var eccentricity: Double
    public var meanAnomaly: Double
    public var inclination: Double
    public var LAN: Double
    public var argumentOfPeriapsis: Double
    
    //-MARK: Epoch
    //Epoch for KSP has been zero for all orbits
    //Generally, this probably needs to be specified
    public var epoch: Double = 0
    func timeFromEpoch(for time: Double) -> Double {
        return time - epoch
    }
    
    public var centralBody: CelestialBody
    
    //-MARK: Orbit Descriptions
    /// Whether or not the orbit is going to clear the body's radius
    public var isSubOrbital: Bool {
        guard isHyperbolic else { return false }
        return self.semiMajorAxis < self.centralBody.radius
    }
    
    /// If the orbit is hyperbolic
    public var isHyperbolic: Bool {
        return self.eccentricity >= 1.0
    }
    
    //-MARK: Initialization
    
    /// Initializes an `Orbit` with the prescribed Keplerian orbital elements
    ///
    /// - Parameters:
    ///   - semiMajorAxis: The semi-major axis of the orbit
    ///   - eccentricity: The eccentricity of the orbit
    ///   - meanAnomaly: The mean anomaly of the orbot
    ///   - inclination: The inclination of the orbit
    ///   - LAN: The longitude of ascending node of the orbit
    ///   - argumentOfPeriapsis: The argument of periapsis of the orbit
    public init(semiMajorAxis: Double, eccentricity: Double, meanAnomaly: Double, inclination: Double, LAN: Double, argumentOfPeriapsis: Double, centralBody: CelestialBody) {
        self.semiMajorAxis = semiMajorAxis
        self.eccentricity = eccentricity
        self.meanAnomaly = meanAnomaly <= 0 ? meanAnomaly + 2.0 * Double.pi : meanAnomaly
        self.inclination = inclination
        self.LAN = LAN
        self.argumentOfPeriapsis = argumentOfPeriapsis
        
        self.centralBody = centralBody
    }
    
    /// Initializes an `Orbit` with the prescribed Keplerian orbital elements using measurements
    ///
    /// - Parameters:
    ///   - semiMajorAxis: The semi-major axis of the orbit expressed as a measurement
    ///   - eccentricity: The eccentricity of the orbit expressed as a measurement
    ///   - meanAnomaly: The mean anomaly of the orbot expressed as a measurement
    ///   - inclination: The inclination of the orbit expressed as a measurement
    ///   - LAN: The longitude of ascending node of the orbit expressed as a measurement
    ///   - argumentOfPeriapsis: The argument of periapsis of the orbit
    public convenience init(semiMajorAxis: Measurement<UnitLength>, eccentricity: Double, meanAnomaly: Measurement<UnitAngle>, inclination: Measurement<UnitAngle>, LAN: Measurement<UnitAngle>, argumentOfPeriapsis: Measurement<UnitAngle>, centralBody: CelestialBody) {
        self.init(semiMajorAxis: semiMajorAxis.converted(to: .meters).value,
                  eccentricity: eccentricity,
                  meanAnomaly: meanAnomaly.converted(to: .radians).value,
                  inclination: inclination.converted(to: .radians).value,
                  LAN: LAN.converted(to: .radians).value,
                  argumentOfPeriapsis: argumentOfPeriapsis.converted(to: .radians).value,
                  centralBody: centralBody)
    }
    
    //-MARK: Alias functions
    
    /// An alias for the eccentricity of the orbit
    private var e: Double {
        return eccentricity
    }
    
    /// An alias for the semi-major axis
    private var a: Double {
        return semiMajorAxis
    }
    
    /// The gravitational parameter of the central body
    private var mu: Double {
        return centralBody.gravitationalParameter
    }
    
    /// An alias for the total energy of the orbit
    private var E: Double {
        return totalEnergy
    }
    
    private var p: Double {
        return parameter
    }
    
    //-MARK: Calculated parameters
    
    /// The parameter, aka semi-latus rectum
    var parameter: Double {
        return a * (1 - e * e)
    }
    
    /// The total energy of the orbit
    var totalEnergy: Double {
        return -mu/(2 * a)
    }
    
    /// The maximum distance from the center of mass of the body that is being orbited
    public var apoapsis: Double {
        return a * (1 + e)
    }
    
    /// The minimum distance from the center of mass of the body that is being orbited
    public var periapsis: Double {
        return a * (1 - e)
    }
    
    /// The minimum altitude from the average surface of the body being orbited
    public var periapsisAltitude: Double {
        return periapsis - centralBody.radius
    }
    
    /// The maximum altitude from the average surface of the body being orbited
    public var apoapsisAltitude: Double {
        return apoapsis - centralBody.radius
    }
    
    /// The altitude at a given time for a given orbit
    public func altitude(atTime time: Double = 0) -> Double {
        return altitude(atTimeFromEpoch: timeFromEpoch(for: time))
    }
    
    /// The altitude at a given time since the epoch of the orbit
    private func altitude(atTimeFromEpoch time: Double = 0) -> Double {
        return radius(atTimeFromEpoch: time) - centralBody.radius
    }

    /// The semi-minor axis of the orbit
    public var semiMinorAxis: Double {
        return semiMajorAxis * sqrt(1.0 - eccentricity * eccentricity)
    }
    
    /// The mean motion of the orbig
    public var meanMotion: Double {
        if semiMajorAxis < 0 {
            return sqrt(mu / pow(-semiMajorAxis, 3))
        } else {
            return sqrt(mu / pow(semiMajorAxis, 3))
        }
    }
    
    public var period: Double {
        let mu = centralBody.gravitationalParameter
        return 2.0 * Double.pi * sqrt(pow(semiMajorAxis, 3) / mu)
    }
    
    public func meanAnomaly(atTime time: Double = 0) -> Double {
        return meanAnomaly(atTimeFromEpoch: timeFromEpoch(for: time))
    }
    
    private func meanAnomaly(atTimeFromEpoch time: Double = 0, in range: Range<Double> = 0..<Double.twoPi) -> Double {
        var meanAnomalyAtTime = self.meanMotion * time + meanAnomaly;
        while meanAnomalyAtTime >= range.upperBound {
            meanAnomalyAtTime -= 2.0 * .pi
        }
        
        while meanAnomalyAtTime < range.lowerBound {
            meanAnomalyAtTime += 2.0 * .pi
        }
        
        return meanAnomalyAtTime
    }
        
    public func hyperbolicAnomaly(fromMeanAnomaly M: Double) -> Double {
        let tolerance = 0.0005
        
        //need to solve M = e * sinh(F) - F
        // start with an initial value of M
        var F: Double = M
        var delta = Double.greatestFiniteMagnitude
        
        while abs(delta) > tolerance {
            let numerator = M - e * sinh(F) + F
            let denominator = e * cosh(F) - 1
            let f1 = F + numerator / denominator
            
            delta = f1 - F
            F = f1
        }
        
        return F
    }
    
    
    public func eccentricAnomaly(fromMeanAnomaly meanAnomaly: Double) -> Double {
        let tolerance = 0.0005
        let eccentricity = self.eccentricity
        
        //need to solve M=E-e\sin E for E
        
        let M = meanAnomaly
        
        var E: Double
        if (eccentricity < 0.8 ) {
            E = meanAnomaly
        } else {
            E = Double.pi
        }
        
        var delta = Double.greatestFiniteMagnitude
        
        while abs(delta) > tolerance {
            delta = (M - E + eccentricity * sin(E)) / (1 - eccentricity * cos(E))
            E = E + delta
        }
        
        return E
    }
    
    public func eccentricAnomaly(atTime time: Double = 0) -> Double {
        return eccentricAnomaly(atTimeFromEpoch: timeFromEpoch(for: time))
    }
    
    private func eccentricAnomaly(atTimeFromEpoch time: Double = 0) -> Double {
        let meanAnomaly = self.meanAnomaly(atTimeFromEpoch: time)
        let eccentricAnomaly = self.eccentricAnomaly(fromMeanAnomaly: meanAnomaly)
        
        return eccentricAnomaly
    }
    
    public func trueAnomaly(atTime time: Double = 0) -> Double {
        trueAnomaly(atTimeFromEpoch: timeFromEpoch(for: time))
    }
    
    private func trueAnomaly(atTimeFromEpoch time: Double = 0) -> Double {
        if self.isHyperbolic {
            let meanAnomaly = self.meanAnomaly(atTimeFromEpoch: time, in: -Double.greatestFiniteMagnitude..<Double.greatestFiniteMagnitude)
            let F = self.hyperbolicAnomaly(fromMeanAnomaly: meanAnomaly)
            return 2 * atan(sqrt((e + 1) / (e - 1)) * tanh(F / 2))
        }
        
        let E = self.eccentricAnomaly(atTimeFromEpoch: time)
        
        /*
        // This is nice in an ideal case, but more research has shown that it has a big problem as e -> 1
        // The denominator becomes zero, and that will trap
        // It's a pretty basic method, and there is a better (?) method implemented below
        let numerator = cos(E) - eccentricity
        let denominator = 1 - eccentricity * cos(E)

        return acos(numerator/denominator)
        */
        
        //https://en.wikipedia.org/wiki/True_anomaly
        let beta = e / (1 + sqrt(1 - e * e))
        return E + 2 * atan(beta * sin(E) / (1 - beta * cos(E)))
         
        /*
         Instead, the code from Matt's library (and Wikipedia):
         
         sin_nu = (np.sin(E)*np.sqrt(1.-e**2.))/(1.-e*np.cos(E))
         cos_nu = (np.cos(E)-e)/(1.-e*np.cos(E))
         nu = np.arctan2(sin_nu,cos_nu)
        
        let sin_nu = (sin(E) * sqrt(1 - e * e)) / (1 - e * cos(E))
        let cos_nu = (cos(E) - e) / (1 - e * cos(E))
        var nu = atan2(sin_nu, cos_nu)
        
        while nu < 0 {
            nu += 2 * .pi
        }
        
        while nu >= .pi * 2 {
            nu -= 2 * .pi
        }
        
        return nu
        */
    }
    
    public func radius(atTime time: Double = 0) -> Double {
        return radius(atTimeFromEpoch: timeFromEpoch(for: time))
    }
    
    private func radius(atTimeFromEpoch time: Double = 0) -> Double {
        let trueAnomaly = self.trueAnomaly(atTimeFromEpoch: time)
        
        if self.isHyperbolic {
            //return a * (1 - e * cosh(trueAnomaly))
            let num = a * (1 - e * e)
            let den = 1 + e * cos(trueAnomaly)
            return num / den
        }
        
        return p / (1 + e * cos(trueAnomaly))
    }
    
    public func orbitalSpeed(atTime time: Double = 0) -> Double {
        return orbitalSpeed(atTimeFromEpoch: timeFromEpoch(for: time))
    }
    
    private func orbitalSpeed(atTimeFromEpoch time: Double = 0) -> Double {
        let mu = centralBody.gravitationalParameter
        let radius = self.radius(atTimeFromEpoch: time)
        return sqrt(mu * (2.0 / radius - 1 / self.semiMajorAxis))
    }
    
    public func timeToApoapsis(atTime time: Double = 0) -> Double {
        return timeToApoapsis(atTimeFromEpoch: timeFromEpoch(for: time))
    }
    
    private func timeToApoapsis(atTimeFromEpoch time: Double = 0) -> Double {
        let effectiveMeanAnomaly = self.meanAnomaly(atTimeFromEpoch: time)
        var target = Double.pi
        if effectiveMeanAnomaly > target {
            target += 2.0 * Double.pi
        }
        
        let meanAnomalyFractionRemaining = (target - effectiveMeanAnomaly) / (2.0 * Double.pi)
        let period = self.period
        
        return period * meanAnomalyFractionRemaining
    }
    
    public func timeToPeriapsis(atTime time: Double = 0) -> Double {
        return timeToPeriapsis(atTimeFromEpoch: timeFromEpoch(for: time))
    }
    
    private func timeToPeriapsis(atTimeFromEpoch time: Double = 0) -> Double {
        let effectiveMeanAnomaly = self.meanAnomaly(atTimeFromEpoch: time)
        
        if self.isHyperbolic {
            let meanAnomalyAtPeriapsis = 0.0
            let time = (effectiveMeanAnomaly - meanAnomalyAtPeriapsis) / meanMotion
            return time
        }
        
        let meanAnomalyFractionRemaining = (2.0 * Double.pi - effectiveMeanAnomaly) / (2.0 * Double.pi)
        let period = self.period
        
        return period * meanAnomalyFractionRemaining
    }
    
    public func timeSincePeriapsis(atTime time: Double) -> Double {
        let ttpe = self.timeToPeriapsis(atTime: time)
        return self.period - ttpe
    }
    
    public var specificAngularMomentum: Double {
        let mu = centralBody.gravitationalParameter
        let a = semiMajorAxis
        let e = eccentricity
        
        return sqrt(mu * a * (1 - e * e))
    }
    
    public func cartesianPosition(atTime time: Double = 0) -> Vector3D {
        return cartesianPosition(atTimeFromEpoch: timeFromEpoch(for: time))
    }
    
    private func cartesianPosition(atTimeFromEpoch time: Double = 0) -> Vector3D {
        func rad(from deg: Double) -> Double {
            return Double.pi * deg / 180
        }
        
        let r = radius(atTimeFromEpoch: time)
        let nu = trueAnomaly(atTimeFromEpoch: time)
        let x = r * (cos(LAN) * cos(nu + argumentOfPeriapsis) - sin(LAN) * sin(nu + argumentOfPeriapsis) * cos(inclination))
        let y = r * (sin(LAN) * cos(nu + argumentOfPeriapsis) + cos(LAN) * sin(nu + argumentOfPeriapsis) * cos(inclination))
        let z = r * sin(inclination) * sin(nu + argumentOfPeriapsis)
        return Vector3D(x: x, y: y, z: z)
    }
    
    public func cartesian(atTime time: Double = 0) -> (position: Vector3D, velocity: Vector3D) {
        return cartesian(atTimeFromEpoch: timeFromEpoch(for: time))
    }
    
    private func cartesian(atTimeFromEpoch time: Double = 0) -> (position: Vector3D, velocity: Vector3D) {
        func rad(from deg: Double) -> Double {
            return Double.pi * deg / 180
        }
        
        //def coe2rv(a, e, i, OMEGA, omega, nu, **kwargs):
        //#accept kwarg input for mu. If none is used, use that of Earth
        //#with a satellite of negligably small mass in km^3/s^2
        //try:
        //    mu = kwargs['mu']
        //except:
        //    mu = 398600.4415

        let i = inclination
        let OMEGA = LAN
        let omega = argumentOfPeriapsis
        let nu = trueAnomaly(atTime: timeFromEpoch(for: time))
        let e = eccentricity
        let mu = self.centralBody.gravitationalParameter
        
        //#define semipameter
        //p = a*(1-e**2)
        let p = semiMajorAxis * ( 1 - e * e )

        //#define r in the perifocal plane
        
        //r_pqw = vstack([
        //    ((p*cos(nu))/(1+e*cos(nu))),
        //    ((p*sin(nu))/(1+e*cos(nu))),
        //    zeros(length)
        //    ])
        let r_pqw = Vector3D(x: p * cos(nu) / (1 + e * cos(nu)),
                             y: p * sin(nu) / (1 + e * cos(nu)),
                             z: 0)
        let r_pqw_simd = r_pqw.simd

        //#define v in the perifocal plane
        //v_pqw = vstack([
        //    (-sqrt(mu/p)*sin(nu)),
        //    (sqrt(mu/p)*(e+cos(nu))),
        //    zeros(length)
        //   ])
        let v_pqw = Vector3D(x: -sqrt(mu/p)*sin(nu),
                             y: sqrt(mu/p)*(e+cos(nu)),
                             z: 0)
        let v_pqw_simd = v_pqw.simd
        
        func rz(_ angle: Double) -> simd_double3x3 {
            simd_double3x3([
                simd_double3( cos(angle), sin(angle), 0),
                simd_double3(-sin(angle), cos(angle), 0),
                simd_double3( 0,            0,          1)
            ])
        }
        
        func rx(_ angle: Double) -> simd_double3x3 {
            simd_double3x3([
                simd_double3(1, 0, 0),
                simd_double3(0, cos(angle), sin(angle)),
                simd_double3(0, -sin(angle), cos(angle))
            ])
        }
        
        let rot1 = rz(-OMEGA)
        let rot2 = rx(-i)
        let rot3 = rz(-omega)
        
        let pqw2ijk = rot3 * rot2 * rot1
        let r_ijk = r_pqw_simd * pqw2ijk
        let v_ijk = v_pqw_simd * pqw2ijk 
                
        //pqw2ijk = rz(-OMEGA).dot(rx(-i).dot(rz(-omega)))
        //r_ijk = pqw2ijk.dot(r_pqw)
        //v_ijk = pqw2ijk.dot(v_pqw)
        //X_ijk = vstack([r_ijk,v_ijk]).T[0]
        
        return (position: Vector3D(from: r_ijk), velocity: Vector3D(from: v_ijk))
    }
    
}

extension Orbit: Equatable {
    public static func == (lhs: Orbit, rhs: Orbit) -> Bool {
        if lhs.semiMajorAxis != rhs.semiMajorAxis { return false }
        if lhs.eccentricity != rhs.eccentricity { return false }
        if lhs.meanAnomaly != rhs.meanAnomaly { return false }
        if lhs.inclination != rhs.inclination { return false }
        if lhs.LAN != rhs.LAN { return false }
        if lhs.argumentOfPeriapsis != rhs.argumentOfPeriapsis { return false }
        if lhs.epoch != rhs.epoch { return false }
        
        if lhs.centralBody != rhs.centralBody { return false }
        
        return true
    }
}
