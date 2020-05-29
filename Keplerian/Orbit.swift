//
//  Orbit.swift
//  KSPCockpitPanel
//
//  Created by Mike Muszynski on 7/22/17.
//  Copyright © 2017 Mike Muszynski. All rights reserved.
//

import Foundation
import simd

public class Orbit {
    //The Keplerian Orbital elements
    public var semiMajorAxis: Double
    public var eccentricity: Double
    public var meanAnomaly: Double
    public var inclination: Double
    public var LAN: Double
    public var argumentOfPeriapsis: Double
    
    //Epoch for KSP has been zero for all orbits
    //Generally, this probably needs to be specified
    public var epoch: Double = 0
    func timeFromEpoch(for time: Double) -> Double {
        return time - epoch
    }
    
    public var centralBody: CelestialBody
    
    public var isSubOrbital: Bool {
        return self.semiMajorAxis < self.centralBody.radius
    }
    
    public var isEscapeOrbit: Bool {
        return self.eccentricity >= 1.0
    }
    
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
    
    public convenience init(semiMajorAxis: Measurement<UnitLength>, eccentricity: Double, meanAnomaly: Measurement<UnitAngle>, inclination: Measurement<UnitAngle>, LAN: Measurement<UnitAngle>, argumentOfPeriapsis: Measurement<UnitAngle>, centralBody: CelestialBody) {
        
        
        self.init(semiMajorAxis: semiMajorAxis.converted(to: .meters).value,
                  eccentricity: eccentricity,
                  meanAnomaly: meanAnomaly.converted(to: .radians).value,
                  inclination: inclination.converted(to: .radians).value,
                  LAN: LAN.converted(to: .degrees).value,
                  argumentOfPeriapsis: argumentOfPeriapsis.converted(to: .radians).value,
                  centralBody: centralBody)
    }
    
    /// The maximum distance from the center of mass of the body that is being orbited
    public var apoapsis: Double {
        return (1 + eccentricity) * semiMajorAxis
    }
    
    /// The minimum distance from the center of mass of the body that is being orbited
    public var periapsis: Double {
        return (1 - eccentricity) * semiMajorAxis
    }
    
    public var periapsisAltitude: Double {
        return periapsis - centralBody.radius
    }
    
    public var apoapsisAltitude: Double {
        return apoapsis - centralBody.radius
    }
    
    public func altitude(atTime time: Double = 0) -> Double {
        return altitude(atTimeFromEpoch: timeFromEpoch(for: time))
    }
    
    private func altitude(atTimeFromEpoch time: Double = 0) -> Double {
        return radius(atTimeFromEpoch: time) - centralBody.radius
    }
    
    public var semiMinorAxis: Double {
        return semiMajorAxis * sqrt(1.0 - eccentricity * eccentricity)
    }
    
    public var meanMotion: Double {
        let gm = centralBody.gravitationalParameter
        return sqrt(gm / pow(semiMajorAxis, 3))
    }
    
    public var period: Double {
        let gm = centralBody.gravitationalParameter
        return 2.0 * Double.pi * sqrt(pow(semiMajorAxis, 3) / gm)
    }
    
    public func meanAnomaly(atTime time: Double = 0) -> Double {
        return meanAnomaly(atTimeFromEpoch: timeFromEpoch(for: time))
    }
    
    private func meanAnomaly(atTimeFromEpoch time: Double = 0) -> Double {
        var meanAnomalyAtTime = self.meanMotion * time + meanAnomaly;
        while meanAnomalyAtTime > (2.0 * Double.pi) {
            meanAnomalyAtTime -= 2.0 * Double.pi
        }
        return meanAnomalyAtTime
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
        let E = self.eccentricAnomaly(atTimeFromEpoch: time)
        let ecc = eccentricity

        let numerator = cos(E) - ecc
        let denominator = 1 - ecc * cos(E)
        
        return acos(numerator/denominator)
    }
    
    public func radius(atTime time: Double = 0) -> Double {
        return radius(atTimeFromEpoch: timeFromEpoch(for: time))
    }
    
    private func radius(atTimeFromEpoch time: Double = 0) -> Double {
        let trueAnomaly = self.trueAnomaly(atTimeFromEpoch: time)
        return semiMajorAxis * (1 - eccentricity * eccentricity) / (1 + eccentricity * cos(trueAnomaly))
    }
    
    public func orbitalSpeed(atTime time: Double = 0) -> Double {
        return orbitalSpeed(atTimeFromEpoch: timeFromEpoch(for: time))
    }
    
    private func orbitalSpeed(atTimeFromEpoch time: Double = 0) -> Double {
        let gm = centralBody.gravitationalParameter
        let radius = self.radius(atTimeFromEpoch: time)
        return sqrt(gm * (2.0 / radius - 1 / self.semiMajorAxis))
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
        
        let meanAnomalyFractionRemaining = (2.0 * Double.pi - effectiveMeanAnomaly) / (2.0 * Double.pi)
        let period = self.period
        
        return period * meanAnomalyFractionRemaining
    }
    
    public var specificAngularMomentum: Double {
        let gm = centralBody.gravitationalParameter
        let a = semiMajorAxis
        let e = eccentricity
        
        return sqrt(gm * a * (1 - e * e))
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
        let x = r * (cos(rad(from: LAN)) * cos(nu + argumentOfPeriapsis) - sin(rad(from: LAN)) * sin(nu + argumentOfPeriapsis) * cos(inclination))
        let y = r * (sin(rad(from: LAN)) * cos(nu + argumentOfPeriapsis) + cos(rad(from: LAN)) * sin(nu + argumentOfPeriapsis) * cos(inclination))
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
        let OMEGA = rad(from: LAN)
        let omega = argumentOfPeriapsis
        let nu = trueAnomaly(atTimeFromEpoch: time)
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
        
        let pqw2ijk = rz(-OMEGA) * rx(-i) * rz(-omega)
        let r_ijk = pqw2ijk * r_pqw_simd
        let v_ijk = pqw2ijk * v_pqw_simd
                
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
