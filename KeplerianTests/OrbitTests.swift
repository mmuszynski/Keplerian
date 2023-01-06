//
//  OrbitTests.swift
//  KeplerianTests
//
//  Created by Mike Muszynski on 5/29/20.
//  Copyright © 2020 Mike Muszynski. All rights reserved.
//

import XCTest
@testable import Keplerian

class OrbitTests: XCTestCase {
    
    func testInitializeOrbitByMeasurement() {
        let sol = CelestialBody(gravitationalParameter: 1.32712440042 * 1e20, radius: 696340, mass: 0)

        let mercuryOrbit = Orbit(semiMajorAxis: (57.91 * 10e6).km,
                                 eccentricity: 0.2056630,
                                 meanAnomaly: 174.796.degrees,
                                 inclination: 7.005.degrees,
                                 LAN: 48.331.degrees,
                                 argumentOfPeriapsis: 29.124.degrees,
                                 centralBody: sol)
        
        let equalOrbit = Orbit(semiMajorAxis: 57.91 * 10e9,
                                 eccentricity: 0.2056630,
                                 meanAnomaly: 3.050765719,
                                 inclination: 0.12226031,
                                 LAN: 48.331 * .pi / 180,
                                 argumentOfPeriapsis: 0.508309691,
                                 centralBody: sol)
        
        XCTAssertEqual(mercuryOrbit.semiMajorAxis, equalOrbit.semiMajorAxis, accuracy: 1e-6)
        XCTAssertEqual(mercuryOrbit.eccentricity, equalOrbit.eccentricity, accuracy: 1e-6)
        XCTAssertEqual(mercuryOrbit.meanAnomaly, equalOrbit.meanAnomaly, accuracy: 1e-6)
        XCTAssertEqual(mercuryOrbit.inclination, equalOrbit.inclination, accuracy: 1e-6)
        XCTAssertEqual(mercuryOrbit.LAN, equalOrbit.LAN, accuracy: 1e-6)
        XCTAssertEqual(mercuryOrbit.argumentOfPeriapsis, equalOrbit.argumentOfPeriapsis, accuracy: 1e-6)
    }
    
    func testTrivialPlanetPositions() {
        let sun = CelestialBody(gravitationalParameter: 1000, radius: 1000, mass: 1000)
        let orbit = Orbit(semiMajorAxis: 2000, eccentricity: 0, meanAnomaly: 0, inclination: 0, LAN: 0, argumentOfPeriapsis: 0, centralBody: sun)
        
        let initialPosition = Vector3D(x: 2000, y: 0, z: 0)
        XCTAssertEqual(initialPosition, orbit.cartesian(atTime: 0).position)
        
        let eccentricOrbit = Orbit(semiMajorAxis: 2000, eccentricity: 0.5, meanAnomaly: 0, inclination: 0, LAN: 0, argumentOfPeriapsis: 0, centralBody: sun)
        XCTAssertEqual(initialPosition, orbit.cartesian(atTime: 0).position)
        
        let expectedPosition = Vector3D(x: -3000, y: 0, z: 0)
        let diff = expectedPosition - eccentricOrbit.cartesian(atTime: 8885.765876316733).position
        XCTAssertEqual(diff.magnitude, 0, accuracy: 1e-6)
    }
    
    func testKerbinOrbitFromMattLibrary() {
        let orbit = Orbit.kerbin
        
        print("testing time since periapsis: \(orbit.timeSincePeriapsis(atTime: 0))")
        XCTAssertEqual(orbit.trueAnomaly(atTime: 0), 179.90874767107852.deg2rad, accuracy: 1e-6)
        XCTAssertEqual(orbit.meanAnomaly(atTime: 0), 179.90874767107852.deg2rad, accuracy: 1e-6)
        XCTAssertEqual(orbit.eccentricAnomaly(atTime: 0), 179.90874767107852.deg2rad, accuracy: 1e-6)
        
        let atPeri = -4599439.396167472
        print("testing time since periapsis: \(orbit.timeSincePeriapsis(atTime: atPeri))")
        XCTAssertEqual(orbit.trueAnomaly(atTime: atPeri), 0, accuracy: 1e-6)
        XCTAssertEqual(orbit.meanAnomaly(atTime: atPeri), 0, accuracy: 1e-6)
        XCTAssertEqual(orbit.eccentricAnomaly(atTime: atPeri), 0, accuracy: 1e-6)
        
        let atRandom = 6123648.12356
        print("testing time since periapsis: \(orbit.timeSincePeriapsis(atTime: atRandom))")
        XCTAssertEqual(orbit.trueAnomaly(atTime: atRandom), 59.437475021205209.deg2rad, accuracy: 1e-6)
        XCTAssertEqual(orbit.meanAnomaly(atTime: atRandom), 59.437475021205209.deg2rad, accuracy: 1e-6)
        XCTAssertEqual(orbit.eccentricAnomaly(atTime: atRandom), 59.437475021205209.deg2rad, accuracy: 1e-6)
        
        let kspTime = KSPDate(year: 1, day: 2, hour: 0, minute: 0, second: 0)
        print("testing time since periapsis: \(orbit.timeSincePeriapsis(atTime: atRandom))")
        print(orbit.anom(atTime: kspTime.timeIntervalSinceReferenceDate))
        let E_m_nu_t = [180.7536395404949, 180.7536395404949, 180.7536395404949, 4621039.396167472]
        XCTAssertEqual(orbit.eccentricAnomaly(atTime: kspTime.timeIntervalSinceReferenceDate).rad2deg, E_m_nu_t[0], accuracy: 1e-6)
        XCTAssertEqual(orbit.meanAnomaly(atTime: kspTime.timeIntervalSinceReferenceDate).rad2deg, E_m_nu_t[1], accuracy: 1e-6)
        XCTAssertEqual(orbit.trueAnomaly(atTime: kspTime.timeIntervalSinceReferenceDate).rad2deg, E_m_nu_t[2], accuracy: 1e-6)
        XCTAssertEqual(orbit.timeSincePeriapsis(atTime: kspTime.timeIntervalSinceReferenceDate), E_m_nu_t[3], accuracy: 1e-6)
    }
    
    
    func testOrbitEquality() {
        XCTAssertEqual(Orbit.kerbin, Orbit.kerbin)
        XCTAssertNotEqual(Orbit.kerbin, Orbit.duna)
    }
    
    func testBodyEquality() {
        XCTAssertEqual(CelestialBody.kerbin, .kerbin)
        XCTAssertNotEqual(CelestialBody.kerbin, .duna)
        XCTAssertNotEqual(CelestialBody.kerbol, .kerbin)
    }
    
    func testHyperbolicOrbit() {
        let orbit = Orbit(semiMajorAxis: (-700000.0).meters, eccentricity: 1.25, meanAnomaly: 0.0.rad, inclination: 0.0.rad, LAN: 0.0.rad, argumentOfPeriapsis: 0.0.rad, centralBody: .kerbin)
        
        XCTAssertGreaterThan(orbit.eccentricity, 1)
        XCTAssertLessThan(orbit.semiMajorAxis, 0)
        XCTAssertEqual(orbit.periapsis, orbit.semiMajorAxis * (1 - orbit.eccentricity))
        XCTAssertEqual(orbit.parameter, orbit.semiMajorAxis * (1 - orbit.eccentricity * orbit.eccentricity))
        XCTAssertEqual(orbit.totalEnergy, -orbit.centralBody.gravitationalParameter / (2 * orbit.semiMajorAxis))
        
        //http://control.asu.edu/Classes/MAE462/462Lecture05.pdf
        //Given time, t, solve for Hyperbolic Mean Anomaly
        let M = orbit.hyperbolicMeanAnomaly(atTime: 0)
        let H = orbit.hyperbolicAnomaly(fromMeanAnomaly: M)
        
        
        XCTFail("unimplemented")
        print(orbit.radius())
        print(orbit.periapsis)
    }
    
    func testKSPOrbit() {
        let orbitURL = Bundle.module.url(forResource: "realKSPOrbit", withExtension: "json")!
        
        do {
            let orbitData = try Data(contentsOf: orbitURL)
            let decoder = JSONDecoder()
            let orbit = try decoder.decode(Orbit.self, from: orbitData)
            XCTAssertEqual(orbit.altitude(), 80.0)
        } catch {
            XCTFail("Unhandled error \(error)")
        }
    }
    
    func testKSPOrbit2() {
        let orbitURL = Bundle.module.url(forResource: "129975m-orbit", withExtension: "json")!
        
        do {
            let orbitData = try Data(contentsOf: orbitURL)
            let decoder = JSONDecoder()
            let orbit = try decoder.decode(Orbit.self, from: orbitData)
            XCTAssertEqual(orbit.altitude(), 129975)
        } catch {
            XCTFail("Unhandled error \(error)")
        }
    }
    
    func testOrbitLoad() {
        let _ = Orbit.exampleEscapeOrbit
    }
        
}

private extension Orbit {
    func anom(atTime time: Double) -> String {
        return "anomalies('t', \(self.timeSincePeriapsis(atTime: time)), \(self.semiMajorAxis), \(self.eccentricity), mu=\(self.centralBody.gravitationalParameter))"
    }
    
    func coe(atTime time: Double) -> String {
        return "\(semiMajorAxis), \(eccentricity), \(inclination), \(LAN), \(argumentOfPeriapsis), \(trueAnomaly(atTime: time) * 180 / .pi), mu=\(centralBody.gravitationalParameter)"
    }
}

private extension Double {
    var deg2rad: Double {
        return self * .pi / 180
    }
    var rad2deg: Double {
        return self * 180 / .pi
    }
}
