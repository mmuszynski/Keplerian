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
        let sol = CelestialBody(gravitationalParameter: 1.32712440042 * 1e20, radius: 696340)

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
                                 LAN: 48.331,
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
        let sun = CelestialBody(gravitationalParameter: 1000, radius: 1000)
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
        
        print("testing time since periapsis: \(orbit.timeSincePeriapsis(forTime: 0))")
        XCTAssertEqual(orbit.trueAnomaly(atTime: 0), 179.90874767107852.deg2rad, accuracy: 1e-6)
        XCTAssertEqual(orbit.meanAnomaly(atTime: 0), 179.90874767107852.deg2rad, accuracy: 1e-6)
        XCTAssertEqual(orbit.eccentricAnomaly(atTime: 0), 179.90874767107852.deg2rad, accuracy: 1e-6)
        
        let atPeri = -4599439.396167472
        print("testing time since periapsis: \(orbit.timeSincePeriapsis(forTime: atPeri))")
        XCTAssertEqual(orbit.trueAnomaly(atTime: atPeri), 0, accuracy: 1e-6)
        XCTAssertEqual(orbit.meanAnomaly(atTime: atPeri), 0, accuracy: 1e-6)
        XCTAssertEqual(orbit.eccentricAnomaly(atTime: atPeri), 0, accuracy: 1e-6)
        
        let atRandom = 6123648.12356
        print("testing time since periapsis: \(orbit.timeSincePeriapsis(forTime: atRandom))")
        XCTAssertEqual(orbit.trueAnomaly(atTime: atRandom), 59.437475021205209.deg2rad, accuracy: 1e-6)
        XCTAssertEqual(orbit.meanAnomaly(atTime: atRandom), 59.437475021205209.deg2rad, accuracy: 1e-6)
        XCTAssertEqual(orbit.eccentricAnomaly(atTime: atRandom), 59.437475021205209.deg2rad, accuracy: 1e-6)
    }
}

//matt's library expects time sine peri
private extension Orbit {
    func timeSincePeriapsis(forTime time: Double) -> Double {
        let ttpe = self.timeToPeriapsis(atTime: time)
        return self.period - ttpe
    }
}

private extension Double {
    var deg2rad: Double {
        return self * .pi / 180
    }
}
