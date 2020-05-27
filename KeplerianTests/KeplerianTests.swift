//
//  KeplerianTests.swift
//  KeplerianTests
//
//  Created by Mike Muszynski on 8/7/17.
//  Copyright © 2017 Mike Muszynski. All rights reserved.
//

import XCTest
@testable import Keplerian

extension Double {
    static var J2000: Double = 2451545.0
}

class KeplerianTests: XCTestCase {
    
    func testMeasurementConversion() {
        let thousand = 1.km
        let alsoThousand = 1000.meters
        
        XCTAssertEqual(thousand, alsoThousand)
        XCTAssertEqual(thousand.value, 1.0)
        XCTAssertEqual(alsoThousand.value, 1000.0)
        
        XCTAssertEqual(Measurement(value: 1, unit: UnitLength.kilometers).converted(to: UnitLength.meters).value, 1000.0)
    }
    
    func testInitializeByMeasurement() {
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
        XCTAssertEqual(initialPosition, orbit.cartesianPosition(atTime: 0))
        
        let eccentricOrbit = Orbit(semiMajorAxis: 2000, eccentricity: 0.5, meanAnomaly: 0, inclination: 0, LAN: 0, argumentOfPeriapsis: 0, centralBody: sun)
        XCTAssertEqual(initialPosition, orbit.cartesianPosition(atTime: 0))
        
        let expectedPosition = Vector3D(x: -3000, y: 0, z: 0)
        let diff = expectedPosition - eccentricOrbit.cartesianPosition(atTime: 8885.765876316733)
        XCTAssertEqual(diff.magnitude, 0, accuracy: 1e-6)
    }
    
    func testTrivialLambertSolver() {
        let sun = CelestialBody(gravitationalParameter: 1000, radius: 1000)
        let orbit = Orbit(semiMajorAxis: 2000, eccentricity: 0, meanAnomaly: 0, inclination: 0, LAN: 0, argumentOfPeriapsis: 0, centralBody: sun)
        
        let initial = orbit.cartesianPosition(atTime: 0)
        let final = orbit.cartesianPosition(atTime: orbit.period / 4.0)

        let solver = LambertSolver(position1: initial, position2: final, dt: 1000)
        let solutions = solver.solve()
        print(solutions)
    }
    
    func testLambertSolver() {
        let position1 = Vector3D(x: 147084764.907217, y: -32521189.6497507, z: 467.190091409394)
        let position2 = Vector3D(x: -88002509.1583767, y: -62680223.1330849, z: 4220331.52492018)
        
        let solver = LambertSolver(position1: position1, position2: position2, dt: (2455610-2455450)*86400, mu: 1.32712440018e11)
        let solutions = solver.solve()
        
        guard !solutions.isEmpty else {
            XCTFail("No solutions found")
            return
        }
        
        XCTAssertEqual(solutions[0].0.x, 4.65144349746008, accuracy: 1e-5)
        XCTAssertEqual(solutions[0].0.y, 26.0824144093203, accuracy: 1e-5)
        XCTAssertEqual(solutions[0].0.z, -1.39306043231699, accuracy: 1e-5)
        XCTAssertEqual(solutions[0].1.x, 16.7926204519414, accuracy: 1e-5)
        XCTAssertEqual(solutions[0].1.y, -33.3516748429805, accuracy: 1e-5)
        XCTAssertEqual(solutions[0].1.z, 1.52302150358741, accuracy: 1e-5)
    }
    
    func testLambertSolver2() {
        let position1 = Vector3D(x: 170145121.321308, y: -117637192.836034, z: -6642044.2724648)
        let position2 = Vector3D(x: -803451694.669228, y: 121525767.116065, z: 17465211.7766441)
        
        let solver = LambertSolver(position1: position1, position2: position2, dt: (2457500-2456300)*86400, mu: 1.32712440018e11)
        let solutions = solver.solve()
                
        XCTAssertEqual(solutions[0].0.x, 13.7407773577481)
        XCTAssertEqual(solutions[0].0.y, 28.8309931231422)
        XCTAssertEqual(solutions[0].0.z, 0.691285008034955)
        XCTAssertEqual(solutions[0].1.x, -0.883933068957334)
        XCTAssertEqual(solutions[0].1.y, -7.98362701426338)
        XCTAssertEqual(solutions[0].1.z, -0.240770597841448)
    }
    
//    func testNonTrivialLambertSolver() {
//        let sun = CelestialBody(gravitationalParameter: 1.32712440018e11, radius: 6390000)
//        let earthOrbit = Orbit(semiMajorAxis: 146.9e6.km,
//                               eccentricity: 0.0167086,
//                               meanAnomaly: 358.617.deg,
//                               inclination: 0.deg,
//                               LAN: 174.9.deg,
//                               argumentOfPeriapsis: 288.1.deg,
//                               centralBody: sun)
//        earthOrbit.epoch = .J2000
//
//        let venusOrbit = Orbit(semiMajorAxis: 10820800.0.km,
//                               eccentricity: 0.006772,
//                               meanAnomaly: 50.115.deg,
//                               inclination: 3.39458.deg,
//                               LAN: 76.680.deg,
//                               argumentOfPeriapsis: 54.884.deg,
//                               centralBody: sun)
//        venusOrbit.epoch = .J2000
//
//        let position1 = earthOrbit.cartesianPosition(atTime: 2455450)
//        let position2 = venusOrbit.cartesianPosition(atTime: 2455610)
//
//        let expectedV1 = Vector3D(x: 4.65144349746008, y: 26.0824144093203, z: -1.39306043231699)
//        let expectedV2 = Vector3D(x: 16.7926204519414, y: -33.3516748429805, z: 1.52302150358741)
//
//        let solver = LambertSolver(position1: position1, position2: position2, dt: 2455610-2455450, mu: sun.gravitationalParameter)
//        let solutions = solver.solve()
//        print(solutions)
//
//        print(position1)
//        print(position2)
//    }
    
}
