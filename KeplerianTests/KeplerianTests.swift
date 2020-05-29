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
    
    func testTrivialLambertSolver() {
        let sun = CelestialBody(gravitationalParameter: 1000, radius: 1000)
        let orbit = Orbit(semiMajorAxis: 2000, eccentricity: 0, meanAnomaly: 0, inclination: 0, LAN: 0, argumentOfPeriapsis: 0, centralBody: sun)
        
        let initial = orbit.cartesianPosition(atTime: 0)
        let final = orbit.cartesianPosition(atTime: orbit.period / 4.0)

        let solver = LambertSolver(position1: initial, position2: final, dt: 1000)
        let solutions = solver.solve()
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
                
        XCTAssertEqual(solutions[0].0.x, 13.7407773577481, accuracy: 1e-5)
        XCTAssertEqual(solutions[0].0.y, 28.8309931231422, accuracy: 1e-5)
        XCTAssertEqual(solutions[0].0.z, 0.691285008034955, accuracy: 1e-5)
        XCTAssertEqual(solutions[0].1.x, -0.883933068957334, accuracy: 1e-5)
        XCTAssertEqual(solutions[0].1.y, -7.98362701426338, accuracy: 1e-5)
        XCTAssertEqual(solutions[0].1.z, -0.240770597841448, accuracy: 1e-5)
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
    
    func testKerbalPorkchop() {
        let kerbinToDuna = Porkchop(from: .kerbin, to: .duna, departureWindow: (0, 400))
        measure {
            try! kerbinToDuna.solve()
        }
    }
    
    func testKerbalPorckchopAccuracy() {
        let kerbinToDuna = Porkchop(from: .kerbin, to: .duna, departureWindow: (100, 400))
        try! kerbinToDuna.solve()
        print(kerbinToDuna.solutionSpace.sorted(by: { $0.dV < $1.dV })[0...10] )
    }
}
