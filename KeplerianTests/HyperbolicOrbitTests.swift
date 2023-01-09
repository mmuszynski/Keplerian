//
//  OrbitTests.swift
//  KeplerianTests
//
//  Created by Mike Muszynski on 5/29/20.
//  Copyright © 2020 Mike Muszynski. All rights reserved.
//

import XCTest
@testable import Keplerian

class HyperbolicOrbitTests: XCTestCase {
    
    func orbit(named filename: String) throws -> Orbit {
        let orbitURL = Bundle.module.url(forResource: filename, withExtension: "json")!
        let orbitData = try Data(contentsOf: orbitURL)
        let decoder = JSONDecoder()
        let orbit = try decoder.decode(Orbit.self, from: orbitData)
        return orbit
    }
    
    func testRealEscapeOrbit() throws {
        let orbit = try orbit(named: "munEscapeAltitude213772160")
        XCTAssertEqual(orbit.altitude(atTime: orbit.epoch), 2137721.60, accuracy: 1.0)
    }

    func testAnotherEscapeOrbit() throws {
        /*
         This orbit was collected on 1/6/22 on the way to the Mun.
         KSP reports the True Anomaly (through Kerbal Engineer) as -127.18985 deg
         Mean Anomaly -135.50482 deg
         Altitude 2119000 (2119.0km)
         */
        
        let reportedTrueAnomaly = -127.18985.deg2rad
        let reportedMeanAnomaly = -135.50482.deg2rad
                
        let orbit = try orbit(named: "munEscape20230106")
        
        print(orbit.timeToPeriapsis(atTime: orbit.epoch))
        print(orbit.hyperbolicAnomaly(fromMeanAnomaly: orbit.meanAnomaly(atTime: 4107753.75)))
        
        XCTAssertEqual(orbit.trueAnomaly(atTime: 4107753.75), reportedTrueAnomaly, accuracy: 0.0001)
        XCTAssertEqual(orbit.altitude(atTime: 4107753.75), 2119000, accuracy: 1000)
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
