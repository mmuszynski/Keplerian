//
//  ManeuverTests.swift
//  
//
//  Created by Mike Muszynski on 12/15/22.
//

import XCTest
@testable import Keplerian

final class ManeuverTests: XCTestCase {
    
    /// Tests whether the prograde direction is correct for an orbit
    func testProgradeDirection() {
        let orbit = Orbit(semiMajorAxis: (700000.0).meters, eccentricity: 0.05, meanAnomaly: 0.0.rad, inclination: 0.0.rad, LAN: 0.0.rad, argumentOfPeriapsis: 0.0.rad, centralBody: .kerbin)
        print(orbit.period)
        print(orbit.epoch)
        
        let maneuver = Orbit.Maneuver(progradeDV: 100.0)
    }

}
