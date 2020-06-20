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
    
    func testAlexMoon() {
        /*
         https://alexmoon.github.io/ksp/#/Kerbin/0/Duna/0/false/ballistic/false/1/1
         Github user alexmoon has a detailed chart for this stuff.
         Can I match up the DV with what he suggests is possible for the trip?
         
         Let's use the lambert solver and see whether it matches up
         */
        
        //Year 1, day 231 at 0:14:24
        let departrureDate = KSPDate(year: 1, day: 231, hour: 0, minute: 14, second: 24)
        
        //Year 2, day 74 at 4:57:36
        let arrivalDate = KSPDate(year: 2, day: 74, hour: 4, minute: 57, second: 36)
        
        
        let solver = try! LambertSolver(orbit1: .kerbin, orbit2: .duna, departureTime: departrureDate, travelTime: arrivalDate.timeIntervalSinceReferenceDate - departrureDate.timeIntervalSinceReferenceDate)
        
        let (departureVelocity, arrivalVelocity, _) = solver.solve()[0]
        let vinf_departure = departureVelocity - Orbit.kerbin.cartesian(atTime: departrureDate.timeIntervalSinceReferenceDate).velocity
        let vinf_arrival = arrivalVelocity - Orbit.duna.cartesian(atTime: arrivalDate.timeIntervalSinceReferenceDate).velocity
        
        print((vinf_arrival + vinf_departure).magnitude)
    }
    
//    func testKerbalPorckchopAccuracy() {
//        let window = (KSPDate(year: 1, day: 100),
//                      KSPDate(year: 1, day: 300))
//        
//        let kerbinToDuna = Porkchop(from: .kerbin, to: .duna, departureWindow: window)
//        try! kerbinToDuna.solve()
//        print(kerbinToDuna.solutionSpace.sorted(by: { $0.dV < $1.dV }).first)
//    }
//    
//    func testKerbalPorckchopAccuracySecondWindow() {
//        let window = (KSPDate(year: 1, day: 900),
//                      KSPDate(year: 1, day: 1200))
//        
//        let kerbinToDuna = Porkchop(from: .kerbin, to: .duna, departureWindow: window)
//        try! kerbinToDuna.solve()
//        print(kerbinToDuna.solutionSpace.sorted(by: { $0.dV < $1.dV }).first)
//    }
//    
    func testPorkchopToEve() {
        let window = (KSPDate(year: 1, day: 0), KSPDate(year: 1, day: 852))
        let chop = Porkchop(from: .kerbin, to: .eve, departureWindow: window)
        try! chop.solve()
    }
}
