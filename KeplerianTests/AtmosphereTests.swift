//
//  AtmosphereTests.swift
//  
//
//  Created by Mike Muszynski on 1/5/23.
//

import XCTest
@testable import Keplerian

final class AtmosphereTests: XCTestCase {

    func testKerbinAtmosphere() throws {
        let kerbin = CelestialBody.kerbin
        XCTAssertEqual(kerbin.safeAltitude, 70000)
        //print(kerbin.atmosphere!.pressure(at: 70000) / kerbin.atmosphere!.pressure(at: 0))
    }

}
