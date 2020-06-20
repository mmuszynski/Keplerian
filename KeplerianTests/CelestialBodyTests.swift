//
//  CelestialBodyTests.swift
//  KeplerianTests
//
//  Created by Mike Muszynski on 6/8/20.
//  Copyright © 2020 Mike Muszynski. All rights reserved.
//

import XCTest
@testable import Keplerian

class CelestialBodyTests: XCTestCase {

    func testSOI() {
        let kerbin = CelestialBody.kerbin
        XCTAssertEqual(kerbin.sphereOfInfluence!, 84159286, accuracy: 1)
    }

}
