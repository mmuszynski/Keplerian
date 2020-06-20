//
//  KSPDateTests.swift
//  KeplerianTests
//
//  Created by Mike Muszynski on 5/31/20.
//  Copyright © 2020 Mike Muszynski. All rights reserved.
//

import XCTest
@testable import Keplerian

class KSPDateTests: XCTestCase {
    
    func testTime() {
        let zero = KSPDate(year: 1, day: 1, hour: 0, minute: 0, second: 0)
        XCTAssertEqual(zero.timeIntervalSinceReferenceDate, 0)
        
        let day2 = KSPDate(year: 1, day: 2, hour: 0, minute: 0, second: 0)
        XCTAssertEqual(day2.timeIntervalSinceReferenceDate, 21600)
        
        let year2 = KSPDate(year: 2, day: 1)
        XCTAssertEqual(year2.timeIntervalSinceReferenceDate, 3600 * 6 * 426)
    }
    
    func testInitFromInterval() {
        let interval: TimeInterval = 3600 * 6 * 426
        let intervalYear = KSPDate(timeIntervalSinceReferenceDate: interval)
        let year2 = KSPDate(year: 2, day: 1, hour: 0, minute: 0, second: 0)
        
        XCTAssertEqual(intervalYear, year2)
    }
    
    func testComparable() {
        let first = KSPDate(year: 1, day: 1, hour: 0, minute: 0, second: 0)
        let second = KSPDate(year: 2, day: 1, hour: 0, minute: 0, second: 0)
        
        XCTAssertTrue(first < second)
    }
    
    func testEquality() {
        let first = KSPDate(year: 1, day: 1, hour: 0, minute: 0, second: 0)
        let alsoFirst = KSPDate(year: 1, day: 1, hour: 0, minute: 0, second: 0)
        
        XCTAssertEqual(first, alsoFirst)
    }
    
    func testAddition() {
        let first = KSPDate(year: 1, day: 1, hour: 0, minute: 0, second: 0)
        let second = KSPDate(year: 2, day: 1, hour: 0, minute: 0, second: 0)
        let kspYear = 1.kerbalYear
        
        XCTAssertEqual(second, first + kspYear)
    }
    
    func testHourRollover() {
        let first = KSPDate(year: 1, day: 1, hour: 0, minute: 60, second: 0)
        let second = KSPDate(year: 1, day: 1, hour: 1, minute: 0, second: 0)
        
        XCTAssertEqual(first.timeIntervalSinceReferenceDate, 60 * 60)
        XCTAssertEqual(second.timeIntervalSinceReferenceDate, 60 * 60)
        XCTAssertEqual(first, second)
    }
    
    func testComponentParts() {
        let date = KSPDate(year: 40, day: 93, hour: 4, minute: 39, second: 12)
        XCTAssertEqual(date.year.value, 40)
        XCTAssertEqual(date.day.value, 93)
        XCTAssertEqual(date.hour.value, 4)
        XCTAssertEqual(date.minute.value, 39)
        XCTAssertEqual(date.second.value, 12)
    }
    
    func testDescription() {
        let date = KSPDate(year: 1, day: 20, hour: 3, minute: 14, second: 15)
        XCTAssertEqual(date.description, "KSPDate: Year 1, Day 20 03:14:15")
    }

}
