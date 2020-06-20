//
//  KerbalTime.swift
//  Keplerian
//
//  Created by Mike Muszynski on 5/30/20.
//  Copyright © 2020 Mike Muszynski. All rights reserved.
//

import Foundation

extension UnitDuration {
    static let kerbalYears = UnitDuration(symbol: "kspyear", converter: UnitConverterLinear(coefficient: 21600 * 426))
    static let kerbalDays = UnitDuration(symbol: "kspday", converter: UnitConverterLinear(coefficient: 21600))
}

/// A structure to encapsulate the Date format that KSP Uses
///
/// Since the date structure in Kerbal Space Program uses a non-standard definition for days and years, it became necessary to create a new data structure to facilitate conversions to and from KSP's version of UT.
public struct KSPDate {
    public var timeIntervalSinceReferenceDate: TimeInterval
    
    private var secondsInYear: Double = 1.kerbalYear.converted(to: .seconds).value
    private var secondsInDay: Double = 1.kerbalDay.converted(to: .seconds).value
    private var secondsInHour: Double = 1.hour.converted(to: .seconds).value
    private var secondsInMinute: Double = 60
    
    var year: Measurement<UnitDuration> {
        let years = floor(timeIntervalSinceReferenceDate / secondsInYear)
        return Measurement(value: years + 1, unit: .kerbalYears)
    }
    var day: Measurement<UnitDuration> {
        var remainder = timeIntervalSinceReferenceDate.remainder(dividingBy: secondsInYear)
        if remainder < 0 { remainder += secondsInYear }
        
        let days = floor(remainder / secondsInDay)
        return Measurement(value: days + 1, unit: .kerbalDays)
    }
    var hour: Measurement<UnitDuration> {
        var remainder = timeIntervalSinceReferenceDate.remainder(dividingBy: secondsInDay)
        if remainder < 0 { remainder += secondsInDay }
        
        let min = floor(remainder / secondsInHour)
        return Measurement(value: min, unit: .minutes)
    }
    var minute: Measurement<UnitDuration> {
        var remainder = timeIntervalSinceReferenceDate.remainder(dividingBy: secondsInHour)
        if remainder < 0 { remainder += secondsInHour }
        
        let min = floor(remainder / 60)
        return Measurement(value: min, unit: .minutes)
    }
    var second: Measurement<UnitDuration> {
        var secs = remainder(timeIntervalSinceReferenceDate, 60)
        if secs < 0 { secs += 60}
        
        return Measurement(value: secs, unit: .seconds)
    }
    
    public init(year: Int, day: Int, hour: Int = 0, minute: Int = 0, second: Int = 0) {
        var accum = TimeInterval(second)
        accum += secondsInMinute * TimeInterval(minute)
        accum += secondsInHour * TimeInterval(hour)
        accum += secondsInDay * TimeInterval(day - 1)
        accum += secondsInYear * TimeInterval(year - 1)
        self.timeIntervalSinceReferenceDate = accum
    }
    
    public init(timeIntervalSinceReferenceDate seconds: TimeInterval) {
        self.timeIntervalSinceReferenceDate = seconds
    }
}

extension KSPDate: CustomStringConvertible {
    public var description: String {
        let formatter = NumberFormatter()
        formatter.format = "00"
        return "KSPDate: Year \(Int(self.year.value)), Day \(Int(self.day.value)) " +
                "\(formatter.string(from: NSNumber(value: self.hour.value))!):" +
                "\(formatter.string(from: NSNumber(value: self.minute.value))!):" +
                "\(formatter.string(from: NSNumber(value: self.second.value))!)"
    }
}

// MARK: Date Arithmetic

/// Adds a time interval to a `KSPDate`
/// - Parameters:
///   - lhs: The original date
///   - rhs: A number of seconds to add
/// - Returns: A `KSPDate` later than the original date by a given `TimeInterval`
public func +(lhs: KSPDate, rhs: TimeInterval) -> KSPDate {
    return KSPDate(timeIntervalSinceReferenceDate: lhs.timeIntervalSinceReferenceDate + rhs)
}

/// Adds a `Measurement` of `UnitDuration` to a `KSPDate`
/// - Parameters:
///   - lhs: The original date
///   - rhs: A measurement of duration
/// - Returns: A `KSPDate` later than the original date by a `Measurement` of  `UnitDuration`
///
/// - Note: Any `UnitDuration` is acceptable, as it will first be converted to seconds and then added as a `TimeInterval`
public func +(lhs: KSPDate, rhs: Measurement<UnitDuration>) -> KSPDate {
    return lhs + TimeInterval(rhs.converted(to: .seconds).value)
}

/// Subtracts a time interval from a `KSPDate`
/// - Parameters:
///   - lhs: The original date
///   - rhs: A number of seconds to subtract
/// - Returns: A `KSPDate` earlier than the original date by a given `TimeInterval`
public func -(lhs: KSPDate, rhs: TimeInterval) -> KSPDate {
    return lhs + -rhs
}

/// Subtracts a `Measurement` of `UnitDuration` from a `KSPDate`
/// - Parameters:
///   - lhs: The original date
///   - rhs: A measurement of duration
/// - Returns: A `KSPDate` earlier than the original date by a `Measurement` of  `UnitDuration`
///
/// - Note: Any `UnitDuration` is acceptable, as it will first be converted to seconds and then subtracted as a `TimeInterval`
public func -(lhs: KSPDate, rhs: Measurement<UnitDuration>) -> KSPDate {
    return lhs - TimeInterval(rhs.converted(to: .seconds).value)
}

extension KSPDate: Comparable {
    public static func < (lhs: KSPDate, rhs: KSPDate) -> Bool {
        return lhs.timeIntervalSinceReferenceDate < rhs.timeIntervalSinceReferenceDate
    }
}

extension BinaryInteger {
    public var kerbalYear: Measurement<UnitDuration> {
        return Measurement(value: Double(self), unit: .kerbalYears)
    }
    
    public var kerbalYears: Measurement<UnitDuration> {
        return self.kerbalYear
    }
    
    public var kerbalDay: Measurement<UnitDuration> {
        return Measurement(value: Double(self), unit: .kerbalDays)
    }
    
    public var kerbalDays: Measurement<UnitDuration> {
        return self.kerbalDay
    }
    
    public var hour: Measurement<UnitDuration> {
        return Measurement(value: Double(self), unit: .hours)
    }
    
    public var hours: Measurement<UnitDuration> {
        return self.hour
    }
}

extension BinaryFloatingPoint {
    public var kerbalYear: Measurement<UnitDuration> {
        return Measurement(value: Double(self), unit: .kerbalYears)
    }
    
    public var kerbalYears: Measurement<UnitDuration> {
        return self.kerbalYear
    }
    
    public var kerbalDay: Measurement<UnitDuration> {
        return Measurement(value: Double(self), unit: .kerbalDays)
    }
    
    public var kerbalDays: Measurement<UnitDuration> {
        return self.kerbalDay
    }
    
    public var hour: Measurement<UnitDuration> {
        return Measurement(value: Double(self), unit: .hours)
    }
    
    public var hours: Measurement<UnitDuration> {
        return self.hour
    }
}
