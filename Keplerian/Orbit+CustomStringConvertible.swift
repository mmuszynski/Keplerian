//
//  Orbit+CustomStringConvertible.swift
//  KSPTelemetry
//
//  Created by Mike Muszynski on 5/24/20.
//  Copyright © 2020 Mike Muszynski. All rights reserved.
//

import Foundation

extension Orbit: CustomDebugStringConvertible {
    public var debugDescription: String {
        let sma = String(describing: semiMajorAxis)
        let ecc = String(describing: eccentricity)
        let meA = String(describing: meanAnomaly)
        let inc = String(describing: inclination)
        let lan = String(describing: LAN)
        let arg = String(describing: argumentOfPeriapsis)
        
        var string = "semiMajorAxis: " + sma + "\r"
        string += "eccentricity: " + ecc + "\r"
        string += "meanAnomaly: " + meA + "\r"
        string += "inclination: " + inc + "\r"
        string += "LAN: " + lan + "\r"
        string += "argumentOfPeriapsis: " + arg + "\r"
        string += "bodyRadius: " + String(describing: centralBody.radius) + "\r"
        string += "bodyGravParameter: " + String(describing: centralBody.gravitationalParameter)

        return string
    }
}
