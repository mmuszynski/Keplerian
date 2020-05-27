//
//  Orbit+TestValues.swift
//  KSPTelemetry
//
//  Created by Mike Muszynski on 5/24/20.
//  Copyright © 2020 Mike Muszynski. All rights reserved.
//

import Foundation

extension Orbit {
    public class var testEscapeOrbit: Orbit {
        let mun = CelestialBody(gravitationalParameter: 65138397184.0, radius: 2e6)
        let orbit = Orbit(semiMajorAxis: -944427.0625, eccentricity: 1.22502624988556, meanAnomaly: 0.196598216891289, inclination: 0.129754841327667, LAN: 162.380401611328, argumentOfPeriapsis: 356.085906982422, centralBody: mun)
        return orbit
    }
    
    public class var exampleCircularOrbit: Orbit {
        let earth = CelestialBody(gravitationalParameter: 3.986e14, radius: 6371000)
        let orbit = Orbit(semiMajorAxis: 35901000, eccentricity: 0.0, meanAnomaly: 0.196598216891289, inclination: 0.129754841327667, LAN: 30, argumentOfPeriapsis: 30, centralBody: earth)
        return orbit
    }
    
    public class var exampleEllipticalOrbit: Orbit {
        let earth = CelestialBody(gravitationalParameter: 3.986e14, radius: 6371000)
        let orbit = Orbit(semiMajorAxis: 35901000, eccentricity: 0.5, meanAnomaly: 0.196598216891289, inclination: 0.129754841327667, LAN: 30, argumentOfPeriapsis: 30, centralBody: earth)
        return orbit
    }
    
    public class var exampleSuborbitalOrbit: Orbit {
        let earth = CelestialBody(gravitationalParameter: 3.986e14, radius: 6371000)
        let orbit = Orbit(semiMajorAxis: 6371000 * 0.5 + 50, eccentricity: 0.9999, meanAnomaly: 0.0, inclination: 0.129754841327667, LAN: 30, argumentOfPeriapsis: 30, centralBody: earth)
        return orbit
    }
    
    public class var kerbin: Orbit {
        Orbit(semiMajorAxis: 13599840256 , eccentricity: 0, meanAnomaly: 3.14, inclination: 0, LAN: 0, argumentOfPeriapsis: 0, centralBody: .kerbol)
    }
    
    public class var duna: Orbit {
        Orbit(semiMajorAxis: 20726155264, eccentricity: 0.051, meanAnomaly: 3.14, inclination: 0.06, LAN: 135.5, argumentOfPeriapsis: 0, centralBody: .kerbol)
    }
}
