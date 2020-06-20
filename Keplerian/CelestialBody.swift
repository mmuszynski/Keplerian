//
//  CelestialBody.swift
//  KSPCockpitPanel
//
//  Created by Mike Muszynski on 7/22/17.
//  Copyright © 2017 Mike Muszynski. All rights reserved.
//

import Foundation

public struct CelestialBody: Equatable, Codable {
    public var gravitationalParameter: Double
    public var radius: Double
    public var atmosphereAltitude: Double = 70000
    public var orbit: Orbit?
    public var mass: Double
    
    public var sphereOfInfluence: Double? {
        guard let orbit = orbit else { return nil }
        let a = orbit.semiMajorAxis
        let m = self.mass
        let M = orbit.centralBody.mass
        
        return a * pow(m / M, 2 / 5)
    }
    
    public var safeApoapsis: Double {
        return self.radius + self.atmosphereAltitude
    }
    
    public var safeAltitude: Double {
        return atmosphereAltitude
    }
    
    public init(gravitationalParameter: Double, radius: Double, mass: Double, orbit: Orbit? = nil) {
        self.gravitationalParameter = gravitationalParameter
        self.radius = radius
        self.orbit = orbit
        self.mass = mass
    }
}
