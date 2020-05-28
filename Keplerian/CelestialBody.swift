//
//  CelestialBody.swift
//  KSPCockpitPanel
//
//  Created by Mike Muszynski on 7/22/17.
//  Copyright © 2017 Mike Muszynski. All rights reserved.
//

import Foundation

public struct CelestialBody: Equatable {
    public var gravitationalParameter: Double
    public var radius: Double
    public var atmosphereAltitude: Double = 70000
    public var orbit: Orbit?
    
    public var safeApoapsis: Double {
        return self.radius + self.atmosphereAltitude
    }
    
    public var safeAltitude: Double {
        return atmosphereAltitude
    }
    
    public init(gravitationalParameter: Double, radius: Double, orbit: Orbit? = nil) {
        self.gravitationalParameter = gravitationalParameter
        self.radius = radius
        self.orbit = orbit
    }
}

extension CelestialBody {
    static var Kerbol = CelestialBody(gravitationalParameter: 1.1723328e18, radius: 261600000.0)
    static var Kerbin = CelestialBody(gravitationalParameter: 3.5316000e12, radius: 600000.0, orbit: .kerbin)
    static var Duna = CelestialBody(gravitationalParameter: 3.0136321e11, radius: 320000.0, orbit: .duna)
}
