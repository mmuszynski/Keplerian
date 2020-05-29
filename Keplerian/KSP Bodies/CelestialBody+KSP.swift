//
//  CelestialBody+KSP.swift
//  Keplerian
//
//  Created by Mike Muszynski on 5/29/20.
//  Copyright © 2020 Mike Muszynski. All rights reserved.
//

import Foundation

extension CelestialBody {
    static var kerbol = CelestialBody(gravitationalParameter: 1.1723328e18, radius: 2.616e9)
    static var kerbin = CelestialBody(gravitationalParameter: 3.5316000e12, radius: 6e6, orbit: .kerbin)
    static var duna = CelestialBody(gravitationalParameter: 3.0136321e11, radius: 3.2e6, orbit: .duna)
    static var mun = CelestialBody(gravitationalParameter: 6.5138398e10, radius: 2e6, orbit: .mun)
    static var minmus = CelestialBody(gravitationalParameter: 1.7658000e9, radius: 6e5, orbit: .minmus)
}
