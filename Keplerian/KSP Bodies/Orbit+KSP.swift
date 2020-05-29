//
//  Orbit+KSP.swift
//  Keplerian
//
//  Created by Mike Muszynski on 5/29/20.
//  Copyright © 2020 Mike Muszynski. All rights reserved.
//

import Foundation

private extension Double {
    var deg2rad: Double {
        return self * .pi / 180
    }
}

extension Orbit {
    public static var kerbin = Orbit(semiMajorAxis: 13599840256,
                                     eccentricity: 0,
                                     meanAnomaly: 3.14,
                                     inclination: 0,
                                     LAN: 0,
                                     argumentOfPeriapsis: 0,
                                     centralBody: .kerbol)
    
    public static var duna = Orbit(semiMajorAxis: 20726155264,
                                   eccentricity: 0.051,
                                   meanAnomaly: 3.14,
                                   inclination: 0.06,
                                   LAN: 135.5,
                                   argumentOfPeriapsis: 0,
                                   centralBody: .kerbol)
                                   
    public static var mun = Orbit(semiMajorAxis: 1.2e7,
                                  eccentricity: 0,
                                  meanAnomaly: 1.7,
                                  inclination: 0,
                                  LAN: 0,
                                  argumentOfPeriapsis: 0,
                                  centralBody: .kerbin)
    
    public static var minmus = Orbit(semiMajorAxis: 4.7e7.meters,
                                     eccentricity: 0.0,
                                     meanAnomaly: 0.9.radians,
                                     inclination: 6.0.degrees,
                                     LAN: 78.0.degrees,
                                     argumentOfPeriapsis: 38.0.degrees,
                                     centralBody: .kerbin)
}
