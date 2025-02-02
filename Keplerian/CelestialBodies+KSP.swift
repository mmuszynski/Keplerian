//
//  CelestialBodies+Universe.swift
//  Keplerian
//
//  Created by Mike Muszynski on 2/2/25.
//

fileprivate extension Double {
    var deg2rad: Double {
        return self * .pi / 180
    }
}

extension CelestialBody {
    public static let kerbol = CelestialBody(gravitationalParameter: 1.1723328e18, radius: 261600000, mass: 1.7565459e28)
    public static let moho = CelestialBody(gravitationalParameter: 1.6860938e11, radius: 250000, mass: 2.5263314e21, orbit: .moho)
    public static let eve = CelestialBody(gravitationalParameter: 8.1717302e12, radius: 700000, mass: 1.2243980e23, orbit: .eve)
    public static let kerbin = {
        var k = CelestialBody(gravitationalParameter: 3.5316000e12, radius: 6e5, mass: 5.2915158e22, orbit: .kerbin)
        k.atmosphereAltitude = 70000
        return k
    }()
    public static let mun = CelestialBody(gravitationalParameter: 6.5138398e10, radius: 2e5, mass: 9.7599066e20, orbit: .mun)
    public static let minmus = CelestialBody(gravitationalParameter: 1.7658000e9, radius: 6e4, mass: 2.6457580e19, orbit: .minmus)
    public static let duna = {
        var d = CelestialBody(gravitationalParameter: 3.0136321e11, radius: 320000, mass: 4.5154270e21, orbit: .duna)
        d.atmosphereAltitude = 50000
        return d
    }()
    public static let ike = CelestialBody(gravitationalParameter: 1.8568369e10, radius: 130000, mass: 2.7821615e20, orbit: .ike)
    public static let jool = CelestialBody(gravitationalParameter: 2.8252800e14, radius: 6000000, mass: 4.2332127e24, orbit: .jool)
    public static let dres = CelestialBody(gravitationalParameter: 2.1484489e10, radius: 138000, mass: 3.2190937e20, orbit: .dres)
    public static let eeloo = CelestialBody(gravitationalParameter: 7.4410815e10, radius: 210000, mass: 1.1149224e21, orbit: .eeloo)
    
    public static let allKSPBodies: [CelestialBody] = [
        .kerbol, .moho, .eve, .kerbin, .mun, .minmus, .duna, .jool, . dres, .eeloo, .ike
    ]
}

extension Orbit {
    public static let moho = Orbit(semiMajorAxis: 5263138304.0.meters,
                                   eccentricity: 0.2,
                                   meanAnomaly: 3.14.rad,
                                   inclination: 7.0.deg,
                                   LAN: 70.0.deg,
                                   argumentOfPeriapsis: 15.0.deg,
                                   centralBody: .kerbol)
    
    public static let eve = Orbit(semiMajorAxis: 9832684544.meters,
                                  eccentricity: 0.01,
                                  meanAnomaly: 3.14.radians,
                                  inclination: 2.1.degrees,
                                  LAN: 15.degrees,
                                  argumentOfPeriapsis: 0.degrees,
                                  centralBody: .kerbol)

    public static let kerbin = Orbit(semiMajorAxis: 13599840256,
                                     eccentricity: 0,
                                     meanAnomaly: 3.14,
                                     inclination: 0,
                                     LAN: 0,
                                     argumentOfPeriapsis: 0,
                                     centralBody: .kerbol)
    
    public static let mun = Orbit(semiMajorAxis: 1.2e7,
                                  eccentricity: 0,
                                  meanAnomaly: 1.7,
                                  inclination: 0,
                                  LAN: 0,
                                  argumentOfPeriapsis: 0,
                                  centralBody: .kerbin)
    
    public static let minmus = Orbit(semiMajorAxis: 4.7e7.meters,
                                     eccentricity: 0.0,
                                     meanAnomaly: 0.9.radians,
                                     inclination: 6.0.degrees,
                                     LAN: 78.0.degrees,
                                     argumentOfPeriapsis: 38.0.degrees,
                                     centralBody: .kerbin)
    
    public static let duna = Orbit(semiMajorAxis: 20726155264.meters,
                                   eccentricity: 0.051,
                                   meanAnomaly: 3.14.radians,
                                   inclination: 0.06.degrees,
                                   LAN: 135.5.degrees,
                                   argumentOfPeriapsis: 0.degrees,
                                   centralBody: .kerbol)
    
    public static let jool = Orbit(semiMajorAxis: 68773560320.0.meters,
                                   eccentricity: 0.05,
                                   meanAnomaly: 1.304.radians,
                                   inclination: 0.1.degrees,
                                   LAN: 52.0.degrees,
                                   argumentOfPeriapsis: 0.0.degrees,
                                   centralBody: .kerbol)
    
    public static let dres = Orbit(semiMajorAxis: 40839348203.0.meters,
                                   eccentricity: 0.145,
                                   meanAnomaly: 5.radians,
                                   inclination: 3.14.degrees,
                                   LAN: 280.degrees,
                                   argumentOfPeriapsis: 90.degrees,
                                   centralBody: .kerbol)
    
    public static let eeloo = Orbit(semiMajorAxis: 90118820000.0.meters,
                                    eccentricity: 0.26,
                                    meanAnomaly: 6.15.radians,
                                    inclination: 3.14.degrees,
                                    LAN: 50.degrees,
                                    argumentOfPeriapsis: 260.degrees,
                                    centralBody: .kerbol)
    
    public static let ike = Orbit(semiMajorAxis: 3200000.meters,
                                  eccentricity: 0.03,
                                  meanAnomaly: 1.7.radians,
                                  inclination: 0.2.degrees,
                                  LAN: 0.degrees,
                                  argumentOfPeriapsis: 0.degrees,
                                  centralBody: .duna)
}
