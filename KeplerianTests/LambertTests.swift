//
//  LambertTests.swift
//  KeplerianTests
//
//  Created by Mike Muszynski on 5/30/20.
//  Copyright © 2020 Mike Muszynski. All rights reserved.
//

import XCTest
@testable import Keplerian

class LambertTests: XCTestCase {

    func testLambertSolver() {
        let position1 = Vector3D(x: 147084764.907217, y: -32521189.6497507, z: 467.190091409394)
        let position2 = Vector3D(x: -88002509.1583767, y: -62680223.1330849, z: 4220331.52492018)
        
        let solver = LambertSolver(position1: position1, position2: position2, dt: (2455610-2455450)*86400, mu: 1.32712440018e11)
        let solutions = solver.solve()
        
        guard !solutions.isEmpty else {
            XCTFail("No solutions found")
            return
        }
        
        XCTAssertEqual(solutions[0].0.x, 4.65144349746008, accuracy: 1e-5)
        XCTAssertEqual(solutions[0].0.y, 26.0824144093203, accuracy: 1e-5)
        XCTAssertEqual(solutions[0].0.z, -1.39306043231699, accuracy: 1e-5)
        XCTAssertEqual(solutions[0].1.x, 16.7926204519414, accuracy: 1e-5)
        XCTAssertEqual(solutions[0].1.y, -33.3516748429805, accuracy: 1e-5)
        XCTAssertEqual(solutions[0].1.z, 1.52302150358741, accuracy: 1e-5)
    }
    
    func testLambertSolver2() {
        let position1 = Vector3D(x: 170145121.321308, y: -117637192.836034, z: -6642044.2724648)
        let position2 = Vector3D(x: -803451694.669228, y: 121525767.116065, z: 17465211.7766441)
        
        let solver = LambertSolver(position1: position1, position2: position2, dt: (2457500-2456300)*86400, mu: 1.32712440018e11)
        let solutions = solver.solve()
                
        XCTAssertEqual(solutions[0].0.x, 13.7407773577481, accuracy: 1e-5)
        XCTAssertEqual(solutions[0].0.y, 28.8309931231422, accuracy: 1e-5)
        XCTAssertEqual(solutions[0].0.z, 0.691285008034955, accuracy: 1e-5)
        XCTAssertEqual(solutions[0].1.x, -0.883933068957334, accuracy: 1e-5)
        XCTAssertEqual(solutions[0].1.y, -7.98362701426338, accuracy: 1e-5)
        XCTAssertEqual(solutions[0].1.z, -0.240770597841448, accuracy: 1e-5)
    }
    
    func testIJKToDuna() {
        let departure = Orbit.kerbin
        let destination = Orbit.duna
        
        let departureTime = KSPDate(year: 1, day: 200, hour: 0, minute: 0, second: 0)
        let arrivalTime = KSPDate(year: 1, day: 300, hour: 0, minute: 0, second: 0)
        
        print(arrivalTime.timeIntervalSinceReferenceDate - departureTime.timeIntervalSinceReferenceDate)
        print(departure.coe2rv(atTime: departureTime.timeIntervalSinceReferenceDate))
        print(destination.coe2rv(atTime: arrivalTime.timeIntervalSinceReferenceDate))
        
        let departureIJK = departure.cartesian(atTime: Double(departureTime.timeIntervalSinceReferenceDate))
        let destinationIJK = destination.cartesian(atTime: Double(arrivalTime.timeIntervalSinceReferenceDate))
        
        let departureIJKTest = [ 13304731552.328235626220703,
                                 -2817760335.656604290008545,
                                 0.000000000000000,
                                 1923.662144854920371,
                                 9083.032403714601969,
                                 0.000000000000000]
        XCTAssertEqual(departureIJK.position.x, departureIJKTest[0], accuracy: 1e-5)
        XCTAssertEqual(departureIJK.position.y, departureIJKTest[1], accuracy: 1e-5)
        XCTAssertEqual(departureIJK.position.z, departureIJKTest[2], accuracy: 1e-5)
        XCTAssertEqual(departureIJK.velocity.x, departureIJKTest[3], accuracy: 1e-5)
        XCTAssertEqual(departureIJK.velocity.y, departureIJKTest[4], accuracy: 1e-5)
        XCTAssertEqual(departureIJK.velocity.z, departureIJKTest[5], accuracy: 1e-5)
        
        let destinationIJKTest = [
            1639398431.6417618,
        19951441866.660767,
        -16105322.423889095,
        -7774.5314008089863,
        342.78087632163533,
        5.450406033953969
        ]
        XCTAssertEqual(destinationIJK.position.x, destinationIJKTest[0], accuracy: 1e-5)
        XCTAssertEqual(destinationIJK.position.y, destinationIJKTest[1], accuracy: 1e-5)
        XCTAssertEqual(destinationIJK.position.z, destinationIJKTest[2], accuracy: 1e-5)
        XCTAssertEqual(destinationIJK.velocity.x, destinationIJKTest[3], accuracy: 1e-5)
        XCTAssertEqual(destinationIJK.velocity.y, destinationIJKTest[4], accuracy: 1e-5)
        XCTAssertEqual(destinationIJK.velocity.z, destinationIJKTest[5], accuracy: 1e-5)
    }
    
    func testIJKEve() {
        let departure = Orbit.kerbin
        let destination = Orbit.eve
        
        let departureTime = KSPDate(year: 1, day: 200, hour: 0, minute: 0, second: 0)
        let arrivalTime = KSPDate(year: 1, day: 300, hour: 0, minute: 0, second: 0)
        
        print(departure.timeSincePeriapsis(forTime: departureTime.timeIntervalSinceReferenceDate))
        print(departure.anom(atTime: departureTime.timeIntervalSinceReferenceDate))
        print(departure.coe(atTime: departureTime.timeIntervalSinceReferenceDate))
        
        print(departure.coe2rv(atTime: departureTime.timeIntervalSinceReferenceDate))
        print(destination.coe2rv(atTime: arrivalTime.timeIntervalSinceReferenceDate))
        
        let departureIJK = departure.cartesian(atTime: Double(departureTime.timeIntervalSinceReferenceDate))
        let destinationIJK = destination.cartesian(atTime: Double(arrivalTime.timeIntervalSinceReferenceDate))
        
        let departureIJKTest = [ 13304731552.328235626220703,
                                 -2817760335.656604290008545,
                                 0.000000000000000,
                                 1923.662144854920371,
                                 9083.032403714601969,
                                 0.000000000000000]
        XCTAssertEqual(departureIJK.position.x, departureIJKTest[0], accuracy: 1e-5)
        XCTAssertEqual(departureIJK.position.y, departureIJKTest[1], accuracy: 1e-5)
        XCTAssertEqual(departureIJK.position.z, departureIJKTest[2], accuracy: 1e-5)
        XCTAssertEqual(departureIJK.velocity.x, departureIJKTest[3], accuracy: 1e-5)
        XCTAssertEqual(departureIJK.velocity.y, departureIJKTest[4], accuracy: 1e-5)
        XCTAssertEqual(departureIJK.velocity.z, departureIJKTest[5], accuracy: 1e-5)
        
        let destinationIJKTest = [  -4190408087.8222322,
                                    -8960007201.2478027,
                                    -277584613.58439338,
                                    9863.4821821708629,
                                    -4512.7719574339681,
                                    -253.44640143248461]
        XCTAssertEqual(destinationIJK.position.x, destinationIJKTest[0], accuracy: 1e-5)
        XCTAssertEqual(destinationIJK.position.y, destinationIJKTest[1], accuracy: 1e-5)
        XCTAssertEqual(destinationIJK.position.z, destinationIJKTest[2], accuracy: 1e-5)
        XCTAssertEqual(destinationIJK.velocity.x, destinationIJKTest[3], accuracy: 1e-5)
        XCTAssertEqual(destinationIJK.velocity.y, destinationIJKTest[4], accuracy: 1e-5)
        XCTAssertEqual(destinationIJK.velocity.z, destinationIJKTest[5], accuracy: 1e-5)
    }
    
    func testLambertKerbinToEve() {
        //Year 3, day 401 at 2:09:36
        let departureTime = KSPDate(year: 3, day: 401, hour: 2, minute: 9, second: 36)
        //Year 4, day 162 at 5:45:36
        let arrivalTime = KSPDate(year: 4, day: 162, hour: 5, minute: 45, second: 36)
        
        let travelTime = arrivalTime.timeIntervalSinceReferenceDate - departureTime.timeIntervalSinceReferenceDate
        
        print(Orbit.kerbin.coe2rv(atTime: departureTime.timeIntervalSinceReferenceDate))
        print(Orbit.eve.coe2rv(atTime: arrivalTime.timeIntervalSinceReferenceDate))
        
        print(departureTime.timeIntervalSinceReferenceDate)
        print(arrivalTime.timeIntervalSinceReferenceDate)
        print(travelTime)
        
        let solver = try! LambertSolver(orbit1: .kerbin, orbit2: .eve, departureTime: departureTime, travelTime: travelTime)
        let solutions = solver.solve()
        
        let v0 = [-3084.93199034, -7812.84494769, 10.28033528]
        let vf = [-1705.22437416,  11666.6233287, -12.45954879]
            
        XCTAssertEqual(solutions[0].0.x, v0[0], accuracy: 1e-5)
        XCTAssertEqual(solutions[0].0.y, v0[1], accuracy: 1e-5)
        XCTAssertEqual(solutions[0].0.z, v0[2], accuracy: 1e-5)
        XCTAssertEqual(solutions[0].1.x, vf[0], accuracy: 1e-5)
        XCTAssertEqual(solutions[0].1.y, vf[1], accuracy: 1e-5)
        XCTAssertEqual(solutions[0].1.z, vf[2], accuracy: 1e-5)
        
        let vinf_departure = solutions[0].0 - Orbit.kerbin.cartesian(atTime: departureTime.timeIntervalSinceReferenceDate).velocity
        let vinf_arrival = solutions[0].1 - Orbit.eve.cartesian(atTime: arrivalTime.timeIntervalSinceReferenceDate).velocity
        
        let vinf_arrive = [  940.42330767,  967.38853121, -416.52301516]
        let vinf_depart = [  390.46934141,  796.65789962,   10.28033528]
        
        XCTAssertEqual(vinf_departure.x, vinf_depart[0], accuracy: 1e-5)
        XCTAssertEqual(vinf_departure.y, vinf_depart[1], accuracy: 1e-5)
        XCTAssertEqual(vinf_departure.z, vinf_depart[2], accuracy: 1e-5)
        
        XCTAssertEqual(vinf_arrival.x, vinf_arrive[0], accuracy: 1e-5)
        XCTAssertEqual(vinf_arrival.y, vinf_arrive[1], accuracy: 1e-5)
        XCTAssertEqual(vinf_arrival.z, vinf_arrive[2], accuracy: 1e-5)
        
        let kerbinParkingOrbit = Orbit(semiMajorAxis: CelestialBody.kerbin.radius + 100000, eccentricity: 0, meanAnomaly: 0, inclination: 0, LAN: 0, argumentOfPeriapsis: 0, centralBody: .kerbin)
        let eveParkingOrbit = Orbit(semiMajorAxis: CelestialBody.eve.radius + 100000, eccentricity: 0, meanAnomaly: 0, inclination: 0, LAN: 0, argumentOfPeriapsis: 0, centralBody: .eve)

        func circularToEscapeDV(from body: CelestialBody, v0: Double, vsoi: Double, relativeInclination: Double) -> Double {
            let mu = body.gravitationalParameter
            let rsoi = body.sphereOfInfluence!
            let v1 = sqrt(vsoi * vsoi + 2 * v0 * v0 - 2 * mu / rsoi)
            
            let r0 = mu / (v0 * v0)
            let e = r0 * v1 * v1 / mu - 1
            let ap = r0 * (1 + e)
            
            if relativeInclination == 0 {
                return v1 - v0
            } else {
                return sqrt(v0 * v0 + v1 * v1 - 2 * v0 * v1 * cos(relativeInclination))
            }
        }
        
        func circularToEscapeDV2(from body: CelestialBody, orbitRadius: Double, vinf: Vector3D) -> Double {
            let mu = body.gravitationalParameter
            let rsoi = body.sphereOfInfluence!
            
            let v0 = sqrt(vinf * vinf + (2 * mu / orbitRadius))
            let dv = v0 - sqrt(mu / orbitRadius)
            return dv
        }
        
        func insertionToCircularDV(to body: CelestialBody, vsoi: Double, v0: Double) -> Double {
            let mu = body.gravitationalParameter
            let rsoi = body.sphereOfInfluence!
            return sqrt(vsoi * vsoi + 2 * v0 * v0 - 2 * mu / rsoi) - v0
        }
        
        func dv(from parkingOrbit: Orbit, vinf: Vector3D) -> Double{
            let mu = parkingOrbit.centralBody.gravitationalParameter
            let soi = parkingOrbit.centralBody.sphereOfInfluence!
            
            let a_hyp = 1 / (2 / soi - (vinf * vinf) / mu)
            let vp = sqrt(mu * (2 / parkingOrbit.radius() - 1 / a_hyp))
            let dv = vp - sqrt(mu / 700000)
            return dv
        }
        
        //let dv_departure = dv(from: kerbinParkingOrbit, vinf: vinf_departure)
        let dv_departure = circularToEscapeDV2(from: .kerbin, orbitRadius: 700000, vinf: vinf_departure)
        let dv_arrival = dv(from: eveParkingOrbit, vinf: vinf_arrival)
        
        print(dv_departure)
        
        let ejectionDV = circularToEscapeDV(from: .kerbin, v0: kerbinParkingOrbit.orbitalSpeed(), vsoi: vinf_departure.magnitude, relativeInclination: Orbit.eve.inclination - Orbit.kerbin.inclination)
        let arrivalDV = insertionToCircularDV(to: .duna, vsoi: vinf_arrival.magnitude, v0: eveParkingOrbit.orbitalSpeed())
        
        print(ejectionDV + arrivalDV)
    }

}

private extension Orbit {
    func anom(atTime time: Double) -> String {
        return "anomalies('t', \(self.timeSincePeriapsis(forTime: time)), \(self.semiMajorAxis), \(self.eccentricity), mu=\(self.centralBody.gravitationalParameter)"
    }
    
    func coe(atTime time: Double) -> String {
        return "\(semiMajorAxis), \(eccentricity), \(inclination.rad2deg), \(LAN.rad2deg), \(argumentOfPeriapsis), \(trueAnomaly(atTime: time) * 180 / .pi), mu=\(centralBody.gravitationalParameter)"
    }
    
    func coe2rv(atTime time: Double) -> String {
        return "coe2rv(" + coe(atTime: time) + ")"
    }
}

//matt's library expects time sine peri
private extension Orbit {
    func timeSincePeriapsis(forTime time: Double) -> Double {
        let ttpe = self.timeToPeriapsis(atTime: time)
        return self.period - ttpe
    }
}

private extension Double {
    var deg2rad: Double {
        return self * .pi / 180
    }
    var rad2deg: Double {
        return self * 180 / .pi
    }
}
