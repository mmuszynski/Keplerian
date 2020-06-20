//
//  Porkchop.swift
//  Keplerian
//
//  Created by Mike Muszynski on 5/27/20.
//  Copyright © 2020 Mike Muszynski. All rights reserved.
//

import Foundation

public class Porkchop {
    public typealias Solution = (departureDate: KSPDate, travelTime: Double, departureDV: Double, arrivalDV: Double)
    
    enum SetupError: Error {
        case invalid
    }
    
    let initialBody: CelestialBody
    let finalBody: CelestialBody
    
    let initialParkingOrbit: Double = 100000
    let finalParkingOrbit: Double = 1000000
        
    var earliestDeparture: KSPDate
    var latestDeparture: KSPDate
    public var departureTimeStep: TimeInterval = 1.kerbalDay.converted(to: .seconds).value
    
    public var minimumTravelTime: TimeInterval = 100.kerbalDays.converted(to: .seconds).value
    public var maximumTravelTime: TimeInterval = 450.kerbalDays.converted(to: .seconds).value
    public var travelTimeStep: TimeInterval = 1.kerbalDay.converted(to: .seconds).value
    
    public var solutionSize: (dateSteps: Int, timeSteps: Int) = (0,0)
    
    public var solutionSpace: [Solution] = []
    var sortedSolutionSpace: [Solution] = []
    public var minimumSolution: Solution? { sortedSolutionSpace.first }
    
    public var deltaVRange: ClosedRange<Double>? {
        guard let lowest = sortedSolutionSpace.first, let highest = sortedSolutionSpace.last else { return nil }
        return (lowest.departureDV + lowest.arrivalDV)...(highest.departureDV + highest.arrivalDV)
    }
    
    public var travelTimeRange: ClosedRange<Double> {
        return minimumTravelTime...maximumTravelTime
    }
    
    public var departureTimeRange: ClosedRange<KSPDate> {
        return earliestDeparture...latestDeparture
    }
    
    public init(from body1: CelestialBody, to body2: CelestialBody, departureWindow: (KSPDate, KSPDate)) {
        self.initialBody = body1
        self.finalBody = body2
        self.earliestDeparture = departureWindow.0
        self.latestDeparture = departureWindow.1
    }
    
    public func solve() throws {
        solutionSpace = []
        solutionSize = (0,0)
        
        guard let orbit1 = initialBody.orbit else { throw SetupError.invalid }
        guard let orbit2 = finalBody.orbit else { throw SetupError.invalid }
        guard orbit1.centralBody == orbit2.centralBody else { throw SetupError.invalid }
        
        let departureRange = stride(from: earliestDeparture.timeIntervalSinceReferenceDate,
                                    to: latestDeparture.timeIntervalSinceReferenceDate,
                                    by: departureTimeStep)
        
        let travelTimeRange = stride(from: minimumTravelTime,
                                     to: maximumTravelTime,
                                     by: travelTimeStep)
        
        for departureTime in departureRange {
            solutionSize.dateSteps += 1
            solutionSize.timeSteps = 0
            
            for travelTime in travelTimeRange {
                solutionSize.timeSteps += 1
                
                let departure_ijk = orbit1.cartesian(atTime: departureTime)
                let arrival_ijk = orbit2.cartesian(atTime: departureTime + travelTime)
                                
                let solver = LambertSolver(position1: departure_ijk.position,
                                           position2: arrival_ijk.position,
                                           dt: travelTime,
                                           mu: orbit1.centralBody.gravitationalParameter)
                
                let (departureVelocity, arrivalVelocity, _) = solver.solve()[0]
                let vinf_departure = departureVelocity - departure_ijk.velocity
                let vinf_arrival = arrivalVelocity - arrival_ijk.velocity
                                
                func dv(fromCircularOrbitOfAltitude altitude: Double, around body: CelestialBody, vinf: Vector3D) -> Double{
                    let radius = altitude + body.radius
                    let mu = body.gravitationalParameter
                    let soi = body.sphereOfInfluence!
                    
                    let a_hyp = 1 / (2 / soi - (vinf * vinf) / mu)
                    let vp = sqrt(mu * (2 / radius - 1 / a_hyp))
                    let dv = vp - sqrt(mu / radius)
                    return dv
                }
                                                
                func circularToEscapeDV(from body: CelestialBody, v0: Double, vsoi: Double, relativeInclination: Double) -> Double {
                    let mu = body.gravitationalParameter
                    let rsoi = body.sphereOfInfluence!
                    let v1 = sqrt(vsoi * vsoi + 2 * v0 * v0 - 2 * mu / rsoi)
                    
                    let r0 = mu / (v0 * v0)
                    let e = r0 * v1 * v1 / mu - 1
                    let _ = r0 * (1 + e)
                    
                    if relativeInclination == 0 {
                        return v1 - v0
                    } else {
                        return sqrt(v0 * v0 + v1 * v1 - 2 * v0 * v1 * cos(relativeInclination))
                    }
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
                
                func dvFromBraeuning(vinf: Vector3D, from body: CelestialBody, orbitAltitude altitude: Double) -> Double {
                    let mu = body.gravitationalParameter
                    let r = body.radius + altitude
                    
                    let v0 = sqrt(vinf * vinf + 2 * mu / r)
                    let orbitalV = sqrt(mu / (r))
                    return v0 - orbitalV
                }
                
                let parkingOrbit1 = Orbit(semiMajorAxis: initialBody.radius + 100000,
                                          eccentricity: 0,
                                          meanAnomaly: 0,
                                          inclination: 0,
                                          LAN: 0,
                                          argumentOfPeriapsis: 0,
                                          centralBody: initialBody)
                
                let parkingOrbit2 = Orbit(semiMajorAxis: finalBody.radius + 100000,
                                          eccentricity: 0,
                                          meanAnomaly: 0,
                                          inclination: 0,
                                          LAN: 0,
                                          argumentOfPeriapsis: 0,
                                          centralBody: finalBody)
                
                let ejectionDV = circularToEscapeDV(from: initialBody,
                                                    v0: parkingOrbit1.orbitalSpeed(),
                                                    vsoi: vinf_departure.magnitude,
                                                    relativeInclination: finalBody.orbit!.inclination - initialBody.orbit!.inclination)
//                let ejectionDV = dvFromBraeuning(vinf: vinf_departure, from: .kerbin, orbitAltitude: 100000)
                let arrivalDV = insertionToCircularDV(to: finalBody,
                                                      vsoi: vinf_arrival.magnitude,
                                                      v0: parkingOrbit2.orbitalSpeed())
                                
                let solution = (departureDate: KSPDate(timeIntervalSinceReferenceDate: departureTime), travelTime: travelTime, departureDV: ejectionDV, arrivalDV: arrivalDV)
                self.solutionSpace.append(solution)
            }
        }
        
        self.sortedSolutionSpace = solutionSpace.sorted {
            ($0.departureDV + $0.arrivalDV) < ($1.departureDV + $1.arrivalDV)
        }
    }
    
    
}
