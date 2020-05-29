//
//  Porkchop.swift
//  Keplerian
//
//  Created by Mike Muszynski on 5/27/20.
//  Copyright © 2020 Mike Muszynski. All rights reserved.
//

import Foundation

class Porkchop {
    
    enum SetupError: Error {
        case invalid
    }
    
    let initialBody: CelestialBody
    let finalBody: CelestialBody
    
    let departureBody = CelestialBody.kerbin
    
    var earliestDeparture: Double
    var latestDeparture: Double
    var minimumTravelTime: Double = 150
    var maximumTravelTime: Double = 450
    
    var solutionSpace: [(departureTime: Double, travelTime: Double, dV: Double)] = []
    var minimumSolution: (departureTime: Double, travelTime: Double, dV: Double)? {
        return solutionSpace.sorted {
            $0.dV < $1.dV
        }.first
    }
    
    init(from body1: CelestialBody, to body2: CelestialBody, departureWindow: (Double, Double)) {
        self.initialBody = body1
        self.finalBody = body2
        self.earliestDeparture = departureWindow.0
        self.latestDeparture = departureWindow.1
    }
    
    func solve() throws {
        solutionSpace = []
        guard let orbit1 = initialBody.orbit else { throw SetupError.invalid }
        guard let orbit2 = finalBody.orbit else { throw SetupError.invalid }
        guard orbit1.centralBody == orbit2.centralBody else { throw SetupError.invalid }
        
        for departure in stride(from: earliestDeparture, to: latestDeparture, by: 1) {
            for travelTime in stride(from: minimumTravelTime, to: maximumTravelTime, by: 1) {
                let departure_sec = departure * 21600
                let travelTime_sec = travelTime * 21600
                
                let departure_ijk = orbit1.cartesian(atTime: departure_sec)
                let arrival_ijk = orbit2.cartesian(atTime: departure_sec + travelTime_sec)
                                
                let solver = LambertSolver(position1: departure_ijk.position,
                                           position2: arrival_ijk.position,
                                           dt: travelTime_sec,
                                           mu: orbit1.centralBody.gravitationalParameter)
                
                let (departureVelocity, arrivalVelocity, _) = solver.solve()[0]
                let vinf_departure = departureVelocity - departure_ijk.velocity
                let vinf_arrival = arrivalVelocity - arrival_ijk.velocity
                
                let totalDV = vinf_departure.magnitude + vinf_arrival.magnitude
                let solution = (departureTime: departure, travelTime: travelTime, dV: totalDV)
                self.solutionSpace.append(solution)
            }
        }
    }
    
    
}
