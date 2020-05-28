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
    
    var earliestDeparture: Double
    var minimumTravelTime: Double = 1
    var maximumTravelTime: Double = 400
    
    init(from body1: CelestialBody, to body2: CelestialBody, leavingAfter time: Double) {
        self.initialBody = body1
        self.finalBody = body2
        self.earliestDeparture = time
    }
    
    func solve() throws {
        guard let orbit1 = initialBody.orbit else { throw SetupError.invalid }
        guard let orbit2 = initialBody.orbit else { throw SetupError.invalid }
        guard orbit1.centralBody == orbit2.centralBody else { throw SetupError.invalid }
        
        for departure in stride(from: earliestDeparture, to: earliestDeparture + 220, by: 1) {
            for travelTime in stride(from: minimumTravelTime, to: maximumTravelTime, by: 1) {
                let solver = LambertSolver(position1: orbit1.cartesianPosition(atTime: departure),
                                           position2: orbit2.cartesianPosition(atTime: departure),
                                           dt: minimumTravelTime * 21600,
                                           mu: orbit1.centralBody.gravitationalParameter)
                
                solver.solve()
            }
        }
    }
    
    
}
