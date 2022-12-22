//
//  File.swift
//  
//
//  Created by Mike Muszynski on 12/15/22.
//

import Foundation

extension Orbit {
    enum Direction {
        case prograde, retrograde
        case normal, antinormal
        case radialIn, radialOut
        
        /// Returns a unit vector in the cartesian direction sepcified
        /// - Parameter orbit: The orbit used to calculate the cartesian direction
        /// - Returns: A unit vector in the direction specified
        func cartesianDirection(for orbit: Orbit) -> Vector3D {
            .zero.normalized
        }
    }
}
