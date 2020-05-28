//
//  VectorMath.swift
//  KSPCockpitPanel
//
//  Created by Mike Muszynski on 7/24/17.
//  Copyright © 2017 Mike Muszynski. All rights reserved.
//

import Foundation
import simd

struct Vector2D: Equatable {
    static func ==(lhs: Vector2D, rhs: Vector2D) -> Bool {
        return lhs.vector3D == rhs.vector3D
    }
    
    var x: Double
    var y: Double
    
    static var zero: Vector2D {
        return Vector2D(x: 0, y: 0)
    }
    
    var vector3D: Vector3D {
        return Vector3D(x: x, y: y, z: 0)
    }
}

public struct Vector3D: Equatable {
    public static func ==(lhs: Vector3D, rhs: Vector3D) -> Bool {
        return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z
    }
    
    var x: Double
    var y: Double
    var z: Double
    
    init(x: Double, y: Double, z: Double) {
        self.x = x
        self.y = y
        self.z = z
    }

    init (from: simd_double3) {
        self = Vector3D(x: from.x, y: from.y, z: from.z)
    }
    
    static var zero: Vector3D {
        return Vector3D(x: 0, y: 0, z: 0)
    }
    
    var norm2: Double {
        return magnitude
    }
    
    var magnitude: Double {
        return sqrt(x * x + y * y + z * z)
    }
    
    mutating func normalize() {
        let mag = self.magnitude
        guard mag > 0 else {
            self = .zero
            return
        }
        x /= mag
        y /= mag
        z /= mag
    }
    
    var normalized: Vector3D {
        var new = self
        new.normalize()
        return new
    }
    
    var simd: simd_double3 {
        return simd_double3(x: x, y: y, z: z)
    }
}

func +(left: Vector3D, right: Vector3D) -> Vector3D {
    return Vector3D(x: left.x + right.x, y: left.y + right.y, z: left.z + right.z)
}

prefix func -(input: Vector3D) -> Vector3D {
    return Vector3D(x: -input.x, y: -input.y, z: -input.z)
}

func -(left: Vector3D, right: Vector3D) -> Vector3D {
    return left + -right
}

func dotProduct(left: Vector3D, right: Vector3D) -> Double {
    return left.x * right.x + left.y * right.y + left.z * right.z
}

func *(lhs: Vector3D, rhs: Vector3D) -> Double {
    return dotProduct(left: lhs, right: rhs)
}

func *(lhs: Vector3D, rhs: Double) -> Vector3D {
    return Vector3D(x: lhs.x * rhs, y: lhs.y * rhs, z: lhs.z * rhs)
}

func *(lhs: Double, rhs: Vector3D) -> Vector3D {
    return rhs * lhs
}

func /(lhs: Vector3D, rhs: Double) -> Vector3D {
    return Vector3D(x: lhs.x / rhs, y: lhs.y / rhs, z: lhs.z / rhs)
}

func crossProduct(left: Vector3D, right: Vector3D) -> Vector3D {
    let x = left.y * right.z - left.z * right.y
    let y = left.z * right.x - left.x * right.z
    let z = left.x * right.y - left.y * right.x
    return Vector3D(x: x, y: y, z: z)
}

infix operator *&: MultiplicationPrecedence
func *&(lhs: Vector3D, rhs: Vector3D) -> Vector3D {
    return crossProduct(left: lhs, right: rhs)
}

