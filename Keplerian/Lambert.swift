//
//  Lambert.swift
//  KSPTelemetry
//
//  Created by Mike Muszynski on 5/24/20.
//  Copyright © 2020 Mike Muszynski. All rights reserved.
//

import Foundation

typealias LambertSolution = (Vector3D, Vector3D, Double)

class LambertSolver {
    
    var mu: Double
    var pos1: Vector3D
    var pos2: Vector3D
    var dt: Double
    var maxRevs: Int = 0
    var prograde = true
    
    init(position1: Vector3D, position2: Vector3D, dt: Double, mu: Double = CelestialBody.kerbol.gravitationalParameter) {
        self.pos1 = position1
        self.pos2 = position2
        self.dt = dt
        self.mu = mu
    }
    
    func solve() ->  [LambertSolution] {
        return attempt1()
    }
    
    /*
    private func attempt3() -> [LambertSolution] {
        //Izzo's version apparently
        let tolerance = 1e-8
        let secsInDay: Double = 86400
        
        //Non-dimensionalize
        let r1 = position1.magnitude
        let r1vec = position1 / r1
        let r2vec = position2 / r1
        let V = sqrt(mu / r1)
        let T = r1 / V
        var tf = dt * secsInDay
        
        //Relevant geometry parameters (dimensionless)
        let mr2vec = r2vec.magnitude
        var dth = acos(max(-1, min(1, (r1vec * r2vec) / mr2vec)))
        
        //Decide whether to use the left or the right branch
        let leftbranch = maxRevs >= 0
        let longway = tf >= 0
        
        let m = abs(maxRevs)
        tf = abs(tf)
        let logt = log(tf)
        
        if longway {
            dth = 2 * .pi - dth
        }
        
        //derived quantities
        let c = sqrt(1 + mr2vec * mr2vec - 2 * mr2vec * cos(dth))
        let s = (1 + mr2vec + c) / 2
        let a_min = s / 2
        let Lambda = sqrt(mr2vec) * cos(dth / 2) / s
        
        //I think this is what he is going for
        let cross = crossProduct(left: r1vec, right: r2vec)
        let mcr = cross.magnitude
        let nrmunit = cross / mcr
        
        //Initial values
        var inn1: Double
        var inn2: Double
        var x1: Double
        var x2: Double
        var y1: Double
        var y2: Double
        
        if m == 0 {
            inn1 = -0.5223
            inn2 = 0.5223
            x1 = log(1 + inn1)
            x2 = log(1 + inn2)
        } else {
            if !leftbranch {
                inn1 = -0.5234
                inn2 = -0.2234
            } else {
                inn1 = +0.7234
                inn2 = +0.5234
            }
            
            x1 = tan(inn1 * .pi / 2)
            x2 = tan(inn2 * .pi / 2)
        }
        
        //initial formation
        //xx = [inn1, inn2]
        //aa = [a1, a2]
        
        let a1 = a_min / (1 - inn1 * inn1)
        let a2 = a_min / (1 - inn2 * inn2)
        
        let beta1 = 2 * asin(sqrt(s-c)/2/a1) * (longway ? 1 : -1)
        let beta2 = 2 * asin(sqrt(s-c)/2/a2) * (longway ? 1 : -1)

        let alfa1 = 2 * acos(max(-1, min(1, inn1)))
        let alfa2 = 2 * acos(max(-1, min(1, inn2)))
        
        let alfa1sin = alfa1 - sin(alfa1)
        let beta1sin = beta1 - sin(beta1)
        let yy1 = a1 * sqrt(a1) * ((alfa1sin - beta1sin) + 2 * Double.pi * Double(m))
        
        let alfa2sin = alfa2 - sin(alfa2)
        let beta2sin = beta2 - sin(beta2)
        let yy2 = a2 * sqrt(a2) * ((alfa2sin - beta2sin) + 2 * Double.pi * Double(m))
        
        //initial estimates for y
        if m == 0 {
            y1 = log(yy1) - logt
            y2 = log(yy2) - logt
        } else {
            y1 = yy1 - tf
            y2 = yy2 - tf
        }
        
        //solve for x
        var error = Double.greatestFiniteMagnitude
        var iterations = 0
        var xnew: Double = 0
        var ynew: Double = 0
        var x: Double
        var a: Double
        var alfa: Double
        var beta: Double
        var tof: Double
        
        while error > tolerance {
            iterations += 1
            xnew = (x1 * y2 - y1 * x2) / (y2 - y1)
            if m == 0 {
                x = exp(xnew) - 1
            } else {
                x = atan(xnew) * 2 / .pi
            }
            
            a = a_min / (1 - x * x)
            if x < 1 {
                beta = 2 * asin(sqrt((s-c)/2/a)) * (longway ? 1 : -1)
                alfa = 2 * acos( max(-1, min(1, x)))
            } else {
                alfa = 2 * acosh(x)
                beta = 2 * asinh(sqrt((s-c)/(-2 * a))) * (longway ? 1 : -1)
            }
            
            if a > 0 {
                tof = a * sqrt(a) * ((alfa - sin(alfa)) - (beta - sin(beta)) + 2 * .pi * Double(m))
            } else {
                tof = -a * sqrt(-a) * ((sinh(alfa) - alfa) - (sinh(beta) - beta))
            }
            
            if m == 0 {
                ynew = log(tof) - logt
            } else {
                ynew = tof - tf
            }
            
            x1 = x2
            x2 = xnew
            
            y1 = y2
            y2 = ynew
            
            error = abs(x1-xnew)
            if iterations > 15 {
                //here you need to solve with a different lambert targeter, and that code is just thrown there below
                break
            }
        }
        
        /*
         % If the Newton-Raphson scheme failed, try to solve the problem
            % with the other Lambert targeter.
            if bad
                % NOTE: use the original, UN-normalized quantities
                [V1, V2, extremal_distances, exitflag] = ...
                    lambert_LancasterBlanchard(r1vec*r1, r2vec*r1, longway*tf*T, leftbranch*m, muC);
                return
            end
         */
        
        if m == 0 {
            x = exp(xnew) - 1
        } else {
            x = atan(xnew) * 2 / .pi
        }
        
        /*
         The solution has been evaluated in terms of log(x+1) or tan(x*pi/2), we
         now need the conic. As for transfer angles near to pi the Lagrange-
         coefficients technique goes singular (dg approaches a zero/zero that is
         numerically bad) we here use a different technique for those cases. When
         the transfer angle is exactly equal to pi, then the ih unit vector is not
         determined. The remaining equations, though, are still valid.
         */
        
        //semimajor axis
        a = a_min / (1 - x * x)

        //psi
        var psi, eta2, eta: Double
        if x < 1 {
            beta = 2 * asin(sqrt((s-c)/2/a)) * (longway ? 1 : -1)
            alfa = 2 * acos( max(-1, min(1, x)))
            psi  = (alfa-beta)/2;
            eta2 = 2*a*sin(psi)*sin(psi)/s;
            eta  = sqrt(eta2);
        } else {
            beta = 2*asinh(sqrt((c-s)/2/a)) * (longway ? 1 : -1)
            alfa = 2*acosh(x);
            psi  = (alfa-beta)/2;
            eta2 = -2*a*sinh(psi)*sinh(psi)/s;
            eta  = sqrt(eta2);
        }
        
        // unit of the normalized normal vector
        let ih = nrmunit * (longway ? 1 : -1)
        
        // unit vector for normalized [r2vec]
        let r2n = r2vec/mr2vec;
        
        let crsprd1 = crossProduct(left: ih, right: r1vec)
        let crsprd2 = crossProduct(left: ih, right: r2n)
        
        // radial and tangential directions for departure velocity
        let Vr1 = 1/eta/sqrt(a_min) * (2*Lambda*a_min - Lambda - x*eta);
        let Vt1 = sqrt(mr2vec/a_min/eta2 * sin(dth/2) * sin(dth/2));
        
        // radial and tangential directions for arrival velocity
        let Vt2 = Vt1/mr2vec;
        let Vr2 = (Vt1 - Vt2)/tan(dth/2) - Vr1;

        // terminal velocities
        let V1 = (Vr1*r1vec + Vt1*crsprd1)*V;
        let V2 = (Vr2*r2n + Vt2*crsprd2)*V;
        
        return [(V1, V2, Double(m))]
    }

    private func attempt2() ->  [(Vector3D, Vector3D, Double)] {
        var solutions = [(Vector3D, Vector3D, Double)]()
        
        let r0 = position1.magnitude
        let rf = position2.magnitude
        
        let cosdnu = (position1 * position2) / (r0 * rf)
        let dnu = acos(cosdnu)
        
        let a = sqrt(r0 * rf * (1 + cosdnu))
        
        if a == 0 && dnu == 0 { return [] }
        
        var c2 = 1.0 / 2.0
        var c3 = 1.0 / 6.0
        
        var xupper = 4.0 * .pi * .pi
        var xlower = -4.0 * .pi
        var x = xupper + xlower / 2
        
        var y: Double = 0
        
        var calculatedDt = 0.0
        while abs(calculatedDt - dt) > 1e-6 {
            while r0 * rf + a * (x * c3 - 1) / sqrt(c2) < 0 {
                x += 0.1
            }
            
            y = r0 * rf + a * (x * c3 - 1) / sqrt(c2)
            let chi = sqrt(y / c2)
            calculatedDt = (chi * chi * chi * c3 + a * sqrt(y)) / sqrt(mu)
            
            if calculatedDt <= dt {
                xlower = x
            } else {
                xupper = x
            }
            
            x = (xupper + xlower) / 2
            
            if x > 1e-6 {
                c2 = (1.0 - cos(sqrt(x))) / x
                c3 = (sqrt(x) - sin(sqrt(x))) / sqrt(x * x * x)
            } else if x < -1e-6 {
                c2 = (1.0 - cosh(sqrt(-x))) / x
                c3 = (sinh(sqrt(-x)) - sqrt(-x)) / sqrt(-x * -x * -x)
            } else {
                c2 = 1/2
                c3 = 1/6
            }
        }
        
        let f = 1 - y / r0
        let gdot = 1 - y / rf
        let g = a * sqrt(y / mu)
        
        let v0 = (position2 - f * position1) * (1 / g)
        let vf = (((gdot * position2) - position1)) * (1 / g)
        
        return [(v0, vf, 0)]
    }
     */
    
    private func attempt1() ->  [(Vector3D, Vector3D, Double)] {
        func acoth(_ x: Double) -> Double {
            0.5 * log((x + 1) / (x - 1))
        }
        
        func acot(_ x: Double) -> Double {
            0.5 * Double.pi - atan(x)
        }
        
        func relativeError(_ a: Double, _ b: Double) -> Double {
            return abs(1.0 - a / b)
        }
        
        var solutions = [(Vector3D, Vector3D, Double)]()
        
        let r1 = pos1.magnitude
        let r2 = pos2.magnitude
        
        let deltaPos = pos2 - pos1
        let c = deltaPos.magnitude
        let m = r1 + r2 + c
        let n = r1 + r2 - c
        
        var transferAngle = acos((pos1 * pos2) / (r1 * r2))
        if (pos1.x * pos2.y - pos1.y * pos2.x) * (prograde ? 1 : -1) < 0 {
            transferAngle = Double.pi * 2.0 - transferAngle
        }
        
        var angleParameter = sqrt(n / m)
        if transferAngle > .pi {
            angleParameter = -angleParameter
        }
        
        let normalizedTime = 4 * dt / sqrt(mu / (m * m * m))
        let parabolicNormalizedTime = 2 / 3 * (1 - angleParameter * angleParameter * angleParameter)
        
        let sqrtMu = sqrt(mu)
        let invSqrtM = 1 / sqrt(m)
        let invSqrtN = 1 / sqrt(n)
        
        func solution(x: Double, y: Double, N: Double) -> (Vector3D, Vector3D, Double) {
            let vc = sqrtMu * (y * invSqrtN + x * invSqrtM)
            let vr = sqrtMu * (y * invSqrtN - x * invSqrtM)
            let ec = deltaPos * (vc / c)
            let v1 = ec + (pos1 * (vr / r1))
            let v2 = ec - (pos2 * (vr / r2))
            
            let Nout = N * Double.pi * 2 + transferAngle
            return (v1, v2, Nout)
        }
        
        func fy(_ x: Double) -> Double {
            sqrt(1 - angleParameter * angleParameter * (1 - x * x)) * (angleParameter < 0 ? -1 : 1)
        }
        
        func ftau(_ x: Double) -> Double {
            if x == 1.0 { return parabolicNormalizedTime - normalizedTime }
            
            let y = fy(x)
            if x > 1 {
                let g = sqrt(x * x - 1)
                let h = sqrt(y * y - 1)
                return (-acoth(x / g) + acoth(y / h) + x * g - y * h) / (g * g * g) - normalizedTime
            } else {
                let g = sqrt(1 - x * x)
                let h = sqrt(1 - y * y)
                return (acot(x / g) - atan(h / y) - x * g + y * h) / (g * g * g) - normalizedTime
            }
        }
        
        //Solution spaces
        if relativeError(normalizedTime, parabolicNormalizedTime) < 1e-6 {
            //Unique, parabolic solution
            
            let x = 1.0
            let y = Double(angleParameter < 0 ? -1 : 1)
            
            solutions.append(solution(x: x, y: y, N: 0))
            
        } else if normalizedTime < parabolicNormalizedTime {
            //Unique, hyperbolic solution
            
            var x1 = 1.0
            var x2 = 1.0
            
            while !(ftau(x2) < 0.0) {
                x1 = x2
                x2 *= 2.0
            }
            
            let x = convergingSolver(x1...x2, tolerance: 1e-12, ftau)
            solutions.append(solution(x: x, y: fy(x), N: 0))
            
        } else {
            //Elliptical solutions (potentially multiple)
            
            //let maxRevs = min(self.maxRevs, Int(floor(normalizedTime / .pi)))
            var minimumEnergyNormalizedTime = acos(angleParameter) + angleParameter * sqrt(1 - angleParameter * angleParameter)
            
            for N in 0...maxRevs {
                if N > 0 && N == maxRevs {
                    func phix(_ x: Double) -> Double {
                        let g = sqrt(1 - x * x)
                        return acot(x / g) - (2 + x * x) * g / (3 * x)
                    }
                    
                    func phiy(_ y: Double) -> Double {
                        let h = sqrt(1 - y * y)
                        return atan(h / y) - (2 + y * y) * h / (3 * y)
                    }
                    
                    var xMT: Double
                    var minimumNormalizedTime: Double
                    //Find the minimum normalized time an N revolution trajectory will take
                    if angleParameter == 1 {
                        xMT = 0
                        minimumNormalizedTime = minimumEnergyNormalizedTime
                    } else if angleParameter == 0 {
                        xMT = convergingSolver(0...1, tolerance: 1e-12) { phix($0) + Double(N) * .pi }
                        minimumNormalizedTime = 2 / (3 * xMT)
                    } else {
                        xMT = convergingSolver(0...1, tolerance: 1e-12) { phix($0) - phiy(fy($0)) + Double(N) * .pi }
                        let fyxMT = fy(xMT)
                        let absfyxMT = abs(fyxMT)
                        
                        minimumNormalizedTime = 2 / 3 * (1 / xMT - angleParameter * angleParameter * angleParameter) / absfyxMT
                    }
                    
                    if relativeError(normalizedTime, minimumNormalizedTime) < 1e-12 {
                        solutions.append(solution(x: xMT, y: fy(xMT), N: n+1 * .pi * 2 - transferAngle))
                        break
                    } else if normalizedTime < minimumNormalizedTime {
                        break
                    } else if normalizedTime < minimumEnergyNormalizedTime {
                        let x = convergingSolver(0...xMT, tolerance: 1e-12, ftau)
                        solutions.append(solution(x: x, y: fy(x), N: Double(N)))
                        let x2 = convergingSolver(xMT...1, tolerance: 1e-12, ftau)
                        solutions.append(solution(x: x2, y: fy(x2), N: Double(N)))
                        break
                    }
                }
                
                if relativeError(normalizedTime, minimumEnergyNormalizedTime) < 1e-6 {
                    solutions.append(solution(x: 0, y: fy(0), N: Double(N)))
                    if N > 0 {
                        let x = convergingSolver(1e-6...1.0, tolerance: 1e-12, ftau)
                        solutions.append(solution(x: x, y: fy(x), N: Double(N)))
                    }
                } else {
                    if N > 0 || normalizedTime > minimumEnergyNormalizedTime {
                        let x = convergingSolver(-1...0, tolerance: 1e-12, ftau)
                        solutions.append(solution(x: x, y: fy(x), N: Double(N)))
                    }
                    
                    if N > 0 || normalizedTime < minimumEnergyNormalizedTime {
                        let x = convergingSolver(0...1.0, tolerance: 1e-12, ftau)
                        solutions.append(solution(x: x, y: fy(x), N: Double(N)))
                    }
                
                    minimumEnergyNormalizedTime += .pi
                }
            }
        }
        
        return solutions
    }
    
}

func convergingSolver(_ range: ClosedRange<Double>, tolerance: Double, _ function: (Double)->Double) -> Double {
    if range.upperBound - range.lowerBound < tolerance { return range.upperBound + range.lowerBound / 2 }
    
    let distance = range.upperBound - range.lowerBound
    let lower = range.lowerBound...range.lowerBound + distance / 2
    let upper = lower.upperBound...range.upperBound
    
    if function(lower.upperBound) < 0 {
        return convergingSolver(upper, tolerance: tolerance, function)
    } else {
        return convergingSolver(lower, tolerance: tolerance, function)
    }
    
}
