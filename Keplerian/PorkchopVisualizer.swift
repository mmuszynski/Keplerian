//
//  PorkchopVisualizer.swift
//  Keplerian
//
//  Created by Mike Muszynski on 6/5/20.
//  Copyright © 2020 Mike Muszynski. All rights reserved.
//

#if os(macOS)
import Cocoa
import Accelerate

public class PorkchopView: NSView {
    
    public var porkchop: Porkchop?
    var gradient: NSGradient {
        NSGradient(colorsAndLocations: (NSColor(calibratedRed: 0, green: 0, blue: 0.5, alpha: 1.0), 0),
                   (NSColor(calibratedRed: 0.5, green: 0.5, blue: 1, alpha: 1.0), 0.005),
                   (.cyan, 0.05),
                   (.green, 0.3),
                   (.yellow, 0.6),
                   (.orange, 0.8),
                   (.red, 0.9))!
    }
    var solutionSpace = [Porkchop.Solution]() {
        didSet {
            self.setNeedsDisplay(self.bounds)
        }
    }
    
    public override func draw(_ dirtyRect: NSRect) {
        drawPorkchop()
    }
    
    func drawPorkchop(accelerated: Bool = false) {
        guard let porkchop = porkchop else { return }
        let departureRange = porkchop.departureTimeRange
        let xScale = departureRange.upperBound.timeIntervalSinceReferenceDate - departureRange.lowerBound.timeIntervalSinceReferenceDate
        let yScale = porkchop.travelTimeRange.upperBound - porkchop.travelTimeRange.lowerBound

        let xTranslate = departureRange.lowerBound.timeIntervalSinceReferenceDate
        let yTranslate = porkchop.travelTimeRange.lowerBound
        
        let translate = CGAffineTransform(translationX: CGFloat(-xTranslate),
                                          y: CGFloat(-yTranslate))
        let scale = CGAffineTransform(scaleX: self.bounds.width / CGFloat(xScale),
                                      y: self.bounds.height / CGFloat(yScale))
        
        if accelerated {
//            let xNum = porkchop.solutionSize.dateSteps
//            let yNum = porkchop.solutionSize.timeSteps
//            
//            let xRange = porkchop.departureTimeRange.lowerBound.timeIntervalSinceReferenceDate...porkchop.departureTimeRange.upperBound.timeIntervalSinceReferenceDate
//            
//            var xValues = [Double]()
//            var yValues = [Double]()
//            var zValues = [Double]()
//            porkchop.solutionSpace.forEach {
//                xValues.append($0.departureDate.timeIntervalSinceReferenceDate)
//                yValues.append($0.travelTime)
//                zValues.append($0.dV)
//            }
//            
//            let translatedXValues = vDSP.add(-xTranslate, xValues)
//            let translatedYValues = vDSP.add(-yTranslate, yValues)
//            
//            let scaledXValues = vDSP.multiply(xScale, translatedXValues)
//            let scaledYValues = vDSP.multiply(yScale, translatedYValues)
//            
//            let translatedZValues = vDSP.add(-porkchop.deltaVRange!.lowerBound, zValues)
//            let scaledZValues = vDSP.divide(translatedZValues, porkchop.deltaVRange!.upperBound - porkchop.deltaVRange!.lowerBound)
//
//            for i in 0..<xValues.count {
//                let (x, y, z) = (CGFloat(scaledXValues[i]), CGFloat(scaledYValues[i]), CGFloat(scaledZValues[i]))
//                
//                let color = gradient.interpolatedColor(atLocation: z)
//                let box = NSRect(x: x,
//                                 y: y,
//                                 width: CGFloat(xScale),
//                                 height: CGFloat(yScale))
//                
//                color.setFill()
//                box.fill()
//            }
        } else {
            func dvToFraction(_ dv: Double) -> CGFloat {
                guard let range = porkchop.deltaVRange else { return 1.0 }
                return CGFloat((dv - range.lowerBound) / (range.upperBound - range.lowerBound))
            }
            
            for solution in porkchop.solutionSpace {
                var box = NSRect(x: solution.departureDate.timeIntervalSinceReferenceDate,
                                 y: solution.travelTime,
                                 width: porkchop.departureTimeStep,
                                 height: porkchop.travelTimeStep)
                box = box.applying(translate)
                box = box.applying(scale)
                
                gradient.interpolatedColor(atLocation: dvToFraction(solution.departureDV + solution.arrivalDV)).setFill()
                box.fill()
            }
        }
    }
}
#endif
