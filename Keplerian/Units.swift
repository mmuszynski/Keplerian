//
//  Units.swift
//  KSPTelemetry
//
//  Created by Mike Muszynski on 5/25/20.
//  Copyright © 2020 Mike Muszynski. All rights reserved.
//

import Foundation

extension Double {
    var meters: Measurement<UnitLength> {
        Measurement(value: self, unit: UnitLength.meters)
    }
    
    var m: Measurement<UnitLength> {
        meters
    }

    var kilometers: Measurement<UnitLength> {
        Measurement(value: self, unit: UnitLength.kilometers)
    }
    
    var km: Measurement<UnitLength> {
        kilometers
    }
    
    var radians: Measurement<UnitAngle> {
        Measurement(value: self, unit: UnitAngle.radians)
    }
    
    var rad: Measurement<UnitAngle> {
        radians
    }
    
    var degrees: Measurement<UnitAngle> {
        Measurement(value: self, unit: UnitAngle.degrees)
    }
    
    var deg: Measurement<UnitAngle> {
        degrees
    }
}
