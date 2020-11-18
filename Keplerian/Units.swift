//
//  Units.swift
//  KSPTelemetry
//
//  Created by Mike Muszynski on 5/25/20.
//  Copyright © 2020 Mike Muszynski. All rights reserved.
//

import Foundation

extension Double {
    public var meters: Measurement<UnitLength> {
        Measurement(value: self, unit: UnitLength.meters)
    }
    
    public var m: Measurement<UnitLength> {
        meters
    }

    public var kilometers: Measurement<UnitLength> {
        Measurement(value: self, unit: UnitLength.kilometers)
    }
    
    public var km: Measurement<UnitLength> {
        kilometers
    }
    
    public var radians: Measurement<UnitAngle> {
        Measurement(value: self, unit: UnitAngle.radians)
    }
    
    public var rad: Measurement<UnitAngle> {
        radians
    }
    
    public var degrees: Measurement<UnitAngle> {
        Measurement(value: self, unit: UnitAngle.degrees)
    }
    
    public var deg: Measurement<UnitAngle> {
        degrees
    }
}
