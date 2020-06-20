import Cocoa
import Keplerian
import PlaygroundSupport

let porkchop = Porkchop(from: .kerbin, to: .eve, departureWindow: (KSPDate(year: 2, day: 100), KSPDate(year: 2, day: 200)))
porkchop.minimumTravelTime = 100.kerbalDays.converted(to: .seconds).value
porkchop.maximumTravelTime = 200.kerbalDays.converted(to: .seconds).value
porkchop.travelTimeStep = 1.kerbalDays.converted(to: .seconds).value
porkchop.departureTimeStep = 1.kerbalDays.converted(to: .seconds).value
try! porkchop.solve()

print(porkchop.minimumSolution)
print(porkchop.minimumSolution!.arrivalDV + porkchop.minimumSolution!.departureDV)
porkchop.minimumSolution?.departureDate.timeIntervalSinceReferenceDate
KSPDate(timeIntervalSinceReferenceDate: porkchop.minimumSolution!.travelTime)
KSPDate(year: 2, day: 160).timeIntervalSinceReferenceDate

let view = PorkchopView(frame: NSRect(x: 0, y: 0, width: 500, height: 500))
view.porkchop = porkchop

PlaygroundPage.current.liveView = view
