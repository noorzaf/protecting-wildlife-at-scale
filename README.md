# Protecting Wildlife at Scale

**A mathematical model for deploying rangers, drones and cameras across Etosha National Park, Namibia.**

📄 **[Read the full paper (PDF)](paper/Protecting_Wildlife_at_Scale.pdf)**

---

## The problem

Etosha National Park covers 22,935 km² of savanna, salt pan and grassland, but has only 295 Ministry staff to protect it from poaching, wildfire and wildlife straying outside the park. The core challenge isn't just too few resources. It's where to put the ones you have.

## The approach

1. **Divide the park into zones.** A 20 × 25 grid (500 zones, about 46 km² each). Each zone gets an importance weight from four factors:

   | Factor | Weight |
   |---|---|
   | Species risk (e.g. black rhino habitat) | 50% |
   | Accessibility to poachers | 30% |
   | Wildlife movement probability | 15% |
   | Environmental hazard | 5% |

2. **Place patrols to cover the most important ground.** An integer programming model assigns 100 active rangers per shift to maximise the total risk-weighted area covered, concentrating effort on high-priority zones instead of spreading it evenly.

3. **Model wildfire as the top threat.** In 2025 a single fire burned 38% of Etosha in one week. Fire ignitions are modelled as a seasonal Poisson process calibrated on 24 years of MODIS satellite data, and fire spread uses Rothermel's surface fire spread model. A fleet of 12–18 thermal-imaging drones flying from six solar-powered docking stations at existing rest camps cuts the need for ground fire patrols by about 68–70%.

4. **Warn rangers before animals leave the park.** A wildlife exit model turns each animal's GPS distance to the boundary and speed towards it into a risk score, and triggers an alert when the score crosses a threshold.

5. **Adapt by season.** In the wet season floods make the Etosha Pan impassable by vehicle, so drones take over coverage there. In the dry season wildlife and poaching both concentrate around the park's 86 waterholes, where fixed cameras reduce the ranger hours needed.

## Key results

- **50 deployable rangers** reach about **87%** of maximum risk-weighted coverage; **80 rangers** reach full coverage.
- Only about **53%** of the 295 staff are patrol-ready at any time, so this sets a practical operating minimum.
- Because conservation funding is partly donor-dependent, the model recommends **one-off technology investment** (drones and waterhole cameras) over a large, recurring expansion of staff.
- The framework transfers to other parks by changing only the inputs. The paper applies it to **Kui Buri National Park (Thailand)** and **Yellowstone (USA)**.

## The code

[`model/etosha_model.cpp`](model/etosha_model.cpp) implements every part of the model and checks each one against its expected range: the zone grid, reachability, greedy patrol optimisation, the seasonal Poisson fire model, the Fire Risk Index, the wildlife exit alerts and the dry/wet season comparison.

```bash
cd model
g++ -std=c++17 -O2 etosha_model.cpp -o etosha_model
./etosha_model
```

Output:

```
============================================================
Zone grid created (500 zones)             [PASS]   500 zones
All h_j in [0, 1]                         [PASS]   bounds satisfied
Reachability matrix built                 [PASS]   avg zones per patrol = 8
Patrol model converged                    [PASS]   500/500 zones covered
Weighted coverage >= 80%                  [PASS]   100.00%
Sum x_i = 100 (resource constraint)       [PASS]   deployed = 100
Annual fire probability computed          [PASS]   0.8244
Sep-Oct fire probability ~32-40%          [PASS]   38.97%
FRI computed for all zones                [PASS]   mean = 0.5861
High-risk zones (FRI > 0.7) identified    [PASS]   32 zones
RI and AS in valid range [0, 1]           [PASS]   all 5 test animals verified
Alert triggers for D <= 200 m, high Vn    [PASS]   Lion D: AS = 0.4938 (alert issued)
Dry season full coverage achieved         [PASS]   100.00%
Wet season coverage maintained > 85%      [PASS]   100.00%
============================================================
Total: 14  Passed: 14  Failed: 0
All model components verified successfully.
```

The zone values in the code are synthetic but reproducible (seeded), built to mirror Etosha's layout: waterholes, park borders and the central pan. They test that the model behaves correctly, not that they match field data exactly.

## Repository layout

```
├── paper/Protecting_Wildlife_at_Scale.pdf   full paper (36 pages)
├── model/etosha_model.cpp                   C++ implementation and verification
└── README.md
```
