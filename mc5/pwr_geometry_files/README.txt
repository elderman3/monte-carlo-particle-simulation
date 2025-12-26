PWR geometry + materials pack (for the provided geometry/material parser)

Units
- Geometry dimensions are in cm (consistent with number densities in atoms/(barn·cm) and microscopic cross sections in barns).

Physical template
- Generic PWR “mini-core” for testing: 3×3 assemblies, each 17×17 pins, active height 366.0 cm.
- Guide tube pattern: indices {2,5,8,11,14}×{2,5,8,11,14} (25 positions): center is an instrument tube (modeled as water-filled guide tube), remaining 24 are guide tubes.
- One control assembly in core center. Corner assemblies are reflector (all water pins). Remaining are fuel assemblies.
- Simple reactor vessel: cylindrical steel wall (inner radius 200.0 cm, outer radius 220.0 cm), height 450.0 cm.
- Coolant region: cylinder inside vessel minus the core cuboid volume.

Versions
1) Normal (logical grouping)
   - vessel → core → assemblies → pins
   - Assemblies directly contain 289 pin subuniverses each.
   - Top-level file: geometry/geometry_normal.txt (also copied to geometry/geometry.txt)

2) Oct8 (space subdivision, ≤8 children/geometry entries per universe)
   - Same physical layout, but assemblies and core are subdivided with a quadtree.
   - Every universe in geometry/oct8 satisfies (nSubUniverses + nGeometries) ≤ 8.
   - Top-level file: geometry/geometry_oct8.txt

3) Aintree2 (binary space subdivision, ≤2 children/geometry entries per universe)
   - Same physical layout, but:
     - vessel is split into two levels (outer steel wall + inner coolant/core),
     - core and assemblies are binary-split until leaves contain ≤2 items,
     - pin cells are decomposed into nested cylindrical sub-universes so each universe also satisfies the ≤2 rule.
   - Every universe in geometry/aintree2 satisfies (nSubUniverses + nGeometries) ≤ 2.
   - Top-level file: geometry/geometry_aintree2.txt

Switching between versions
- The code reads geometry/geometry.txt by default.
- To test another version, replace geometry/geometry.txt with one of:
  - geometry/geometry_normal.txt
  - geometry/geometry_oct8.txt
  - geometry/geometry_aintree2.txt

Materials
- material/water_borated.txt  : H1, O16, B10, B11 (1000 ppm B by weight, 0.70 g/cc)
- material/fuel_uo2_4p.txt    : U235/U238/O16 (4% enrichment, 10.5 g/cc UO2)
- material/zircaloy.txt       : Zr90 (6.55 g/cc)
- material/helium_gap.txt     : He4 (very low number density)
- material/b4c.txt            : B10/B11/C12 (2.52 g/cc)
- material/steel.txt          : Fe56 only (simplified steel)

Data
- data/*.dat are simplified parametric cross sections (not full ENDF). They are intended to validate parsing and geometry traversal, and to provide “physically plausible” reaction channels (scatter/capture/fission) for performance testing.
