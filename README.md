# Flight Path Project

## TODO
- [ ] Implement Erosion-based Terrain Generation as described in the paper
  - [ ] Generate Base Terrain
    - [ ] Simulate 1/f noise by implement either Spectral Sythesis Method or Midpoint Displacement Method
    - [ ] Implement Voronoi diagrams to create distinct features
    - [ ] Combine Voronoi Diagram with 1/f noise.
    - [ ] Apply perturbation filter to 'crumple' straight lines from the voronoi diagram
  - [ ] Implement Erosion Algorithms
    - [ ] Thermal Erosion
    - [ ] Hydraulic Erosian
    - [ ] The Proposed Algorithm in the paper
- [ ] Render the Height Map

- [ ] Create Overhead Minimap 
- [ ] Create First-person 'plane view'

- [ ] Calculate optimal path through the terrain
- [ ] Smoothly fly through the terrain
