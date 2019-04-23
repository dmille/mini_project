# Flight Path Project
Project for Geometric Modelling at Saarland University

## TODO

### Camera Views
- [x] Render the Height Map

- [x] Create Overhead Minimap 
- [x] Create First-person 'plane view'

### Shortest Path
- [x] Calculate optimal path through the terrain
- [x] Smoothly fly through the terrain

### Terrain Generation
- [x] Implement Erosion-based Terrain Generation as described in the paper
  - [x] Generate Base Terrain
    - [ ] Simulate 1/f noise by implement either Spectral Sythesis Method or Midpoint Displacement Method
    - [ ] Implement Voronoi diagrams to create distinct features
    - [ ] Combine Voronoi Diagram with 1/f noise.
    - [ ] Apply perturbation filter to 'crumple' straight lines from the voronoi diagram
  - [ ] Implement Erosion Algorithms
    - [ ] Thermal Erosion
    - [ ] Hydraulic Erosian
    - [ ] The Proposed Algorithm in the paper
