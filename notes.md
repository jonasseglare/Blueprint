# Notes/logbook

## Log 2022-04-18: Next steps

Task breakdown, cutting plans:

* ~~Add an AnnotationData class with extra information for each annotation, e.g. drill radius~~
* ~~Fix the cutting plan so that the beam dir points in the y-direction~~
* ~~Given a list of cutting plans and a beam length, optimize their order. Have a flag whether they are flippable or not.
  ALGORITHM: Start with longest cutting plan. Loop over remaining ones and pick the greatest next one that fits inside.
  Try flipping it if that makes it overall shorter.~~
* ~~Given a list of simple cutting plans, make a function pack that translates them so that they don't intersect.~~
* ~~Add flag whether a cutting plan can be mirrored or not.~~
* ~~Compress the graphics to be rendered~~
* ~~Fit the layout inside the viewable area and render.~~
  - ~~Collapse annotations that are close to each other.~~
  - ~~Generate unique annotation labels, as needed.~~
  - ~~Render the beam label.~~
* ~~A BeamGroup <: AbstractBeamGroup: A set of beams (maybe just a vector, or something that can be expanded to a vector)~~
  ~~that are rendered to either~~
  - ~~Plans~~
  - ~~A 3d model.~~
  ~~Every beam gets a unique number (index in the expanded vector). This is assigned at the beginning.~~
  ~~When rendering plans, we~~
  1. ~~Group beams by similar characteristics. (same face and beam dims)~~
  2. ~~For each group, optimize layout.~~
  3. ~~Render the plans and the tables.~~
* ~~Label every beam in a plan~~
* ~~Label every annotation in a plan~~
* ~~Render the beam boundaries themselves~~ <-- Skips this, just ads clutter.
* ~~Given a beamgroup for the entire design, generate beam cutting plans for common beams~~
* ~~Write a function to check if a beam exists at all (if the intersections is the empty set).~~
* ~~Return a table with the meaning of different annotations in the figure.~~
* ~~Render full report in HTML~~
* ~~Render full report in Markdown~~
* ~~Render a 3d model for an AbstractComponent in STL format (for online viewing)~~
* ~~Line projection of polyhedron faces in a small depth span, in order to make diagrams~~

[About coordinates](https://juliagraphics.github.io/Luxor.jl/stable/explanation/basics/)
[How the @png macro works](https://juliagraphics.github.io/Luxor.jl/stable/tutorial/basictutorial/#What-you-need)

### Cad formats

* [STL](https://en.wikipedia.org/wiki/STL_(file_format)):
  - Human readable
  - Triangles
* [STEP](https://en.wikipedia.org/wiki/ISO_10303-21)
  - Human readable
  - Instance names
  - Many entities


### Continuation

More features to make it useful.

* ~~Set memberships for beams, components, etc: Render one mesh for each set.~~
* Beam array
* Drilling routines for groups
* Miter cut
* Write function that, given a beam and a direction, it returns the plane key of the beam with the normal that points in that direction. Good for selecting a side of the beam without knowing its exact key.
* ~~Write functions to work with groups of beams: translation, etc. Maybe a an abstract class.~~
