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
* Render a *sequence* of cutting plans instead of just one. Make sure that the output image has adequate size.
* Return a table with the meaning of different annotations in the figure.
* Adjust the scale at which a cutting plan is rendered based on the *minimum* bbox side.

[About coordinates](https://juliagraphics.github.io/Luxor.jl/stable/explanation/basics/)
[How the @png macro works](https://juliagraphics.github.io/Luxor.jl/stable/tutorial/basictutorial/#What-you-need)

### Continuation

* Write function that, given a beam and a direction, it returns the plane key of the beam with the normal that points in that direction. Good for selecting a side of the beam without knowing its exact key.
* Write a function to check if a beam exists at all (if the intersections is the empty set).
* Write functions to work with *groups of beams*: translation, etc. Maybe a an abstract class. <-- SPECIFY THIS.
* Render HTML report.
