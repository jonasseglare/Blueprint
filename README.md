# Blueprint

## How to work with the project (personal notes)

### Shortcuts

Start new project:
```
julia
]
generate Blueprint
```
### Install package
Type this:
```
]
add Setfield
Ctrl + C
```

### Unit testing

Call `make test`. See [this Stack Overflow question](https://stackoverflow.com/q/70772117/1008794) for more information.


About dependencies, inheriting in test environment:
https://discourse.julialang.org/t/inheriting-package-dependencies-in-test-environment/68505/2

### In Emacs

`C-c C-p` : start REPL
`C-c C-a` : activate this project
`C-c C-b` : load module

### Module reloading

[Revise based workflows](https://docs.julialang.org/en/v1/manual/workflow-tips/#Revise-based-workflows)

## Objective: A CAD-modeling library for woodworking

This plan conforms with Roamdap.

### Task: Export to 3D

Parts:

* Create a sample drawing
* Implement colors (RgbColor) for BeamSpec
* Mesh from a beam, from polyhedron
* Export mesh, visualize
* Merge

#### Done

See commit: [Add mesh export](https://github.com/jonasseglare/Blueprint/commit/aaf4c05bc53a8932252f28686abfc22ed58b5685)

### Task: Render diagram

Steps:

* Givet en beam och en sida, exportera CuttingPlan
* Use Cairo/Luxor to 
* merge

### Task: Design insulation boards

* Make a new module
* Make the design
* Render 3D
* Render diagram
* Merge


