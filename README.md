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

## Objective: A CAD-modeling library for woodworking [rdap: 33% done]

This plan conforms with Roamdap.
An Objective is done when all subtasks and subtasks are done, thereby acting like a grouping.
A Task is done when it has at least one subsection marked done.
A Task can be marked closed if it is no longer pursued. It cannot at the same time be done.
A Task is open as long as it is not closed or done.

### Task: Export to 3D [rdap: 100% done]

Parts:

* Create a sample drawing
* Implement colors (RgbColor) for BeamSpec
* Mesh from a beam, from polyhedron
* Export mesh, visualize
* Merge

#### Done

See commit: [Add mesh export](https://github.com/jonasseglare/Blueprint/commit/aaf4c05bc53a8932252f28686abfc22ed58b5685)

### Task: Organize and render diagram [rdap: 0% done]

Steps:

* Givet en beam och en sida, exportera CuttingPlan
  - Borrhål
  - Hur man ska såga (konturen)
* Givet flera cutting plans för samma beam-typ, optimera dessa givet brädans längd. Det ger oss en sammansatt cutting plan.
* Rendera diagram med hjälp av Luxor eller Cairo, till PDF kanske.
* Merge

### Task: Exportera tabell med mått

Detta är en fortsättning på den andra uppgiften.

* Givet en cutting plan, mata ut 

### Task: Design insulation boards [rdap: 0% done]

* Make a new module
* Make the design
* Render 3D
* Render diagram
* Merge

### Roamdap summary

Todo
