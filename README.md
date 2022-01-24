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
