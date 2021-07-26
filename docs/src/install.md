# Installation

`ArrayTools` can be installed as any other [offical Julia
packages](https://pkg.julialang.org/):

```julia
… pkg> add ArrayTools
```

where `… pkg>` represents the package manager prompt (the ellipsis `…` denote
your current environment).  To start Julia's package manager, launch Julia and,
at the [REPL of
Julia](https://docs.julialang.org/en/stable/manual/interacting-with-julia/),
hit the `]` key; you should get the above `… pkg>` prompt.  To revert to
Julia's REPL, just hit the `Backspace` key at the `… pkg>` prompt.

To check whether the `ArrayTools` package works correctly, type:

```julia
… pkg> test ArrayTools
```

Later, to update to the last version (and run tests), you can type:

```julia
… pkg> update ArrayTools
… pkg> test ArrayTools
```

If something goes wrong, it may be because you already have an old version of
`ArrayTools`.  Uninstall `ArrayTools` as follows:

```julia
… pkg> rm ArrayTools
… pkg> gc
… pkg> add https://github.com/emmt/ArrayTools.jl
```

before re-installing.
