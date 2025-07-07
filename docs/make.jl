using yk
using Documenter

DocMeta.setdocmeta!(yk, :DocTestSetup, :(using yk); recursive=true)

makedocs(;
    modules=[yk],
    authors="YuttariKanata <yuttarikanata@gmail.com> and contributors",
    sitename="yk.jl",
    format=Documenter.HTML(;
        canonical="https://YuttariKanata.github.io/yk.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/YuttariKanata/yk.jl",
    devbranch="master",
)
