# Introduction

**valency-anndata** (pronounced "valencian data") is an experimental toolkit for exploration of vote matrices of "valence data" (_participants x statements_), also known as _polislike_ data.

We aim to grow and support a community of data scientists, statisticians, software developers, UI designers, and deliberative process facilitators who wish to democratize insights from valence data, and strengthen democratic self-governance.

## Attribution

We build on the following best-practices and standards:

- [anndata][], a tool for working with 2D annotated data matrices,
- [scanpy][], a toolkit for exploring gene expression matrices (_cells x genes_), and
- [the scverse][], the community that has grown around the above Python tooling.

We are inspired by the solutions to comparable challenges that are being resolved in the single-cell -omics community, we aspire to bring comparable rigour, collaboration, and reproducibility to the social sciences and its corresponding valence data.

## Commands

* `mkdocs new [dir-name]` - Create a new project.
* `mkdocs serve` - Start the live-reloading docs server.
* `mkdocs build` - Build the documentation site.
* `mkdocs -h` - Print help message and exit.

## Project layout

    mkdocs.yml    # The configuration file.
    docs/
        index.md  # The documentation homepage.
        ...       # Other markdown pages, images and other files.

<!-- Links -->
   [anndata]: https://scanpy.readthedocs.io/en/stable/
   [scanpy]: https://scanpy.readthedocs.io/en/latest/
   [the scverse]: https://scverse.org/