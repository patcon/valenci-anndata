# Technical Onboarding

There is a lot of context on this project, for which we will eventually have better documentation. For now, here are some curated outside resources to understand some of the context and goals of people who work on polislike tools.

!!! warning

    🚧 _This page is a work-in-progress._ 🚧

## Polislike Pipelines

- `#todo` visually document general pipeline ([advanced pipeline](https://github.com/polis-community/red-dwarf/issues/53#issuecomment-2988756562))
- `#todo` contextual explanation of pipeline

## Algorithms

Polislike pipelines require 2 algorithms, and a 3rd optional one:

1. an **imputation** algorithm, for filling in best-guesses for all missing votes (required for next step)
2. a **dimensional reduction** algorithm, for projecting high-dimensional participants into low-dimensional 2D maps
3. a **clustering** algorithm (optional), for giving each participant a group label, with which to run statistics about statement votes
    - technically, a clustering algorithm is also required, and it's just optionally a _computational_ algorithm. But humans can use a subjective mental "algorithm" to paint the labels, e.g. "this looks like a blob on the map". In order to get aggregate stats, we need at least one label, so we can calculate "label #1 vs everyone else".

### Blank-filling (aka imputation)

`#todo`

### Dimensional Reduction

- Video: [Bluffer's Guide to Dimensional Reduction](https://www.youtube.com/watch?v=9iol3Lk6kyU) (2019)
    - What: A very accessible talk by the creator of a famous dimensional reduction algorithm, summarizing the landscape of algorithms.
    - Why you should watch: You wish you had more intuitions around how these algorithms work, and you perhaps don't feel "mathematical" enough to understand most of what's been published.
- `#todo` local vs global structure

### Clustering

`#todo`

### Further Algorithms

- comment routing algorithm
- group-informed consensus algorithm
- consensus statement selection algorithm
- group-representative statement selection algorithm

## Single-Cell "Omics"

- `#todo` compare/contrast problem domains
    - single-cell gene expression matrices
        - cell x gene matrix. response = gene expression.
        - non-valenced (only +), variable magnitude (0 and up).
        - No `NaN` responses, but many zeros are "missing", without clear distinction.
    - single-person vote reaction matrices
        - participant x statement matrix. response = vote reaction.
        - binary valenced (- or +), unit magnitide (0 or 1).
        - Response are zeros and "missing" (`NaN`), with clear distinction.

## Glossary

- vote data. valence data. polislike data.
    - data vs matrix?
- perspective map. perspective landscape.
    - constrast to "opinion map".
- participants. rows. observations. samples.
- statements/comments. columns. variables. features. stimuli?
- reactions. votes.
- imputation. blank-filling. missing. NaN.
- reducer. dimensional reduction.
- clustering. community detection.
    - labels.
- scikit learn.
    - estimator.
- anndata.
- single-cell omics. transcriptomics. gene expression. transcripts.
    - scverse.
        - scanpy.
- red-dwarf
- machine learning vs large-language models
- others?
    - agora
    - polisNL. partici.app.
    - Polis EU
    - DigiFinland. Voxit.