# Changelog

All notable changes to this project will be documented in this file.

## [1.8.1] - 2024-08-05

### Miscellaneous Tasks

- Update dependencies

## [1.8.0] - 2024-07-22

### Bug Fixes

- Effective leave gathering for grafted subtrees
- ELC computation when grafting into an empty branch

### Documentation

- Add README

### Features

- Add the `--synteny-threshold` option

### Miscellaneous Tasks

- Release sylvanite version 1.8.0

## [1.7.0] - 2024-01-26

### Bug Fixes

- Drop heuristic for subtree grafting

### Miscellaneous Tasks

- Clippy
- Release sylvanite version 1.7.0

## [1.6.2] - 2023-11-29

### Features

- Add `--dups-from-sequence`

### Miscellaneous Tasks

- Update syntesuite
- Release sylvanite version 1.6.2

## [1.6.1] - 2023-07-16

### Miscellaneous Tasks

- Update dependencies
- Release sylvanite version 1.6.1

## [1.6.0] - 2023-07-03

### Bug Fixes

- Do not use non-initialized cache
- Use a different ELC computation scheme

### Features

- Add a `--tracked` debug flag

### Miscellaneous Tasks

- Update dependencies
- Release sylvanite version 1.6.0

### Hack

- Give more leeway on the threshold

## [1.5.14] - 2023-06-14

### Miscellaneous Tasks

- Downgrade dependency for Guix
- Release sylvanite version 1.5.14

## [1.5.13] - 2023-06-14

### Miscellaneous Tasks

- Update dependencies
- Release sylvanite version 1.5.13

## [1.5.12] - 2023-06-14

### Miscellaneous Tasks

- Add git-cliff to dependencies
- Update dependencies
- Release sylvanite version 1.5.11
- Release sylvanite version 1.5.12

## [1.5.10] - 2023-04-27

### Features

- Indicate the grafting criterion used

### Miscellaneous Tasks

- Release sylvanite version 1.5.10

## [1.5.9] - 2023-03-30

### Miscellaneous Tasks

- Release sylvanite version 1.5.9

### Performance

- Pre-prune the tree before grafting subtrees

## [1.5.8] - 2023-03-27

### Bug Fixes

- Use extended as clusters in pathologic cases

### Miscellaneous Tasks

- Release sylvanite version 1.5.8

## [1.5.7] - 2023-03-27

### Bug Fixes

- Remove debug prints

### Miscellaneous Tasks

- Release sylvanite version 1.5.7

## [1.5.6] - 2023-03-27

### Bug Fixes

- Metagenes

### Miscellaneous Tasks

- Update dependencies
- Release sylvanite version 1.5.6

### Performance

- Use as many iterators as possible; cache descendants nodes
- Implement clutchable leaves caching

## [1.5.5] - 2023-03-16

### Bug Fixes

- UNKNWN species for fanned out tandems

### Miscellaneous Tasks

- Release sylvanite version 1.5.5

## [1.5.4] - 2023-03-16

### Miscellaneous Tasks

- Downgrade dependencies for Guix
- Release sylvanite version 1.5.4

## [1.5.3] - 2023-03-16

### Features

- Make tandem compression/fanout available under --merge-tandems

### Miscellaneous Tasks

- Update dependencies
- Release sylvanite version 1.5.3

### Performance

- Use smallvec to greatly spare allocations

### Refactor

- Rename proteins to genes

## [1.5.2] - 2023-01-30

### Bug Fixes

- Do not output internal data

### Miscellaneous Tasks

- Release sylvanite version 1.5.2

## [1.5.1] - 2023-01-23

### Miscellaneous Tasks

- Update dependencies
- Release sylvanite version 1.5.1

## [1.5.0] - 2023-01-23

### Miscellaneous Tasks

- Remove unused dependency
- Release sylvanite version 1.5.0

## [1.4.2] - 2023-01-23

### Miscellaneous Tasks

- Update dependencies
- Release sylvanite version 1.4.2

## [1.4.1] - 2023-01-18

### Miscellaneous Tasks

- Update syntesuite
- Release sylvanite version 1.4.1

## [1.4.0] - 2023-01-18

### Miscellaneous Tasks

- Release sylvanite version 1.4.0

### Refactor

- Extract syntesuite

## [1.3.0] - 2023-01-18

### Documentation

- Document CLI options

### Features

- Better message error in case of missing files
- Add `build-database`

### Miscellaneous Tasks

- Sylva should not handle species tree conventions
- Update dependencies
- Stderrlog -> buche
- Add error types
- Update dependencies
- Standardize errors
- Update newick
- Release sylvanite version 1.3.0

### Refactor

- Simplify database handling
- Use clap_verbose

## [1.2.1] - 2023-01-04

### Bug Fixes

- Out of bounds

### Miscellaneous Tasks

- Gitignore++
- Release sylvanite version 1.2.1

## [1.2.0] - 2023-01-02

### Bug Fixes

- Use correct names for trees

### Miscellaneous Tasks

- Release sylvanite version 1.2.0

## [1.1.0] - 2023-01-01

### Features

- Introduce a partial-caching of the DB
- Create output dirs if they do not exist

### Miscellaneous Tasks

- Clippy
- Dead code pruning
- Release sylvanite version 1.1.0

## [1.0.1] - 2022-12-30

### Miscellaneous Tasks

- Downgrade rusqlite for Guix packaging
- Release sylvanite version 1.0.1

## [1.0.0] - 2022-12-30

### Bug Fixes

- Duplication growth CS computation

### Features

- Better computation of synteny when assembling the final tree
- Satellites and solos are now injected as cardinal one subtrees
- Grafted subtrees are first ordered by decreasing topo depth
- Add git-cliff & cargo-release

### Miscellaneous Tasks

- Renaming
- Release sylvanite version 1.0.0

<!-- generated by git-cliff -->
