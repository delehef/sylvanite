# Sylvanite

Sylvanite builds putative gene trees for a given list of gene families from both
syntenic and and sequence information, aiming at leveraging either when it is at
its most favorable timescale. Indeed, if syntenic information retains more
information on large chronological scales where sequence information tends to
become blurry, sequence information offers better granularity at smaller
timescale where syntenic information low resoltution fails to provide sufficient
resolution.

## Usage

Sylvanite work process requires the following data sources:
  - a species tree covering all the concerned species;
  - one or more gene families, _i.e._ a set of genes ID known to be descending
    of the same common ancestral gene;
  - pairwise genetic distance informations between all the genes of the same
    family, for instance obtained through family-wise alignment, _e.g._ with
    Diamond, a scoring scheme being then applied;
  - sytenic information for the involved genes, typically encoded as GFF files
    for the concerned species.

From there, a three steps pipeline is applied.


#### Examples

The following examples assume a folder structure akin to this one:

```
workdir
|
+-- db.sqlite        # generated by sylvanite build-database
|
+-- species_tree.nhx # species tree, newick format
|
+-- divergences/     # sequence pairwise distance matrices
|
+-- families/        # folder containing the gene families
|
+-- gff/             # GFF files for all the species
|
+-- syntenies/       # syntenic pairwise similarity matrices
|
+-- trees/           # final output
```

### Database Creation

Firstly, a SQLite database containing encoding all required syntenic information
in an easy-to-process way must be created. For this, the `sylvanite
build-database` command is used.

It requires the gene families for which the trees must be computed and the GFF
files for all the concerned species. Gene families must be provided in text
format, with one gene ID per line. GFF file must follow the GFF3 format.

Gene families files may be provided either individually, or through a folder
containing all of them, which prove much less cumbersome as soon as more than a
couple families are being processed.

The method to extract gene IDs and species names from GFF files default to
regexps that matches the metadata of the GFF files provided by ENSEMBL; however,
they can be altered with respectively `--id-pattern` and `--species-pattern` to
accomodate for GFF files following other conventions.

#### Example

This invokation would create, in `workdir/db.sqlite`, the database from GFFs
files obtained from ENSEMBL stored in `workdir/gffs/` and the families files
grouped in `workdir/families/`.

```
sylvanite build-database --gffs workdir/gff/ --families workdir/families/ --out workdir/db.sqlite"
```

#### `sylvanite build-database` reference

```
USAGE:
    sylvanite build-database [OPTIONS] --families <FAMILIES> --out <OUTFILE> <--gffs <GFFS>|--emf <EMF>>

OPTIONS:
        --families <FAMILIES>
            the files and/or directories containing the gene families to process

        --gffs <GFFS>
            the directory containing the GFF3 files to process; those can be gzipped

    -h, --help
            Print help information

        --id-pattern <ID_PATTERN>
            regex to extract feature name from GFF ID field; must contain a named capture group `id`
            [default: gene:(?P<id>.*)]

        --id-type <ID_TYPE>
            the features to extract from GFF files [default: gene]

    -o, --out <OUTFILE>
            where to write the database

    -q, --quiet
            Less output per occurrence

        --species-pattern <SPECIES_PATTERN>
            regex to extract species name from GFF file name; must contain a named capture group
            `species` [default: (?P<species>.*).gff3]

    -v, --verbose
            More output per occurrence
```

### Syntenic Alignment & Scoring

From the previosuly created database, syntenic landscapes of the genes in each
family are aligned in a pairwise fashion, then scored to generate distance
matrices akin the the ones coming from sequences alignment.

This requires the gene families files, as well as the syntenic database computed
in the previous step.

#### Example

This invokation would, from the previously computed database in
`workdir/db.sqlite` and the gene families files gathered in `workdir/families/`,
compute all the pairwise syntenic matrice distances for each family and store
them in the `wokrdir/syntenies` folder.

```
find workdir/families/ -type f | parallel -j2 sylvanite align --database workdir/db.sqlite -o workdir/syntenies {}
```

#### `sylvanite align` reference

```
sylvanite-align
Create intra-family syntenic distance matrices from the provided syntenics database

USAGE:
    sylvanite align [OPTIONS] --database <DATABASE> <INFILES>...

ARGS:
    <INFILES>...    one or more gene families for which to build the syntenic distance matrices

OPTIONS:
    -b, --bar                    if set, display a progress bar
    -D, --database <DATABASE>    the path to the genomes database; to be built with `build-database`
    -h, --help                   Print help information
    -o, --outdir <OUTDIR>        where to store the computed matrices
    -q, --quiet                  Less output per occurrence
    -v, --verbose                More output per occurrence
```

### Tree Reconstruction

Finally, from the provided distances matrices for sequence and synteny, the
species tree and the database, the `build-trees` subcommand is invoked to
generate the trees. It requires access to the database created in the first
step, the syntenic distance matrices, the sequence distance matrices, the
species tree and the gene families files.

#### Example

This would combine all the data previosuly computed to generate the

```
find workdir/families/ -type f | parallel -j2 sylvanite sylva --database workdir/db.sqlite --divergences workdir/divergences/ --syntenies workdir/syntenies -S workdir/species_tree.nhx --outdir workdir/trees/ {}"
```

#### `sylvanite sylva` reference

```
sylvanite-sylva
Create gene family trees from gene families, syntenic & sequence distance matrices, and syntenic
database

USAGE:
    sylvanite sylva [OPTIONS] --database <DATABASE> --species-tree <SPECIES_TREE> --syntenies <SYNTENIES> --divergences <DIVERGENCES> <INFILES>...

ARGS:
    <INFILES>...    the gene family files to build trees from

OPTIONS:
    -d, --divergences <DIVERGENCES>
            where to find the sequence distance matrices; can be a file or a path

    -D, --database <DATABASE>
            the path to the genomes database; to be built with `build-database`

        --dups-from-sequence
            if set, use sequence information rather than synteny sequence to resolve duplication
            arms

    -h, --help
            Print help information

        --merge-tandems
            if set, merge tandem genes together and span them post-hoc

        --no-overwrite
            if set, do not overwrite already existing files

    -o, --outdir <OUTDIR>
            where to write the created trees

    -q, --quiet
            Less output per occurrence

    -s, --syntenies <SYNTENIES>
            where to find the syntenic distance matrices; can be a file or a path

    -S, --species-tree <SPECIES_TREE>
            the species tree to use

        --synteny-threshold <SYNTENY_THRESHOLD>
            if set, override the automated synteny threshold choice

        --timings <TIMINGS>
            if set, where to write the computation time statistics

        --tracked <TRACKED>
            an optional set of comma-separated genes to track (for debug purposes)

    -v, --verbose
            More output per occurrence
```