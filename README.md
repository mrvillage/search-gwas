# search-gwas

A simple CLI tool to download and search the [GWAS Catalog](https://www.ebi.ac.uk/gwas/).

## Installation

```bash
cargo install search-gwas
```

## Usage

```bash
search-gwas trait hypothyroidism -g COL5A2,TSHR
# OR
search-gwas trait hypothyroidism -g COL5A2 -g TSHR
```

### Additional options
`-a` show full association data
`-l` show PubMed links instead of IDs
`-c` output CSV data
