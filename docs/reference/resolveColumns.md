# Resolve column names from a prefix pattern or explicit vector

Helper to standardise column selection across functions. If `prefix` is
a single string it is used as a `grepl` pattern against `df_names`; if
it is a character vector of length \> 1 it is used directly as column
names.

## Usage

``` r
resolveColumns(df_names, prefix)
```

## Arguments

- df_names:

  Character vector of column names (i.e. `names(df)`).

- prefix:

  Either a single character string pattern or a character vector of
  explicit column names.

## Value

A character vector of resolved column names.
