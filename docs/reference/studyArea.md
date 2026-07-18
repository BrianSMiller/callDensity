# Calculate study area for call density estimate.

Estimate study area for call density given a radius from the hydrophone.
If only a radius is provided, then the area=pi\*radius^2 is returned.
Optionally, if a matrix of truncation distances are provided, then the
study area is calculated as the sum of the areas of each radial transect
(assuming N transects will have the same 360/N degrees of angular
coverage around the circle).

## Usage

``` r
studyArea(w, truncationDistance = w)
```

## Arguments

- w:

  - Radius of study area

- truncationDistance:

  - Matrix of truncation distances `[1 x Nr]`. Nr is the number of
    radial transects in the simulation (i.e. the same as the number of
    TL profiles).

## Value

Area of the simulation
