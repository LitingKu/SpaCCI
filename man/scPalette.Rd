% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visiualization.R
\name{scPalette}
\alias{scPalette}
\title{Generate a Color Palette}
\usage{
scPalette(n)
}
\arguments{
\item{n}{An integer specifying the number of colors needed.}
}
\value{
A character vector of colors in hexadecimal format.
}
\description{
This function generates a color palette. It selects colors from a predefined color space, and if more colors are needed than are available in the predefined set, it generates a palette using color interpolation.
}
\examples{
\donttest{
# Generate a palette with 5 colors
palette <- scPalette(5)
print(palette)

# Generate a palette with 30 colors
large_palette <- scPalette(30)
print(large_palette)
}
}
