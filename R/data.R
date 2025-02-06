#' Example Antigenic Mapping Data
#'
#' HI titers of Influenza antigens and antisera published in Smith et al., 2004 were 
#' used to find the antigenic relationships and coordinates of the antigens.
#' It can be used for mapping. The data captures how different influenza virus strains 
#' (antigens) react with antisera from infected individuals.
#'
#' @format A data frame with 285 rows and 11 variables:
#' \describe{
#'   \item{V1}{First dimension coordinate from 5D mapping}
#'   \item{V2}{Second dimension coordinate from 5D mapping}
#'   \item{V3}{Third dimension coordinate from 5D mapping}
#'   \item{V4}{Fourth dimension coordinate from 5D mapping}
#'   \item{V5}{Fifth dimension coordinate from 5D mapping}
#'   \item{name}{Strain identifier}
#'   \item{antigen}{Logical; TRUE if point represents an antigen}
#'   \item{antiserum}{Logical; TRUE if point represents an antiserum} 
#'   \item{cluster}{Factor indicating antigenic cluster assignment (A/H3N2 1968-2003)}
#'   \item{color}{Color assignment for visualization}
#'   \item{year}{Year of strain isolation}
#' }
#' @source Arhami and Rohani 2025
#'         \doi{}
"example_positions"

#' H3N2 Influenza HI Assay Data from Smith et al. 2004
#'
#' Hemagglutination inhibition (HI) assay data for influenza A/H3N2 viruses
#' spanning 35 years of evolution.
#'
#' @format A data frame with the following variables:
#' \describe{
#'   \item{virusStrain}{Character. Virus strain identifier}
#'   \item{serumStrain}{Character. Antiserum strain identifier}
#'   \item{titer}{Numeric. HI assay titer value}
#'   \item{virusYear}{Numeric. Year virus was isolated}
#'   \item{serumYear}{Numeric. Year serum was collected}
#'   \item{cluster}{Factor. Antigenic cluster assignment}
#'   \item{color}{Character. Color code for visualization}
#' }
#' @source Smith et al. (2004) Science, 305(5682), 371-376.
"h3n2_data"

#' HIV Virus Metadata
#'
#' Reference information for HIV virus strains used in neutralization assays.
#'
#' @format A data frame with the following variables:
#' \describe{
#'   \item{Virus.name}{Character. Virus strain identifier}
#'   \item{Country}{Character. Country of origin}
#'   \item{Subtype}{Character. HIV subtype}
#'   \item{Year}{Numeric. Year of isolation}
#' }
#' @source Los Alamos HIV Database (https://www.hiv.lanl.gov/)
"hiv_viruses"

#' HIV Neutralization Assay Data
#'
#' IC50 neutralization measurements between HIV viruses and antibodies.
#'
#' @format A data frame with the following variables:
#' \describe{
#'   \item{Antibody}{Character. Antibody identifier}
#'   \item{Virus}{Character. Virus strain identifier}
#'   \item{IC50}{Numeric. IC50 neutralization value}
#' }
#' @source Los Alamos HIV Database (https://www.hiv.lanl.gov/)
"hiv_titers"