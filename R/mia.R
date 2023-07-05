#' \code{mia} Package.
#'
#' \code{mia} implements tools for microbiome analysis based on the
#' \code{SummarizedExperiment}, \code{SingleCellExperiment} and
#' \code{TreeSummarizedExperiment} infrastructure. Data wrangling and analysis
#' in the context of taxonomic data is the main scope. Additional functions for
#' common task are implemented such as community indices calculation and
#' summarization.
#'
#' @name mia-package
#' @docType package
#' @seealso \link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment} class
NULL

#' @import methods
#' @import TreeSummarizedExperiment
#' @import DelayedArray
#' @import scater
#' @importFrom dplyr %>%
#' @importFrom rlang sym :=
NULL

#' @title mia datasets
#'
#' @description
#' These datasets are conversions of the \pkg{phyloseq} datasets
#' \code{GlobalPatterns}, \code{enterotype}, and \code{esophagus} into the
#' \link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}
#' data container.
#'
#' \code{dmn_se} contains an example \code{SummarizedExperiment} derived
#' from data in the \pkg{DirichletMultinomial} package. See
#' \code{?calculateDMN} for more details.
#'

#' Global patterns of 16S rRNA diversity at a depth of millions of sequences per sample (2011)
#'
#' This work compared the microbial communities from 25 environmental samples and three known ``mock communities''
#' at a an averag depth of 3.1 million reads per sample.
#' Authors reproduced diversity patterns seen in many other 
#' published studies, while invesitigating technical issues/bias by 
#' applying the same techniques to simulated microbial communities of known
#' composition. Many thanks to J. Gregory Caporaso for providing the OTU-clustered data files
#' for inclusion in the \pkg{phyloseq} package, from which this data \code{TreeSummarizedExperiment}
#' version was then converted.
#'
#' @name mia-datasets
#' @aliases GlobalPatterns
#' @keywords data
#' @usage data(GlobalPatterns)
#' @author Caporaso, J. G., et al.
#' @docType data
#' @references
#' Caporaso, J. G., et al. (2011). 
#' Global patterns of 16S rRNA diversity at a depth of millions of sequences per sample.
#' PNAS, 108, 4516-4522.
#' \url{http://www.pnas.org/content/108/suppl.1/4516.short}
"GlobalPatterns"

#' Enterotype data
#'
#' The enterotype data of the human gut microbiome (Arumugam et al. 2011) includes
#' taxonomic profiling for 280 fecal samples from 22 subjects based on shotgun DNA sequencing.
#' The authors claimed that the data naturally clumps into three community-level clusters, or
#' ``enterotypes'', that are not immediately explained by sequencing technology or demographic 
#' features of the subjects. A later addendum (2014) the authors stated that enterotypes
#' "should not be seen as discrete clusters, but as a way of stratifying samples to reduce complexity."
#'
#' @name mia-datasets
#' @aliases enterotype
#' @keywords data
#' @usage data(enterotype)
#' @author Arumugam, M., Raes, J., et al.
#' @docType data
#' @references
#' Arumugam, M., et al. (2011). Enterotypes of the human gut microbiome.
#' Nature, 473(7346), 174-180.
#' \url{http://www.nature.com/doifinder/10.1038/nature09944}
#' Supplemental information includes subject data. 
#' OTU-clustered data was initially downloaded from the publicly-accessible:
#' \url{http://www.bork.embl.de/Docu/Arumugam_et_al_2011/downloads.html}
#'
#' Arumugam, M., et al. (2014). Addendum: Enterotypes of the human gut microbiome.
#' Nature 506, 516 (2014). \url{https://doi.org/10.1038/nature13075}
"enterotype"

#' Small example dataset from a human esophageal community
#' 
#' The esophagus data set from Pei et al. (2004)
#' includes 3 samples from 3 human adults based on biopsies analysed with 16S rDNA PCR.
#' The 16S rRNA sequence processing has been provided in the mothur wiki
#' at the link below. 
#'
#' @name mia-datasets
#' @aliases esophagus
#' @keywords data
#' @usage data(esophagus)
#' @author Pei et al. \email{zhiheng.pei@@med.nyu.edu}.
#' @docType data
#' @references 
#' Pei, Z., Bini, E. J., Yang, L., Zhou, M., Francois, F., & Blaser, M. J. (2004). 
#' Bacterial biota in the human distal esophagus.
#' Proceedings of the National Academy of Sciences of the United States of America, 101(12), 4250-4255.
#' \url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC384727}
#'
#' McMurdie, J. & Holmes, S. (2013) \emph{phyloseq}: An R Package for reproducible interactive analysis
#' and graphics of microbiome census data. PLoS ONE. 8(4):e61217.
#' \url{https://doi.org/10.1371/journal.pone.0061217}
#'
#' Mothur-processed files and the sequence data can be downloaded at:
#' \url{http://www.mothur.org/wiki/Esophageal_community_analysis}
"esophagus"
#'

#' Twins data set for Dirichlet Multinomial Mixtures (DMM)
#'
#' This data set from Turnbaugh et al. (2009) was used to introduce
#' Dirichlet Multinomial Mixtures (DMM) for microbiota stratification by
#' Holmes et al. (2012).
#' 
#' @name mia-datasets
#' @aliases dmn_se, twins
#' @keywords data
#' @usage data(dmn_se)
#' @author Turnbaugh, PJ et al.
#' @references
#' Holmes I, Harris K, Quince C (2012).
#' Dirichlet Multinomial Mixtures: Generative Models for Microbial Metagenomics.
#' PLoS ONE 7(2): e30126. \url{https://doi.org/10.1371/journal.pone.0030126}
#'
#' Turnbaugh PJ, Hamady M, Yatsunenko T, Cantarel BL, Duncan A, et al. (2009).
#' A core gut microbiome in obese and lean twins. Nature 457: 480–484. 
#' \url{https://doi.org/10.1038/nature07540}
"dmn_se"

#' peerj13075
#' 
#'
#' PeerJ data by Potbhare et al. (2022) includes skin microbial profiles of 58 volunteers with multiple factors. 
#' 16S r-RNA sequencing of V3-V4 regions was done to generate millions of read using illumina platform.
#' A standard bioinformatic and statistical analysis done to explore skin bacterial diversity and its association with age, diet, geographical locations.
#' The authors investigated significant association of skin microbiota with individual’s geographical location.
#'
#' @name mia-datasets
#' @aliases peerj13075
#' @keywords data
#' @usage data(peerj13075)
#' @author Potbhare, R., et al.
#' @docType data
#' @references
#' Potbhare, R., RaviKumar, A., Munukka, E., Lahti, L., & Ashma, R. (2022). 
#' Skin microbiota diversity among genetically unrelated individuals of Indian origin. 
#' PeerJ, 10, e13075.
#' \url{https://peerj.com/articles/13075/}
#' Supplemental information includes OTU table and taxonomy table and publicly-accessible from: 
#' \url{https://www.doi.org/10.7717/peerj.13075/supp-1}
#' \url{https://www.doi.org/10.7717/peerj.13075/supp-2}
"peerj13075"

#' Multiomics dataset from a rat experiment studying effect of fat and prebiotics in diet
#' 
#' The HintikkaXO dataset contains high-throughput profiling data from 40 rat 
#' samples, including 39 biomarkers, 38 metabolites (NMR), and 12706 OTUs from 
#' 318 species, measured from Cecum. This is diet comparison study with High/Low 
#' fat diet and xylo-oligosaccaride supplementation. Column metadata is common 
#' for all experiments (microbiota, metabolites, biomarkers) and includes the 
#' following fields:
#' 
#' \itemize{
#'   \item{Sample: Sample ID (character)}
#'   \item{Rat: Rat ID (factor)}
#'   \item{Site: Site of measurement ("Cecum"); single value}
#'   \item{Diet: Diet group (factor; combination of the Fat and XOS fields)}
#'   \item{Fat: Fat in Diet (factor; Low/High)}
#'   \item{XOS: XOS Diet Supplement (numeric; 0/1)}
#' }
#' 
#' Row metadata of the microbiota data contains taxonomic information on the 
#' Phylum, Class, Order, Family, Genus, Species, and OTU levels.
#' 
#' Biomarker data contains 39 biomarkers.
#' 
#' Metabolite data contains 38 NMR metabolites.
#' 
#' @name mia-datasets
#' @aliases HintikkaXOData
#' @keywords data
#' @usage data(HintikkaXOData)
#' @author Leo Lahti et al.
#' @docType data
#' @references
#' Hintikka L et al. (2021): Xylo-oligosaccharides in prevention of hepatic 
#' steatosis and adipose tissue inflammation: associating taxonomic and 
#' metabolomic patterns in fecal microbiotas with biclustering. International 
#' Journal of Environmental Research and Public Health 18(8):4049 
#' \url{https://doi.org/10.3390/ijerph18084049}
#' 
"HintikkaXOData"

#' Tengeler2020
#' 
#' Tengeler data by Tengeler et al. (2022) includes gut microbiota profiles of 27 persons with ADHD. 
#' A standard bioinformatic and statistical analysis done to  demonstrate that altered microbial 
#' composition could be a driver of altered brain structure and function and concomitant changes in the animals’ behavior.
#' They investigated this by colonizing young, male, germ-free C57BL/6JOlaHsd mice with microbiota from individuals with and without ADHD.
#'
#' @name mia-datasets
#' @aliases Tengeler2020
#' @keywords data
#' @usage data(Tengeler2020)
#' @author A.C. Tengeler, et al.
#' @docType data
#' @references
#' Tengeler, A.C., Dam, S.A., Wiesmann, M. et al. 
#' Gut microbiota from persons with attention-deficit/hyperactivity disorder affects the brain in mice. 
#' Microbiome 8, 44 (2020). 
#' \url{https://doi.org/10.1186/s40168-020-00816-x}
#' Supplemental information includes Home-cage activity, methods, results and imaging parameters and publicly-accessible from: 
#' \url{https://static-content.springer.com/esm/art%3A10.1186%2Fs40168-020-00816-x/MediaObjects/40168_2020_816_MOESM1_ESM.docx}
#' \url{https://static-content.springer.com/esm/art%3A10.1186%2Fs40168-020-00816-x/MediaObjects/40168_2020_816_MOESM2_ESM.docx}
#' \url{https://static-content.springer.com/esm/art%3A10.1186%2Fs40168-020-00816-x/MediaObjects/40168_2020_816_MOESM3_ESM.docx}
#' 
"Tengeler2020"
