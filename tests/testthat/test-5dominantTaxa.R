context("dominantTaxa")

test_that("dominantTaxa", {

    test_dominantTaxa <- function(tse){
        # Counts table should not be changed
        expect_equal(assays(dominantTaxa(tse))$counts, assays(tse)$counts)

        # Should get same result
        expect_equal(colData(dominantTaxa(tse, name="dominant"))$dominant,
                     c("36155",  "256977",
                       "71074",  "331820", "331820", "98605",  "484436", "64396",  "360229", "114821",
                       "279599", "329744", "549656", "549656", "549656", "557211", "534609", "12812",
                       "87194",  "319044", "319044", "331820", "189047", "550960", "550960", "550960"))
        expect_equal(colData(dominantTaxa(tse, name="dominant"))$dominant, 
                     colData(addDominantTaxa(tse, name="dominant"))$dominant)

        expect_equal(colData(dominantTaxa(tse, rank = "Order"))$dominant_taxa,
                     c("Solibacterales",
                       "Spartobacteriales", "Sphingobacteriales", "Bacteroidales", "Bacteroidales", "Lactobacillales", "Pseudomonadales",
                       "Lactobacillales", "Neisseriales", "Clostridiales", "Nostocales", "Actinomycetales", "Stramenopiles",
                       "Stramenopiles", "Stramenopiles", "Synechococcales", "Flavobacteriales", "Flavobacteriales",
                       "Chromatiales", "Desulfobacterales", "Chromatiales", "Clostridiales", "Clostridiales", "Clostridiales",
                       "Clostridiales", "Clostridiales"))
        expect_equal(colData(dominantTaxa(tse, rank = "Order"))$dominant_taxa, 
                     colData(addDominantTaxa(tse, rank = "Order"))$dominant_taxa)

        expect_equal(colData(dominantTaxa(tse, rank = "Class"))$dominant_taxa,
                     c("Solibacteres",
                       "Spartobacteria", "Alphaproteobacteria", "Bacteroidia", "Bacteroidia", "Bacilli", "Gammaproteobacteria",
                       "Bacilli", "Betaproteobacteria", "Clostridia", "Nostocophycideae", "Actinobacteria", "Chloroplast",
                       "Chloroplast", "Chloroplast", "Alphaproteobacteria", "Flavobacteria", "Chloroplast", "Gammaproteobacteria",
                       "Gammaproteobacteria", "Deltaproteobacteria", "Clostridia", "Clostridia", "Clostridia", "Clostridia", "Clostridia"))
        expect_equal(colData(dominantTaxa(tse, rank = "Class"))$dominant_taxa, 
                     colData(addDominantTaxa(tse, rank = "Class"))$dominant_taxa)

        expect_equal(colData(dominantTaxa(tse, rank = "Kingdom"))$dominant_taxa,
                     c("Bacteria",
                       "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria",
                       "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria",
                       "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria"))
        expect_equal(colData(dominantTaxa(tse, rank = "Kingdom"))$dominant_taxa, 
                     colData(addDominantTaxa(tse, rank = "Kingdom"))$dominant_taxa)

    }

    # TSE object
    data(GlobalPatterns)
    test_dominantTaxa(GlobalPatterns)
    tse <- GlobalPatterns
    assay(tse,"counts") <- DelayedArray(assay(tse,"counts"))
    test_dominantTaxa(tse)

    # TSE object
    tse <- microbiomeDataSets::dietswap()

    expect_equal(head(colData(dominantTaxa(tse, rank = "Genus"))$dominant_taxa, 20),
                 c("Bacteroides vulgatus et rel.",
                   "Prevotella melaninogenica et rel.", "Prevotella melaninogenica et rel.", "Oscillospira guillermondii et rel.",
                   "Prevotella melaninogenica et rel.", "Prevotella melaninogenica et rel.", "Clostridium cellulosi et rel.",
                   "Prevotella melaninogenica et rel.", "Streptococcus bovis et rel.", "Prevotella melaninogenica et rel.",
                   "Prevotella melaninogenica et rel.", "Prevotella melaninogenica et rel.", "Oscillospira guillermondii et rel.",
                   "Oscillospira guillermondii et rel.", "Prevotella melaninogenica et rel.", "Prevotella melaninogenica et rel.",
                   "Prevotella melaninogenica et rel.", "Oscillospira guillermondii et rel.", "Oscillospira guillermondii et rel.",
                   "Prevotella melaninogenica et rel." ))
    expect_equal(colData(dominantTaxa(tse, rank = "Genus"))$dominant_taxa, 
                 colData(addDominantTaxa(tse, rank = "Genus"))$dominant_taxa)

    expect_equal(head(colData(dominantTaxa(tse, rank = "Family"))$dominant_taxa, 20),
                 c("Bacteroidetes",
                   "Clostridium cluster IV", "Bacteroidetes", "Clostridium cluster IV", "Bacteroidetes", "Bacteroidetes",
                   "Clostridium cluster IV", "Bacteroidetes", "Clostridium cluster IV", "Clostridium cluster IV", "Bacteroidetes",
                   "Bacteroidetes", "Clostridium cluster IV", "Clostridium cluster IV", "Bacteroidetes", "Clostridium cluster IV",
                   "Bacteroidetes", "Clostridium cluster IV", "Clostridium cluster IV", "Bacteroidetes"))
    expect_equal(colData(dominantTaxa(tse, rank = "Family"))$dominant_taxa, 
                 colData(addDominantTaxa(tse, rank = "Family"))$dominant_taxa)

    expect_equal(head(colData(dominantTaxa(tse, rank = "Phylum"))$dominant_taxa, 20),
                 c("Bacteroidetes",
                   "Firmicutes", "Bacteroidetes", "Firmicutes", "Bacteroidetes", "Bacteroidetes", "Firmicutes", "Bacteroidetes",
                   "Firmicutes", "Firmicutes", "Bacteroidetes", "Firmicutes", "Firmicutes", "Firmicutes", "Bacteroidetes",
                   "Firmicutes", "Bacteroidetes", "Firmicutes", "Firmicutes", "Bacteroidetes"))
    expect_equal(colData(dominantTaxa(tse, rank = "Phylum"))$dominant_taxa, 
                 colData(addDominantTaxa(tse, rank = "Phylum"))$dominant_taxa)
})


test_that("summarizeDominantTaxa", {

    test_summarizeDominantTaxa <- function(tse){
        expect_equal(summarizeDominantTaxa(tse)$dominant_taxa,
                     c("331820",
                       "549656", "550960", "319044", "114821", "12812",  "189047", "256977",
                       "279599", "329744", "360229", "36155", "484436", "534609", "557211", "64396",
                       "71074", "87194", "98605" ))

        expect_equal(summarizeDominantTaxa(tse, rank = "Kingdom", name = "dominant_kingdom")$dominant_kingdom,
                     c("Bacteria"))

        expect_equal(summarizeDominantTaxa(tse, rank = "Order", name = "dominant_kingdom")$rel.freq,
                     c(23.1, 11.5, 7.7, 7.7, 7.7, 7.7, 3.8, 3.8, 3.8, 3.8, 3.8, 3.8, 3.8, 3.8, 3.8))

        expect_equal(summarizeDominantTaxa(tse, rank = "Class", name = "dominant_kingdom")$rel.freq.pct,
                     c("23%", "15%", "12%", "8%", "8%", "8%", "4%", "4%", "4%", "4%", "4%", "4%", "4%"))
    }

    # TSE object
    data(GlobalPatterns)
    test_summarizeDominantTaxa(GlobalPatterns)
    tse <- GlobalPatterns
    assay(tse,"counts") <- DelayedArray(assay(tse,"counts"))
    test_summarizeDominantTaxa(tse)

    # TSE object
    tse <- microbiomeDataSets::dietswap()

    expect_equal(summarizeDominantTaxa(tse, rank = "Genus")$n,
                 c(104, 57, 43, 9, 2, 2, 2, 1, 1, 1))

    expect_equal(summarizeDominantTaxa(tse, group = "nationality")$nationality,
                 as.factor(c("AFR", "AAM", "AAM", "AFR", "AAM", "AAM", "AAM", "AAM", "AFR",
                   "AAM", "AAM", "AAM", "AFR", "AFR")))

    expect_equal(summarizeDominantTaxa(tse, rank = "Family", group = "sex")$dominant_taxa,
                 c("Bacteroidetes", "Bacteroidetes", "Clostridium cluster IV", "Clostridium cluster IV",
                   "Clostridium cluster XIVa", "Verrucomicrobia"))

    expect_equal(summarizeDominantTaxa(tse, rank = "Phylum", group = "timepoint", name = "test")$n,
                 c(24, 23, 20, 19, 19, 19, 19, 18, 17, 16, 15, 13))

    expect_equal(summarizeDominantTaxa(tse, rank = "Genus", group = "timepoint.within.group", name = "testing")$rel.freq,
                 c(46.4, 47.3, 27.3, 24.1, 20.5, 18.2, 4.5, 3.6, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9))

    expect_equal(summarizeDominantTaxa(tse, group = "bmi_group")$rel.freq.pct,
                 c("51%", "66%", "39%", "26%", "28%", "19%", "20%", "20%", "7%", "5%", "4%", "3%", "3%", "3%", "2%", "1%", "1%", "1%", "1%"))

})

