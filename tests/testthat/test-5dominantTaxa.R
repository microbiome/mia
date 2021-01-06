context("dominantTaxa")

test_that("dominantTaxa", {

    # TSE object
    data(GlobalPatterns)
    tse <- GlobalPatterns

    # Counts table should not be changed
    testthat::expect_equal(assays(mia::dominantTaxa(tse))$counts, assays(tse)$counts)

    # Should get same result
    testthat::expect_equal(colData(mia::dominantTaxa(tse, name="dominant"))$dominant, c("36155",  "256977",
    "71074",  "331820", "331820", "98605",  "484436", "64396",  "360229", "114821",
    "279599", "329744", "549656", "549656", "549656", "557211", "534609", "12812",
    "87194",  "319044", "319044", "331820", "189047", "550960", "550960", "550960"))

    testthat::expect_equal(colData(mia::dominantTaxa(tse, rank = "Order"))$dominant_taxa, c("Solibacterales",
    "Spartobacteriales", "Sphingobacteriales", "Bacteroidales", "Bacteroidales", "Lactobacillales", "Pseudomonadales",
    "Lactobacillales", "Neisseriales", "Clostridiales", "Nostocales", "Actinomycetales", "Stramenopiles",
    "Stramenopiles", "Stramenopiles", "Synechococcales", "Flavobacteriales", "Flavobacteriales",
    "Chromatiales", "Desulfobacterales", "Chromatiales", "Clostridiales", "Clostridiales", "Clostridiales",
    "Clostridiales", "Clostridiales"))

    testthat::expect_equal(colData(mia::dominantTaxa(tse, rank = "Class"))$dominant_taxa, c("Solibacteres",
    "Spartobacteria", "Alphaproteobacteria", "Bacteroidia", "Bacteroidia", "Bacilli", "Gammaproteobacteria",
    "Bacilli", "Betaproteobacteria", "Clostridia", "Nostocophycideae", "Actinobacteria", "Chloroplast",
    "Chloroplast", "Chloroplast", "Alphaproteobacteria", "Flavobacteria", "Chloroplast", "Gammaproteobacteria",
    "Gammaproteobacteria", "Deltaproteobacteria", "Clostridia", "Clostridia", "Clostridia", "Clostridia", "Clostridia"))

    testthat::expect_equal(colData(mia::dominantTaxa(tse, rank = "Kingdom"))$dominant_taxa, c("Bacteria",
    "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria",
    "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria",
    "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria"))

    # TSE object
    tse <- microbiomeDataSets::dietswap()

    testthat::expect_equal(head(colData(mia::dominantTaxa(tse, rank = "Genus"))$dominant_taxa, 20), c("Bacteroides vulgatus et rel.",
    "Prevotella melaninogenica et rel.", "Prevotella melaninogenica et rel.", "Oscillospira guillermondii et rel.",
    "Prevotella melaninogenica et rel.", "Prevotella melaninogenica et rel.", "Clostridium cellulosi et rel.",
    "Prevotella melaninogenica et rel.", "Streptococcus bovis et rel.", "Prevotella melaninogenica et rel.",
    "Prevotella melaninogenica et rel.", "Prevotella melaninogenica et rel.", "Oscillospira guillermondii et rel.",
    "Oscillospira guillermondii et rel.", "Prevotella melaninogenica et rel.", "Prevotella melaninogenica et rel.",
    "Prevotella melaninogenica et rel.", "Oscillospira guillermondii et rel.", "Oscillospira guillermondii et rel.",
    "Prevotella melaninogenica et rel." ))

    testthat::expect_equal(head(colData(mia::dominantTaxa(tse, rank = "Family"))$dominant_taxa, 20), c("Bacteroidetes",
    "Clostridium cluster IV", "Bacteroidetes", "Clostridium cluster IV", "Bacteroidetes", "Bacteroidetes",
    "Clostridium cluster IV", "Bacteroidetes", "Clostridium cluster IV", "Clostridium cluster IV", "Bacteroidetes",
    "Bacteroidetes", "Clostridium cluster IV", "Clostridium cluster IV", "Bacteroidetes", "Clostridium cluster IV",
    "Bacteroidetes", "Clostridium cluster IV", "Clostridium cluster IV", "Bacteroidetes"))

    testthat::expect_equal(head(colData(mia::dominantTaxa(tse, rank = "Phylum"))$dominant_taxa, 20), c("Bacteroidetes",
    "Firmicutes", "Bacteroidetes", "Firmicutes", "Bacteroidetes", "Bacteroidetes", "Firmicutes", "Bacteroidetes",
    "Firmicutes", "Firmicutes", "Bacteroidetes", "Firmicutes", "Firmicutes", "Firmicutes", "Bacteroidetes",
    "Firmicutes", "Bacteroidetes", "Firmicutes", "Firmicutes", "Bacteroidetes"))

})


test_that("getDominantTaxa", {

    # TSE object
    data(GlobalPatterns)
    tse <- GlobalPatterns

    testthat::expect_equal(mia::getDominantTaxa(tse)$dominant_taxa, c("331820",
    "549656", "550960", "319044", "114821", "12812",  "189047", "256977",
    "279599", "329744", "360229", "36155", "484436", "534609", "557211", "64396",
    "71074", "87194", "98605" ))

    testthat::expect_equal(mia::getDominantTaxa(tse, rank = "Kingdom", name = "dominant_kingdom")$dominant_kingdom,
    c("Bacteria"))

    testthat::expect_equal(mia::getDominantTaxa(tse, rank = "Order", name = "dominant_kingdom")$rel.freq,
                 c(23.1, 11.5, 7.7, 7.7, 7.7, 7.7, 3.8, 3.8, 3.8, 3.8, 3.8, 3.8, 3.8, 3.8, 3.8))

    testthat::expect_equal(mia::getDominantTaxa(tse, rank = "Class", name = "dominant_kingdom")$rel.freq.pct,
                 c("23%", "15%", "12%", "8%", "8%", "8%", "4%", "4%", "4%", "4%", "4%", "4%", "4%"))

    # TSE object
    tse <- microbiomeDataSets::dietswap()

    testthat::expect_equal(mia::getDominantTaxa(tse, rank = "Genus")$n,
                 c(104, 57, 43, 9, 2, 2, 2, 1, 1, 1))

    testthat::expect_equal(mia::getDominantTaxa(tse, group = "nationality")$nationality,
                 as.factor(c("AFR", "AAM", "AAM", "AFR", "AAM", "AAM", "AAM", "AAM", "AFR",
                   "AAM", "AAM", "AAM", "AFR", "AFR")))

    testthat::expect_equal(mia::getDominantTaxa(tse, rank = "Family", group = "sex")$dominant_taxa,
                 c("Bacteroidetes", "Bacteroidetes", "Clostridium cluster IV", "Clostridium cluster IV",
                   "Clostridium cluster XIVa", "Verrucomicrobia"))

    testthat::expect_equal(mia::getDominantTaxa(tse, rank = "Phylum", group = "timepoint", name = "test")$n,
                 c(24, 23, 20, 19, 19, 19, 19, 18, 17, 16, 15, 13))

    testthat::expect_equal(mia::getDominantTaxa(tse, rank = "Genus", group = "timepoint.within.group", name = "testing")$rel.freq,
                 c(46.4, 47.3, 27.3, 24.1, 20.5, 18.2, 4.5, 3.6, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9))

    testthat::expect_equal(mia::getDominantTaxa(tse, group = "bmi_group")$rel.freq.pct,
                 c("51%", "66%", "39%", "26%", "28%", "19%", "20%", "20%", "7%", "5%", "4%", "3%", "3%", "3%", "2%", "1%", "1%", "1%", "1%"))

})

