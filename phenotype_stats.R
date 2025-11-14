# phenotype count table
phenotypes <- c("Cardiac", "Swim_Blader", "Cranio_facial",
                "Otolith", "Tail", "Eye")

totals <- c(ctrl = 44, capgb = 68, nefla = 72, rdh5 = 72, myh7 = 25)

counts <- matrix(c(
  # ctrl  capgb  nefla  rdh5  myh7
  3,     7,     8,    15,    1,   # Cardiac
  1,    19,    20,    31,   10,   # Swim Bladder
  0,     1,    12,    11,    1,   # Cranio-facial
  0,     0,     0,     1,    1,   # Otolith
  0,     3,     5,     2,    0,   # Tail
  0,     1,     0,     0,    0    # Eye
), nrow = 6, byrow = TRUE,
dimnames = list(phenotypes, names(totals)))

# Fisher test comparing each mutant vs WT
run_fisher <- function(pheno, mutant) {
  
  wt_yes  <- counts[pheno, "ctrl"]
  wt_no   <- totals["ctrl"] - wt_yes
  
  mut_yes <- counts[pheno, mutant]
  mut_no  <- totals[mutant] - mut_yes
  
  tbl <- matrix(c(wt_yes, wt_no, mut_yes, mut_no),
                nrow = 2, byrow = TRUE,
                dimnames = list(
                  c("WT","Mutant"),
                  c("Yes","No")
                ))
  
  out <- fisher.test(tbl)
  
  data.frame(
    phenotype = pheno,
    mutant = mutant,
    wt_yes = wt_yes,
    wt_no = wt_no,
    mut_yes = mut_yes,
    mut_no = mut_no,
    p_value = out$p.value,
    odds_ratio = out$estimate
  )
}

mutants <- c("capgb", "nefla", "rdh5", "myh7")

# apply to all phenotype × mutant combinations:
results <- do.call(rbind,
                   lapply(phenotypes, function(p)
                     do.call(rbind,
                             lapply(mutants, function(m)
                               run_fisher(p, m)
                             ))))

# Multiple testing correction (FDR)
results$FDR <- p.adjust(results$p_value, method = "BH")

# Add significance symbols
results$signif <- cut(results$FDR,
                      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                      labels = c("***", "**", "*", ""))

results
