#!/usr/bin/env Rscript
# Level-3 (Q, and the default QCDPB betas): sesame vs R, given the same mask.
#
#   Rscript tests/compare_qcdpb.R <platform> <prefix> <c_qmask.txt> <c_betas.txt>
#
# R's own qualityMask cannot be the reference: on HM450/EPIC it still reads
# the older KYCG tracks, a different definition from the store's M_* tracks,
# and on EPICv2/MSA its mask object is an older build (NUMERICS.md, "Mask
# lineage"). So R is handed the CLI's Q set and must agree on two things:
#
#   (1) the mask: every probe the CLI masks exists in R's SigDF, and R's
#       qualityMask(mask = <CLI set>) masks exactly that set -- exact; and
#       R's own qualityMask differs from it by a pinned count (zero on all
#       four once R reads sesameData's 20261008 masks);
#   (2) the betas of the default pipeline, R prepSesame("QCDPB") with that
#       mask vs the CLI's QCDPB: the residual is D/P/B annotation lineage, and
#       it is PINNED per platform, not toleranced. The count of NA mismatches
#       must equal today's, the probes off by more than 1e-3 must be exactly
#       today's, and every other probe must stay under today's ceiling. A
#       change in the first two, either way, fails: re-measure and re-pin on
#       purpose.

suppressMessages(library(sesame))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4)
    stop("usage: compare_qcdpb.R <platform> <prefix> <c_qmask.txt> <c_betas.txt>")
platform <- args[1]; prefix <- args[2]

## Pins, measured 2026-10-12 on the InfiniumAnnotation v8.3 store, sesame
## 1.31.5, sesameData 1.29.10. namis = probes NA on one side only; big =
## probes with |diff| > 1e-3; ceiling = max |diff| over all the rest.
## v8.3 calls 63 probes Infinium-I that sesameData's older manifest calls
## Infinium-II (HM450 30, EPIC 32, MSA 1; DIVERGENCES.md C:pneg-design-type):
## the big sets below are exactly those probes, plus cg09334382 on EPIC.
pins <- list(
    HM450  = list(namis = 28, ceiling = 5.3e-4, big =
        c("cg02909136", "cg05618386", "cg06775939", "cg07434271",
          "cg09443652", "cg09823095", "cg11900657", "cg18477191",
          "cg20249919", "cg24746100", "cg25109431", "cg25967516",
          "cg27207470", "cg27505684")),
    ## cg09334382: Infinium-II in sesameData, Infinium-I in v8 (corrected
    ## Illumina manifest) -- the design-type case run_dyebiasL.sh describes
    EPIC   = list(namis = 92, ceiling = 9.9e-4, big =
        c("cg02909136", "cg05618386", "cg06775939", "cg07434271",
          "cg09334382", "cg09823095", "cg11900657", "cg15499467",
          "cg18477191", "cg20249919", "cg22409342", "cg24746100",
          "cg25109431", "cg25967516", "cg27207470", "cg27505684")),
    EPICv2 = list(namis = 0, big = character(0), ceiling = 1.5e-5),
    MSA    = list(namis = 20, big = "cg12810120_TC11", ceiling = 2.1e-4))

## R's OWN qualityMask vs the CLI's, pinned as c(C-only, R-only). Exact
## (0/0) on all four once R reads sesameData's KYCG.<PLAT>.Mask.20261008,
## which are built from the same v8.2 .cm (sesame >= 1.31.6, sesameData >=
## 1.31.1). Before that, EPICv2/MSA read older objects (pins below) and
## HM450/EPIC the older KYCG tracks, which are not compared.
new_masks <- sprintf("KYCG.%s.Mask.20261008", platform) %in%
    sesameData::sesameDataList()$Title
native_pins <- if (new_masks) {
    list(HM450 = c(0, 0), EPIC = c(0, 0), EPICv2 = c(0, 0), MSA = c(0, 0))
} else {
    list(EPICv2 = c(166, 1079), MSA = c(76, 9))
}

## With sesameData's v8.2-built objects (the KYCG masks above and the
## <PLAT>.address rebuilt from the store ordering), R reads exactly what the
## CLI reads and the betas agree to ~1e-15 on all four (measured 2026-10-08
## on a private install): no NA mismatch, nothing over 1e-12.
if (new_masks) pins <- lapply(pins, function(p)
    list(namis = 0, big = character(0), ceiling = 1e-12))
pin <- pins[[platform]]
if (is.null(pin)) stop(sprintf("compare_qcdpb.R has no pin for %s", platform))

cmask <- readLines(args[3])
cb <- read.table(args[4], colClasses = c("character", "numeric"))
cB <- setNames(cb$V2, cb$V1)
tag <- sprintf("%s/%s", platform, basename(prefix))

sdf <- suppressWarnings(readIDATpair(prefix, platform = platform))
fail <- character(0)

## (1) the mask
absent <- setdiff(cmask, sdf$Probe_ID)
if (length(absent))
    fail <- c(fail, sprintf("%d masked probes not in R's SigDF (e.g. %s)",
                            length(absent), absent[1]))
sq <- qualityMask(sdf, mask = cmask)
rmask <- sq$Probe_ID[sq$mask]
if (!setequal(rmask, cmask))
    fail <- c(fail, sprintf("mask differs: %d R-only, %d C-only",
                            length(setdiff(rmask, cmask)),
                            length(setdiff(cmask, rmask))))

## (1b) R's own mask, where R reads the same tracks
if (!is.null(native_pins[[platform]])) {
    nq <- qualityMask(sdf)
    nmask <- nq$Probe_ID[nq$mask]
    got <- c(length(setdiff(cmask, nmask)), length(setdiff(nmask, cmask)))
    cat(sprintf("%s native R qualityMask: %d probes, C-only %d, R-only %d\n",
                tag, length(nmask), got[1], got[2]))
    if (any(got != native_pins[[platform]]))
        fail <- c(fail, sprintf("native mask C-only/R-only %d/%d, pinned %d/%d",
                                got[1], got[2], native_pins[[platform]][1],
                                native_pins[[platform]][2]))
} else {
    cat(sprintf("%s native R qualityMask: not compared, sesame %s reads the %s\n",
                tag, packageVersion("sesame"), "older KYCG tracks"))
}

## (2) the betas, the rest of the default pipeline after that Q
rB <- getBetas(prepSesame(sq, "CDPB"))
ids <- intersect(names(rB), names(cB))
namis <- ids[xor(is.na(rB[ids]), is.na(cB[ids]))]
d <- abs(rB[ids] - cB[ids]); d <- d[!is.na(d)]
big <- sort(names(d)[d > 1e-3])
rest <- if (any(d <= 1e-3)) max(d[d <= 1e-3]) else 0

cat(sprintf(paste("%s QCDPB: Q %d probes (= R), n=%d, NA-mism=%d,",
                  ">1e-3: %d [%s], max of rest %.3e\n"),
            tag, length(cmask), length(ids), length(namis), length(big),
            paste(big, collapse = " "), rest))

if (is.na(pin$namis)) {
    fail <- c(fail, "no pin yet -- measured values above")
} else {
    if (length(namis) != pin$namis)
        fail <- c(fail, sprintf("NA mismatches %d, pinned %d",
                                length(namis), pin$namis))
    if (!identical(big, sort(pin$big)))
        fail <- c(fail, sprintf(">1e-3 set [%s], pinned [%s]",
                                paste(big, collapse = " "),
                                paste(pin$big, collapse = " ")))
    if (rest > pin$ceiling)
        fail <- c(fail, sprintf("max of rest %.3e over the pinned %.3e",
                                rest, pin$ceiling))
}

if (length(fail)) {
    cat(sprintf("FAIL %s: %s\n", tag, paste(fail, collapse = "; ")))
    quit(status = 1)
}
cat(sprintf("ok   %s: mask identical, betas at the pinned residual\n", tag))
