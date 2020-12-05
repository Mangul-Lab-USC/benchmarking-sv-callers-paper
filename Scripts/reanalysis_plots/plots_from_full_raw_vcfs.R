library(StructuralVariantAnnotation)
library(tidyverse)
theme_set(theme_bw())
setwd("~/benchmarking-sv-callers-paper") # set to project root
samples = c("AKR_J", "A_J", "BALB_cJ", "C3H_HeJ", "CBA_J", "DBA_2J", "LP_J")

# Load truth VCF
truth_bpgr = list()
for (sample in samples) {
	gr = breakpointRanges(readVcf(paste0("Data/gold_standard/mouse_vcf/", sample, "_reference.vcf")))
	gr$caller = "truth"
	gr$true_length = gr$svLen
	gr$tp = TRUE
	gr$tpdup = FALSE
	gr$sample = sample
	truth_bpgr[[sample]] = gr
}

# Load caller VCFs
# Note: vcf.gz and .bcf must be converted to vcf using gunzip or bcftools view .bcf > .vcf
caller_bpgr = lapply(samples, function(x) list())
names(caller_bpgr) = samples
for (sample in samples) {
	write(paste("Loading", sample), stderr())
	# Load GRIDSS
	gr = breakpointRanges(readVcf(paste0("Data/raw_data/mouse/raw_vcf/gridss/", sample, "/", sample, ".chr19.full.1/gridss.annotated.vcf")))
	gr$sample = sample
	gr$caller = "gridss"
	caller_bpgr[[sample]][["gridss"]] = gr
	
	if (sample != "AKR_J") { # AKR_J is a 0 byte file
		# Load lumpy
		# lumpy filenames are inconsistent do we need to pattern match
		gr = breakpointRanges(readVcf(list.files("Data/raw_data/mouse/raw_vcf/lumpy", full.names=TRUE, pattern=paste0("lumpexpress.", sample, ".*sorted.bam.modified.vcf$"))))
		gr$sample = sample
		gr$caller = "lumpy"
		caller_bpgr[[sample]][["lumpy"]] = gr
	}
	
	# Load delly
	gr = breakpointRanges(readVcf(paste0("Data/raw_data/mouse/raw_vcf/delly.mouse/full/delly.", sample, "_chr19_sorted.bam.output.bcf.vcf")))
	gr$sample = sample
	gr$caller = "delly"
	caller_bpgr[[sample]][["delly"]] = gr
	
	# Load breakdancer
	gr = breakpointRanges(readVcf(paste0("Data/raw_data/mouse/raw_vcf/breakdancer.mouse/full/", sample, "_chr19_sorted.bam.vcf")))
	gr$sample = sample
	gr$caller = "breakdancer"
	caller_bpgr[[sample]][["breakdancer"]] = gr
	
	# Load sniffles
	# TODO: Add support for sniffles' TRA records (do not conform to the VCF specifications) to the StructuralVariantAnnotation package
	# for this analysis, we can just delete the TRA records from the VCF since we're ignoring them anyway
	#gr = breakpointRanges(readVcf(paste0("Data/raw_data/mouse/raw_vcf/sniffels/full/", sample, "_chr19_sorted.sniffles_output.vcf")))
	#gr$sample = sample
	#gr$caller = "sniffles"
	#caller_bpgr[[sample]][["sniffles"]] = gr
}

annotate_calls = function(truth_gr, caller_gr, maxerrorbp=100) {
	caller_gr$QUAL[is.na(caller_gr$QUAL)] = 0
	hits = findBreakpointOverlaps(caller_gr, truth_gr, maxgap=maxerrorbp, restrictMarginToSizeMultiple=0.5)
	# Only the best QUAL TP is considered a TP, the remaining are considered FP
	hits = as.data.frame(hits) %>%
		mutate(QUAL=caller_gr$QUAL[queryHits]) %>%
		group_by(subjectHits) %>%
		arrange(desc(QUAL)) %>%
		mutate(isBestHit=row_number() == 1) %>%
		ungroup()
	besthits = hits %>% filter(isBestHit)
	otherhits = hits %>% filter(!isBestHit)
	caller_gr$tp = FALSE
	caller_gr$tp[besthits$queryHits] = TRUE
	caller_gr$tpdup = FALSE
	caller_gr$tpdup[otherhits$queryHits] = TRUE
	caller_gr$true_length = NA_integer_
	caller_gr$true_length[hits$queryHits] = truth_gr[hits$subjectHits]$svLen
	caller_gr$svtype = simpleEventType(caller_gr)
	return(caller_gr)
}
# pair-wise correlation
for (sample in samples) {
	for (caller1 in names(caller_bpgr[[sample]])) {
		for (caller2 in names(caller_bpgr[[sample]])) {
			if (!(sample == "AKR_J" & "lumpy" %in% c(caller1, caller2))) {
				mcols(caller_bpgr[[sample]][[caller1]])[paste0(caller2, "_hits")] = countBreakpointOverlaps(caller_bpgr[[sample]][[caller1]], caller_bpgr[[sample]][[caller2]], maxgap=100, restrictMarginToSizeMultiple=0.5)
			}
		}
	}
}

# truth annotation
for (sample in samples) {
	for (caller in names(caller_bpgr[[sample]])) {
		caller_bpgr[[sample]][[caller]] = annotate_calls(truth_bpgr[[sample]], caller_bpgr[[sample]][[caller]])
	}
}
allgr = c(unlist(GRangesList(lapply(caller_bpgr, function(x) unlist(GRangesList(x))))), unlist(GRangesList(truth_bpgr)))
deldf = as.data.frame(allgr) %>%
	filter(svtype == "DEL" & strand=="+") %>% # just take the lower breakend
	mutate(svLen=abs(svLen)) %>% # most callers don't follow the specs: DELs should have -ve length
	mutate(svLen > 50) %>%
	mutate(sizebin=ifelse(svLen < 100, "<100", ifelse(svLen <= 1000, "100-999", "1000+")))

ggplot(deldf %>%
			 	filter(tp & caller != "truth") %>%
			 	group_by(caller, sizebin) %>%
			 	mutate(mean_error=mean(abs(svLen) - abs(true_length)))) +
	aes(x=abs(svLen) - abs(true_length)) +
	facet_grid(sizebin ~ caller, scales="free") +
	geom_histogram() +
	geom_vline(aes(xintercept=mean_error), colour="blue", linetype="dotted") +
	scale_x_continuous(limits=c(-100, 100)) +
	labs(title="Alternative Figure 2: binned event size error distribution", x="error (bp)")
ggsave("alternate_figure2_size_binnedhistogram.pdf")

ggplot(deldf %>%
			 	filter(tp & caller != "truth") %>%
			 	group_by(caller) %>%
			 	mutate(mean_error=mean(abs(svLen) - abs(true_length)))) +
	aes(y=abs(svLen) - abs(true_length), x=caller) +
	geom_violin() +
	labs(title="Alternative Figure 2: event size error distribution", y="error (bp)")
ggsave("alternate_figure2_voilin.pdf")

rocby = function(df, ...) {
	groupingCols <- quos(...)
	df %>%
		replace_na(list(QUAL=0)) %>%
		mutate(QUAL=round(QUAL)) %>%
		group_by(caller, !!!groupingCols) %>%
		dplyr::select(QUAL, tp) %>%
		group_by(caller, !!!groupingCols, QUAL) %>%
		summarise(tp=sum(tp), ncalls=n()) %>%
		group_by(caller, !!!groupingCols) %>%
		arrange(desc(QUAL)) %>%
		mutate(cumtp=cumsum(tp), cumncalls=cumsum(ncalls)) %>%
		group_by(!!!groupingCols) %>%
		mutate(total_tp=max(cumtp)) %>% # truth set is in here so it'll be the max
		ungroup() %>%
		mutate(
			sensitivity=cumtp/total_tp,
			precision=cumtp/cumncalls)
}

overall_roc_df = bind_rows(
	rocby(deldf) %>% mutate(callset="All calls"),
	rocby(deldf %>% filter(FILTER %in% c("PASS", "."))) %>% mutate(callset="PASS only"))

ggplot(overall_roc_df) +
	aes(x=sensitivity, y=precision, color=caller, linetype=callset) +
	geom_line() +
	geom_point(data=overall_roc_df %>% group_by(caller, callset) %>% top_n(1, sensitivity), aes(shape=callset)) +
	labs(title="Precision recall with QUAL and FILTER")
ggsave("roc_overall.pdf")

ggplot(rocby(deldf, sample)) +
	aes(x=sensitivity, y=precision, color=caller) +
	geom_point() +
	facet_wrap(~ sample) +
	labs(title="Precision recall by sample")
ggsave("roc_by_sample.pdf")

ggplot(rocby(deldf, sizebin)) +
	aes(x=sensitivity, y=precision, color=caller) +
	geom_point() +
	facet_wrap(~ sizebin) +
	labs(title="Precision recall by size bin")
ggsave("roc_by_sizebin.pdf")


deldf %>%
	group_by(caller, sizebin) %>%
	summarise(n=n(), delly_hits=sum(delly_hits > 0, na.rm=TRUE))






























