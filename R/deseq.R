# Summarize DESeq2Results with most extreme transcript
# @param res  a \code{DESeq2Results} object
# @param fannot  feature annotation in same order as \code{res}
deseq_summarize_extreme <- function(res, fannot) {
	s <- data.frame(
		gene_name = fannot$gene_name,
		transcript = fannot$ensembl_transcript,
		stat = res$stat
	);

	s.max <- group_by(s, gene_name) %>%
		summarize(
			idx = which.max(abs(stat)),
			stat = stat[idx],
			transcript = transcript[idx]
		) %>% ungroup();

	res.d <- data.frame(
		transcript = fannot$ensembl_transcript,
		base_mean = res$baseMean,
		delta = res$log2FoldChange,
		p = res$pvalue,
		q = res$padj
	);

	left_join(s.max, res.d)
}

deseq_annotate_significant <- function(des, fdr.cut, delta.cut, n=NULL) {
	des <- mutate(des,
		keep = q < fdr.cut & abs(delta) > delta.cut,
	);

	if (!is.null(n)) {
		top.genes <- with(filter(des, keep), gene_name[order(q)][1:n]);
		label <- with(des, ifelse(gene_name %in% top.genes, as.character(gene_name), NA));
	} else {
		label <- des$gene_name;
	}

	mutate(des,
		label = label,
		direction = factor(
			ifelse(keep, ifelse(delta < 0, "Down", "Up"), "NS"),
			c("Down", "Up", "NS")
		)
	)
}

deseq_volcano_plot <- function(des, fdr.cut=0.01, delta.cut=1, n=50) {
	des.a <- deseq_annotate_significant(des, fdr.cut, delta.cut, n);

	ggplot(des.a,
		aes(
			x=delta, y=q, colour=direction,
			alpha=ifelse(keep, 0.7, 0.2),
			label=label
		)
	) + theme_clean() +
		geom_vline(xintercept=0) +
		geom_vline(xintercept=c(delta.cut, -delta.cut), linetype=3, colour="grey60") +
		geom_hline(yintercept=fdr.cut, linetype=3, colour="grey60") +
		geom_point(show.legend=FALSE) +
		geom_label_repel(show.legend=FALSE) +
		scale_y_continuous(trans=revlog_trans(10), sec.axis = dup_axis(name=NULL)) +
		scale_colour_manual(values=diff.group.cols) +
		xlab("log2 fold change") + ylab("false discovery rate")
}

deseq_plot_stat <- function(d) {
	ggdotchart(
		d,
		x = "gene_name", 
		y = "stat",
		color = "direction",
		palette = diff.group.cols,
		add = "segments",
		sorting = "descending",
		rotate = TRUE
	) + geom_hline(yintercept=0, color="grey90") +
		xlab("") + ylab("Wald statistic") + guides(color=FALSE)
}

