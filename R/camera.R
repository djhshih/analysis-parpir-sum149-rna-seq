library(limma)

camera_single <- function(s, gsets=NULL, index=NULL) {
	if (is.null(index)) {
		index <- lapply(gsets,
			function(gset) {
				names(s) %in% gset
			}
		);
	}

	d <- cameraPR(s, index, sort=FALSE);
	d$delta <- unlist(lapply(index,
		function(i) {
			mean(s[i]) - mean(s[!i])
		}
	))
	d$z1 <- unlist(lapply(index,
		function(i) {
			mean(s[i])
		}
	));

	d
}

camera_batch <- function(y, gsets) {
	index <- lapply(gsets,
		function(gset) {
			rownames(y) %in% gset
		}
	);

	res <- lapply(1:ncol(y), function(j) {
		camera_single(y[,j], index=index)
	});
	names(res) <- colnames(y);

	res
}

camera_transform <- function(y, gsets) {
	index <- lapply(gsets,
		function(gset) {
			rownames(y) %in% gset
		}
	);

	res <- apply(y, 2, function(s) {
		d <- cameraPR(s, index, sort=FALSE);
		z <- log10(d$FDR);

		# flip the sign of enriched results
		idx <- d$Direction == "Up";
		z[idx] <- -z[idx]

		z
	});
	rownames(res) <- names(index);

	res
}

cam_volcano_plot <- function(cam, fdr.cut=0.01, z.cut=0.5, label.size=4, rename_gset=identity) {
	if (is.null(cam$gset)) {
		cam$gset = rownames(cam);
	}

	cam <- mutate(cam,
		keep = FDR < fdr.cut & abs(z1) > z.cut,
		group = factor(
			ifelse(keep, ifelse(delta < 0, "Down", "Up"), "NS"),
			c("Down", "Up", "NS")
		),
		gset = rename_gset(gset)
	);

	ggplot(cam, aes(x=z1, y=FDR, alpha=keep, colour=group, label=ifelse(keep, gset, NA))) +
		theme_clean() +
		geom_vline(xintercept=0) +
		geom_vline(xintercept=c(z.cut, -z.cut), linetype=3, colour="grey60") +
		geom_hline(yintercept=fdr.cut, linetype=3, colour="grey60") +
		geom_point(show.legend=FALSE) +
		geom_label_repel(show.legend=FALSE, nudge_y=0.1, nudge_x=0.1, size=label.size) +
		scale_y_continuous(trans=revlog_trans(10), sec.axis = dup_axis(name=NULL)) +
		scale_colour_manual(values=diff.group.cols) +
		xlab("mean expression difference") + ylab("false discovery rate")
}

plot_gene_set_density <- function(y, genes=NULL) {
	if (is.null(genes)) {
		y.sel <- y;
	} else {
		y.sel <- y[rownames(y) %in% genes, ];
	}
	y.sel.m <- melt(y.sel, varnames=c("gene", "comparison"));

	ggplot(y.sel.m, aes(x=value, fill=comparison)) + theme_clean() +
		geom_vline(xintercept=0, colour="grey60", linetype=2) +
		geom_density(alpha=0.5) +
		scale_fill_manual(values=clone.cols) +
		facet_wrap(~ comparison, ncol=1) +
		guides(fill=FALSE) +
		xlab("standardized difference vs. parental")
}

