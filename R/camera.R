library(limma)
library(ggplot2)
library(ggrepel)

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

cam_volcano_plot <- function(cam, fdr.cut=0.01, z.cut=0.5, nudge=0.1, label.size=4, down.sample=NA, rename_gset=identity, show.legend=FALSE, two.sided=TRUE, ylim=NULL) {
	if (is.null(cam$gset)) {
		cam$gset = rownames(cam);
	}

	cam <- mutate(cam,
		gset = rename_gset(gset)
	);

	if (is.null(cam$keep)) {
		cam <- mutate(cam, keep = FDR < fdr.cut & abs(z1) > z.cut);
	}

	if (is.null(cam$group)) {
		cam <- mutate(cam,
			group = factor(
				ifelse(keep, ifelse(z1 < 0, "Down", "Up"), "NS"),
				c("Down", "Up", "NS")
			)
		);
	}

	if (!is.na(down.sample)) {
		# hide a set of insigificant data points
		idx <- which(cam$FDR >= fdr.cut);
		# higher FDR => higher chance of being omitted
		omit <- sample(idx, floor(down.sample * length(idx)), prob = cam$FDR[idx]);
		cam$FDR[omit] <- NA;
	}

	g <- ggplot(cam,
		aes(
			x=z1, y=FDR, colour=group,
			alpha = ifelse(keep, 0.7, 0.2),
			label = ifelse(keep, gset, NA)
		)
	) + theme_clean();

	if (two.sided) {
		g <- g + geom_vline(xintercept=0) + 
			scale_y_continuous(trans=revlog_trans(10), limits=ylim, sec.axis = dup_axis(name=NULL));
	} else {
		g <- g + scale_y_continuous(trans=revlog_trans(10), limits=ylim);
	}

	g <- g + geom_vline(xintercept=c(z.cut, -z.cut), linetype=3, colour="grey60") +
		geom_hline(yintercept=fdr.cut, linetype=3, colour="grey60") +
		geom_point(show.legend=show.legend) +
		geom_label_repel(show.legend=FALSE, nudge_y=nudge, nudge_x=nudge, size=label.size) +
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

