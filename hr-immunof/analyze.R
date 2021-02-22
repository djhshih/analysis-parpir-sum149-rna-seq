library(io)
library(dplyr)
library(ggplot2)
library(binom)
library(ggsci)


shift <- function(x, n=1) {
		if (n == 0) {
					x
	} else {
				c(tail(x, -n), head(x, n))
		}
}

clones <- c("P", "RC1", "RC2", "RC3", "RC4", "RC5", "RC6", "RC7");
clone.cols <- shift(pal_d3()(length(clones)), -1);
names(clone.cols) <- clones;

sample.cols <- c("#444444FF", clone.cols);
names(sample.cols) <- c("MCF10A", paste0("SUM149-", clones));


out.fname <- filename("hr-immunof");
xs <- qread("data", pattern="\\.csv");

ss <- strsplit(names(xs), "_");
samples <- toupper(unlist(lapply(ss, function(s) s[2])));
targets <- sub("GAMMA", "g", toupper(unlist(lapply(ss, function(s) s[3]))));

d <- do.call(rbind, mapply(
	function(i, sample_id, tt) {
		ss <- strsplit(xs[[i]]$filename, "-");
		mutate(xs[[i]],
			sample = sample_id,
			target = tt,
			channel = unlist(lapply(ss, function(s) sub("C", "", s[1]))),
			image = unlist(lapply(ss, function(s) sub("Image", "", s[2]))),
		) %>% select(-filename)
	},
	names(xs),
	samples,
	targets,
	SIMPLIFY=FALSE
));
rownames(d) <- NULL;

assign_contrast <- function(samples) {
	ifelse(
		grepl("SUM149", samples),
		ifelse(samples == "SUM149-P", 0, 1),
		NA
	)
}

d <- mutate(d, contrast = assign_contrast(sample));

theme_clean_x45 <- function() {
	theme_bw() +
	theme(
		panel.grid = element_blank(),
		strip.background = element_blank(),
		axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
	)
}


target.vals <- unique(targets);
names(target.vals) <- target.vals;
hs <- lapply(target.vals,
	function(tt) {
		wilcox.test(
			n_nuclear_blobs ~ contrast,
			filter(d, target == tt)
		)
	}
);

hs.t <- lapply(target.vals,
	function(tt) {
		t.test(
			n_nuclear_blobs ~ contrast,
			filter(d, target == tt)
		)
	}
);

ps <- unlist(lapply(hs, function(h) h$p.value));

qdraw(
	ggplot(d, aes(x = sample, y = n_nuclear_blobs)) +
		theme_clean_x45() +
		geom_boxplot() +
		geom_jitter(width=0.2, alpha=0.2, aes(colour = sample))  +
		scale_colour_manual(values=sample.cols) + guides(colour=FALSE) +
		scale_y_log10() +
		facet_wrap(
			~ factor(
				target,
				levels = names(ps),
				labels = sprintf("%s (p = %s)", names(ps),
				unlist(lapply(ps, format, digits=2)))
			),
			ncol=1
		) +
		ylab("# nuclear foci") + xlab("")
	,
	width = 2.5, height = 6,
	file = insert(out.fname, tag="n-foci", ext="pdf")
)


n.cut <- 10;

pos.d0 <- group_by(d, sample, target) %>%
	summarize(
		positivity = sum(n_nuclear_blobs >= n.cut),
		n_cells = n(),
	) %>%
	mutate(target = paste0(target, "+"));

filter(d, sample == "SUM149-P", image == 18, cell == 1)

res.d <- group_by(d, sample, image, cell) %>%
	summarize(
		positive =
			n_nuclear_blobs[target == "g-H2AX"] >= n.cut &
			n_nuclear_blobs[target == "RAD51"] < n.cut
	);

res.pos.d <- group_by(res.d, sample) %>% 
	summarize(
		target = "g-H2AX+ RAD51-",
		positivity = sum(positive),
		n_cells = n()
	);

pos.d <- rbind(
	pos.d0,
	res.pos.d
);


cf <- with(pos.d, binom.confint(positivity, n_cells, method="exact"));
pos.d <- cbind(pos.d, select(cf, mean, lower, upper));

dd <- filter(pos.d) %>%
	mutate(contrast = assign_contrast(sample)) %>%
	filter(!is.na(contrast)) %>%
	group_by(contrast, target) %>%
	summarize(x1 = sum(positivity), x0 = sum(n_cells) - x1) %>%
	ungroup();

# compare parental vs. resistant clones in terms of the fraction of
# foci-positive nuclei
foci.vals <- c("g-H2AX+", "RAD51+", "g-H2AX+ RAD51-");
names(foci.vals) <- foci.vals;
hs.pos <- lapply(foci.vals,
	function(tt) {
		ct <- filter(dd, target == tt) %>% select(x0, x1) %>% as.matrix();
		fisher.test(ct)
	}
);

ps.pos <- lapply(hs.pos, function(h) h$p.value);

qdraw(
	ggplot(pos.d, aes(x = sample, y = mean, ymin = lower, ymax = upper, fill = sample)) +
		theme_clean_x45() +
		geom_col() +
		geom_errorbar(width=0.5) +
		scale_fill_manual(values=sample.cols) + guides(fill=FALSE) +
		facet_wrap(
			~ factor(
				target,
				levels = names(ps.pos),
				labels = sprintf("%s (p = %s)", names(ps.pos),
				unlist(lapply(ps.pos, format, digits=2)))
			),
			ncol=1
		) +
		ylim(0, 1) +
		ylab("fraction of nuclei") + xlab("")
	,
	width = 2.5, height = 6,
	file = insert(out.fname, tag="fraction", ext="pdf")
)

