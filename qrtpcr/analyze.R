library(io)
library(dplyr)
library(ggplot2)
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


ref.target <- "18S";
ref.sample <- "SUM149-P";

# alpha level for confidence interval
alpha <- 0.2;     
z <- qnorm(1 - alpha/2)


# input data
x1 <- qread("data/qrtpcr_2021-01-20.csv");
x2 <- qread("data/qrtpcr_2021-01-25.csv");
out.fname <- filename("qrtpcr", tag="sum149");

# normalize experiments toward the last experiment

m1 <- x1 %>% select(target, sample, cq) %>% 
	group_by(target, sample) %>%
	summarize(cq=mean(cq));
	
m2 <- x2 %>% select(target, sample, cq) %>% 
	group_by(target, sample) %>%
	summarize(cq=mean(cq));

mc <- rbind(m1, m2) %>% group_by(target, sample) %>%
	summarize(dcq = cq - last(cq)) %>%
	filter(dcq > 0)

mct <- mc %>% group_by(target) %>% summarize(dcq_mean = mean(dcq), dcq_sd = sd(dcq));
#mcs <- mc %>% group_by(sample) %>% summarize(dcq_mean = mean(dcq), dcq_sd = sd(dcq));

summary(mc$dcq)

# substract background difference from experiment 1
x1b <- left_join(x1, mct) %>% mutate(cq = cq - dcq_mean) %>% 
	select(target, sample, cq);

# combined experiments
xcb <- rbind(x1b, select(x2, target, sample, cq)) %>% 
	group_by(target, sample) %>%
	filter(!is.na(cq))


# apply reference target and sample normalization by the delta-delta-cq method

# compute reference target mean csq
xcrt <- filter(xcb, target == ref.target) %>% 
	summarize(cq_ref = mean(cq)) %>% ungroup();

# normalize against reference target
xcn <- left_join(xcb, select(xcrt, -target), by="sample") %>%
	mutate(dcq = cq - cq_ref);

xcrs <- filter(xcn, sample == ref.sample) %>%
	summarize(dcq_ref = mean(dcq)) %>% ungroup();

xcnn <- left_join(xcn, select(xcrs, -sample), by="target") %>%
	mutate(ddcq = dcq - dcq_ref);

qwrite(xcnn, insert(out.fname, tag="ddcq", ext="tsv"));


s0 <- summarize(xcnn, ddcq_mean=mean(ddcq), ddcq_se=sd(ddcq)) %>%
	mutate(
		y = 2^(-ddcq_mean),
		ymin = 2^(-ddcq_mean - z*ddcq_se),
		ymax = 2^(-ddcq_mean + z*ddcq_se)
	);

qwrite(s0, insert(out.fname, tag="summary", ext="tsv"));


# compare all samples against the ref sample

targets <- unique(xcnn$target);
names(targets) <- targets;

# test results
hs <- lapply(
	targets,
	function(tt) {
		d <- filter(xcnn, target == tt);
		with(d, t.test(cq[sample != ref.sample], cq[sample == ref.sample]))
	}
);

# p values
ps <- unlist(lapply(hs, function(h) h$p.value));

# summary differential expression data
s <- filter(s0, target != ref.target) %>%
	mutate(sample = sub("SUM149-", "", sample)) %>%
	mutate(target = factor(target,
		levels=names(ps),
		labels=sprintf("%s (p = %s)", names(ps),
			# format number individually to avoid all numbers being shown in
			# scientific notation
			unlist(lapply(ps, format, digits=2)))
	));

# plot differential expression

integer_breaks <- function(lim) seq(floor(lim[1]), ceiling(lim[2]));
integer_limits <- function(lim) c(floor(lim[1]), ceiling(lim[2]));

qdraw(
	ggplot(s, aes(
		x=reorder(sample, desc(sample)), fill=sample,
		y=y, ymin=ymin, ymax=ymax
	)) +
		geom_hline(aes(yintercept=1), colour="grey80") +
		geom_col() + geom_errorbar(width=0.5) +
		theme_bw() + coord_flip() +
		facet_wrap(~ target, scale="free_x") +
		xlab("") + ylab("relative expression by qRT-PCR") +
		scale_fill_manual(values=clone.cols) +
		scale_y_continuous(breaks=integer_breaks, limits=integer_limits) +
		guides(fill = FALSE) +
		theme(
			panel.grid = element_blank(),
			strip.background = element_blank()
		)
	,
	file = insert(out.fname, ext="pdf"),
	width = 5, height = 2.5
)

