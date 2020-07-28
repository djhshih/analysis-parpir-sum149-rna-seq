library(io)

# Given model
# E[y] = a x + b z + c x:z + d
# 
#

# interaction effect
interact <- qread("parpi-resist_deseq-stat_treatment-clone-interaction_interactions.mtx");

# treatment effect
treatment <- qread("parpi-resist_treatment-clone-interaction_treatment.rnk");

treatment <-  treatment[match(rownames(interact), names(treatment))];

stopifnot(all(names(treatment) == rownames(interact), na.rm=TRUE))

clone.treatment <- interact + treatment;

out.fname <- filename("parpi-resist", tag=c("deseq-stat", "treatment-clone-interaction"), ext="mtx");

qwrite(clone.treatment, insert(out.fname, c("clone", "treated-vs-untreated")));


# clone effect
clone <- qread("parpi-resist_deseq-stat_treatment-clone-interaction_clones-vs-parental.mtx");

stopifnot(all(rownames(interact) == rownames(clone)))

treatment.clone <- interact + clone;

qwrite(treatment.clone, insert(out.fname, c("treatment", "clones-vs-parental")));

