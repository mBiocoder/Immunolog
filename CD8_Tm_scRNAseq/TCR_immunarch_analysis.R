library(dplyr)
library(immunarch)

file_path = "./raw_data/immunarch/"
immdata_10x <- repLoad(file_path)
names(immdata_10x)


#Get the most abundant clonotypes
top(immdata_10x$data[[1]])

#Filter functional / non-functional / in-frame / out-of-frame clonotypes
coding(immdata_10x$data[[1]]) #list of data frames with coding sequences
noncoding(immdata_10x$data[[1]]) #list of data frames with non-coding sequences
nrow(inframes(immdata_10x$data[[1]])) #number of filtered sequences is 1389
nrow(outofframes(immdata_10x$data[[1]])) #out-of-frame clonotypes is 0

#Basic analysis
exp_vol <- repExplore(immdata_10x$data, .method = "volume")
p1 <- vis(exp_vol, .by = c("Condition"), .meta = immdata_10x$meta)
p1
fixVis(p1)


exp_len <- repExplore(immdata_10x$data, .method = "len", .col = "aa")
exp_cnt <- repExplore(immdata_10x$data, .method = "count")

p2 <- vis(exp_len)
p3 <- vis(exp_cnt)

p2
p3


imm_pr <- repClonality(immdata_10x$data, .method = "clonal.prop")
imm_pr

imm_top <- repClonality(immdata_10x$data, .method = "top", .head = c(100, 1000, 3000))
imm_top

imm_rare <- repClonality(immdata_10x$data, .method = "rare", .bound = c(1, 3, 10, 30, 100))
imm_rare

imm_hom <- repClonality(immdata_10x$data,
                        .method = "homeo",
                        .clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1)
)
imm_hom


vis(imm_top) + vis(imm_top, .by = "Condition", .meta = immdata_10x$meta)
vis(imm_rare) + vis(imm_rare, .by = "Condition", .meta = immdata_10x$meta)
vis(imm_hom) + vis(imm_hom, .by = c("Condition"), .meta = immdata_10x$meta)


#Repertoire overlap
imm_ov1 <- repOverlap(immdata_10x$data, .method = "public", .verbose = F)
imm_ov2 <- repOverlap(immdata_10x$data, .method = "morisita", .verbose = F)

p1 <- vis(imm_ov1)
p2 <- vis(imm_ov2, .text.size = 2)

p1 + p2

imm_ov1 <- repOverlap(immdata_10x$data, .method = "overlap", .verbose = F)
imm_ov2 <- repOverlap(immdata_10x$data, .method = "jaccard", .verbose = F)

p1 <- vis(imm_ov1)
p2 <- vis(imm_ov2, .text.size = 2)

p1 + p2

imm_ov1 <- repOverlap(immdata_10x$data, .method = "tversky", .verbose = F)
imm_ov2 <- repOverlap(immdata_10x$data, .method = "cosine", .verbose = F)

p1 <- vis(imm_ov1)
p2 <- vis(imm_ov2, .text.size = 2)

p1 + p2




# Repertoire diversity

# Chao1 diversity measure
div_chao <- repDiversity(immdata_10x$data, "chao1")

# Hill numbers
div_hill <- repDiversity(immdata_10x$data, "hill")

# D50
div_d50 <- repDiversity(immdata_10x$data, "d50")

# Ecological diversity measure
div_div <- repDiversity(immdata_10x$data, "div")

# Gini Simpson index
gini_simp <- repDiversity(immdata_10x$data, "gini.simp")

# Gini index
gini <- repDiversity(immdata_10x$data, "gini")

# Inverse simpson index
inv_simp <- repDiversity(immdata_10x$data, "inv.simp")



p3 <- vis(div_chao, .by = c("Condition"), .meta = immdata_10x$meta)
p4 <- vis(div_d50, .by = "Condition", .meta = immdata_10x$meta)
p6 <- vis(div_div)
p7 <- vis(gini_simp, .by = "Condition", .meta = immdata_10x$meta)
p8 <- vis(gini, .by = "Condition", .meta = immdata_10x$meta)
p9 <- vis(inv_simp, .by = "Condition", .meta = immdata_10x$meta)

p3 + p4
p6 + p7
p9

# Hill numbers
repDiversity(
  .data = immdata_10x$data, .method = "hill", .max.q = 6,
  .min.q = 1, .do.norm = NA, .laplace = 0
) %>% vis()

#Rarefaction
imm_raref <- repDiversity(immdata_10x$data, "raref", .verbose = F)

p1 <- vis(imm_raref)
p1
