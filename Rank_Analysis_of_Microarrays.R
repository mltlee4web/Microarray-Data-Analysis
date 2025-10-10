##################################
#
# RANK ANALYSIS OF MICROARRAYS (RAM)
#
# The following two examples were given in in Lee, Gray, Bjorkbacka, Freeman.
# Generalized Rank Tests for Replicated Microarray Data, 
# Statistical Applications in Genetics and Molecular Biology,(2005),
# Vol. 4: No. 1, Article 3. http://www.bepress.com/sagmb/vol4/iss1/art3
#   
# This program was written by Dr. Robert J. Gray.
# User documentation was provided by Dr. Weiliang Qiu.
# ##################################

# For Unix version of 'permax',
# 'permax' should be installed in the home directory.

library(permax, lib.loc="~/.")
options(width=90)

##########################
# Example 1: One Treatment Applied to 4 Wild-type and 4 Knock-out mice
# (Two arrays made for each mouse)
#########################
# Read input data
# The data file 'data4ranksum.dat' has 3339 rows and 17 columns.
# The 1st row contains variable name. The 1st column contains gene id.
# 'data1' is a 3338 by 16 data frame.
# 'header=T" indates that the first row contains variable names.
# The rows correspond to genes that have complete observations from all 8 mice.
# The columns correspond to mice and their replicates.

data1 <- read.table("data4ranksum.dat", header=T, row.names=1)

# Wild type mouse 1:4
# 'grp' indicates the column numbers in 'data1' corresponding to group 1
# (i.e. wild type mice 1, 2, 3, 4).
# The 1st replicates of these 4 mice correspond to columns 1, 2, 3, 4 of
# 'data1'; The 2nd replicates correspond to columns 9, 10, 11, 12 of 'data1'.
# 1:4 means 1, 2, 3, 4
# 9:12 means 9, 10, 11, 12
# 'grp' is (1,2,3,4,9,10,11,12)

grp <- c(1:4,9:12) 

# cluster=cy3,cy5 from each mouse
# 'clust' is a vector of cluster membership indicators for the columns 
# of 'data1'
# That is, columns 1 and 9 are from cluster 1 
#          columns 2 and 10 are from cluster 2 
#          ................................
#          columns 8 and 16 are from cluster 8 
# 'clust' is (1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8)

clust <- c(1:8,1:8)

# call the function permax()
#
# for each gene (each with two replicates), calculate the generalized ranksum 
# statistic by using all possible combinations of permuations.
#
# mouse  W1  W2  W3  W4  K1  K2  K3  K4
#repeat1 X11 X21 X31 X41 Y11 Y21 Y31 Y41
#repeat2 X12 X22 X32 X42 Y12 Y22 Y32 Y42
#
# Wi means wild-tye mouse i, i=1,2,3,4
# Kj means knock-out mouse j, j=1,2,3,4
#
# e.g. rank(X11) is among all 16 observations: X11, X21, X31, X41, Y11, Y21,
# Y31, Y41, X12, X22, X32, X42, Y12, Y22, Y32, Y42.
#
# We can randomly assign 4 mice as wild-type 
# and the remaining 4 mice as knock-out.
# There are a total of 70 (8 choose 4) possible permutations.
#
# The 16 columns of 'data1' correspond to
# X11 X21 X31 X41 Y11 Y21 Y31 Y41 X12 X22 X32 X42 Y12 Y22 Y32 Y42
#
# nperm=0 means all possible combinations are enumerated
#
# nl: A lower tail critical value is determined using the 'nl'th
#      largest value of the test statistics.  The number of genes
#      with statistics <= this critical value is determined for each
#      combination.  
# nr: An upper tail critical value is determined using the
#     'n-nr+1'st largest value of the test statistics, where 'n' is
#     the number of genes.  The number of genes with statistics >=
#     this critical value is determined for each combination.  
#
# 'cluster' is specified, 'stratify=FALSE', 'permute.cluster=TRUE'
# it is assumed that group or treatment is defined at the cluster level
# (that is, all members of each cluster are in the same group). In this case
# statistics are computed the same as in the independence case, but only
# whole clusters are permuted to preserve the within cluster dependence.

res1 <- permax(data=data1,ig1=grp,nperm=0, ranks=TRUE,cluster=clust,
               stratify=FALSE, permute.cluster=TRUE, nl=2,nr=5,expord=TRUE)

summary(res1)

# Plots statistics vs expected order statisics computed under the
#   permutation distribution for output from permax
# del: if 'del>0', then lines with slope=1 and intercepts '+del' and
#   '-del' will be added to the plot

plot.expord(x=res1,del=5)

######
# no ranks ('rank=FALSE'), i.e. use T-statistic
#####

res11 <- permax(data=data1,ig1=grp,nperm=0, ranks=FALSE,cluster=clust,
                permute.cluster=TRUE,nl=6,nr=23,expord=TRUE)
plot.expord(x=res11)


####################################
# Example 2: Two Treatments Applied to the Same Set of Mice
####################################

# read data
# The data file 'data4signedrank.dat' has 1324 rows and 25 columns.
# The 1st column is gene id.
# 'data2' is a 1324 by 24 data frame
# The rows correspond to genes.
# The columns correspond to intensity ratios of 6 mice and their replicates
# 'header=T" means the first row is variable names.
# 'row.names=1' means the first column of the file 'data4signedrank.dat' is
# row names.

data2<-read.table("data4signedrank.dat", header=T, row.names=1)

# 'clust12' is a vector of cluster membership indicators for the columns 
# of 'data2'
# That is, columns 1 and 7 are from cluster 1
#          columns 2 and 8 are from cluster 2
#          ................................
#          columns 6 and 12 are from cluster 6
#          columns 13 and 19 are from cluster 7
#          columns 14 and 20 are from cluster 8
#          ................................
#          columns 18 and 24 are from cluster 12
# 'clust' is (1,2,3,4,5,6,1,2,3,4,5,6,7,8,9,10,11,12,7,8,9,10,11,12)

clust12 <- c(1:6,1:6,7:12,7:12)

options(width=95)

# 'grp' indicates the column numbers in 'data2' corresponding to group 1
# (i.e. mice with treatment A).
# The 1st replicates of these 6 mice correspond to columns 1,...,6 of
# 'data2'; The 2nd replicates correspond to columns 13,...,18 of 'data2'.
# 'grp' is (1,2,3,4,5,6,13,14,15,16,17,18)

grp <- c(1:6,13:18) 


###
# RS1: The Generalized Ranksum Test (All 24 Ratios Ranked Together)
###
#        m1  m2  m3  m4  m5  m6  m1  m2  m3  m4  m5  m6
#repeat1 X11 X21 X31 X41 X51 X61 Y11 Y21 Y31 Y41 Y51 Y61
#repeat2 X12 X22 X32 X42 X52 X62 Y12 Y22 Y32 Y42 Y52 Y62
#
# We rank Xij, Yij, i=1,...,6, j=1,...,6 together (24 ranks)
# e.g. rank(X11) is the rank among X11, X21, X31,X41,X51, X61, Y11, Y21,
# Y31, Y41, Y51, Y61, X12, X22, X32, X42, X52, X62, Y12, Y22, Y32, Y42,
# Y52, Y62.
#
# For each pair (X11, Y11), (X21, Y21), (X31, Y31), (X41, Y41), (X51, Y51),
# (X61, Y61), (X12, Y12), (X22, Y22), (X32, Y32), (X42, Y42), (X52, Y52),
# (X62, Y62), we have two possible permutations.
# e.g.
# Permuation1: X11 is the score from treatment A, 
#              Y11 is the score from treatment B
# Permuation2: Y11 is the score from treatment A, 
#              X11 is the score from treatment B
# There are totally 2^12=4096 possible permutations.
# cluster on mouse-cy combination, using ranks
#
# 
# The 24 columns of 'data2' correspond to
# X11 X21 X31 X41 X51 X61 Y11 Y21 Y31 Y41 Y51 Y61 X12 X22 X32 X42 X52 X62 
# Y12 Y22 Y32 Y42 Y52 Y62
# 
# 'cluster' is specified, 'stratify=FALSE', 'permute.cluster=FALSE'
# statistics are computed as in the independence case, but columns are only
# permuted within clusters, since under the null it is assumed that
# columns from the same cluster are exchangeable, while columns from
# different clusters may not be.

res2.RS1 <- permax(data=data2,ig1=grp,nperm=0,min.np=0,logs=TRUE,
                   ranks=TRUE, cluster=clust12,stratify=FALSE,signed.rank=FALSE,
                   nl=38,nr=18, expord=TRUE, permute.cluster=FALSE)
summary(object=res2.RS1,nl=38,nr=18)


#####
# ST: Signed Rank Test for Blocked Comparison of Two Treatments
#####
#         mouse1     mouse2     mouse3     mouse4     mouse5      mouse6
#repeat1 (X11, Y11) (X21, Y21) (X31, Y31) (X41, Y41) (X51, Y51) (X61, Y61)
#repeat2 (X12, Y12) (X22, Y22) (X32, Y32) (X42, Y42) (X52, Y52) (X62, Y62)
#
# sr1=sign(Y11-X11)*rank(|Y11-X11|), sr2=sign(Y21-X21)*rank(|Y21-X21|), ...
# sr12=sign(|Y62-X62|)*rank(|Y62-X62|)
# 
# Analysis is based on the signed ranks sr1, sr2, ..., sr12
#
# For each pair (Xij, Xij), there are two possible permutations.
# Permuation1: Xij is the score for treament A, Yij is the score for teatment B
# Permuation2: Yij is the score for treament A, Xij is the score for teatment B
# There are totally 2^12=4096 possible permutations.
# signed-rank test 
#
#
# 'cluster' is specified, 'stratify=FALSE', 'permute.cluster=FALSE'
# statistics are computed as in the independence case, but columns are only
# permuted within clusters, since under the null it is assumed that
# columns from the same cluster are exchangeable, while columns from
# different clusters may not be.
#
# 'signed.rank=TRUE' indicates to use signed rank test

res2.ST <- permax(data=data2,ig1=grp,nperm=0,min.np=0,logs=TRUE,
                  ranks=TRUE, cluster=clust12,stratify=FALSE,signed.rank=TRUE,
                  nl=34,nr=21, expord=TRUE, permute.cluster=FALSE)
summary(object=res2.ST, nl=34, nr=21)

#####
# RS2: The Generalized Ranksum Test (Separate Ranking of Two Sets of 12
#      Intensity Ratios)
#####
#        m1  m2  m3  m4  m5  m6  m1  m2  m3  m4  m5  m6
#repeat1 X11 X21 X31 X41 X51 X61 Y11 Y21 Y31 Y41 Y51 Y61
#repeat2 X12 X22 X32 X42 X52 X62 Y12 Y22 Y32 Y42 Y52 Y62
#
# We rank Xi1, Yj1, i=1,...,6, j=1,...,6 together (12 ranks)
# e.g. rank(X11) is the rank among X11, X21, X31,X41,X51, X61, Y11, Y21,
# Y31, Y41, Y51, Y61.
#
# We rank Xi2, Yi2, i=1,...,6, j=1,...,6 together (12 ranks)
# e.g. rank(X12) is the rank among X12, X22, X32,X42,X52, X62, Y12, Y22,
# Y32, Y42, Y52, Y62.
#
# For each pair (X11, Y11), (X21, Y21), (X31, Y31), (X41, Y41), (X51, Y51),
# (X61, Y61), (X12, Y12), (X22, Y22), (X32, Y32), (X42, Y42), (X52, Y52),
# (X62, Y62), we have two possible permutations.
# e.g.
# Permuation1: X11 is the score from treatment A, 
#              Y11 is the score from treatment B
# Permuation2: Y11 is the score from treatment A, 
#              X11 is the score from treatment B
# There are a total of 2^12=4096 possible permutations.
#
# cluster on mouse-cy combination, separately ranking cy3 and cy5 values
# 'data2[,1:12]' is the first 12 columns of the data frame 'data2'
# 'data2[,13:24]' is a matrix consists of the columns 13, 14,...,24 of the
# data frame 'data2'
#
# 'data3' is a 1324 by 24 matrix. rows correspond to genes.
# Elements of 'data3' are ranks.
# e.g. 'data3[1,2]' is the rank of X21 for gene 1 among the 12 scores 
# X11, X21, X31, X41, X51, X61, Y11, Y21, Y31, Y41, Y51, Y61.
# 'data[3, 14]' is the rank of X12 for gene 3 among the 12 genes
# X12, X22, X32, X42, X52, X62, Y12, Y22, Y32, Y42, Y52, Y62.
#
# The command 'apply(data2[,1:12], 1, rank)' means that we rank elements of
# each row of the matrix 'data2[,1:12]'.
# The function 't(A)' means we get the transpose of the matrix 'A'.
# The command 'apply(data2[,1:12], 1, rank)' will produce 12 by 1324 matrix
# instead of 1324 by 12 matrix.

data3 <- cbind(t(apply(data2[,1:12],1,rank)),t(apply(data2[,13:24],1,rank)))

dimnames(data3) <- dimnames(data2)

# 'cluster' is specified, 'stratify=FALSE', 'permute.cluster=FALSE'
# statistics are computed as in the independence case, but columns are only
# permuted within clusters, since under the null hypothesis it is assumed that
# columns from the same cluster are exchangeable, while columns from
# different clusters may not be.

res2.RS2 <- permax(data=data3,ig1=grp,nperm=0,min.np=0,logs=FALSE,
                   ranks=TRUE, cluster=clust12,stratify=FALSE,
                   signed.rank=FALSE,nl=33,nr=18, expord=TRUE,
                   permute.cluster=FALSE)

summary(object=res2.RS2,nl=33,nr=18)