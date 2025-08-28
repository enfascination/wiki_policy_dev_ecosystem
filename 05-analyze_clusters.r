library(assertthat)

### case of k=2
if (FALSE) {
  print( "case of k=2")
  load("c2_results_002.Rdata")
  assignments <- c()
  lengths <- c()
  lengths_med <- c()
  lengths_min <- c()
  lengths_max <- c()
  lengths_2 <- c()
  lengths_2rat <- c()
  for (seq in seq_of_seqs_c2) {
    assignments <-  c(assignments, table(seq$assignment))
    c1 <- subset(seq, assignment == 1)
    c2 <- subset(seq, assignment == 2)
    lengths <- c(lengths, round(mean(c1$len_seq)), round(mean(c2$len_seq)))
    lengths_med <- c(lengths_med, median(c1$len_seq), median(c2$len_seq))
    lengths_min <- c(lengths_min, min(c1$len_seq), min(c2$len_seq))
    lengths_max <- c(lengths_max, max(c1$len_seq), max(c2$len_seq))
    lengths_2 <- c(lengths_2, sum(c1$len_seq == 2), sum(c2$len_seq == 2))
    lengths_2rat <- c(lengths_2rat, sum(c1$len_seq == 2)/nrow(c1), sum(c2$len_seq == 2)/nrow(c2))
  }
  tbl <- table(assignments)
  tbl <- tbl[as.integer(names(tbl)) <= nrow(seq)/2]
  print("distribution of the size of the small cluster")
  print(tbl)
  assert_that(100 == sum(tbl))
  plot(tbl)
  print("distribution of the sizes of sequences in the clusters (mean and median)")
  print(table(lengths))
  assert_that(200 == sum(table(lengths)))
  print(table(lengths_med))
  assert_that(200 == sum(table(lengths_med)))
  print("distributions of the min and maximum lengths of sequences in the clusters)")
  print(table(lengths_min))
  print(table(lengths_max))
  print("amounts of len=2 sequences in clusters (count and ratio)")
  print(table(lengths_2))
  plot(table(lengths_2))
  print(table(signif(lengths_2rat, 2)))
  plot(table(signif(lengths_2rat, 2)))
}

### case of k=3
if (FALSE) {
  print( "case of k=3")
  load("c3_results_003.Rdata")
  assignments1 <- c()
  assignments2 <- c()
  assignments3 <- c()
  lengths1 <- c()
  lengths2 <- c()
  lengths3 <- c()
  lengths_med1 <- c()
  lengths_med2 <- c()
  lengths_med3 <- c()
  for (seq in seq_of_seqs_c3) {
    b1 <- subset(seq, assignment == 1)
    b2 <- subset(seq, assignment == 2)
    b3 <- subset(seq, assignment == 3)
    c_largest <- which.max(c(nrow(b1), nrow(b2), nrow(b3)))
    if (c_largest == 1) {
      c1 <- b1
    } else if (c_largest == 2) {
      c1 <- b2
    } else if (c_largest == 3) {
      c1 <- b3
    }
    c_smallest <- which.min(c(nrow(b1), nrow(b2), nrow(b3)))
    assert_that( c_smallest != c_largest)
    if (c_smallest == 1) {
      c3 <- b1
    } else if (c_smallest == 2) {
      c3 <- b2
    } else if (c_smallest == 3) {
      c3 <- b3
    }
    c_middle <- setdiff(1:3, c(c_largest, c_smallest))
    if (c_middle == 1) {
      c2 <- b1
    } else if (c_middle == 2) {
      c2 <- b2
    } else if (c_middle == 3) {
      c2 <- b3
    }
    assert_that( 3 == length(unique( c(c_largest, c_smallest, c_middle))))
    assignments1 <-  c(assignments1, table(c1$assignment))
    assignments2 <-  c(assignments2, table(c2$assignment))
    assignments3 <-  c(assignments3, table(c3$assignment))
    lengths1 <- c(lengths1, round(mean(c1$len_seq)))
    lengths2 <- c(lengths2, round(mean(c2$len_seq)))
    lengths3 <- c(lengths3, round(mean(c3$len_seq)))
    lengths_med1 <- c(lengths_med1, median(c1$len_seq))
    lengths_med2 <- c(lengths_med2, median(c2$len_seq))
    lengths_med3 <- c(lengths_med3, median(c3$len_seq))
  }
  tbl1 <- table(assignments1)
  tbl2 <- table(assignments2)
  tbl3 <- table(assignments3)
  #tbl <- tbl[as.integer(names(tbl)) <= nrow(seq)/3]
  print("distribution of the size of the small cluster")
  print(tbl3)
  print(sum(tbl3))
  #assert_that(100 == sum(tbl))
  plot(tbl3)
  print("distribution of the sizes of sequences in the small cluster (mean)")
  print(table(lengths3))
  plot(table(lengths3))
  #assert_that(300 == sum(table(lengths)))
  #print(table(lengths_med))
  #assert_that(300 == sum(table(lengths_med)))
  print("")
  print("distribution of the size of the middle cluster")
  print(tbl2)
  print(sum(tbl2))
  plot(tbl2)
  print("distribution of the sizes of sequences in the middle cluster (mean)")
  print(table(lengths2))
  plot(table(lengths2))
  print("")
  print("distribution of the size of the large cluster")
  print(tbl1)
  print(sum(tbl1))
  plot(tbl1)
  print("distribution of the sizes of sequences in the large cluster (mean)")
  print(table(lengths1))
  plot(table(lengths1))
}

### case of k=4
if (TRUE) {
  print( "case of k=4")
  load("c4_results_001.Rdata")
  assignments1 <- c()
  assignments2 <- c()
  assignments3 <- c()
  assignments4 <- c()
  lengths1 <- c()
  lengths2 <- c()
  lengths3 <- c()
  lengths4 <- c()
  for (i in 1:length(seq_of_seqs_c4)) {
    print(i)
    seq <- seq_of_seqs_c4[[i]]
    b1 <- subset(seq, assignment == 1)
    b2 <- subset(seq, assignment == 2)
    b3 <- subset(seq, assignment == 3)
    b4 <- subset(seq, assignment == 4)
    c_idxs <- sort(c(nrow(b1), nrow(b2), nrow(b3), nrow(b4)), index.return=TRUE)$ix
    cc <- list(b1, b2, b3, b4)[c_idxs]
    assert_that(nrow(cc[[1]]) <=  nrow(cc[[2]]))
    assert_that(nrow(cc[[2]]) <=  nrow(cc[[3]]))
    assert_that(nrow(cc[[3]]) <=  nrow(cc[[4]]))
    c1 <- cc[[1]]
    c2 <- cc[[2]]
    c3 <- cc[[3]]
    c4 <- cc[[4]]
    assignments1 <-  c(assignments1, table(c1$assignment))
    assignments2 <-  c(assignments2, table(c2$assignment))
    assignments3 <-  c(assignments3, table(c3$assignment))
    assignments4 <-  c(assignments4, table(c4$assignment))
    lengths1 <- c(lengths1, round(mean(c1$len_seq)))
    lengths2 <- c(lengths2, round(mean(c2$len_seq)))
    lengths3 <- c(lengths3, round(mean(c3$len_seq)))
    lengths4 <- c(lengths4, round(mean(c4$len_seq)))
  }
  tbl1 <- table(assignments1)
  tbl2 <- table(assignments2)
  tbl3 <- table(assignments3)
  tbl4 <- table(assignments4)
  #tbl <- tbl[as.integer(names(tbl)) <= nrow(seq)/3]
  print("distribution of the size of the small cluster")
  print(tbl4)
  print(sum(tbl4))
  plot(tbl4)
  print("distribution of the sizes of sequences in the small cluster (mean)")
  print(table(lengths4))
  plot(table(lengths4))
  print("")
  print("distribution of the size of the middle small cluster")
  print(tbl3)
  print(sum(tbl3))
  plot(tbl3)
  print("distribution of the sizes of sequences in the middle small cluster (mean)")
  print(table(lengths3))
  plot(table(lengths3))
  print("")
  print("distribution of the size of the middle large cluster")
  print(tbl2)
  print(sum(tbl2))
  plot(tbl2)
  print("distribution of the sizes of sequences in the middle large cluster (mean)")
  print(table(lengths2))
  plot(table(lengths2))
  print("")
  print("distribution of the size of the large cluster")
  print(tbl1)
  print(sum(tbl1))
  plot(tbl1)
  print("distribution of the sizes of sequences in the large cluster (mean)")
  print(table(lengths1))
  plot(table(lengths1))
}

### reporting draft
if (FALSE) {
print("The procedure uses a simple hill climbing procedures to finds the best clustering from an initial random allocation of ...")
print("By this methodology we make the statistical argument for a single noisy path, rather than multiple paths.  Specifically, under this methodology, for values of k in {1,4}, we find that the statistical fit of k=1 (one cluster; no clustering) is significantly better than any clustering k>1, and that for all other clusterings investigated, all smaller values of k provide dramatically better statistical fits than any larger value.")
print("A descriptive analysis of the cluster contents reinforces our statistical conclusions.At k=2, all best clusterings of each run of the algorithm were different from each other, providing weak evidence for k=1. Stronger evidence came from what these distinct clusterings had in common: they all resulted in a main large cluster of about 200 sequences that includes all long sequences (all sequences longer than 5 or 6) and a smaller 'degenerate' or miscellaneous cluster of the poorest fit sequences. Most of the optimization process consisted of drawing long sequences from one of the randomly generated clusters to the other, until a small number of only short poorly-fit sequences remained in the smaller cluster.  Importantly, in all runs of the algorithm this smaller miscellaneious cluster was 75% 'sequences' of length 1.  Sequences of length 1 resulted from represent the 46 small wikis that have not yet developed past having a single shared policy. While it is theoretically possible that this second cluster of younger smaller wikis together follow their own policy development process distinct from the larger more developed wikis (an interpretation consistent with multiple viable policy development trajectories), the statistical comparison suggests that this pattern is an artifact of the optimization algorithm we developed, which proceeds from its random initialization towards putting all sequences into a single cluster, longer sequences first, until no single sequence from the miscellaneous cluster  incrementally improves the quality of the clustering.")
print("The evidence for one cluster continues when we interrogate the contents of the clusters for k=3 and k=4. The k=3 and k=4 results share with k=2 a small miscellaneous cluster composed mostly of single-token sequences, with the remaining clusters determined by splitting the larger cluster into 2 or 3 components.  But again, different runs of the algorithm divide the large cluster up differently depending on initial conditions, suggesting that there is no clear or consistent basis for justifying multiple paths. We established that these were inconsistent by measuring how often every pair of longer sequences ends up in the same cluster with each new run of the algorithm.  A clustering over these pairwise patterns in clusterings shows that independent runs of the algorithm at k=3 and k=4 fail to converge independently on similar partitionings, corroborating our conclusion that the ideal k=1, and therefore the population of wikis follows only one single developmental path. ")
      print("With this evidence at k=1-4 that smaller values of k are (astronomically) better fits, that values of k>1 put all short sequences together, and that values of k>2 fail to find consistent partitionings of the long sequences (independent runs of the algorithm produce very different partitionings), we concluded that we would not find any convincing evidence for k>4 paths outperforming the k=1 partitioning")
      print("Hence our ultimate conclusion, that the language editions of Wikipedia share a single developmental sequence as opposed to multiple alterantive developmental sequences.")
}


