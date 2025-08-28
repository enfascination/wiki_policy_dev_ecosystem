# We've already got a strong hunch that there's only a single generating process underlying all sequences, but we don't know the right statistical tool for disqualifying the hypotheses of >1 path. Well, we have one idea: we could take a model fitting approach, define n generating processes (for the n hypotheses that 1--n paths exist), calculate likelihoods, optimize parameters, and generate test statistics.

#I have  independent observations of 200 systems as they progress through up to 60 states (most have only progressed through 10). I want to test the competing hypotheses that their progress is driven by 1 versus 2 versus n markov chains ("generating processes").  I'm conceptualizing a process as a markov model with a strong disposition to go through stages linearly: at any stage a system can transition to any stage but depending on its generating process there is a single unique stage that it tends to transition through, until it has been through all 60. Under that framework it is very straightforward to use data to propose the optimal generating process under n=1: just produce an empirical transition matrix and that's a candidate optimal generating process with a calculable log likelihood.

#To produce two I can imagine a two step optimization process that (somehow) divides the 200 systems into 2 groups that are mostly likely to be driven by distinct generating processes, generates those processes, calculates their log likelihoods of producing their own group of systems, then returns to tweak the assignment of systems to those groups (somehow), recalculates the two generating processes, recalculates the log likelihood, and so on, until likelihood is optimized.

#And to produce n generating processes is just like producing two.

#With this there is an optimal LL associated with each value of n from 1 to 60 (won't theorize more groups than states) and the minimum of those is the basis for a statistical claim for the most likely number of underlying generating processes.

#My main question: does this exist already?  Keep in mind that most work on markov-y sequence analysis, including HMMs, assumes a smaller number of states, a larger number of sequences, longer sequences in which states can be revisited, and complete sequences (many of our systems are just partially on their developmental path through the 60 states). So I'm concerned that it doesn't exist even though it seems like it should.

#My secondary question: I'm stuck on implementing this myself because I'm stuck on an algorithm for proposing changes in group assignment membership that could be plugged into an optimizer like R's `optim()`. Should I just evaluate all 200 possible incremental changes in membership (for the two group case) and pick the one that best minimizes LL, or find the member of 200 whose sequence is most unlike its group mates (how)

# [X] 1) load data
# [X] 2) sort data
# [X] 3) compare two sequences
# [X] 4) generate random seuqnecs andcompare them
# 5) get optim minimally working for anything
# 6) build a sequence tweaker
# 7) plug it all together
# 1) get a pile of sequences
#    (function: take number of culstuers and assign sequences to them evenly)
# 2) generate a 60x60 transition matrix from that pile 
#    * (function: model from pile of sequences)
#      * model from empty pile is uniform
#      * diagonals of these models are 0
#      * not symmetric
#      * unpopulated cells are 0.000001
#      * has an initial state - that you can't transition to
# 3) DONE: calculate the log likelihood of all sequences in the pile from their own transition matrix
#    * df = (60^2) - 60 - 1 for each transition matrix/cluster/model, where 60 is the number of rows/columns in that matrix
#      * be careful if you have an initial state: df = 60^2 - 60 - 1 + 59 because of the probabilities of transitioning from initial state into anything else
#    * (function: likelihood of sequence under model)
#    * (function: likelihood of all sequences from each likelihood)
#    * DONE: (function: df from set of models)
#    * (function: probability from likelihood and df)
# 4) DONE: do it for the whole pile: this is the likelihood of one cluster, to beat
#    * already compared to null, which is astronomical
# 5) now create random pairs of piles of sequences
#    * random selection
# 6) for each sequence in each pile, ask if it's more likely under its own jackknifed transition matrix or under all other piles' transition matrices 
#    * (function: removing a sequence from its pile)
# 7) take the sequence with the biggest improvement under another pile and move it to that other pile, breaking ties by favoring moves from largest to smallest pile
#    * test if ties occur at all
#    * if no sequence in the large pile improves, repeat for all sseq
#    * (function: adding a sequence to a new pile)
# 8) stopping rules: 
#    1) no sequence is better in another pile
#    2) one pile is empty
#    3) detect or prevent thrashing.  jackknifing may be enough to prevent thrashing
#    4) report which stopping rule got engaged
# 9) after convergence calculate loglikelihood and penalized likelihood
# 10) repeat a few times from randomly generated piles to see variety in likelihoods, or their distribution

#install.packages(c("tidyverse", "bit", "stringdist", "assertthat"))
library(tidyverse)
library(readr)
library(assertthat)
library(stringdist)
library(nnet)

INCLUDE_TESTS = FALSE
### helpers
seqdist <- function(sequence, template, codes) {
   lseq = strsplit(sequence, '')[[1]]
   ltemplate = strsplit(template, '')[[1]]
   distance = 0
   for (i in seq(str_length(lseq))) {
     c <- lseq[i]
     distance = distance + chardist(c, ltemplate[i], codes)
   }
   return( distance )
}
chardist <- function(a,b,template) {
  return( abs( which(b == template) - which(a == template)))
}

assign_clusters_0 <- function(seqs, n_clusters) {
   seqs$cluster <- sample.int(n_clusters, nrow(seqs), replace=T)
   return( seqs ) 
}
assign_clusters <- function(n, n_clusters) {
   assignments <- sample.int(n_clusters, n, replace=T)
   return( assignments ) 
}

make_transition_matrix_template <- function(codes, constant=0.000001) {
   n = length(codes)
   #names = strsplit(codes, '')[[1]]
   names = codes
   matrix_template <- matrix(rep(constant,n*n), n,n, dimnames=list(names, names)) - diag(constant,n)
   if (TRUE) {
      matrix_template[,1] <- 0 # first column is zeros since you can't transition to initial state, only from it
   }
   return(matrix_template)
}
get_trans_df <- function(n) {
   if (TRUE) {
      return( (n^2) - n - 1 - (n-1) ) #substract first column, which is all zeros, since you can't transition into the opening state '-'
   } else {
      return( (n^2) - n - 1 )
   }
}

make_transition_matrix_null <- function(codes) {
   n = length(codes)
   #names = strsplit(codes, '')[[1]]
   names = codes
   constant = 1 / (get_trans_df(n) + 1)
   matrix_template <- matrix(rep(constant,n*n), n,n, dimnames=list(names, names)) - diag(constant,n)
   if (TRUE) {
      matrix_template[,1] <- 0 # first column is zeros since you can't transition to initial state, only from it
   }
   return(matrix_template)
}

make_transition_matrix <- function(seqs, codes) {
#    * (function: model from pile of sequences)
#      * model from empty pile is uniform
#      * diagonals of these models are 0
#      * not symmetric
#      * unpopulated cells are 0.000001
   trans <- make_transition_matrix_template(codes)
   n_transitions <- 0
   sequence_lengths <- rep(0, nrow(seqs))
   n_sequences <- 0
   #print(seqs)
   #print("SIZEBEFORE ITERATING")
   #print(nrow(seqs))
   for (i in 1:nrow(seqs)) { # for every sequence
      #print(i)
      #print(seqs[i,])
      sequence <- seqs[i,"sequence"][[1]]
      if (str_length(sequence) > 1) {
         l_sequence <- strsplit(sequence, '')[[1]]
         sequence_lengths[i] <- length(l_sequence)
         n_sequences <- n_sequences + 1
         #print(seq(length(l_sequence)-1))
         #print(dim(trans))
         for (j in seq(length(l_sequence)-1)) { # for every pair of entries
            the_row <- l_sequence[[j]]
            the_col <- l_sequence[[j+1]]
            #print(c(the_row,the_col))
            #print(rownames(trans))
            trans[the_row, the_col] = trans[the_row, the_col] + 1
            n_transitions = n_transitions + 1
         }
      }
   }
   # normalize
   trans <- trans / sum(trans)
   # check
   if (INCLUDE_TESTS) {
      assert_that(n_transitions == sum(sequence_lengths) - n_sequences, msg="in make_transition_matrix, 1")
      assert_that(round(sum(trans)) == 1, msg="in make_transition_matrix, 2")
      assert_that(nrow(trans) == ncol(trans), msg="in make_transition_matrix, 3")
      for (i in nrow(trans)) {
         assert_that(trans[i,i] == 0, msg="in make_transition_matrix, 4")
      }
      assert_that(all(trans[,0] == 0), msg="in make_transition_matrix, 5")
   }
   return(trans)
}

trans_loglik <- function(sequence, trans) {
   s = sequence
   loglik_trans = rep(0, length(s)-1)
   #print(sequence)
   #print(trans)
   #print(rownames(trans))
   #print(colnames(trans))
   #print(dim(trans))
   for (i in 1:(length(s)-1)) {
      assert_that(trans[s[i], s[i+1]] != 0, msg="in trans_loglik")
      loglik_trans[i] = log(trans[s[i], s[i+1]])
   }
   return(sum(loglik_trans))
}

trans_loglik_one_cluster <- function(seqs, codes, null_model=F) {
   loglik = 0
   for (i in 1:nrow(seqs)) {
      sequence = seqs[i, "sequence"][[1]]
      #print(sequence)
      sequence <- strsplit(sequence, '')[[1]]
      #print(sequence)
      if (!null_model) {
         trans <- make_transition_matrix(seqs[-i,], codes)
      } else {
         #trans <- make_transition_matrix(seqs[-i,], codes)
         trans <- make_transition_matrix_null(codes)
      }
      loglik = loglik + trans_loglik(sequence, trans)
   }
   return( loglik) 
}

# full data frame, and stats
policy_adoption <- read_delim("generated_data/policy_page_links-post03.csv", col_names=T, delim=',')
dim(policy_adoption)
str(policy_adoption)
spec(policy_adoption)

# focus down to relevant data, still one row per policy per lang
flows <- arrange(policy_adoption, lang, order_within_wiki) %>% select(lang, order_within_wiki, policy_name_en)
# policies, one row per policy, with appearance frequencies, and one-letter codes ("symbol") for string-based edit distances
policy_names <- policy_adoption %>% count(policy_name_en, sort=T)
if (TRUE) {
   codes <- c('-',letters, LETTERS, c(0,1,2,3,4,5,6,7,8,9))
   codes <- codes[1:(nrow(policy_names)+1)]
   policy_names <- cbind(policy_names, symbol=codes[-1])
} else {
   codes <- c(letters, LETTERS, c(0,1,2,3,4,5,6,7,8,9))
   codes <- codes[1:nrow(policy_names)]
   policy_names <- cbind(policy_names, symbol=codes)
}
str_codes <- paste0(codes, collapse='')
# enrich flows with policy symbols
flows <- flows %>% left_join(policy_names[,c(1,3)], by=join_by(policy_name_en))

# table of policy adoption sequences represented as string. one row per language edition
seqs <- flows %>% group_by(lang) %>% mutate(sequence = paste0(symbol, collapse = "")) 
seqs <- seqs %>% filter(order_within_wiki==0) %>% select("lang","sequence")
if (TRUE) {
   seqs$sequence <- paste('-', seqs$sequence, sep='')
}
seqs <- seqs %>% mutate(len_seq = str_length(sequence)) %>% arrange(desc(len_seq))

# tests
assert_that(1==1, msg="dummy")
assert_that(length(unique(policy_adoption$policy_name_en)) == nrow(policy_names), msg="dummy")
assert_that(length(unique(policy_adoption$lang)) == nrow(seqs), msg="dummy")

assert_that( chardist('-', 'a', codes) == 1, msg="dummy")
assert_that( chardist('a', 'a', codes) == 0, msg="dummy")
assert_that( chardist('a', 'b', codes) == 1, msg="dummy")
assert_that( chardist('b', 'a', codes) == 1, msg="dummy")
assert_that( chardist('1', '2', codes) == 1, msg="dummy")
assert_that( chardist('a', 'A', codes) == 26, msg="dummy")
assert_that( chardist('Z', 'z', codes) == 26, msg="dummy")
assert_that( chardist('0', '7', codes) == 7, msg="dummy")
assert_that( chardist('-', codes[length(codes)], codes) == (length(codes) - 1) , msg="dummy")
assert_that( chardist('a', codes[length(codes)], codes) == (length(codes) - 1 - 1) , msg="dummy")

assert_that( seqdist('-', 'a', codes) == 1 , msg="dummy")
assert_that( seqdist('a', 'a', codes) == 0 , msg="dummy")
assert_that( seqdist('-', str_codes, codes) == 0 , msg="dummy")
assert_that( seqdist('a', str_codes, codes) == 1 , msg="dummy")
assert_that( seqdist(str_codes, str_codes, codes) == 0 , msg="dummy")
assert_that( seqdist('b', str_codes, codes) == 2 , msg="dummy")
assert_that( seqdist('-ab', str_codes, codes) == 0 , msg="dummy")
assert_that( seqdist('-ba', str_codes, codes) == 2 , msg="dummy")
assert_that( seqdist('-bb', str_codes, codes) == 1 , msg="dummy")
assert_that( seqdist('-ca', str_codes, codes) == 3 , msg="dummy")
assert_that( seqdist('-z', str_codes, codes) == 26-1 , msg="dummy")
assert_that( seqdist('-ABC', str_codes, codes) == 26*3 , msg="dummy")

assert_that( round(sum(make_transition_matrix_null(codes[1:6]))) == 1, msg="dummy")
assert_that( sum(make_transition_matrix_null(codes[1:5])) == 1, msg="dummy")
assert_that( round(sum(make_transition_matrix_null(codes))) == 1, msg="dummy")
if (INCLUDE_TESTS) {
   make_transition_matrix_template(codes[1:3])
   make_transition_matrix_null(codes[1:3])
   make_transition_matrix_null(codes)
   make_transition_matrix(seqs, codes)
   a <- signif(make_transition_matrix(cbind("sequence"=c('-abc', '-bca', '-cab', '-bac', '-acb', '-cba')), codes[1:4]), 5)
   b <- signif(make_transition_matrix_null(codes[1:4]), 5)
   assert_that( all( a == b) , msg="dummy")
}

assert_that(get_trans_df( length(codes)) == (sum(ifelse(make_transition_matrix_null(codes) !=0,1,0)) - 1), msg="dummy")
assert_that(get_trans_df( 3 ) == (sum(ifelse(make_transition_matrix_null(codes[1:3]) !=0,1,0)) - 1), msg="dummy")


if (INCLUDE_TESTS) {
   print("solving one cluster")
   G2_r <- -2*trans_loglik_one_cluster(seqs, codes)
   G2_g <- -2*trans_loglik_one_cluster(seqs, codes, null_model=T)
   k_r = get_trans_df( length(codes) )
   k_g = 1
   n = nrow(seqs)
   BIC = ((G2_r^2 - G2_g^2) - (k_g - k_r) * log(n)) / 2
   BF = exp(BIC)
}

find_two_clusters <- function(seqs, init_assignments) {
   clusters = c(1,2)
   to_break = FALSE
   iterates = 0
   # 5) now create random pairs of piles of sequences
   #    * random selection
   #print(colnames(seqs))
   seqs <- cbind(seqs, assignment=init_assignments)
   #print(colnames(seqs))
   #print("IN f_t_c")
   while (!to_break) {
      print( paste0("ALGO ITERATE ", iterates))
      #print(colnames(seqs))
      #print(dim(seqs))
      # 6) for each sequence in each pile, ask if it's more likely under its own jackknifed transition matrix or under all other piles' transition matrices 
      #    * (function: removing a sequence from its pile)
      seqs$LLdiff <- 0
      seqs$LLdiff_raw <- 0
      seqs$LL_this <- 0
      seqs$LL_other <- 0
      seqs$cl_this <- 0
      seqs$cl_other <- 0
      for (i in 1:nrow(seqs)) {
         #print(colnames(seqs))
         #print(dim(seqs))
         sequence <- seqs[i, "sequence"][[1]]
         s <- strsplit(sequence, '')[[1]]
         cluster <- seqs[i, "assignment"][[1]]
         cluster_other <- clusters[which(clusters != cluster)]
         assert_that( seqs[i, "assignment"][[1]] != cluster_other )
         assert_that( seqs[i, "assignment"][[1]] == cluster )
         assert_that( cluster != cluster_other)
         #print(i)
         #print(dim(seqs))
         #print(dim(seqs[-i,]))
         #print("CHECK")
         #print(-i)
         #print(seqs[-i,])
         #print(subset(seqs[-i,], assignment==cluster))
         #print(rep(0, nrow(subset(seqs[-i,], assignment==cluster))))
         #print("GOING IN, SIZE FIRST")
         #b4 <- dim(seqs)
         #print(dim(seqs))
         clust1 <- seqs[-i,]
         clust1 <- subset(clust1, assignment==cluster)
         clust2 <- subset(seqs, assignment==cluster_other)
         trans1 <- make_transition_matrix(clust1, codes)
         #print("NOW SECOnD")
         #print(c(b4, dim(seqs)))
         trans2 <- make_transition_matrix(clust2, codes)
         ll1 <- trans_loglik(s, trans1) / seqs$len_seq[i]
         ll2 <- trans_loglik(s, trans2) / seqs$len_seq[i]
         seqs[i, 'LLdiff_raw'] <- trans_loglik(s, trans2) - trans_loglik(s, trans1)
         seqs[i, 'LLdiff'] <- ll2 - ll1 
         seqs[i, 'LL_other'] <- ll2 
         seqs[i, 'LL_this'] <- ll1 
         seqs[i, 'cl_this'] <- cluster
         seqs[i, 'cl_other'] <- cluster_other
         #print(colnames(seqs))
         #print(dim(seqs))
         assert_that( seqs[i, "assignment"][[1]] != cluster_other , msg=paste0( ))
         assert_that( seqs[i, "assignment"][[1]] == cluster , msg=paste0( ))
      }
      # which.is.max is from nnet, better than which.max because it breaks ties randomly
      winner <- which.is.max(seqs$LLdiff) 
      #print(winner)
      #print( seqs[winner,])
      winner_LLdiff <- max(seqs$LLdiff) 
      winner_LLdiff2 <- seqs[winner, "LLdiff"][[1]]
      winner_cluster  <- seqs[winner, "assignment"][[1]]
      winner_cluster_other <- clusters[which(clusters != winner_cluster)]
      assert_that(seqs[winner,"LLdiff"] == winner_LLdiff, msg=paste0("ERROR KLJHVBJN", seqs[winner,"LLdiff"], winner_LLdiff))
      assert_that( winner_LLdiff == winner_LLdiff2, msg=paste0("ERROR *YGBJ", winner_LLdiff, winner_LLdiff2, winner))
      assert_that( winner_cluster == seqs[winner, "cl_this"][[1]], msg=paste0( "ERROR KJHGTYHJSSS", winner_cluster , seqs[winner, "cl_this"][[1]]))
      assert_that( winner_cluster != seqs[winner, "cl_other"][[1]], msg=paste0("ERROR JHGHJ", winner_cluster , seqs[winner, "cl_other"][[1]] ))
      # 7) take the sequence with the biggest improvement under another pile and move it to that other pile, breaking ties by favoring moves from largest to smallest pile
      #    * test if ties occur at all
      #    * if no sequence in the large pile improves, repeat for all sseq
      #    * (function: adding a sequence to a new pile)
      # 8) stopping rules: 
      #    1) no sequence is better in another pile
      #    2) one pile is empty
      #    3) detect or prevent thrashing.  jackknifing may be enough to prevent thrashing
      #    4) report which stopping rule got engaged
      if (sum(seqs$LLdiff == winner_LLdiff) > 1) {
         print(paste0("ANNOUNCE LSDAFHLJKDFF: lot of ties:", sum(seqs$LLdiff == winner_LLdiff)))
      }
      if (winner_LLdiff > 0) {
         #print(paste0("ANNOUNCE GHJKJHD: found a sequence to change to ", winner_cluster_other))
         #print(winner_cluster)
         #print(winner_cluster_other)
         #print(seqs[winner,])
         seqs_test <- seqs
         seqs$assignment[winner] <- winner_cluster_other
         print(seqs[winner,])
      } else {
         # no sequence can improve by switching. search has converged. We can break
         print("ANNOUNCE DSKLFJ: no sequence can improve by switching. search has converged.")
         seqs$iterates[winner] <- iterates
         to_break = TRUE
      }
      if ((sum(seqs$LLdiff)==0) | (sum(seqs$LLdiff)==nrow(seqs) )) {
         # one of the the clusters are empty
         print("ANNOUNCE &*UHBD:  one of the the clusters are empty")
         to_break = TRUE
      }
      iterates = iterates + 1
      print( paste0("end of iterate ", iterates))
      print( "cluster sizes:")
      print( table(seqs$assignment))
      print( "sequence sizes by cluster:")
      print( seqs %>% group_by(assignment) %>% summarize(mean(len_seq)))
      if (FALSE) {
         print( "ll's by cluster")
         c1_lls <- trans_loglik_one_cluster(subset(seqs, assignment==1), codes)
         c2_lls <- trans_loglik_one_cluster(subset(seqs, assignment==2), codes)
         print( c(c1_lls, c2_lls))
      }
   }
   return(seqs)
}

find_three_clusters <- function(seqs, init_assignments) {
   clusters = c(1,2,3)
   to_break = FALSE
   iterates = 0
   # 5) now create random pairs of piles of sequences
   #    * random selection
   #print(colnames(seqs))
   seqs <- cbind(seqs, assignment=init_assignments)
   #print(colnames(seqs))
   #print("IN f_t_c")
   while (!to_break) {
      print( paste0("ALGO ITERATE ", iterates))
      #print(colnames(seqs))
      #print(dim(seqs))
      # 6) for each sequence in each pile, ask if it's more likely under its own jackknifed transition matrix or under all other piles' transition matrices 
      #    * (function: removing a sequence from its pile)
      seqs$LLdiff2 <- 0
      seqs$LLdiff3 <- 0
      #seqs$LLdiff2_raw <- 0
      #seqs$LLdiff3_raw <- 0
      seqs$LL_this <- 0
      seqs$LL_other <- 0
      seqs$LL_other_other <- 0
      seqs$cl_this <- 0
      seqs$cl_other <- 0
      seqs$cl_other_other <- 0
      for (i in 1:nrow(seqs)) {
         #print(colnames(seqs))
         #print(dim(seqs))
         sequence <- seqs[i, "sequence"][[1]]
         s <- strsplit(sequence, '')[[1]]
         cluster <- seqs[i, "assignment"][[1]]
         cluster_other <- clusters[which(clusters != cluster)][1]
         cluster_other_other <- clusters[which(clusters != cluster)][2]
         assert_that( seqs[i, "assignment"][[1]] != cluster_other )
         assert_that( seqs[i, "assignment"][[1]] != cluster_other_other )
         assert_that( seqs[i, "assignment"][[1]] == cluster )
         assert_that( cluster != cluster_other)
         assert_that( cluster != cluster_other_other)
         assert_that( cluster_other != cluster_other_other)
         #print(i)
         #print(dim(seqs))
         #print(dim(seqs[-i,]))
         #print("CHECK")
         #print(-i)
         #print(seqs[-i,])
         #print(subset(seqs[-i,], assignment==cluster))
         #print(rep(0, nrow(subset(seqs[-i,], assignment==cluster))))
         #print("GOING IN, SIZE FIRST")
         #b4 <- dim(seqs)
         #print(dim(seqs))
         clust1 <- seqs[-i,]
         clust1 <- subset(clust1, assignment==cluster)
         clust2 <- subset(seqs, assignment==cluster_other)
         clust3 <- subset(seqs, assignment==cluster_other_other)
         trans1 <- make_transition_matrix(clust1, codes)
         #print("NOW SECOnD")
         #print(c(b4, dim(seqs)))
         trans2 <- make_transition_matrix(clust2, codes)
         trans3 <- make_transition_matrix(clust3, codes)
         ll1 <- trans_loglik(s, trans1) / seqs$len_seq[i]
         ll2 <- trans_loglik(s, trans2) / seqs$len_seq[i]
         ll3 <- trans_loglik(s, trans3) / seqs$len_seq[i]
         seqs[i, 'LL_this'] <- ll1 
         seqs[i, 'LL_other'] <- ll2 
         seqs[i, 'LL_other_other'] <- ll3 
         #seqs[i, 'LLdiff2_raw'] <- trans_loglik(s, trans2) - trans_loglik(s, trans1)
         #seqs[i, 'LLdiff3_raw'] <- trans_loglik(s, trans3) - trans_loglik(s, trans1)
         seqs[i, 'LLdiff2'] <- ll2 - ll1 
         seqs[i, 'LLdiff3'] <- ll3 - ll1 
         seqs[i, 'cl_this'] <- cluster
         seqs[i, 'cl_other'] <- cluster_other
         seqs[i, 'cl_other_other'] <- cluster_other_other
         #print(colnames(seqs))
         #print(dim(seqs))
         assert_that( seqs[i, "assignment"][[1]] != cluster_other_other , msg=paste0( ))
         assert_that( seqs[i, "assignment"][[1]] != cluster_other , msg=paste0( ))
         assert_that( seqs[i, "assignment"][[1]] == cluster , msg=paste0( ))
      }
      # which.is.max is from nnet, better than which.max because it breaks ties randomly
      contestant2 <- which.is.max(seqs$LLdiff2) 
      contestant3 <- which.is.max(seqs$LLdiff3) 
      if (max(seqs$LLdiff2) > max(seqs$LLdiff3)) {
         winner <- which.is.max(seqs$LLdiff2) 
         winner_LLdiff <- seqs[winner, "LLdiff2"][[1]]
         winner_cluster_from  <- seqs[winner, "assignment"][[1]]
         winner_cluster_to <- clusters[which(clusters != winner_cluster_from)][1]
         print(paste0("watch diff2:", ll2, ll1, max(seqs$LLdiff2), collapse=' '))
      } else {
         winner <- which.is.max(seqs$LLdiff3) 
         winner_LLdiff <- seqs[winner, "LLdiff3"][[1]]
         winner_cluster_from  <- seqs[winner, "assignment"][[1]]
         winner_cluster_to <- clusters[which(clusters != winner_cluster_from)][2]
         print(paste0("watch diff3:", ll3, ll1, max(seqs$LLdiff3), collapse=' '))
      }
      #print(winner)
      #print( seqs[winner,])
      assert_that( winner_LLdiff %in% c(seqs[winner,"LLdiff2"], seqs[winner,"LLdiff3"]), msg=paste0("ERROR KLJHVBJN", seqs[winner,"LLdiff2"], seqs[winner,"LLdiff3"], winner_LLdiff))
      #assert_that( winner_LLdiff == winner_LLdiff2, msg=paste0("ERROR *YGBJ", winner_LLdiff, winner_LLdiff2, winner))
      assert_that( winner_cluster_from == seqs[winner, "cl_this"][[1]], msg=paste0( "ERROR KJHGTYHJSSS", winner_cluster_from , seqs[winner, "cl_this"][[1]]))
      assert_that( winner_cluster_from != seqs[winner, "cl_other"][[1]], msg=paste0("ERROR JHGHJ", winner_cluster_from , seqs[winner, "cl_other"][[1]] ))
      assert_that( winner_cluster_from != seqs[winner, "cl_other_other"][[1]], msg=paste0("ERROR JHGHJ", winner_cluster_from , seqs[winner, "cl_other_other"][[1]] ))
      # 7) take the sequence with the biggest improvement under another pile and move it to that other pile, breaking ties by favoring moves from largest to smallest pile
      #    * test if ties occur at all
      #    * if no sequence in the large pile improves, repeat for all sseq
      #    * (function: adding a sequence to a new pile)
      # 8) stopping rules: 
      #    1) no sequence is better in another pile
      #    2) one pile is empty
      #    3) detect or prevent thrashing.  jackknifing may be enough to prevent thrashing
      #    4) report which stopping rule got engaged
      if (sum(seqs$LLdiff == winner_LLdiff) > 1) {
         print(paste0("ANNOUNCE LSDAFHLJKDFF: lot of ties:", sum(winner_LLdiff %in% c(seqs$LLdiff2, seqs$LLdiff3 ))))
      }
      if (winner_LLdiff > 0) {
         print(paste0("ANNOUNCE GHJKJHD: found a sequence to change to ", winner_cluster_from, winner_cluster_to, winner, collapse=' '))
         print(as.data.frame(seqs[winner,]))
         seqs_test <- seqs
         seqs$assignment[winner] <- winner_cluster_to
         print(as.data.frame(seqs[winner,]))
      } else {
         # no sequence can improve by switching. search has converged. We can break
         print(paste0("ANNOUNCE DSKLFJ: no sequence can improve by switching. search has converged. ",  iterates, ' ', winner, ' ', dim(seqs), ' ', length(seqs$iterates), ' ', length(seqs$iterates[winner]), length(seqs$iterates), ' ', length(seqs$iterates[winner])))
         #print(dim(seqs))
         #print(winner)
         #print(seqs$iterates)
         #print(seqs$iterates[winner])
         #print(iterates)
         #print(str(winner))
         #print(str(seqs$iterates))
         #print(str(seqs$iterates[winner]))
         #print(str(iterates))
         seqs$iterates[winner] <- iterates
         to_break = TRUE
      }
      if ((sum(seqs$LLdiff2)==0) | (sum(seqs$LLdiff2)==nrow(seqs) ) | (sum(seqs$LLdiff3)==0) | (sum(seqs$LLdiff3)==nrow(seqs) )) {
         # one of the the clusters are empty
         print("ANNOUNCE &*UHBD:  one of the the clusters are empty")
         to_break = TRUE
      }
      iterates <- iterates + 1
      print( paste0("end of iterate ", iterates))
      print( "cluster sizes:")
      print( table(seqs$assignment))
      print( "sequence sizes by cluster:")
      print( seqs %>% group_by(assignment) %>% summarize(mean(len_seq)))
      if (FALSE) {
         print( "ll's by cluster")
         c1_lls <- trans_loglik_one_cluster(subset(seqs, assignment==1), codes)
         c2_lls <- trans_loglik_one_cluster(subset(seqs, assignment==2), codes)
         print( c(c1_lls, c2_lls))
      }
   }
   return(seqs)
}

find_k_clusters <- function(seqs, init_assignments, k) {
   clusters = 1:k
   to_break = FALSE
   iterates = 0
   # 5) now create random pairs of piles of sequences
   #    * random selection
   #print(colnames(seqs))
   seqs <- cbind(seqs, assignment=init_assignments)
   #print(colnames(seqs))
   #print("IN f_t_c")
   while (!to_break) {
      print( paste0(k , '-', "ALGO ITERATE ", iterates))
      #print(colnames(seqs))
      #print(dim(seqs))
      # 6) for each sequence in each pile, ask if it's more likely under its own jackknifed transition matrix or under all other piles' transition matrices 
      #    * (function: removing a sequence from its pile)
      seqs$LLdiff1 <- 0
      diff1_idx <-  which(names(seqs) == 'LLdiff1', arr.ind=TRUE)
      # k minus one for off by one, minus one again because there is one less diff than lls total.
      if (k > 2) {
         diff_idxs <-  diff1_idx:(diff1_idx+k-1-1) 
      } else if (k == 2) {
         diff_idxs <-  c(diff1_idx)
      } else {print("can't do k=1")}
      seqs[,diff_idxs ] <- 0
      names(seqs)[diff_idxs ] <- paste("LLdiff", 1:(k-1), sep='')
      LL_idxs <- (diff1_idx+k-1):(diff1_idx+k+k-1-1)
      seqs[,LL_idxs ] <- 0
      names(seqs)[LL_idxs ] <- paste("LL", 1:k, sep='')
      seqs$cl_this <- 0
      for (i in 1:nrow(seqs)) {
         #print(colnames(seqs))
         #print(dim(seqs))
         sequence <- seqs[i, "sequence"][[1]]
         s <- strsplit(sequence, '')[[1]]
         cluster_i <- seqs[i, "assignment"][[1]]
         cluster_others <- clusters[which(clusters != cluster_i)]
         assert_that( seqs[i, "assignment"][[1]] == cluster_i , msg='tt')
         for (c in 1:length(cluster_others)) {
            #print(seqs[i, "assignment"][[1]])
            #print(cluster_others[c])
            #print(cluster_i)
            #print(length(seqs[i, "assignment"][[1]]))
            #print(length(cluster_others[c]))
            assert_that( seqs[i, "assignment"][[1]] != cluster_others[c], msg='xyz' )
            assert_that( cluster_i != cluster_others[c], msg="hey")
         }
         #print(dim(seqs))
         #print(dim(seqs[-i,]))
         #print("CHECK")
         #print(-i)
         #print(seqs[-i,])
         #print(subset(seqs[-i,], assignment==cluster_i))
         #print(rep(0, nrow(subset(seqs[-i,], assignment==cluster_i))))
         #print("GOING IN, SIZE FIRST")
         #b4 <- dim(seqs)
         #print(dim(seqs))
         the_clusters <- as.list(rep(0,k))
         the_matrices <- as.list(rep(0,k))
         the_clusters[[cluster_i]] <- seqs[-i,]
         the_clusters[[cluster_i]] <- subset(the_clusters[[cluster_i]], assignment==cluster_i)
         for (c in cluster_others) {
            the_clusters[[c]] <- subset(seqs, assignment==c)
         }
         for (c in clusters) {
            the_matrices[[c]] <- make_transition_matrix(the_clusters[[c]], codes)
            an_ll <- trans_loglik(s, the_matrices[[c]]) / seqs$len_seq[i]
            seqs[i, LL_idxs[c]] <- an_ll 
         }
         #print("NOW SECOnD")
         #print(c(b4, dim(seqs)))
         for (c in 1:length(cluster_others)) {
            seqs[i, diff_idxs[c]] <- seqs[i, LL_idxs[cluster_others[c]] ] - seqs[i, LL_idxs[cluster_i] ]
         }
         #seqs[i, 'LLdiff2_raw'] <- trans_loglik(s, trans2) - trans_loglik(s, trans1)
         #seqs[i, 'LLdiff3_raw'] <- trans_loglik(s, trans3) - trans_loglik(s, trans1)
         seqs[i, 'cl_this'] <- cluster_i
         #print(colnames(seqs))
         #print(dim(seqs))
         for (c in cluster_others) {
            assert_that( seqs[i, "assignment"][[1]] != c , msg=paste0('z' ))
            assert_that( seqs[i, "assignment"][[1]] != c , msg=paste0('y' ))
            assert_that( seqs[i, "assignment"][[1]] == cluster_i , msg=paste0('x' ))
         }
      }
      # which.is.max is from nnet, better than which.max because it breaks ties randomly
       
      contestants_idx <- as.data.frame(which(seqs[,diff_idxs] == max(seqs[,diff_idxs]), arr.ind=TRUE))
      winner_idx <- unlist(contestants_idx[sample(nrow(contestants_idx)),][1,])
      winner <- winner_idx[1]
      winner_LLdiff <- seqs[winner_idx[1], diff_idxs[winner_idx[2]]][[1]]
      winner_cluster_from  <- seqs[winner, "assignment"][[1]]
      winner_cluster_others <- clusters[-winner_cluster_from]
      winner_cluster_to <- winner_cluster_others[winner_idx[2]]
      assert_that( signif(winner_LLdiff,3) %in% signif(seqs[winner,diff_idxs], 3), msg=paste0("ERROR KLJHVBJN"))
      #assert_that( winner_LLdiff == winner_LLdiff2, msg=paste0("ERROR *YGBJ", winner_LLdiff, winner_LLdiff2, winner))
      assert_that( winner_cluster_from == seqs[winner, "cl_this"][[1]], msg=paste0( "ERROR KJHGTYHJSSS", winner_cluster_from , seqs[winner, "cl_this"][[1]]))
      assert_that( winner_cluster_from != winner_cluster_to, msg=paste0( "ERROR OIDHJ9_d", winner_cluster_from , winner_cluster_to ))
      # 7) take the sequence with the biggest improvement under another pile and move it to that other pile, breaking ties by favoring moves from largest to smallest pile
      #    * test if ties occur at all
      #    * if no sequence in the large pile improves, repeat for all sseq
      #    * (function: adding a sequence to a new pile)
      # 8) stopping rules: 
      #    1) no sequence is better in another pile
      #    2) one pile is empty
      #    3) detect or prevent thrashing.  jackknifing may be enough to prevent thrashing
      #    4) report which stopping rule got engaged
      if (nrow(contestants_idx) > 1) {
         print(paste0("ANNOUNCE LSDAFHLJKDFF: lot of ties:"))
      }
      if (winner_LLdiff > 0) {
         print(paste0("ANNOUNCE GHJKJHD: found a sequence to change to ", winner_cluster_from, winner_cluster_to, winner, collapse=' '))
         #print(as.data.frame(seqs[winner,]))
         seqs_test <- seqs
         seqs$assignment[winner] <- winner_cluster_to
         #print(as.data.frame(seqs[winner,]))
      } else {
         # no sequence can improve by switching. search has converged. We can break
         print(paste0("ANNOUNCE DSKLFJ: no sequence can improve by switching. search has converged. ",  iterates, ' ', winner, ' ', dim(seqs), ' ', length(seqs$iterates), ' ', length(seqs$iterates[winner]), length(seqs$iterates), ' ', length(seqs$iterates[winner])))
         #print(dim(seqs))
         #print(winner)
         #print(seqs$iterates)
         #print(seqs$iterates[winner])
         #print(iterates)
         #print(str(winner))
         #print(str(seqs$iterates))
         #print(str(seqs$iterates[winner]))
         #print(str(iterates))
         seqs$iterates[winner] <- iterates
         to_break = TRUE
      }
      if ( length(unique( seqs$assignment)) != length(clusters)){
         # one of the the clusters are empty
         print("ANNOUNCE &*UHBD:  one of the the clusters are empty")
         to_break = TRUE
      }
      iterates <- iterates + 1
      print( paste0("end of iterate ", iterates))
      print( "cluster sizes:")
      print( table(seqs$assignment))
      print( "sequence sizes by cluster:")
      print( seqs %>% group_by(assignment) %>% summarize(mean(len_seq)))
      if (FALSE) {
         print( "ll's by cluster")
         for (c in clusters) {
            c_lls <- trans_loglik_one_cluster(subset(seqs, assignment==c), codes)
            print( c_lls)
         }
      }
   }
   return(seqs)
}


# summary: spplit the sequences into two piles, build two transition matrices, run all sequences through both transition matrices (jacknifed to exclude sequence in question "out of sample"), calculate LLs, take the sequence with a better LL in the other cluster, move it, repeat.
# function for LL per transition
# stick find_two_clusters into a loop 
if (FALSE) {
   reps = 100
   seqs$iterates <- 0
   seq_of_seqs_c2 = as.list(rep(0,reps))
   for (i in 1:reps) {
      print(c("REPETITION ", i))
      seq_of_seqs_c2[[i]] <- find_two_clusters(seqs, assign_clusters(nrow(seqs), 2))
      save(seq_of_seqs_c2, file="c2_results_003.Rdata")
   }
} 
if (FALSE) {
   reps = 100
   seqs$iterates <- 0
   seq_of_seqs_c2 = as.list(rep(0,reps))
   seq_of_seqs_c3 = as.list(rep(0,reps))
   for (i in 1:reps) {
      print(c("REPETITION ", i))
      seq_of_seqs_c3[[i]] <- find_three_clusters(seqs, assign_clusters(nrow(seqs), 3))
      save(seq_of_seqs_c3, file="c3_results_003.Rdata")
   }
}
if (TRUE) {
   reps = 100
   seqs$iterates <- 0
   seq_of_seqs_c2 = as.list(rep(0,reps))
   seq_of_seqs_c3 = as.list(rep(0,reps))
   seq_of_seqs_c4 = as.list(rep(0,reps))
   for (i in 1:reps) {
      print(c("REPETITION ", i))
      seq_of_seqs_c4[[i]] <- find_k_clusters(seqs, assign_clusters(nrow(seqs), 4), 4)
      save(seq_of_seqs_c4, file="c4_results_002.Rdata")
   }
}
