context("data fetch")

source("helper-organismdb.R")
source("helper-makeEset.R")

test_that("block data fetch works", {
  gr <- GRanges(seqnames="chr1", ranges=IRanges::IRanges(start=1:10, width=1),
                 seqinfo=Seqinfo(seqnames="chr1",genome="hcb"))
  ms_obj <- epivizrData::register(gr)
  query <- GRanges(seqnames="chr1", ranges=IRanges::IRanges(start=2, end=6))
  res <- ms_obj$get_rows(query, character())
  out <- list(globalStartIndex=2,
              useOffset=FALSE,
              values=list(id=2:6,
                chr=rep("chr1", 5),
                start=2:6,
                end=2:6,
                metadata=NULL))
  expect_equal(res, out)
})

test_that("block fetch works on unsorted data", {
  gr <- GRanges(seqnames="chr1", ranges=IRanges::IRanges(start=10:1, width=1),
                seqinfo=Seqinfo(seqnames="chr1",genome="hcb"))
  ms_obj <- epivizrData::register(gr)

  query <- GRanges(seqnames="chr1", ranges=IRanges::IRanges(start=2, end=6))
  res <- ms_obj$get_rows(query, NULL)
  out <- list(globalStartIndex=2,
              useOffset=FALSE,
              values=list(id=2:6,
                chr=rep("chr1", 5),
                start=2:6,
                end=2:6,
                metadata=(NULL)))
  expect_equal(res, out)
})

test_that("data fetch works on bp data", {
  gr <- GRanges(seqnames="chr1", 
                ranges=IRanges::IRanges(start=seq(1,100,by=5), width=1), 
                score1=seq(1,100,by=5), score2=-seq(1,100,by=5),
                seqinfo=Seqinfo(seqnames=c("chr1","chr2"), genome="hcb"))
  
  ms_obj <- epivizrData::register(gr, type="bp")

  query <- GRanges(seqnames="chr1", ranges=IRanges::IRanges(start=2,end=6))
  res <- ms_obj$get_rows(query, NULL)
  out <- list(globalStartIndex=2,
              useOffset=FALSE,
              values=list(id=list(2),
                chr=list("chr1"),
                start=list(6),
                end=list(6),
                metadata=NULL))

  expect_equal(res, out)
  #print(res);print(out)

  res <- ms_obj$get_values(query, c("score1"))
  out <- list(globalStartIndex=2,
              values=list(6))
  expect_equal(res, out)  
})

test_that("data fetch works on bp data with NAs", {
  # TODO: fix this test
  gr <- GRanges(seqnames="chr1", 
                ranges=IRanges::IRanges(start=seq(1,100,by=5), width=1), 
                score1=seq(1,100,by=5), score2=-seq(1,100,by=5),
                seqinfo=Seqinfo(seqnames=c("chr1","chr2"),genome="hcb"))
  gr$score2[1:10] <- NA

  ms_obj <- epivizrData::register(gr, type="bp")
  query <- GRanges("chr1", IRanges::IRanges(start=2, end=6))
  res <- ms_obj$get_rows(query, NULL)
  out <- list(globalStartIndex=NULL,
             useOffset=FALSE,
             values=list(id=list(),
                         chr=list(),
                         start=list(),
                         end=list(),
                         metadata=NULL))
  expect_equal(res, out)
 
  res <- ms_obj$get_values(query, c("score1"))
  out <- list(globalStartIndex=NULL, values=list())
  expect_equal(res,out)
})


test_that("feature data fetch works", {
  skip_if_not_installed("hgu133plus2.db")
  eset <- make_test_eset()
  ms_obj <- epivizrData::register(eset, columns=c("SAMP_1", "SAMP_2"))
  query <- GRanges(seqnames="chr6", 
                   ranges=IRanges::IRanges(start=30000000,end=40000000))

  olaps <- findOverlaps(query, ms_obj$.object)
  hits <- unique(subjectHits(olaps))
  hits <- seq(min(hits), max(hits))
  tmp <- ms_obj$.object[hits,]

  m <- match(rowRanges(tmp)$PROBEID, featureNames(eset))
  mat <- exprs(eset)[m, c("SAMP_1", "SAMP_2")]

  res <- ms_obj$get_rows(query, c("PROBEID","SYMBOL"))
  
  out <- list(globalStartIndex=min(hits),
              useOffset=FALSE,
              values=list(
                id=hits,
                chr=as.vector(seqnames(tmp)),
                start=start(tmp),
                end=end(tmp),
                metadata=list(PROBEID=rowRanges(tmp)$PROBEID,
                  SYMBOL=rowRanges(tmp)$SYMBOL)
                ))
  expect_equal(res, out)
  #print(res); print(out)
  
  res <- ms_obj$get_values(query, "SAMP_1")
  out <- list(globalStartIndex=min(hits),
              values=unname(mat[,"SAMP_1"]))
  #print(res);print(out)
  expect_equal(res,out)
})

test_that("geneinfo fetch works", {
  skip_if_not_installed("bumphunter")
  gr <- make_test_gene_info()
  msmt <- epivizrData::register(gr, type="gene_info")
  query <- GRanges("chr11", IRanges::IRanges(start=102500000, end=103000000))
  res <- msmt$get_rows(query, c("gene", "exon_starts", "exon_ends"))
  
  msGR <- msmt$.object
  olaps <- findOverlaps(query, msGR)
  hits <- subjectHits(olaps)
  hits <- seq(min(hits), max(hits))
  tmp <- msGR[hits,]
  
  out <- list(globalStartIndex=hits[1],
              useOffset=FALSE,
              values=list(
                id=hits,
                chr=as.vector(seqnames(tmp)),
                start=start(tmp),
                end=end(tmp),
                metadata=list(gene=unname(as.character(tmp$Gene)),
                              exon_starts=unname(lapply(start(tmp$Exons),paste,collapse=",")),
                              exon_ends=unname(lapply(end(tmp$Exons), paste, collapse=","))),
                strand=unname(as.character(strand(tmp)))))
  #print(res); print(out)
  expect_equal(res, out)
})

