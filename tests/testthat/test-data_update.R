context("update measurement")

source("helper-makeEset.R")

test_that("update block works", {
  gr1 <- GRanges(seqnames="chr1", ranges=IRanges::IRanges(start=1:10, width=1))
  gr2 <- GRanges(seqnames="chr12", ranges=IRanges::IRanges(start=1:1000,width=10))
  
  ms_obj <- epivizrData::register(gr1)
	ms_obj$update(gr2)
	expect_identical(as(ms_obj$.object, "GRanges"), gr2)
})


test_that("update bp works", {
  gr1 <- GRanges(seqnames="chr1", 
                 ranges=IRanges::IRanges(start=seq(1,100,by=25), width=1), 
                 score1=rnorm(length(seq(1,100,by=25))),
                 score2=rnorm(length(seq(1,100,by=25))))
  gr2 <- GRanges(seqnames="chr12", 
                 ranges=IRanges::IRanges(start=seq(1,100,by=25), width=1), 
                 score1=rnorm(length(seq(1,100,by=25))),
                 score2=rnorm(length(seq(1,100,by=25))))
  gr3 <- gr1
  gr3$score1 <- NULL

  ms_obj <- epivizrData::register(gr1, type="bp")
	ms_obj$update(gr2)
	expect_identical(as(ms_obj$.object, "GRanges"), gr2)
	
	expect_error(msObj$update(gr3))
}) 

test_that("update feature works", {
  sset <- make_test_SE()
  sset2 <- sset[2:10,]
  sset3 <- sset2[,-1]

  ms_obj <- epivizrData::register(sset, columns=c("A","B"), assay="counts2")
	ms_obj$update(sset2)
	
	expect_identical(as(rowRanges(ms_obj$.object), "GRanges"), rowRanges(sset2))
	expect_identical(assays(ms_obj$.object), assays(sset2))
	expect_identical(colData(ms_obj$.object), colData(sset2))

	expect_error(ms_obj$update(gr3))
}) 

