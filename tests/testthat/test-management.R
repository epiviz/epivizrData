context("measurement management")

source("helper-organismdb.R")
source("helper-makeEset.R")

test_that("add measurement works without datasource_name", {
  server <- epivizrServer::createServer()
  mgr <- createMgr(server)
  
  se <- make_test_SE()
  ms_obj <- mgr$add_measurements(se, columns=c("A", "B"), assay="counts2", send_request=FALSE)
  expect_equal(ms_obj$get_name(), "se")
})

test_that("add measurement works without connection", {
  server <- epivizrServer::createServer()
  mgr <- createMgr(server)
  
  se <- make_test_SE()
  ms_obj <- mgr$add_measurements(se, "example", columns=c("A", "B"), assay="counts2", send_request=FALSE)
  ms_id <- ms_obj$get_id()
  ms_name <- ms_obj$get_source_name()
  
  expect_equal(mgr$num_datasources(), 1)
  expect_false(is.null(mgr$.ms_list[[ms_id]]))

  rngs <- unname(sapply(c("A","B"), function(col) range(pretty(range(assay(se,"counts2")[,col], na.rm=TRUE)))))
  exp_ms <- lapply(c("A","B"), function(col) {
    i <- match(col,c("A","B"))
    list(
         name=col,
         type="feature",
         datasourceGroup=ms_id,
         defaultChartType="ScatterPlot",
         dataprovider=character(),
         annotation=list(Treatment=as.character(colData(se)[i,])),
         minValue=rngs[1,i],
         maxValue=rngs[2,i],
         metadata=c("probeid"),
         id=col,
         datasourceId=ms_id,
         datasourceName=ms_name
         )
  })
  
  ms_record <- mgr$.ms_list[[ms_id]]
  expect_equal(lapply(ms_record$measurements, as.list), exp_ms)
  expect_equal(ms_record$name, "example")
  expect_identical(ms_record$obj, ms_obj)
  expect_false(ms_record$connected)
})

test_that("get_measurements works without connection", {
  skip_if_not_installed("hgu133plus2.db")
  gr1 <- GRanges(seqnames="chr1", ranges=IRanges::IRanges(start=1:10, width=100))
  gr2 <- GRanges(seqnames="chr2", ranges=IRanges::IRanges(start=2:20, width=100))
  gr3 <- GRanges(seqnames="chr1", ranges=IRanges::IRanges(start=seq(1,100,by=25), width=1), 
                 score=rnorm(length(seq(1,100,by=25))))
  eset <- make_test_eset()
  
  server <- epivizrServer::createServer()
  mgr <- epivizrData::createMgr(server)
  
  msObj1 <- mgr$add_measurements(gr1, "dev1", send_request=FALSE); msId1=msObj1$get_id(); msName1=msObj1$get_source_name()
  msObj2 <- mgr$add_measurements(gr2, "dev2", send_request=FALSE); msId2=msObj2$get_id(); msName2=msObj2$get_source_name()
  msObj3 <- mgr$add_measurements(gr3, "dev3", send_request=FALSE, type="bp"); msId3=msObj3$get_id(); msName3=msObj3$get_source_name()
  msObj4 <- mgr$add_measurements(eset, "dev4", send_request=FALSE, columns=c("SAMP_1", "SAMP_2")); msId4=msObj4$get_id(); msName4=msObj4$get_source_name()
    
  rngs3 <- range(pretty(range(mcols(gr3)[,"score"],na.rm=TRUE)))
  rngs4 <- sapply(1:2, function(i) range(pretty(range(exprs(eset)[,paste0("SAMP_",i)],na.rm=TRUE))))
    
  res <- mgr$get_measurements()
    
  expMs <- list(
      
      name=c("dev1", "dev2", "score", paste0("SAMP_",1:2)),
      type=c(rep("range", 2), rep("feature",3)),
      datasourceGroup=c(msId1, msId2, msId3, rep(msId4,2)),
      defaultChartType=c(rep("BlocksTrack", 2), "LineTrack", rep("ScatterPlot",2)),
      annotation=c(rep(list(NULL),3), 
                   lapply(1:2, function(i) 
                     list(a=as.character(pData(eset)[i,1]), 
                          b=as.character(pData(eset)[i,2])))),
      minValue=c(rep(NA,2), rngs3[1], rngs4[1,]),
      maxValue=c(rep(NA,2), rngs3[2], rngs4[2,]),
      metadata=c(list(NULL), list(NULL), list(NULL), lapply(1:2,function(i) c("PROBEID","SYMBOL"))),
      id=c(msId1, msId2, "score", paste0("SAMP_",1:2)),
      datasourceId=c(msId1, msId2, msId3, rep(msId4,2)),
      datasourceName=c(msName1, msName2, msName3, rep(msName4, 2))
    )
    
  for (entry in names(expMs)) {
    expect_equal(res[[entry]], expMs[[entry]], label=entry)
  }
  
  for (id in ls(mgr$.ms_list)) {
    expect_false(mgr$.ms_list[[id]]$connected)
  }
})

test_that("rm_measurements works without connection", {
  gr <- GRanges(seqnames="chr1", ranges=IRanges::IRanges(start=seq(1,100,by=25), width=1),
                score1=rnorm(length(seq(1,100,by=25))),
                score2=rnorm(length(seq(1,100,by=25))))

  server <- epivizrServer::createServer()
  mgr <- createMgr(server)
  
  ms_obj <- mgr$add_measurements(gr, "dev1", send_request=FALSE, type="bp")
  id <- ms_obj$get_id()
  mgr$rm_measurements(id)
  
  expect_equal(mgr$num_datasources(), 0)
  expect_true(is.null(mgr$.ms_list[[id]]))
    
  ms_obj <- mgr$add_measurements(gr, "dev1", send_request=FALSE, type="bp")
  mgr$rm_measurements(ms_obj)
    
  expect_equal(mgr$num_datasources(), 0)
  expect_true(is.null(mgr$.ms_list[[ms_obj$get_id()]]))
})

test_that("rm_allMeasurements works without connection", {
  gr1 <- GRanges(seqnames="chr1", ranges=IRanges::IRanges(start=1:10, width=100))
  gr2 <- GRanges(seqnames="chr2", ranges=IRanges::IRanges(start=2:20, width=100))
  gr3 <- GRanges(seqnames="chr1", ranges=IRanges::IRanges(start=seq(1,100,by=25), width=1), 
                 score=rnorm(length(seq(1,100,by=25))))

  server <- epivizrServer::createServer()
  mgr <- createMgr(server)
  
  ms1 <- mgr$add_measurements(gr1, "dev1", send_request=FALSE); msId1 <- ms1$get_id()
  ms2 <- mgr$add_measurements(gr2, "dev2", send_request=FALSE); msId2 <- ms2$get_id()
  ms3 <- mgr$add_measurements(gr3, "dev3", send_request=FALSE, type="bp"); msId3 <- ms3$get_id()

  expect_equal(mgr$num_datasources(), 3)
  mgr$rm_all_measurements()
  expect_equal(mgr$num_datasources(), 0)
})

test_that("list_measurements works", {
  skip_if_not_installed("hgu133plus2.db")
  gr1 <- GRanges(seqnames="chr1", ranges=IRanges(start=1:10, width=100))
  gr2 <- GRanges(seqnames="chr2", ranges=IRanges(start=2:20, width=100))
  gr3 <- GRanges(seqnames="chr1", ranges=IRanges(start=seq(1,100,by=25), width=1), 
                 score=rnorm(length(seq(1,100,by=25))))
  eset <- make_test_eset()
  
  mgr <- createMgr(epivizrServer::createServer())
  ms1 <- mgr$add_measurements(gr1, "dev1", send_request=FALSE); msId1 <- ms1$get_id()
  ms2 <- mgr$add_measurements(gr2, "dev2", send_request=FALSE); msId2 <- ms2$get_id()
  ms3 <- mgr$add_measurements(gr3, "dev3", send_request=FALSE, type="bp"); msId3 <- ms3$get_id()
  ms4 <- mgr$add_measurements(eset, "dev4", send_request=FALSE, columns=c("SAMP_1", "SAMP_2"))
  msId4 <- ms4$get_id()
    
  ms <- mgr$list_measurements()
  expected_df <- data.frame(id=c(msId1, msId2, msId3, msId4),
                            name=paste0("dev", 1:4),
                            length=c(length(gr1), length(gr2), length(gr3), length(ms4$.object)),
                            connected=rep("", 4),
                            columns=c("","","score",paste0("SAMP_",1:2,collapse=",")),
                            stringsAsFactors=FALSE)
  expect_equal(ms, expected_df)
})
