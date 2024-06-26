source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

## note this testing is not intended to blur into general model processing testing.
## It assumes the model processing is ok and really tests traversal of the graph to get dependencies

## also self=FALSE is not a deep testing need because this is processed in R after the deeper processing (graph traversal) in C++

## first model, with no criss-crossing dependencies.
test_that("getDependencies in first model with no criss-crossing dependencies", {
    c1 <- nimbleCode({
        a1 ~ dnorm(0,1)
        a2 ~ dnorm(0,1)
        b1 ~ dnorm(a1, 1) ## only a1 
        b2 ~ dnorm(a1, sd = a2) ## a1 and a2
        c1 <- a1 + 1 
        d1 ~ dnorm(c1, 1) ## only a1 via c1
        c2 <- a1 + 2
        c3 <- c2^2
        e1 ~ dnorm(c3, 1) ## only a1 via c2 and d1
        c4 <- a2 + 3
        f1 ~ dnorm(c3, sd = c4) ## a1 via c1 and d1; a2 via c3
        g1 ~ dnorm(f1, 1)
    })
    
    m1 <- nimbleModel(c1)
    ans1 <- m1$topologicallySortNodes(c('a1','b1','b2','c1','d1','c2','c3','e1','f1'))
    
    ## basic cases from a1
    expect_identical(m1$getDependencies('a1'), ans1)
    expect_identical(m1$getDependencies(c('a1','a1')), ans1)
    expect_identical(m1$getDependencies(c('a1','c1')), ans1)
    expect_identical(m1$getDependencies(c('c1','a1')), ans1)
    expect_identical(m1$getDependencies(c('a1', 'f1')), c(ans1, 'g1'))
    expect_identical(m1$getDependencies(c('f1', 'a1')), c(ans1, 'g1'))
    expect_identical(m1$getDependencies(c('a1'), downstream=TRUE), c(ans1, 'g1'))
    
    ## basic cases from a1 and a2
    ans2 <- m1$topologicallySortNodes(c('a1','a2','b1','b2','c4','c1','d1','c2','c3','e1','f1'))
    expect_identical(m1$getDependencies(c('a1','a2')), ans2)
    expect_identical(m1$getDependencies(c('a2','a1')), ans2)
    
    ## some omit cases
    ans3 <- m1$topologicallySortNodes(c('a1','b1','b2','c1','d1','c2'))
    expect_identical(m1$getDependencies('a1', omit = 'c3'), ans3)
    
    ans4 <- m1$topologicallySortNodes(c('a1','a2','b1','b2','c1','d1','c2','c3','e1','f1'))
    expect_identical(m1$getDependencies(c('a2','a1'), omit = 'c4', downstream = TRUE), c(ans4, 'g1'))
    
    ## some cases starting from deterministic
    expect_identical(m1$getDependencies('c1'), c('c1','d1'))
    expect_identical(m1$getDependencies('c2'), c('c2','c3','e1','f1'))
    expect_identical(m1$getDependencies(c('c1','c2')), c('c1','c2','d1','c3','e1','f1'))
})
    
################
    ## second model, with some criss-crossing dependencies
test_that("getDependencies in second model with some criss-crossing dependencies", {
    c2 <- nimbleCode({
        a1 <- 1
        b1 ~ dnorm(a1, 1)
        c1 <- b1 + 1
        c2 <- b1 + 2
        d1 ~ dnorm(c1, sd = c2)
        e1 ~ dnorm(d1, 1)
        f1 ~ dnorm(d1, sd = c1)
        g1 ~ dnorm(e1, sd = c1)
        c3 <- c1 + c2
        h1 ~ dnorm(c3, 1)
        h2 ~ dnorm(c3, sd = c1)
    })

    m2 <- nimbleModel(c2)

    ans1 <- m2$topologicallySortNodes(c('a1','b1'))
    expect_identical(m2$getDependencies('a1'), ans1)

    ans2 <- m2$topologicallySortNodes(c('b1','c1','c2','d1','f1','g1','c3','h1','h2'))
    expect_identical(m2$getDependencies('b1'), ans2)
    expect_identical(m2$getDependencies(c('b1','c1')), ans2)
    expect_identical(m2$getDependencies(c('c1','b1')), ans2)
    expect_identical(m2$getDependencies(c('b1','c1','f1')), ans2)
    expect_identical(m2$getDependencies(c('f1','c1','b1')), ans2)
    expect_identical(m2$getDependencies(c('f1','b1','c1')), ans2)

    ans3 <- m2$topologicallySortNodes(c('b1','c1','c2','d1','e1','f1','g1','c3','h1','h2'))
    expect_identical(m2$getDependencies(c('b1','c1','d1')), ans3)
    expect_identical(m2$getDependencies(rev(ans2)), ans3)
    expect_identical(m2$getDependencies(c('b1','c1','e1')), ans3)
    expect_identical(m2$getDependencies(c('e1','c1','b1')), ans3)

    ans4 <- m2$topologicallySortNodes(c('c2','d1','c3','h1','h2'))
    expect_identical(m2$getDependencies(c('c2')), ans4)
    expect_identical(m2$getDependencies(c('c3','c2')), ans4)
    expect_identical(m2$getDependencies(c('c2','h1','c3')), ans4)
})

test_that("basic use of getParents works", {
    code=nimbleCode({
        for(j in 1:J) {
            for(i in 1:I)
                y[j,i] ~ dnorm(theta[j], sigma)
            theta[j] ~ dnorm(mu, sd = tau)
        }
        mu ~ dnorm(mu0,1)
        mu0 <- mu00
        mu00 <- mu000
        sigma ~ dunif(sigma0,1)
        sigma1 ~ dunif(0,1)
        tau ~ dunif(0,1)
    })
    I <- 4
    J <- 3
    constants <- list(I = I, J = J)
    m <- nimbleModel(code,
                     data = list(y = matrix(rnorm(I*J), J, I)),
                     constants = constants)
    expect_identical(m$getParents('mu', stochOnly =  TRUE), character(0))
    expect_identical(m$getParents('mu'), c('mu00', 'mu0'))
    expect_identical(m$getParents('y', stochOnly = TRUE),
                     c('sigma','theta[1]','theta[2]','theta[3]'))
    expect_identical(m$getParents('y'),
                     c('sigma','lifted_d1_over_sqrt_oPsigma_cP', 'theta[1]','theta[2]','theta[3]'))
    expect_identical(m$getParents('y[2, 1:3]', stochOnly = TRUE),
                     c('sigma','theta[2]'))
    expect_identical(m$getParents('theta'),
                     c('tau', 'mu'))
    
})


## First model from test-getDependencies
test_that("getParents works in model with no criss-crossing dependencies", {
  c1 <- nimbleCode({
    a1 ~ dnorm(0,1)
    a2 ~ dnorm(0,1)
    b1 ~ dnorm(a1, 1) ## only a1 
    b2 ~ dnorm(a1, sd = a2) ## a1 and a2
    c1 <- a1 + 1 
    d1 ~ dnorm(c1, 1) ## only a1 via c1
    c2 <- a1 + 2
    c3 <- c2^2
    e1 ~ dnorm(c3, 1) ## only a1 via c2 and d1
    c4 <- a2 + 3
    f1 ~ dnorm(c3, sd = c4) ## a1 via c1 and d1; a2 via c3
    g1 ~ dnorm(f1, 1)
    h1 ~ dnorm(lho, 1) # left-hand-side only
  })
  m1 <- nimbleModel(c1)

  expect_identical(m1$getParents("f1"), c("a1", "a2", "c2", "c4", "c3"))
  expect_identical(m1$getParents("f1", immediateOnly = TRUE), c("c4", "c3"))
  expect_identical(m1$getParents(c("f1", "g1"), stochOnly = TRUE), c("a1", "a2", "f1"))
  expect_identical(m1$getParents(c("f1", "g1"), immediateOnly = TRUE), c("c4", "c3", "f1"))
  expect_identical(m1$getParents(c("g1", "f1"), stochOnly = TRUE), c("a1", "a2", "f1"))
  expect_identical(m1$getParents(c("f1", "g1"), stochOnly = TRUE, self = TRUE), c("a1", "a2", "f1", "g1"))
  expect_identical(m1$getParents(c("c3", "c2", "e1"), stochOnly = FALSE), c("a1", "c2", "c3"))
  expect_identical(m1$getParents(c("c3", "c2", "e1"), immediateOnly = TRUE, stochOnly = FALSE), c("a1", "c2", "c3"))
  expect_identical(m1$getParents("h1", includeRHSonly = TRUE, stochOnly = FALSE), c("lho"))
})

test_that("getParents works in model with some criss-crossing dependencies", {
  ## Second model from test-getDependencies
  c2 <- nimbleCode({
    a1 <- 1
    b1 ~ dnorm(a1, 1)
    c1 <- b1 + 1
    c2 <- b1 + 2
    d1 ~ dnorm(c1, sd = c2)
    e1 ~ dnorm(d1, 1)
    f1 ~ dnorm(d1, sd = c1)
    g1 ~ dnorm(e1, sd = c1)
    c3 <- c1 + c2
    h1 ~ dnorm(c3, 1)
    h2 ~ dnorm(c3, sd = c1)
  })
  m2 <- nimbleModel(c2)

  expect_identical(m2$getParents("h2", stochOnly = TRUE), "b1")
  expect_identical(m2$getParents("h2", stochOnly = FALSE), c("b1", "c1", "c2", "c3"))
  expect_identical(m2$getParents("h2", stochOnly = FALSE, determOnly = TRUE), c("c1", "c2", "c3"))
  expect_identical(m2$getParents("h2", stochOnly = TRUE, self = TRUE), c("b1", "h2"))
  expect_identical(m2$getParents("g1", stochOnly = TRUE), c("b1", "e1"))
  expect_identical(m2$getParents("g1", stochOnly = FALSE), c("b1", "c1", "e1"))
  expect_identical(m2$getParents("b1", stochOnly = TRUE), character())
  expect_identical(m2$getParents("b1", stochOnly = FALSE), c("a1"))
  expect_identical(m2$getParents("f1", stochOnly = TRUE), c("b1", "d1"))
  expect_identical(m2$getParents("f1", stochOnly = FALSE), c("b1", "c1", "d1"))
  expect_identical(m2$getParents(c("e1", "h1"), stochOnly = TRUE), c("b1", "d1"))
  expect_identical(m2$getParents("d1", stochOnly = TRUE), c("b1"))
  expect_identical(m2$getParents("d1", stochOnly = TRUE, upstream = TRUE ), c("b1"))
  expect_identical(m2$getParents("d1", upstream = TRUE, stochOnly = FALSE ),
                   c("a1", "b1","c1","c2"))
  expect_identical(m2$getParents("e1", upstream = TRUE, stochOnly = FALSE ),
                   c("a1", "b1", "c1", "c2", "d1"))
  expect_identical(m2$getParents("e1", omit = "d1", upstream = TRUE, stochOnly = FALSE ),
                   character()) 
  expect_identical(m2$getParents("e1", omit = "d1", stochOnly = FALSE ),
                   character()) 
  expect_identical(m2$getParents("e1", omit = c("c1", "c2"), upstream = TRUE, stochOnly = FALSE ),
                   "d1") 
  m2$setData(list(d1 = 5))
  expect_identical(m2$getParents("g1", upstream = TRUE, includeData = FALSE, stochOnly = TRUE ),
                   c("b1", "e1"))
})

test_that("getParents works in model with LHSinferred (aka split) nodes", {
  c3 <- nimbleCode({
    a[1:3] ~ dmnorm( mu[1:3], cov = cov[1:3, 1:3])
    b[1:3] <- a[1] + c(1, 2, 3)
    sig <- sqrt(var)
    var ~ dunif(0, 1)
    c[1] ~ dnorm(b[1], sd = sig)
    a2 ~ dnorm(0,1)
    sig2 ~ dunif(0,1)
    c[2] ~ dnorm(a[2], sd = sig2)
  })
  m3 <- nimbleModel(c3, inits = list(mu = 1:3, cov = diag(3)))

  expect_identical(m3$getParents("c[1]", stochOnly = TRUE), c("var", "a[1:3]"))
  expect_identical(m3$getParents("c[2]", stochOnly = TRUE), c("sig2", "a[1:3]"))
  expect_identical(m3$getParents("c[1:2]", stochOnly = TRUE), c("var", "sig2", "a[1:3]"))
  expect_identical(m3$getParents("c[1:2]", immediateOnly = TRUE, stochOnly = TRUE), c("sig2", "a[1:3]"))
  expect_identical(m3$getParents("sig", stochOnly = TRUE), c("var"))
  expect_identical(m3$getDependencies("a[1]"), c("a[1:3]", "b[1:3]", "c[1]"))
  expect_identical(m3$getDependencies("a[2]"), c("a[1:3]", "c[2]"))
  expect_identical(m3$getDependencies("b[2]"), c("b[1:3]"))
  expect_identical(m3$getDependencies("b[1]"), c("b[1:3]", "c[1]"))
  expect_identical(m3$getDependencies("b[1]", returnScalarComponents = TRUE), c("b[1]", "b[2]", "b[3]", "c[1]"))
})

test_that('getNodeNames and getDependencies with predictiveNode options', {
    
    code <- nimbleCode({
        a ~ dnorm(0, 1)
        for(i in 1:8) {
            b[i] ~ dnorm(a, 1)
            c[i] ~ dnorm(b[i], 1)
            d[i] ~ dnorm(c[i], 1)
        }
    })

    Rmodel <- nimbleModel(code)
    
    expect_identical(Rmodel$getPredictiveRootNodeIDs(), 1L)
    expect_identical(Rmodel$getPredictiveNodeIDs(), 1:length(Rmodel$getNodeNames()))
    expect_identical(Rmodel$getNodeNames(includePredictive = FALSE), character())
    expect_identical(Rmodel$getNodeNames(predictiveOnly = TRUE), Rmodel$getNodeNames())

    Rmodel$resetData()
    Rmodel$setData(list(b=1:8, c=1:8, d=1:8))
    
    expect_identical(Rmodel$getPredictiveRootNodeIDs(), integer())
    expect_identical(Rmodel$getPredictiveNodeIDs(), integer())
    expect_identical(Rmodel$getNodeNames(includePredictive = FALSE), Rmodel$getNodeNames())
    expect_identical(Rmodel$getNodeNames(predictiveOnly = TRUE), character())

    Rmodel$resetData()
    Rmodel$setData(list(
               b = rep(c(0, NA), each = 4),
               c = rep(c(0, NA, 0, NA), each = 2),
               d = rep(c(0, NA), 4)))

    expect_identical(Rmodel$getPredictiveRootNodeIDs(), as.integer(c(9, 13)))
    expect_identical(Rmodel$getPredictiveNodeIDs(), as.integer(c(9, 13, 17, 19, 21, 23, 25)))
    expect_identical(Rmodel$getNodeNames(includePredictive = FALSE), setdiff(Rmodel$getNodeNames(), Rmodel$modelDef$maps$graphID_2_nodeName[Rmodel$getPredictiveNodeIDs()]))
    expect_identical(Rmodel$getNodeNames(predictiveOnly = TRUE), Rmodel$modelDef$maps$graphID_2_nodeName[Rmodel$getPredictiveNodeIDs()])

    expect_identical(Rmodel$getDependencies('a', includePredictive = FALSE), Rmodel$expandNodeNames(c('a', 'b[1:7]')))
    expect_identical(Rmodel$getDependencies('a', predictiveOnly = TRUE), 'b[8]')

    expect_identical(Rmodel$getDependencies('b', includePredictive = FALSE), Rmodel$expandNodeNames(c('b[1:7]', 'c[1:3]', 'c[5:7]')))
    expect_identical(Rmodel$getDependencies('b', predictiveOnly = TRUE), c('b[8]', 'c[4]', 'c[8]'))

    expect_identical(Rmodel$getDependencies('c', includePredictive = FALSE), Rmodel$expandNodeNames(c('c[1:3]', 'c[5:7]', 'd[1]', 'd[3]', 'd[5]', 'd[7]')))
    expect_identical(Rmodel$getDependencies('c', predictiveOnly = TRUE), c('c[4]', 'c[8]', 'd[2]', 'd[4]', 'd[6]', 'd[8]'))

    expect_identical(Rmodel$getDependencies('d', includePredictive = FALSE), c('d[1]', 'd[3]', 'd[5]', 'd[7]'))
    expect_identical(Rmodel$getDependencies('d', predictiveOnly = TRUE), c('d[2]', 'd[4]', 'd[6]', 'd[8]'))

    Rmodel$resetData()

    expect_identical(Rmodel$getPredictiveRootNodeIDs(), 1L)
    expect_identical(Rmodel$getPredictiveNodeIDs(), 1:length(Rmodel$getNodeNames()))
    expect_identical(Rmodel$getNodeNames(includePredictive = FALSE), character())
    expect_identical(Rmodel$getNodeNames(predictiveOnly = TRUE), Rmodel$getNodeNames())
})

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
