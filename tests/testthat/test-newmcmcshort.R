source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

context("Testing of default MCMC")

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)
nimblePPBranchSamplerSetting <- getNimbleOption('MCMCjointlySamplePredictiveBranches')
nimbleOptions(MCMCjointlySamplePredictiveBranches = FALSE)

## If you do *not* want to write to results files
##    comment out the sink() call below.  And consider setting verbose = FALSE 
## To record a new gold file, nimbleOptions('generateGoldFileForMCMCtesting') should contain the path to the directory where you want to put it
## e.g. nimbleOptions(generateGoldFileForMCMCtesting = getwd())
## Comparison to the gold file won't work until it is installed with the package.

goldFileName <- 'mcmcTestLog_Correct.Rout'
tempFileName <- 'mcmcTestLog.Rout'
generatingGoldFile <- !is.null(nimbleOptions('generateGoldFileForMCMCtesting'))
outputFile <- if(generatingGoldFile) file.path(nimbleOptions('generateGoldFileForMCMCtesting'), goldFileName) else tempFileName

## capture warnings
sink_with_messages(outputFile)

nimbleProgressBarSetting <- nimbleOptions('MCMCprogressBar')
nimbleOptions(MCMCprogressBar = FALSE)

## tests of classic BUGS examples

test_mcmc('blocker', numItsC = 1000, resampleData = TRUE)
# 100% coverage; looks fine

test_mcmc('bones', numItsC = 10000, resampleData = TRUE)
# 100% coverage; looks fine

test_mcmc('dyes', numItsC = 1000, resampleData = TRUE)
# 100% coverage; looks fine

test_mcmc('equiv', numItsC = 1000, resampleData = TRUE)
# looks good
# testing: tau[2]=97.95, 198.8 ; tau[1]=102.2,55
# phi = -.008,.052; pi = -.1805,.052

test_mcmc('line', numItsC = 1000, resampleData = TRUE)
# 100% coverage; looks fine

test_mcmc('oxford', numItsC = 1000, resampleData = TRUE)
# probably ok; seems to overcover for 'b', but 'b' in this
# parameteriz'n is a top-level node and the multiplic'n
# by sigma seems to lead to frequentist overcoverage
# similar results in JAGS

test_mcmc('pump', numItsC = 1000, resampleData = TRUE)
# 100% coverage; looks fine

test_mcmc('rats', numItsC = 1000, resampleData = TRUE)
# 93.8% coverage; looks fine and compares well to JAGS
# however in resampleData, one of the taus wildly misses

test_mcmc('seeds', numItsC = 1000, resampleData = TRUE)
# fine

test_mcmc('dugongs', numItsC = 1000, resampleData = TRUE)
# 100% coverage; looks fine


test_mcmc('epil', model = 'epil2.bug', inits = 'epil-inits.R',
              data = 'epil-data.R', numItsC = 1000, resampleData = TRUE)
# looks ok

test_mcmc('epil', model = 'epil3.bug', inits = 'epil-inits.R',
              data = 'epil-data.R', numItsC = 1000, resampleData = TRUE)
# looks ok

test_mcmc('seeds', model = 'seedsuni.bug', inits = 'seeds-init.R',
              data = 'seeds-data.R', numItsC = 1000, resampleData = TRUE)
# looks fine - intervals for b's seem a bit large but probably ok
# particularly since default seeds.bug seems fine
# results compared to JAGS look fine
test_mcmc('seeds', model = 'seedssig.bug', inits = 'seeds-init.R',
              data = 'seeds-data.R', numItsC = 1000, resampleData = TRUE)
# looks fine - intervals for b's seem a bit large but probably ok

test_mcmc('birats', model = 'birats1.bug', inits = 'birats-inits.R',
              data = 'birats-data.R', numItsC = 1000, resampleData = TRUE)
# seems fine

test_mcmc('birats', model = 'birats3.bug', inits = 'birats-inits.R',
              data = 'birats-data.R', numItsC = 1000, resampleData = TRUE)
# seems fine

test_mcmc('birats', model = 'birats2.bug', inits = 'birats-inits.R',
            data = 'birats-data.R', numItsC = 1000, resampleData = TRUE)
# looks fine now that values() returns in order
# result changes as of v0.4 because in v0.3-1 'omega.beta' was found
# as both topNode and nontopNode and was being simulated into
# incorrectly in resampleData - this affected values further downstream

test_mcmc('ice', model = 'icear.bug', inits = 'ice-inits.R',
          data = 'ice-data.R', numItsC = 1000, resampleData = TRUE,
          knownFailures = list('coverage' = 'KNOWN ISSUE: coverage is low'))
# resampleData gives very large magnitude betas because beta[1],beta[2] are not
# actually topNodes because of (weak) dependence on tau, and
# are simulated from their priors to have large magnitude values

test_that('ice example reworked', {
                                        # rework ice example so that beta[1] and beta[2] will be top nodes
    system.in.dir(paste("sed 's/tau\\*1.0E-6/1.0E-6/g' icear.bug > ", file.path(tempdir(), "icear.bug")), dir = system.file('classic-bugs','vol2','ice', package = 'nimble'))
    test_mcmc(model = file.path(tempdir(), "icear.bug"), inits = system.file('classic-bugs', 'vol2', 'ice','ice-inits.R', package = 'nimble'), data = system.file('classic-bugs', 'vol2', 'ice','ice-data.R', package = 'nimble'), numItsC = 1000, resampleData = TRUE, avoidNestedTest = TRUE)
                                        # looks fine, but alpha and beta values shifted a bit (systematically) relative to JAGS results - on further inspection this is because mixing for this model is poor in both NIMBLE and JAGS - with longer runs they seem to agree (as best as one can tell given the mixing without doing a super long run)
})

test_mcmc('beetles', model = 'beetles-logit.bug', inits = 'beetles-inits.R',
          data = 'beetles-data.R', numItsC = 1000, resampleData = TRUE)
                                        # getting warning; deterministic model node is NA or NaN in model initialization
                                        # weirdness with llike.sat[8] being NaN on init (actually that makes sense), and with weird lifting of RHS of llike.sat

test_that('leuk example setup', {
    writeLines(c("var","Y[N,T],","dN[N,T];"), con = file.path(tempdir(), "leuk.bug")) ## echo doesn't seem to work on Windows
                                        # need nimStep in data block as we no longer have step
    system.in.dir(paste("cat leuk.bug >> ", file.path(tempdir(), "leuk.bug")), dir = system.file('classic-bugs','vol1','leuk',package = 'nimble'))
                                        # need nimStep in data block as we no longer have step
    system.in.dir(paste("sed -i -e 's/step/nimStep/g'", file.path(tempdir(), "leuk.bug")))
    
    test_mcmc(model = file.path(tempdir(), "leuk.bug"), name = 'leuk', inits = system.file('classic-bugs', 'vol1', 'leuk','leuk-init.R', package = 'nimble'), data = system.file('classic-bugs', 'vol1', 'leuk','leuk-data.R', package = 'nimble'), numItsC = 1000,
              results = list(mean = list(beta = 1.58), sd = list(beta = 0.43)),
              resultsTolerance = list(mean = list(beta = 0.02), sd = list(beta = 0.02)), avoidNestedTest = TRUE)
})

test_that('salm example setup', {
    writeLines(paste("var","logx[doses];"), con = file.path(tempdir(), "salm.bug"))
    system.in.dir(paste("cat salm.bug >>", file.path(tempdir(), "salm.bug")), dir = system.file('classic-bugs','vol1','salm', package = 'nimble'))
    test_mcmc(model = file.path(tempdir(), "salm.bug"), name = 'salm', inits = system.file('classic-bugs', 'vol1', 'salm','salm-init.R', package = 'nimble'), data = system.file('classic-bugs', 'vol1', 'salm','salm-data.R', package = 'nimble'), numItsC = 1000, avoidNestedTest = TRUE)
                                        # looks good compared to JAGS
})

test_that('air example setup', {
    file.copy(system.file('classic-bugs','vol2','air','air.bug', package = 'nimble'), file.path(tempdir(), "air.bug"), overwrite=TRUE)
    system.in.dir(paste("sed -i -e 's/mean(X)/mean(X\\[\\])/g'", file.path(tempdir(), "air.bug")))
    test_mcmc(model = file.path(tempdir(), "air.bug"), name = 'air', inits = system.file('classic-bugs', 'vol2', 'air','air-inits.R', package = 'nimble'), data = system.file('classic-bugs', 'vol2', 'air','air-data.R', package = 'nimble'), numItsC = 1000, avoidNestedTest = TRUE)
                                        # theta[2] posterior is a bit off from JAGS - would be worth more investigation
})

test_that('jaw-linear setup', {
    system.in.dir(paste("sed 's/mean(age)/mean(age\\[1:M\\])/g' jaw-linear.bug > ", file.path(tempdir(), "jaw-linear.bug")), dir = system.file('classic-bugs','vol2','jaw', package = 'nimble')) # alternative way to get size info in there
    test_mcmc(model = file.path(tempdir(), "jaw-linear.bug"), name = 'jaw-linear', inits = system.file('classic-bugs', 'vol2', 'jaw','jaw-inits.R', package = 'nimble'), data = system.file('classic-bugs', 'vol2', 'jaw','jaw-data.R', package = 'nimble'), numItsC = 1000, avoidNestedTest = TRUE) # , knownFailures = list('R MCMC' = 'Cholesky of NA matrix fails in R 3.4.2 in calculate(model) of initializeModel() but not in R 3.4.1'))
})
## note R MCMC used to fail when tried to do Cholesky of 0 matrix in 2-point method, but no longer doing multiplicative link for Wishart targets
                                      
test_mcmc('pump',
          resampleData = TRUE,
          results = list(mean = list(
                             "theta[1]" = 0.06,
                             "theta[2]" = 0.10,
                             "theta[9]" = 1.58,
                             "theta[10]" = 1.97,
                             alpha = 0.73,
                             beta = 0.98)),
          resultsTolerance = list(mean = list(
                                      "theta[1]" = 0.01,
                                      "theta[2]" = 0.01,
                                      "theta[9]" = 0.05,
                                      "theta[10]" = 0.05,
                                      alpha = 0.1,
                                      beta = 0.1)))


test_that('gap setup', {
    ## LogProb gap: bug fixed in after v0.3
    ## Problem that occurred in v0.3: because of gap in logProb_a (i.e. logProb_a[2]
    ## is defined but logProb_a[1] is not)
    ## Because logProbs get scrambled, the random walk sampler would always accept,
    ## meaning the sd of proposal steps approaches Inf
    gapCode <- nimbleCode({
	a[1] <- 1
	a[2] ~ dnorm(0,1)
    })
    
    test_mcmc(model = gapCode, seed = 0, numItsC = 100000,
              results = list(mean = list(`a[2]` = 0) ),
              resultsTolerance = list(mean = list(`a[2]` = 0.1)),
              samplers = list(list(type = 'RW', target = 'a[2]')),
              avoidNestedTest = TRUE
              )
})

if(.Platform$OS.type == 'windows') {
    message("Stopping tests now in Windows to avoid crashing until we can unload compiled projects")
    message("To continue testing use 'mcmc2' tests")
    q("no")
}

### Daniel's world's simplest MCMC demo

test_that('very simple example setup', {
    code <- nimbleCode({
        x ~ dnorm(0, 2)
        y ~ dnorm(x+1, 3)
        z ~ dnorm(y+2, 4)
    })
    data = list(y = 3)
    
    test_mcmc(model = code,
              name = 'very simple example',
              data = data,
              resampleData = FALSE,
              results = list(
                  mean = list(x = 6/5, z = 5),
                  sd = list(x = 1/sqrt(5), z = 1/2)),
              resultsTolerance = list(mean = list(x = .1, z = .1),
                                      sd = list(x = .05, z = .05)),
              avoidNestedTest = TRUE)
})

### linear Gaussian state-space model (of length 5)

test_that('linear Gaussian state-space model MCMC works', {
  set.seed(0)
  n <- 5
  a <- 6
  b <- 0.8
  sigmaPN <- 2
  sigmaOE <- 4
  set.seed(0)
  x <- numeric(n)
  y <- numeric(n)
  x[1] <- 1
  y[1] <- rnorm(1, x[1], sigmaOE)
  for(i in 2:n) {
    x[i] <- rnorm(1, a+b*x[i-1], sigmaPN)
    y[i] <- rnorm(1, x[i], sigmaOE)
  }
  code <- nimbleCode({
    a ~ dnorm(0, sd = 10000)
    b ~ dnorm(0, sd = 10000)
    sigmaPN ~ dunif(0, 10000)
    sigmaOE ~ dunif(0, 10000)
    x[1] ~ dnorm(0, sd = 10000)
    y[1] ~ dnorm(x[1], sd = sigmaOE)
    for(t in 2:N) {
      x[t] ~ dnorm(a + b*x[t-1], sd = sigmaPN)
      y[t] ~ dnorm(x[t], sd = sigmaOE)
    }
  })
  constants <- list(N = length(y))
  data <- list(y = y)
  inits <- list(a=6, b=0.8, sigmaOE=4, sigmaPN=2, x=y+rnorm(length(y)))

  test_mcmc(model = code,
            name = 'linear Gaussian state-space model',
            data = c(data, constants),
            inits = inits,
            seed = 123,
            resampleData = FALSE,
            results = list(
              mean = list(a = 37.36, b = -2.26, sigmaPN = 5.27, sigmaOE = 5.08),
              sd = list(a = 30.14, b = 2.60, sigmaPN = 6.54, sigmaOE = 2.96)),
            ## The expected results come from exactly these run conditions,
            ## so they are essentially exact MCMC chain replication results.
            ## For that reason, tolerances are all set to the precision with
            ## which the numbers were entered, to nearest 0.01.
            resultsTolerance = list(mean = list(a = .01, b = .01, sigmaPN=0.01, sigmaOE=0.01),
                                    sd = list(a = .01, b = .01, sigmaPN=0.01, sigmaOE=0.01)),
            avoidNestedTest = TRUE)
})


### basic block sampler example

test_that('basic no-block sampler setup', {
    code <- nimbleCode({
        for(i in 1:3) {
            x[i] ~ dnorm(0, 1)
            y[i] ~ dnorm(x[i], 2)
        }
    })
    data = list(y = -1:1)
    
    test_mcmc(model = code,
              name = 'basic no-block sampler',
              data = data,
              resampleData = FALSE,
              results = list(
                  mean = list(x = c(-2/3,0,2/3)),
                  var = list(x = rep(1/3,3))),
              resultsTolerance = list(mean = list(x = rep(.1,3)),
                                      var = list(x = rep(.05,3))),
              avoidNestedTest = TRUE)
    
    
    test_mcmc(model = code,
              name = 'basic block sampler on scalars',
              data = data, resampleData = FALSE,
              results = list(
                  mean = list(x = c(-2/3,0,2/3)),
                  var = list(x = rep(1/3,3))),
              resultsTolerance = list(mean = list(x = rep(.1,3)),
                                      var = list(x = rep(.05,3))),
              samplers = list(
                  list(type = 'RW_block', target = 'x[1]'),
                  list(type = 'RW_block', target = 'x[2]'),
                  list(type = 'RW_block', target = 'x[3]')
              ), removeAllDefaultSamplers = TRUE, numItsC = 10000,
              avoidNestedTest = TRUE)
    
    test_mcmc(model = code,
              name = 'basic block sampler on vector',
              data = data,
              resampleData = FALSE,
              results = list(
                  mean = list(x = c(-2/3,0,2/3)),
                  var = list(x = rep(1/3,3))),
              resultsTolerance = list(mean = list(x = rep(.1,3)),
                                      var = list(x = rep(.05,3))),
              samplers = list(
                  list(type = 'RW_block', target = 'x', control = list(adaptInterval = 500))
              ), numItsC = 10000, avoidNestedTest = TRUE)  
})

### slice sampler example

test_that('slice sampler example setup', {
    code <- nimbleCode({
        z ~ dnorm(0, 1)
        normal5_10 ~ dnorm(5, sd = 10)
        beta1_1 ~ dbeta(1, 1)
        beta3_5 ~ dbeta(3, 5)
        binom10_p5 ~ dbin(size=10, prob=0.5)
        binom20_p3 ~ dbin(size=20, prob=0.3)
    })
    
    test_mcmc(model = code,
              name = "slice sampler example",
              resampleData = FALSE,
              results = list(
                  mean = list(z = 0, "beta1_1" = 0.5, "beta3_5" = 3/(3+5),
                              "binom10_p5" = 10*.5, "binom20_p3" = 20*.3),
                  sd = list(z = 1, "beta1_1" = sqrt(1/12),
                            "beta3_5" = sqrt(3*5/((3+5)^2*(3+5+1))),
                            "binom10_p5" = sqrt(10*.5*.5),
                            "binom20_p3" = sqrt(20*.3*.7))),
              resultsTolerance = list(
                  mean = list(z = 0.1, "beta1_1" = 0.5, "beta3_5" = .2,
                              "binom10_p5" = .25, "binom20_p3" = .25),
                  sd = list(z = .1, "beta1_1" = .05, "beta3_5" = .03,
                            "binom10_p5" = .2, "binom20_p3" = .25)),
              samplers = list(list(type = 'slice', target = 'z', control = list(adaptInterval = 10)),
                              list(type = 'slice', target = 'normal5_10', control = list(adaptInterval = 10)),
                              list(type = 'slice', target = 'beta1_1', control = list(adaptInterval = 10)),
                              list(type = 'slice', target = 'beta3_5', control = list(adaptInterval = 10)),
                              list(type = 'slice', target = 'binom10_p5', control = list(adaptInterval = 10)),
                              list(type = 'slice', target = 'binom20_p3', control = list(adaptInterval = 10))),
              avoidNestedTest = TRUE)
    })


### elliptical slice sampler 'ess'

test_that('elliptical slice sampler setup', {
    set.seed(0)
    ESScode <- quote({
        x[1:d] ~ dmnorm(mu_x[1:d], prec = prec_x[1:d, 1:d])
        y[1:d] ~ dmnorm(x[1:d], prec = prec_y[1:d, 1:d])
    })
    d <- 3
    mu_x <- rnorm(d)
    temp <- array(rnorm(d^2), c(d,d))
    prec_x <- solve(temp %*% t(temp))
    temp <- array(rnorm(d^2), c(d,d))
    prec_y <- solve(temp %*% t(temp))
    y <- rnorm(d)
    ESSconstants <- list(d = d, mu_x = mu_x, prec_x = prec_x, prec_y = prec_y)
    ESSdata <- list(y = y)
    ESSinits <- list(x = rep(0, d))
    
    test_mcmc(model = ESScode, data = c(ESSconstants, ESSdata), inits = ESSinits,
              name = 'exact values of elliptical slice sampler',
              seed = 0,
              exactSample = list(`x[1]` = c(-0.492880566939352, -0.214539223107114, 1.79345037297218, 1.17324496091208, 2.14095077672555, 1.60417482445964, 1.94196916651627, 2.66737323347255, 2.66744178776022, 0.253966883192744), `x[2]` = c(-0.161210109217102, -0.0726534676226932, 0.338308532423757, -0.823652445515156, -0.344130712698579, -0.132642244861469, -0.0253168895009594, 0.0701624130921676, 0.0796842215444978, -0.66369112443311), `x[3]` = c(0.278627475932455, 0.0661336950029345, 0.407055002920732, 1.98761228946318, 1.0839897275519, 1.00262648370199, 0.459841485268785, 2.59229443025387, 1.83769567435409, 1.92954706515119)),
              samplers = list(list(type = 'ess', target = 'x')))
    
    test_mcmc(model = ESScode, data = c(ESSconstants, ESSdata), inits = ESSinits,
              name = 'results to tolerance of elliptical slice sampler',
              results = list(mean = list(x = c(1.0216463, -0.4007247, 1.1416904))),
              resultsTolerance = list(mean = list(x = c(0.01, 0.01, 0.01))),
              numItsC = 100000,
              samplers = list(list(type = 'ess', target = 'x')), avoidNestedTest = TRUE)
    })


### demo2 of check conjugacy

test_that('beta-binom conjugacy setup', {
    code <- nimbleCode({
        x ~ dbeta(3, 13)
        y[1] ~ dbin(x, 10)
        y[2] ~ dbin(x, 20)
    })
    data = list(y = c(3,4))
    
    test_mcmc(model = code, name = 'check of beta-binom conjugacy', data = data, exactSample = list(x = c(0.195510839527966, 0.332847482503424,0.247768152764931, 0.121748195439553, 0.157842271774841, 0.197566496350904, 0.216991517500577, 0.276609942874852, 0.165733872345582, 0.144695512780252)), seed = 0, avoidNestedTest = TRUE)
})
### checkConjugacy_demo3_run.R - various conjugacies

test_that('various conjugacies setup', {
    code <- nimbleCode({
        x ~ dgamma(1, 1)       # should satisfy 'gamma' conjugacy class
        a  ~ dnorm(0, x)     # should satisfy 'norm' conjugacy class
        a2 ~ dnorm(0, tau = 3*x+0)
        b  ~ dpois(0+5*x)
        b2 ~ dpois(1*x*1)
        c ~ dgamma(1, 7*x*5)
        for(i in 2:3) {
            jTau[i] <- 1
            jNorm[i] ~ dnorm(c * (a+3) - i, var = jTau[i])
            kTauSd[i] <- 2
            kLogNorm[i] ~ dlnorm(0 - a - 6*i, kTauSd[i])
        }
        jNorm[1] <- 0
        kLogNorm[1] <- 0
    })

    sampleVals = list(x = c(3.950556165467749, 1.556947815895538, 1.371834934033851, 2.036442813764752, 2.247416118159410, 2.537131924778210, 2.382184991769738, 2.653737836857812, 2.934255734970981, 3.007873553270551),
                      c = c(0.010341199485849559, 0.010341199485849559, 0.003846483017887228, 0.003846483017887228, 0.003846483017887228, 0.006269117826484087, 0.009183580181658716, 0.009183580181658716, 0.006361841408434201, 0.006361841408434201))
    
    test_mcmc(model = code, name = 'check various conjugacies', exactSample = sampleVals, seed = 0, mcmcControl = list(scale=0.01), avoidNestedTest = TRUE)
    ## with fixing of jNorm[1] and kLogNorm[1] we no longer have: knownFailures = list('R C samples match' = "KNOWN ISSUE: R and C posterior samples are not equal for 'various conjugacies'"))
})

### Weibull-gamma conjugacy
test_that('Weibull-gamma conjugacy setup', {
    y <- 3; depShape <- 2; c <- 2; shape <- 1; rate <- 2
    code <- nimbleCode({
        y ~ dweib(shape = depShape, lambda = c*theta)
        theta ~ dgamma(shape, rate = rate)
    })
    m <- nimbleModel(code, data = list(y = y), inits = list(c = c, theta = 1),
                     constants = list(depShape = depShape, shape = shape, rate = rate))
    conf <- configureMCMC(m)
    samplers <- conf$getSamplers()
    expect_identical(samplers[[1]]$name, 'conjugate_dgamma_dweib_multiplicative',
                                   info = "dweibull-dgamma conjugacy with dependency using lambda not detected")
    mcmc <- buildMCMC(conf)
    comp <- compileNimble(m, mcmc)
    set.seed(0)
    comp$mcmc$run(10)
    smp <- as.matrix(comp$mcmc$mvSamples)

    manualSampler <- function(n, y, depShape, c, shape, rate) {
        out <- rep(0, n)
        shape = shape + 1
        rate = rate + c*y^depShape
        set.seed(0)
        out <- rgamma(n, shape, rate = rate)
        return(out)
    }
    smpMan <- manualSampler(10, y, depShape, c, shape, rate)

    expect_identical(smp[,1], smpMan,
                                   info = "NIMBLE gamma-Weibull conjugate sampler and manual sampler results differ")
})    


sink(NULL)


options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
nimbleOptions(MCMCprogressBar = nimbleProgressBarSetting)
nimbleOptions(MCMCjointlySamplePredictiveBranches = nimblePPBranchSamplerSetting)
