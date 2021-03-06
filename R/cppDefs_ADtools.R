## Convert one symbol object for a C++ var into a symbol object for C++ templated CppAD code
# See symbolTable2templateTypeSymbolTable
cppVarSym2templateTypeCppVarSym <- function(oldSym, addRef = FALSE, clearRef = FALSE, replacementBaseType = 'TYPE_', replacementTemplateArgs = list()) {
    if(oldSym$baseType == 'double') {
        newSym <- cppVarFull(name = oldSym$name, baseType = replacementBaseType, ref = addRef, templateArgs = replacementTemplateArgs)
        return(newSym)
    }

    newSym <- oldSym$copy()
    if(newSym$baseType == 'NimArr') {
        if(newSym$templateArgs[[2]] == 'double') {
            if(length(replacementTemplateArgs)==0)
                newSym$templateArgs[[2]] <- replacementBaseType
            else
                newSym$templateArgs[[2]] <- cppVarFull(name='', baseType = replacementBaseType, templateArgs = replacementTemplateArgs)
            if(clearRef)
                newSym$ref <- FALSE
        }
    }
    newSym
}

## Convert a symbol table for C++ vars into a symbol table for C++ for templated CppAD code
## For CppAD, we wrap C++ code in template<class TYPE_> 
## and replace any double with TYPE_
## This includes NimArr<nDim, double> with NimArr<nDim, TYPE_>
## and similar treatmnt for Eigen templated types.
symbolTable2templateTypeSymbolTable <- function(symTab, addRef = FALSE, clearRef = FALSE, replacementBaseType = 'TYPE_', replacementTemplateArgs = list()) {
    newSymTab <- symbolTable()
    symNames <- symTab$getSymbolNames()
    for(sn in symNames) {
        oldSym <- symTab$getSymbolObject(sn)
        newSym <- cppVarSym2templateTypeCppVarSym(oldSym, addRef = addRef, clearRef = clearRef, replacementBaseType = replacementBaseType, replacementTemplateArgs = replacementTemplateArgs)
        newSymTab$addSymbol(newSym)
    }
    newSymTab
}

## This makes a Cpp function definition object wrapped in template<class TYPE_> and with
## doubles converted to TYPE_s (including in templated use if NimArr and Eigen).
## This is called from an existing version of the cppFunctionDef and returns a separate one
makeTypeTemplateFunction = function(newName, .self) { 
    newCppFunDef <- RCfunctionDef$new(static = TRUE)
    ## use typedefs to change nimble's general typedefs for Eigen locally
    typeDefs <- symbolTable()
    typeDefs$addSymbol(cppVarFull(baseType = "typedef typename EigenTemplateTypes<TYPE_>::typeEigenMapStrd", name = "EigenMapStrd") ) ## these coerces the cppVar system to generate a line of typedef code for us
    typeDefs$addSymbol(cppVarFull(baseType = "typedef typename EigenTemplateTypes<TYPE_>::typeMatrixXd", name = "MatrixXd") )
    newCppFunDef$name <- newName
    newCppFunDef$template <- cppVarFull(name = character(), baseType = 'template', templateArgs = list('class TYPE_'))
    newCppFunDef$args <- symbolTable2templateTypeSymbolTable(.self$args, addRef = TRUE)
    localArgs <- symbolTable2templateTypeSymbolTable(.self$code$objectDefs)
    newCppFunDef$returnType <- cppVarSym2templateTypeCppVarSym(.self$returnType)
    newCppFunDef$code <- cppCodeBlock(code = .self$code$code, objectDefs = localArgs, typeDefs = typeDefs, cppADCode = TRUE,
                                      generatorSymTab = .self$code$objectDefs)
    newCppFunDef
}

recurseSetCppADExprs <- function(code, logicVal = TRUE){
  if(inherits(code, 'exprClass')){
      code$cppADCode <- logicVal
    for(i in seq_along(code$args)){
      recurseSetCppADExprs(code$args[[i]], logicVal)
    }
  }
}



## This makes the function to be called once for CppAD taping
## It sets up AD variables, copies from regular variables into them
## calls the templated version of the member function
## copies the results back out.
## Not that values in the regular variables are not really important during taping.
## Currently those values are intialized to 0.5, which should satisfy needs for (-inf, inf), [0, inf) and [0, 1].
## Ending the tape is not done here.  That is done from the calling function
## (which is in permanent C++, not generated from R)
## We do not assume that in the target function the arguments are independent variables and the
## returned value is the dependent variable.  Those are set by the independentVarNames and dependentVarNames
makeADtapingFunction <- function(newFunName = 'callForADtaping', targetFunDef, ADfunName, independentVarNames, dependentVarNames, isNode, allFunDefs) {
    ## Make new function definition to call for taping (CFT)
    CFT <- RCfunctionDef$new(static = TRUE)
    CFT$returnType <- cppVarFull(baseType = "CppAD::ADFun", templateArgs = list('double'), ptr = 1, name = 'RETURN_TAPE_') ##cppVoid()
    CFT$name <- newFunName
    
    ## args will always be same; these do not depend on case.  actually now these will be empty.
    CFT$args <- symbolTable()
    ## create vector< CppAD::AD<double> > ADindependentVars
    ADindependentVarsSym <- cppVarFull(name = 'ADindependentVars', baseType = 'vector', templateArgs = list( cppVarFull(baseType = 'CppAD::AD', templateArgs = 'double', name = character()) ), ref = FALSE) ## was ref = TRUE if taking as argument
    ## create vector< CppAD::AD<double> ADresponseVars
    ADresponseVarsSym <- cppVarFull(name = 'ADresponseVars', baseType = 'vector', templateArgs = list( cppVarFull(baseType = 'CppAD::AD', templateArgs = 'double', name = character()) ), ref = FALSE) ## ditto
    ## Add them to arguments symbol table ## switch design and make these local
##    CFT$args$addSymbol( ADindependentVarsSym )
##    CFT$args$addSymbol( ADresponseVarsSym )
    ## Make local AD variables for all function inputs and outputs
    ## e.g. if the original targetFun takes NimArr<1, double>, it's templated CppAD version will take NimArr<1, TYPE_>
    ## Next line creates local variables for passing to that templated CppAD version
    localVars <- symbolTable2templateTypeSymbolTable(targetFunDef$args, clearRef = TRUE, replacementBaseType = 'CppAD::AD', replacementTemplateArgs = list('double'))
    if(isNode){
      localVars$removeSymbol('ARG1_INDEXEDNODEINFO__')
      indexNodeInfoSymbol <- symbolInternalType(name = 'ARG1_INDEXEDNODEINFO__', argList = list('indexedNodeInfoClass'))
    }
    
    ## and similar for the return variable
    initADptrCode <- cppLiteral("RETURN_TAPE_ = new CppAD::ADFun<double>;")
    ansSym <- cppVarSym2templateTypeCppVarSym(targetFunDef$returnType, clearRef = TRUE, replacementBaseType = 'CppAD::AD', replacementTemplateArgs = list('double'))
    ansSym$name <- 'ANS_'
    localVars$addSymbol(ansSym)
    symNames <- localVars$getSymbolNames()
    ## set up a set of index variables for copying code, up to six to be arbitrary (allowing up to 6-dimensional nimble objects to be handled)
    indexVarNames <- paste0(letters[9:14],'_')
    for(ivn in indexVarNames)
        localVars$addSymbol( cppVar(name = ivn, baseType = 'int') )    
    
    ## set any sizes, which must be known
    nimbleSymTab <- targetFunDef$RCfunProc$compileInfo$newLocalSymTab

    ## This creates lines like setSize(z, 2 3)
    ## which the C++ output generator turns into something like z.resize(2, 3)
    setSizeLines <- vector('list', length(symNames) + 2) ## extra 2 are for the ADindependentVars and ADresponseVars
    iNextLine <- 1
    
    for(iSym in seq_along(symNames)) {
        thisSymName <- symNames[iSym]
        if(thisSymName == 'ANS_') {
            thisSym <- targetFunDef$RCfunProc$compileInfo$returnSymbol
        } else {
            thisSym <- nimbleSymTab$getSymbolObject(thisSymName)
        }
        if(thisSym$nDim > 0) {
            setSizeCall <- do.call('call',c(list('setSize', quote(as.name(thisSymName))), as.list(thisSym$size))) 
            setSizeLines[[iNextLine]] <- setSizeCall ##RparseTree2ExprClasses(setSizeCall)
            iNextLine <- iNextLine + 1
        } else {
            setSizeLines[[iNextLine]] <- NULL
        }
    }

    localVars$addSymbol( ADindependentVarsSym )
    localVars$addSymbol( ADresponseVarsSym )
    localVars$addSymbol( CFT$returnType )

    ## call CppAD::Independent(ADindependentVars)
    ## This starts CppADs taping system
    CppADindependentCode <- quote(`CppAD::Independent`(ADindependentVars)) ##nimble:::RparseTree2ExprClasses(quote(`CppAD::Independent`(ADindependentVars)))

    ## make copying blocks into independent vars
    ## This looks like e.g.
    ## for(i_ in 1:3) {ADindependentVars[netIncrement_] = x[i]; netIncrement_ <- netIncrement + 1;}
    numIndependentVars <- length(independentVarNames)
    copyIntoIndepVarCode <- vector('list', numIndependentVars+1)
    ## create the netIncrement_ variable and code to initialize it to 1
    localVars$addSymbol( cppVar(name = 'netIncrement_', baseType = 'int') )
    copyIntoIndepVarCode[[1]] <- quote(netIncrement_ <- 1) 
    ## getting the sizes is going to be trickier when an independent var is really an expression, in particular with indexing, like model$x[3]
    ## for now let's assume only cleanly defined vars.
    ## one approach would be intermediate variables
    totalIndependentLength <- 0

    for(ivn in seq_along(independentVarNames)) {
        thisName <- independentVarNames[ivn]
        thisSym <- nimbleSymTab$getSymbolObject(thisName)
        if(thisSym$nDim > 0) {
            thisSizes <- thisSym$size
            sizeList <- lapply(thisSizes, function(x) c(1, x))
            names(sizeList) <- indexVarNames[1:length(sizeList)]
            newRcode <- makeCopyingCodeBlock(as.name(thisName), quote(ADindependentVars), sizeList, indicesRHS = FALSE, incrementIndex = quote(netIncrement_), isNode)
            copyIntoIndepVarCode[[ivn+1]] <- newRcode 
            totalIndependentLength <- totalIndependentLength + prod(thisSizes)
        } else {
            copyIntoIndepVarCode[[ivn+1]] <- substitute({LHS <- ADindependentVars[netIncrement_]; netIncrement_ <- netIncrement_ + 1}, list(LHS = as.name(thisName))) 
            totalIndependentLength <- totalIndependentLength + 1
        }
    }

    ## put dummy values in ADindependentVars
    dummyValueRcode <- substitute(for(III in 1:TOTLENGTH) ADindependentVars[III] = 0.5, list(III = as.name(indexVarNames[1]), TOTLENGTH = totalIndependentLength))
    
    if(isNode){
      dummyIndexNodeInfoCode <- list(cppLiteral('indexedNodeInfo ARG1_INDEXEDNODEINFO__ = generateDummyIndexedNodeInfo();'))
    }
    else   dummyIndexNodeInfoCode <- list()
    ## call the taping function
    TCFcall <- do.call('call', c(list(ADfunName), lapply(targetFunDef$args$getSymbolNames(), as.name)), quote = TRUE)
    tapingCallRCode <- substitute(ANS_ <- TCF, list(TCF = TCFcall))
    
    ## make copying blocks from dependent vars
    numDependentVars <- length(dependentVarNames)
    copyFromDepVarCode <- vector('list', numDependentVars+1)
    copyFromDepVarCode[[1]] <- quote(netIncrement_ <- 1) 
    totalDepLength <- 0;
    for(ivn in seq_along(dependentVarNames)) {
        thisName <- dependentVarNames[ivn]
        if(thisName == 'ANS_') {
            thisSym <- targetFunDef$RCfunProc$compileInfo$returnSymbol
        } else {
            thisSym <- nimbleSymTab$getSymbolObject(thisName)
        }
        if(thisSym$nDim > 0) {
            thisSizes <- thisSym$size
            sizeList <- lapply(thisSizes, function(x) c(1, x))
            names(sizeList) <- indexVarNames[1:length(sizeList)]
            newRcode <- makeCopyingCodeBlock(quote(ADresponseVars), as.name(thisName), sizeList, indicesRHS = TRUE, incrementIndex = quote(netIncrement_))
            copyFromDepVarCode[[ivn+1]] <- newRcode 
            totalDepLength <- totalDepLength + prod(thisSizes)
        } else {
            copyFromDepVarCode[[ivn+1]] <- substitute({ADresponseVars[netIncrement_] <- RHS; netIncrement_ <- netIncrement_ + 1}, list(RHS = as.name(thisName))) 
            totalDepLength <- totalDepLength + 1
        }
    }

    ## Now that we know how big ADindependenVars and ADresponseVars should be, 
    ## we can make two more entries to setSizeCalls for them
    ## Note that code for these will appear above code that uses them.
    setSizeLines[[iNextLine]] <- substitute(cppMemberFunction(resize(ADindependentVars, TIL)), list(TIL = totalIndependentLength))
    iNextLine <- iNextLine + 1
    setSizeLines[[iNextLine]] <- substitute(cppMemberFunction(resize(ADresponseVars, TDL)), list(TDL = totalDepLength))

    ## line to finish taping
    finishTapingCall <- cppLiteral('RETURN_TAPE_->Dependent(ADindependentVars, ADresponseVars);')

    ADoptimizeCalls <- list(cppLiteral("std::cout<<\"size before optimize = \"<< RETURN_TAPE_->size_var() <<\"\\n\";"),
                            cppLiteral("RETURN_TAPE_->optimize();"),
                            cppLiteral("std::cout<<\"size after optimize = \"<< RETURN_TAPE_->size_var() <<\"\\n\";"))

    returnCall <- cppLiteral("return(RETURN_TAPE_);")
    
    ## Finally put together all the code, parse it into the nimble exprClass system,
    ## and add it to the result (CFT)
    allRcode <- do.call('call', c(list('{'), setSizeLines, dummyIndexNodeInfoCode, list(initADptrCode, dummyValueRcode, CppADindependentCode), copyIntoIndepVarCode, list(tapingCallRCode), copyFromDepVarCode, list(finishTapingCall), ADoptimizeCalls, list(returnCall)), quote=TRUE)
    allCode <- RparseTree2ExprClasses(allRcode)
    CFT$code <- cppCodeBlock(code = allCode, objectDefs = localVars)
    CFT
}



makeStaticInitClass <- function(cppDef, derivMethods) {
    cppClass <- cppClassDef(name = paste0('initTape_', cppDef$name), useGenerator = FALSE)
    globalsDef <- cppGlobalObjects(name = paste0('initTapeGlobals_', cppDef$name))
    globalsDef$objectDefs[['staticInitClassObject']] <- cppVarFull(baseType = paste0('initTape_', cppDef$name),
                                                                   name = paste0('initTape_', cppDef$name, '_Object_'))
    initializerCodeList <- list()
    for(derivFun in derivMethods){
        initializerDef <- cppFunctionDef(name = paste0('initTape_', cppDef$name), returnType = emptyTypeInfo())
        ## use of parse instead of substitute is so R CMD check doesn't flag CLASSNAME:: as an unmentioned dependency on package named CLASSNAME
        initializerCodeList <- c(initializerCodeList,
                                 parse(text = paste0("push_back(", cppDef$name, "::allADtapePtrs_, ",
                                                     cppDef$name, "::", paste0(derivFun, "_callForADtaping_"), "() )"))[[1]])
        ## initializerCodeList <- c(initializerCodeList, substitute(push_back(CLASSNAME::allADtapePtrs_, CLASSNAME::ADTAPINGNAME() ), list(CLASSNAME = as.name(cppDef$name),ADTAPINGNAME = as.name(paste0(derivFun, "_callForADtaping_")))))
    }
    initializerCode <- do.call('call', c('{', initializerCodeList), quote = TRUE)
    initializerECcode <- RparseTree2ExprClasses(initializerCode)
    initializerDef$code <- cppCodeBlock(code = initializerECcode, objectDefs = symbolTable())
    cppClass$functionDefs[['initializer']] <- initializerDef
    cppClass$globalObjectsDefs[['globals']] <- globalsDef
    cppClass
}

makeADargumentTransferFunction <- function(newFunName = 'arguments2cppad', targetFunDef, independentVarNames, funIndex = 0, parentsSizeAndDims) {
    ## modeled closely parts of /*  */
    ## needs to set the ADtapePtr to one element of the ADtape
    TF <- RCfunctionDef$new() ## should it be static?
    TF$returnType <- cppVarFull(baseType = 'nimbleCppADinfoClass', ref = TRUE, name = 'RETURN_OBJ')
    TF$name <- newFunName
    localVars <- symbolTable() 
    isNode <- !inherits(parentsSizeAndDims, 'uninitializedField')
    if(!isNode)
      TF$args <- targetFunDef$args
    else{
      TF$args <- symbolTable()
      indexNodeInfoSym <- targetFunDef$args$getSymbolObject('ARG1_INDEXEDNODEINFO__')
      indexNodeInfoSym$name <-'INDEXEDNODEINFO_' ## to conform with original R function indexing
      TF$args$addSymbol(indexNodeInfoSym)   
    }
    
    ## set up index vars (up to 6)
    indexVarNames <- paste0(letters[9:14],'_')
    for(ivn in indexVarNames)
        localVars$addSymbol( cppVar(name = ivn, baseType = 'int') )    
    
    nimbleSymTab <- targetFunDef$RCfunProc$compileInfo$newLocalSymTab

    ## assign tape ptr code
    assignTapePtrCode <- substitute(memberData(ADtapeSetup, ADtape) <- allADtapePtrs_[FUNINDEX], list(FUNINDEX = funIndex)) ## This will have to become a unique index in general. -1 added during output
    
    ## create code to copy from arguments into the independentVars
    numIndependentVars <- length(independentVarNames)
    copyIntoIndepVarCode <- vector('list', numIndependentVars+1)
    ## create the netIncrement_ variable and code to initialize it to 1
    localVars$addSymbol( cppVar(name = 'netIncrement_', baseType = 'int') )
    copyIntoIndepVarCode[[1]] <- quote(netIncrement_ <- 1) 
    totalIndependentLength <- 0
    for(ivn in seq_along(independentVarNames)) {
        thisName <- independentVarNames[ivn]
        thisSym <- nimbleSymTab$getSymbolObject(thisName)
        if(isNode){
          nameSubList <- targetFunDef$RCfunProc$nameSubList
          thisName <- names(nameSubList)[sapply(nameSubList, function(x) return(as.character(x) == thisName))]
          thisModelElementNum <- as.numeric(gsub(".*([0-9]+)$", "\\1", thisName)) ## Extract 1, 2, etc. from end of arg name.
          thisName <- sub("_[0-9]+$","",thisName)
          thisModelName <- paste0('model_', thisName) ## Add model_ at beginning and remove _1, _2, etc. at end of arg name.
          thisSizeAndDims <- parentsSizeAndDims[[thisName]][[thisModelElementNum]]
        }
        if(thisSym$nDim > 0) {
            thisSizes <- thisSym$size
            if(isNode){
              sizeList <- list()
              for(i in 1:length(thisSizeAndDims$lengths)){
                sizeList[[i]] <-  list(thisSizeAndDims$indexExpr[[i]], 
                                   parse(text = paste0(deparse(thisSizeAndDims$indexExpr[[i]]), ' + ', thisSizeAndDims$lengths[i], ' - ', 1))[[1]])
              }
            }
            else{
              sizeList <- lapply(thisSizes, function(x) list(1, x))
            }
            names(sizeList) <- indexVarNames[1:length(sizeList)]
            newRcode <- makeCopyingCodeBlock(quote(memberData(ADtapeSetup, independentVars)), as.name(thisModelName), sizeList, indicesRHS = TRUE, incrementIndex = quote(netIncrement_), isNode)
            copyIntoIndepVarCode[[ivn+1]] <- newRcode 
            totalIndependentLength <- totalIndependentLength + prod(thisSizes)
        } 
        else {
          if(isNode){
            if(length(parentsSizeAndDims[[thisName]][[thisModelElementNum]]$lengths) > 1)
              indexBracketInfo <- paste0('(', 
                                         paste0(sapply(parentsSizeAndDims[[thisName]][[thisModelElementNum]]$indexExpr, deparse), collapse = ', '),
                                        ')')
            else
              indexBracketInfo <- paste0('[', deparse(parentsSizeAndDims[[thisName]][[thisModelElementNum]]$indexExpr[[1]]), ']')
            indexName <- paste0("cppLiteral('(**", thisModelName, ")')", indexBracketInfo)
            RHS <- parse(text = substitute(INDEXNAME, list(INDEXNAME = as.name(indexName))))[[1]]
          }
          else{
            RHS <- as.name(thisName)
          } 
          copyIntoIndepVarCode[[ivn+1]] <- substitute({memberData(ADtapeSetup, independentVars)[netIncrement_] <- RHS; netIncrement_ <- netIncrement_ + 1}, list(RHS = RHS)) 
          totalIndependentLength <- totalIndependentLength + 1
        }
    }
    setSizeLine <- substitute(cppMemberFunction(resize(memberData(ADtapeSetup, independentVars), TIL)), list(TIL = totalIndependentLength))
    returnCall <- cppLiteral("return(ADtapeSetup);")
    
    allRcode <- do.call('call', c(list('{'), list(assignTapePtrCode), list(setSizeLine), copyIntoIndepVarCode, list(returnCall)), quote=TRUE)
    allCode <- RparseTree2ExprClasses(allRcode)
    TF$code <- cppCodeBlock(code = allCode, objectDefs = localVars)
    TF
}

## Generate a block of code for copying to or from CppAD objects, to or from original C++ objects
## On the CppAD side, we are always flattening to 1D.
##
## The code this generates is embedded in the ADtapingFunction made by makeADtapingFunction
##
## Note this does some work similar to BUGScontextClass::embedCodeInForLoop
makeCopyingCodeBlock <- function(LHSvar, RHSvar, indexList, indicesRHS = TRUE, incrementIndex, isNode) {
  indexNames <- names(indexList)
  indexedBracketExpr <- do.call('call', c(list('[', as.name('TO_BE_REPLACED')), lapply(indexNames, as.name)), quote = TRUE)
  if(indicesRHS) {
    if(isNode)   RHS <- eval(substitute(substitute(indexedBracketExpr, list(TO_BE_REPLACED = cppLiteral(paste0('(**', deparse(RHSvar), ')')))), list(indexedBracketExpr = indexedBracketExpr)))
    else RHS <- eval(substitute(substitute(indexedBracketExpr, list(TO_BE_REPLACED = RHSvar)), list(indexedBracketExpr = indexedBracketExpr)))
    LHS <- substitute(A[i], list(A = LHSvar, i = incrementIndex))
  } else {
    LHS <- eval(substitute(substitute(indexedBracketExpr, list(TO_BE_REPLACED = LHSvar)), list(indexedBracketExpr = indexedBracketExpr)))
    RHS <- substitute(A[i], list(A = RHSvar, i = incrementIndex))
  }
  innerCode <- substitute({LHS <- RHS; incrementIndex <- incrementIndex + 1;}, list(LHS = LHS, RHS = RHS, incrementIndex = incrementIndex))
  for(i in length(indexList):1) {
    newForLoop <- substitute(for(NEWINDEX_ in NEWSTART_:NEWEND_) INNERCODE, list(NEWINDEX_ = as.name(indexNames[i]), NEWSTART_ = indexList[[i]][[1]], NEWEND_ = indexList[[i]][[2]], INNERCODE = innerCode))
    innerCode <- newForLoop
  }
  innerCode
}
