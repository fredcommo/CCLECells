# Use cell lines

require(synapseClient)
require(affy)
synapseLogin()

source("/Users/fredcommo/Documents/MyProjects/virtualIC50/analysis/data_utils_2.R")
source("/Users/fredcommo/Documents/MyProjects/virtualIC50/analysis/celline_2_tcga_pipeline.R")

.getDataEntity <- function(synId){
  e <- loadEntity(synId)
  Data <- read.csv(file.path(e$cacheDir, e$files), header = TRUE, sep = ' ',
                   stringsAsFactors = FALSE)
  return(Data)
}


#########################################
#########################################
#########################################
#########################################

ccle <- loadEntity('syn1671195')

setwd('/Users/fredcommo/Documents/MyProjects/CellLines/')

ccleGE <- ccle$objects$ccle_data$ccle_exp
ccleCNV <- ccle$objects$ccle_data$ccle_cnv
ccleMut <- ccle$objects$ccle_data$ccle_mut

ccleIC <- read.table('CCLE_IC50Norm.txt', sep = '\t', stringsAsFactors = FALSE)
ccleA <- read.table('CCLE_ActAreaNorm.txt', sep = '\t', stringsAsFactors = FALSE)
ccleDrug <- read.table('CCLE_DrugList.txt', sep = '\t', stringsAsFactors = FALSE)
ccleCells <- read.table('CCLE_CellsList.txt', sep = '\t', stringsAsFactors = FALSE)

# Rename ccleDrug columns
colnames(ccleDrug) <- c('inSanger', 'drugName', 'brandName', 'Target', 'Mechanism',
                        'Class', 'Phase', 'Provider')

# rename ccleCells columns
colnames(ccleCells) <- c('inSanger', 'GEid', 'cellName', 'alias', 'gender', 'tissue',
                         'histology', 'subType', 'notes','source', 'Expression.arrays',
                         'SNP.arrays', 'Oncomap', 'Hybrid.Capture.Sequencing')

ccleCells$GEid[202:203] <- c('G401_KIDNEY', 'G402_KIDNEY')
ccleCells$GEid[771] <- 'RKN_OVARY'
colnames(ccleGE)[200] <- colnames(ccleCNV)[196] <- 'G402_KIDNEY'

# Rename 'X17AAG' as '_17AAG'
colnames(ccleIC)[1] <- colnames(ccleA)[1] <- ccleDrug$drugName[1] <- '_17AAG'

# Remove TT_OESOPHAGUS from all data sets
ccleCells <- ccleCells [-which(rownames(ccleCells) == 'T0T'),]
ccleGE <- ccleGE[,-which(colnames(ccleGE) == 'TT_OESOPHAGUS')]
ccleCNV <- ccleCNV[,-which(colnames(ccleCNV) == 'TT_OESOPHAGUS')]
ccleMut <- ccleMut[,-which(colnames(ccleMut) == 'TT_OESOPHAGUS')]

# Remove OVCAR3 & SKNBE diplicated
ccleCells <- ccleCells[-c(826, 678),]

ccleCells <- ccleCells[rownames(ccleCells),]
ccleIC <- ccleIC[rownames(ccleIC),colnames(ccleIC)]
ccleA <- ccleA[rownames(ccleA), colnames(ccleA)]
ccleGE <- ccleGE[,order(colnames(ccleGE))]
ccleCNV <- ccleCNV[,order(colnames(ccleCNV))]
ccleMut <- ccleMut[,order(colnames(ccleMut))]
ccleDrug <- ccleDrug[order(ccleDrug$drugName),]

# Match colnames GE & CNV with ccleCells$CCLE.name
# Match ccleIC & ccleA with rownames(ccleCells)

cellNames <- toupper(gsub('\\W', '', rownames(ccleCells)))
ICnames <- toupper(gsub('\\W', '', rownames(ccleIC)))
ICmatch <- lapply(1:nrow(ccleCells), function(i){
  cellName <- cellNames[i]
  matchName <- NA
  cat(i, cellName, '\t')
  if(any(grep(paste0('^', cellName, '$'), ICnames)))
    matchName <- rownames(ccleIC)[grep(paste0('^', cellName, '$'), ICnames)]
  cat(matchName, '\n')
  return (data.frame(cellName = rownames(ccleCells)[i], ICid = matchName))
})
ICmatch <- do.call(rbind, ICmatch)
head(ICmatch, n = 20); dim(ICmatch); sum(!is.na(ICmatch$ICid)); nrow(ccleIC)


# Get cells for which IC50 exists and GE exists
GEnames <- toupper(gsub('_.*', '', colnames(ccleGE)))
cellMatchGE <- lapply(1:length(cellNames), function(i){
  cellName <- as.character(cellNames[i])
  matchName <- NA
  #cat(i, cellName, '\t')
  if(any(grepl(paste0('^',cellName, '$'), GEnames)))
    matchName <- colnames(ccleGE)[grep(paste0('^',cellName, '$'), GEnames)]
  #cat(matchName, '\n')
  return(data.frame(cellName = cellName, GEname = matchName))
})
cellMatchGE <- do.call(rbind, cellMatchGE)
dim(cellMatchGE); head(cellMatchGE)

CNVnames <- toupper(gsub('_.*', '', colnames(ccleCNV)))
cellMatchCNV <- lapply(1:length(cellNames), function(i){
  cellName <- as.character(cellNames[i])
  matchName <- NA
  #cat(i, cellName, '\t')
  if(any(grepl(paste0('^',cellName, '$'), CNVnames)))
    matchName <- colnames(ccleCNV)[grep(paste0('^',cellName, '$'), CNVnames)]
  #cat(matchName, '\n')
  return(data.frame(cellName = cellName, CNVname = matchName))
})
cellMatchCNV <- do.call(rbind, cellMatchCNV)
dim(cellMatchCNV); head(cellMatchCNV)

mutNames <- toupper(gsub('_.*', '', colnames(ccleMut)))
cellMatchMut <- lapply(1:length(cellNames), function(i){
  cellName <- as.character(cellNames[i])
  matchName <- NA
  #cat(i, cellName, '\t')
  if(any(grepl(paste0('^',cellName, '$'), mutNames)))
    matchName <- colnames(ccleMut)[grep(paste0('^',cellName, '$'), mutNames)]
  #cat(matchName, '\n')
  return(data.frame(cellName = cellName, mutName = matchName))
})

cellMatchMut <- do.call(rbind, cellMatchMut)
dim(cellMatchMut); head(cellMatchMut)

# Tranform ccleMut into binary values
mutNames <- rownames(ccleMut)
ccleMut <- lapply(1:nrow(ccleMut), function(i){ifelse(ccleMut[i,]=='0', 0, 1)})
ccleMut <- do.call(rbind, ccleMut)

# Add annotations in ccleCells
fullInfo <- cbind.data.frame(ccleCells[,-c(1:2)], drugTestId = ICmatch$ICid, GEid = cellMatchGE$GEname,
                             CNVid = cellMatchCNV$CNVname, MUTid = cellMatchMut$mutName)

# Check the cell names
all(rownames(ccleIC) == fullInfo$drugTestId[!is.na(fullInfo$drugTestId)])
all(rownames(ccleA) == fullInfo$drugTestId[!is.na(fullInfo$drugTestId)])
all(colnames(ccleGE) == fullInfo$GEid[!is.na(fullInfo$GEid)])
all(colnames(ccleCNV) == fullInfo$CNVid[!is.na(fullInfo$CNVid)])
all(colnames(ccleMut) == fullInfo$MUTid[!is.na(fullInfo$MUTid)])

  
# Add ICids, GE ids, CNV ids and mutation ids to Info
synIds <- data.frame(Data = c('GeneExpr', 'CopyNumber', 'Mutations', 'IC50', 'AUC', 'DrugInfo', 'CellsInfo'),
                     synId = c(NA, NA, NA, NA, NA, NA, NA))

#########################################
# Build a Robject

# Constructor.
setClass('cellLinesObj',
         representation(synId = 'data.frame',						# a data frame containing the synIds!
                        exprSet = 'data.frame',					# a data.frame containing the expression set
                        cnvSet = 'data.frame', 					# a data.frame containing the CNV.
                        mutSet = 'data.frame',   				# a data.frame containing the mutations.
                        IC50 = 'data.frame',				    # a data.frame containing the IC50.
                        AUC = 'data.frame',					    # a data.frame containing the AUC50.
                        drugInfo = 'data.frame',  			# a data.frame containing information about drugs.
                        cellsInfo = 'data.frame')				# # a data.frame containing information about the cell lines.
      )

cellsObj <- function(synId, exprSet, cnvSet, mutSet, IC50, AUC, drugInfo, cellsInfo)
        {
        new('cellLinesObj', synId = synIds, exprSet = exprSet, cnvSet = cnvSet, mutSet = mutSet,
            IC50 = IC50, AUC = AUC, drugInfo = drugInfo, cellsInfo = cellsInfo)
        }

# Accessors
setGeneric("getInfo", function(object, arg = NULL,...) standardGeneric("getInfo"))
setGeneric("getExprs", function(object, arg = NULL,...) standardGeneric("getExprs"))
setGeneric("getCNV", function(object, arg = NULL,...) standardGeneric("getCNV"))
setGeneric("getMut", function(object, arg = NULL,...) standardGeneric("getMut"))
setGeneric("getIC50", function(object, arg = NULL,...) standardGeneric("getIC50"))
setGeneric("getAUC", function(object, arg = NULL,...) standardGeneric("getAUC"))
setGeneric("getDrugs", function(object, arg = NULL,...) standardGeneric("getDrugs"))
setGeneric("getCells", function(object, arg = NULL,...) standardGeneric("getCells"))

setMethod('getInfo', signature = 'cellLinesObj', function(object, arg = NULL,...){out = data.frame(object@synId); return(out)})
setMethod('getExprs', signature = 'cellLinesObj', function(object, arg = NULL,...){return(object@exprSet)})
setMethod('getCNV', signature = 'cellLinesObj', function(object, arg = NULL,...){return(object@cnvSet)})
setMethod('getMut', signature = 'cellLinesObj', function(object, arg = NULL,...){return(object@mutSet)})
setMethod('getIC50', signature = 'cellLinesObj', function(object, arg = NULL,...){return(object@IC50)})
setMethod('getAUC', signature = 'cellLinesObj', function(object, arg = NULL,...){return(object@AUC)})
setMethod('getDrugs', signature = 'cellLinesObj', function(object, arg = NULL,...){return(object@drugInfo)})
setMethod('getCells', signature = 'cellLinesObj', function(object, arg = NULL,...){return(object@cellsInfo)})

# ShowMethod
setGeneric("show", function(object, arg = NULL,...) standardGeneric("show"))
setMethod('show', signature = 'cellLinesObj',
          function(object){
            out <- rbind.data.frame(dim(getExprs(object)), dim(getCNV(object)),
                                    dim(getMut(object)), dim(getIC50(object)),
                                    dim(getAUC(object)), dim(getDrugs(object)),
                                    dim(getCells(object)))
            colnames(out) <- c('nRow', 'nCol')
            out <- cbind.data.frame(getInfo(object), out)
            cat('\nInstance of class', class(object), 'with', length(slotNames(object)), 'slot(s)', '\n\n')
            cat('Use getInfo(object) to get the synIds')
            cat('\n\nAccessor functions:
getExprs(object), getCNV(object), getMut(object), getIC50(object),
getAUC(object), getDrugs(object), getCells(object)\n\n')
            print(out)
            })

ccleCells <- cellsObj(synId = synIds,
                        exprSet = as.data.frame(ccleGE),
                        cnvSet = as.data.frame(ccleCNV),
                        mutSet = as.data.frame(ccleMut),
                        IC50 = as.data.frame(ccleIC),
                        AUC = as.data.frame(ccleA),
                        drugInfo = as.data.frame(ccleDrug),
                        cellsInfo = as.data.frame(fullInfo))
# Store locally
setwd('/Users/fredcommo/Documents/MyProjects/CellLines/CCLECells/')
save(ccleCells, file = 'ccleCells.RData')


#############################################################
#
# Doesn't work yet!
#
#############################################################

# Store in synapse
# demo of how to create a file then use the uploaded file in a wiki
projectId <- "syn1867264"

# Store RData manually
myPath <- paste0(getwd(),'/')
myData <- Data(list(name = "ccleCells.RData", parentId = projectId))
myData <- addFile(myData, 'ccleCells.RData')
myData <- synStore(myData) # then go to synapse and upload the .RData manually

# Store RObject constructor and accessors
myCode <- Code(list(name = "cellLinesObjClass.R", parentId = projectId))
myCode <- addFile(myCode, paste0(myPath, 'cellLinesObjClass.R'))
myCode <- synStore(myCode)

# Store Rcode
myCode <- Code(list(name = "buildCCLEObj.R", parentId = projectId))
myCode <- addFile(myCode, paste0(myPath, 'buildCCLEObj.R'))
myCode <- synStore(myCode)

# Add a provenance
myData <- getEntity('syn1867264')
# used(myData) <- list(list(entity = myCode, wasExecuted=T),
#               list(entity = getEntity('syn427896'), wasExecuted=F),
#               list(entity = getEntity('syn1417761'), wasExecuted=F),
#                    list(entity = getEntity('syn1836914'), wasExecuted=F),
#                    list(entity = getEntity('syn1836934'), wasExecuted=F),
#                    list(entity = getEntity('syn1836925'), wasExecuted=F))
# myData <- synStore(myData)



