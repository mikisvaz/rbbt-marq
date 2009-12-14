library(GEOquery);

# The original version of the function failed if the dataset had extra probe ids not in the platform.
# This version fixes that
"GDS2eSet" <-
  function(GDS,do.log2=FALSE,GPL=NULL,AnnotGPL=TRUE) {
    require(Biobase)
                                        # exclude non-numeric columns
    if(is.null(GPL)) {
      GPL <- getGEO(Meta(GDS)$platform,AnnotGPL=AnnotGPL)
    }
    ord.table <- match(Table(GDS)[,1],Table(GPL)[,1])
    inc.columns <- grep('GSM',colnames(Table(GDS)))
    mat <- suppressWarnings(as.matrix(apply(Table(GDS)[,inc.columns],2,
                                            function(x) {as.numeric(as.character(x))})))
    if(do.log2) {
      expr <- log2(mat)
    } else {
      expr <- mat
    }
    rownames(expr) <- as.character(Table(GDS)$ID_REF)
    tmp <- Columns(GDS)
    rownames(tmp) <- as.character(tmp$sample)
    pheno <- new("AnnotatedDataFrame",data=tmp)
    mabstract=ifelse(is.null(Meta(GDS)$description),"",Meta(GDS)$description)
    mpubmedids=ifelse(is.null(Meta(GDS)$pubmed_id),"",Meta(GDS)$pubmed_id)
    mtitle=ifelse(is.null(Meta(GDS)$title),"",Meta(GDS)$title)
    dt <- Table(GPL)
    rownames(dt) <- as.character(dt$ID)
    featuredata <- new('AnnotatedDataFrame',data=dt[ord.table,],
                       varMetadata=data.frame(Column=Columns(GPL)[,1],
                         labelDescription=Columns(GPL)[,2]))

    # use !is.na(ord.table) to remove extra probe ids in GDS and not in GPL
    eset <- new('ExpressionSet',exprs=expr[!is.na(ord.table),],phenoData=pheno,
                featureData=featuredata[!is.na(ord.table),],
                experimentData=new("MIAME",
                  abstract=mabstract,
                  title=mtitle,
                  pubMedIds=mpubmedids,
                  other=Meta(GDS)))
    return(eset)
  }

