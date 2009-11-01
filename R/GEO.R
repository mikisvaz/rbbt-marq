library(Biobase);
library(GEOquery);

GEO.values <- function(data){
    values <- MA.process(data$m, data$conditions, data$two.channel)

    if (length(values$ratios) == 0){
        return(NULL);
    }else{
        ratios   = as.data.frame(values$ratios);
        t        = as.data.frame(values$t);
        p.values = as.data.frame(values$p.values);


        # Calculate orders from best information
        best = vector();
        names = vector();
        for (name in colnames(ratios)){
            if (sum(colnames(t) == name) > 0){
                best = cbind(best, t[,name]);
                names = c(names, name);
            }else{
                best = cbind(best, ratios[,name]);
                names = c(names, paste(name,'[ratio]', sep=" "));
            }
        }
        rownames(best) <- rownames(ratios)
        orders   = as.data.frame(MA.get_order(best));
        colnames(orders) <- names

        return(list(ratios = ratios, t = t, p.values = p.values, orders = orders));
    }
}



GEO.get <- function(name, cachedir = NULL){
    if (is.null(cachedir)){
        object  <- getGEO(name);
    }else{
        filename = dir(cachedir,pattern= paste(name, 'soft', sep='.'))[1]
        if (is.na(filename)){
            object  <- getGEO(name, destdir=cachedir);
        }else{
            object  <- getGEO(name, filename = paste(cachedir,filename,sep="/"));
        }
    }

    object
}

# Parses phenotype data structure to determine the conditions
GEO.find_conditions <- function(ph) {
    # Find conditions
    labels = varLabels(ph);

    # Removed not usefull
    labels = labels[-grep('sample|description|individual|none',labels)]

    # Find labels with replicates

    labels
}

GEO.GDS.data <- function(name, id.field = NULL, translation.file = NULL, cachedir = NULL){

    gds = GEO.get(name, cachedir);

    gpl_name = Meta(gds)$platform
    gpl = GEO.get(gpl_name, cachedir);
    description = paste(Meta(gds)$title, Meta(gds)$description, sep="\n--\n")

 
    two.channel <- attr(gds,'header')$channel_count == '2'
    value.type  <-  attr(gds,'header')$value_type
    if (length(grep('log',value.type)) == 0){
        do.log2 = TRUE;
    }else{
        do.log2 = FALSE;
    }
    eSet <- GDS2eSet(gds, GPL=gpl, do.log2 = do.log2 );

    ph = phenoData(eSet)
    conditions = sapply(GEO.find_conditions(ph), function(condition){ph[[condition]]});

    m = exprs(eSet)

    trans = NULL
    if (!is.null(id.field)){
        trans = featureData(eSet)[[id.field]];
    }
    if (!is.null(translation.file)){
        trans = scan(file=translation.file,what=character(),sep="\n",quiet=T);
    }

    if (!is.null(trans)){
        m <- MA.translate(m, trans);
    }



    return (list(conditions = conditions, m = m, two.channel = two.channel, description = description))
}


GEO.GDS.process <- function(name, prefix, id.field = NULL, translation.file = NULL,cachedir=NULL){
    tryCatch(
    {
        gds.data = GEO.GDS.data(name, id.field, translation.file, cachedir)
        values = GEO.values(gds.data)
        if (is.null(values)){
            write(file=paste(prefix,'skip',sep="."), "No suitable samples for analysis" );
        }else{
            MA.save(prefix, values$orders, values$ratios, values$t, values$p.values, colnames(values$orders), gds.data$description);
        }
    }
    ,
        error=function(x){ 
             print("Exception caught");
             print(x);
             write(file=paste(prefix,'skip',sep="."), paste("An exception was caught during the analysis.",x,sep="\n") );
        }
    )
}




GEO.GSE.data <- function(gsms, conditions, do.log2 = NULL, translation.file = NULL, use.fields = NULL, cachedir = NULL){

    c = sapply(conditions,function(x){x});
    colnames(c) <- names(conditions);
    conditions = c;

    # Get gene names
    gene.names = vector();
    for(name in gsms){
        gsm = GEO.get(name, cachedir);
        data = attr(dataTable(gsm),'table');
        if (!is.null(use.fields)){
            gsm.names = levels(factor(data[,use.fields[name]]))
        }else{
            gsm.names = levels(factor(data[,1]))
        }
        gene.names <- c(gene.names, gsm.names)
    }

    gene.names= sort(unique((gene.names)));
    gene.names = gene.names[gene.names != ""];

    
    # Merge Arrays
    init <- rep(NA,length(gene.names))
    names(init) <- gene.names;
    m <- NULL;
    for(name in gsms){
        gsm = GEO.get(name, cachedir);
        two.channel <- attr(gsm,'header')$channel_count == '2';
        data = attr(dataTable(gsm),'table')
        values = as.numeric(data[,"VALUE"])

        # Merge duplicates
       

        if (!is.null(use.fields)){
            names = factor(data[,use.fields[name]])
        }else{
            names = factor(data[,1])
        }
        
        values = tapply(values, factor(names), mean, na.rm = T)
        n <- names(values)
        values = as.vector(values)
        names(values) <- n
 
        if (is.null(m)){
            m = init
            m[names(values)] = values
        }else{
            new = init
            new[names(values)] = values;
            m = cbind(m, new)
        }
    }
    colnames(m) <- as.list(gsms);



    if (is.null(m)){
        return(NULL);
    }

    trans = NULL
    if (!is.null(translation.file)){
        trans = read.table(file=translation.file, sep="\t",header=F)[,1];
    }
    if (!is.null(trans)){
        m <- MA.translate(m, trans);
    }

    if(is.null(do.log2)){
      do.log2 = MA.guess.do.log2(m, two.channel);
    }

    if (do.log2){
        m <- log2(m);
    }

    # The log may have created infinite values
    m[is.infinite(m)] = NA

    return (list(conditions = conditions, m = m, two.channel = two.channel))
}

GEO.GSE.process <- function(gsms, conditions, prefix, do.log2 = NULL, translation.file = NULL, use.field = NULL, title = NULL,  description = NULL,cachedir=NULL){
    tryCatch(
    {
        gse.data = GEO.GSE.data(gsms, conditions, do.log2, translation.file, use.field, cachedir)
        values = GEO.values(gse.data)
        full_desc = paste(title, description, sep="\n--\n")
        if (is.null(values)){
            write(file=paste(prefix,'skip',sep="."), "No suitable samples for analysis" );
        }else{
            MA.save(prefix, values$orders, values$ratios, values$t, values$p.values, colnames(values$orders), full_desc );
        }
    }
    ,
        error=function(x){ 
             print("Exception caught");
             print(x);
             write(file=paste(prefix,'skip',sep="."), paste("An exception was caught during the analysis.",x,sep="\n") );
        }
    )

}


GEO.GPL.process <- function(name, prefix, id.field = NULL, cachedir = NULL){
    gpl = GEO.get(name, cachedir);
    if (is.null(id.field)){
        ids <- Table(gpl)$ID;
    }else{
        ids <- Table(gpl)[[id.field]];
    }

    write.table(file=paste(prefix,'codes',sep="."), ids, sep="\t",  row.names=F, col.names=F, quote=F);
}


