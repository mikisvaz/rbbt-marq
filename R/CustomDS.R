CustomDS.process <- function(prefix, cross_platform, conditions, description, two.channel,do.log2 = NULL){
  tryCatch({
    
    # Read data
    genes = scan(file=paste(prefix,'codes',sep = "/"), what=character(),sep="\n",quiet=T);
    m = matrix(scan(file=paste(prefix,'values',sep = "/"), sep="\t",quiet=T),length(genes),byrow=T);
    rownames(m) <- genes;

   
    # Read conditions
    conditions.list = vector();
    names = vector();
    for (file in conditions) {
      filename = paste(prefix,file,sep="/");
      values = scan(file=filename, what=character(), sep = "\n");
      names = c(names, file);
      conditions.list = cbind(conditions.list, values);
    }

    colnames(conditions.list) <- names;




    # Do log2 if needed
    if (is.null(do.log2)){
      do.log2 = MA.guess.do.log2(m, two.channel);
    }
    if (do.log2){
        m <- log2(m);
    }
    m[is.infinite(m)] = NA; # The log may have created infinite values

    # Make translations to genes
    if (cross_platform){
        trans = scan(file=paste(prefix,'cross_platform',sep = "/"),what=character(),sep="\n",quiet=T);
        if (!is.null(trans)){
            m <- MA.translate(m, trans);
        }
    }

    # Once the data is read, change the prefix for cross_platform ids
    if (cross_platform){
      prefix = paste(prefix,'cross_platform', sep="_");
    }

    values <- MA.process(m, conditions.list, two.channel);

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
    rownames(best) <- rownames(ratios);
    orders   = as.data.frame(MA.get_order(best));
    colnames(orders) <- names;


    if (is.null(values)){
        write(file=paste(prefix,'skip',sep="."), "no suitable samples for analysis" );
    }else{
        MA.save(prefix, orders, ratios, t, p.values, colnames(orders), description);
    }
  },
    error=function(x){ 
        print("exception caught");
        print(x);
        write(file=paste(prefix,'skip',sep="."), paste("an exception was caught during the analysis.",x,sep="\n") );
    }
  )
}
