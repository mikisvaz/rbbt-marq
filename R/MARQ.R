library('yaml');

MARQ.config = yaml.load_file('~/.MARQ');


####################################################
# GEO platforms and datasets helper functions

MARQ.GEO.path <- function(dataset, cross_platform = FALSE, datadir = NULL){
    if (is.null(datadir) && exists('MARQ.config')){
       datadir= paste(MARQ.config$datadir, 'GEO', sep="/");
    }

    if (is.null(datadir)){
        print("No datadir specified and no default found (MARQ.config$datadir");
        exit(-1);
    }

    if ( length(grep('_cross_platform', dataset)) == 0 && cross_platform){
        dataset = paste(dataset, '_cross_platform', sep = "");
    
    }
    

    files = Sys.glob(paste(datadir,'*', '*', paste(dataset, 'orders', sep="."), sep="/"));
    
    if (length(files) == 0){
        return(NULL);
    }
    else{
        return(sub('.orders','', files[1]));
    }
}

MARQ.GEO.platform <- function(dataset, datadir = NULL){

    path = MARQ.GEO.path(dataset, datadir = datadir);

    if (is.null(path)){ return(NULL);}

    return(sub(".*(GPL\\d+).*","\\1", path, perl = TRUE));
}

MARQ.GEO.platform.path <- function(platform, datadir = NULL){
    if (is.null(datadir) && exists('MARQ.config')){
       datadir= paste(MARQ.config$datadir, 'GEO', sep="/");
    }

    if (is.null(datadir)){
        print("No datadir specified and no default found (MARQ.config$datadir");
        exit(-1);
    }

    return(paste(datadir, platform, sep="/"));
}

MARQ.GEO.platform.datasets <- function(platform, cross_platform = TRUE, series = TRUE, datadir = NULL){
    if (cross_platform){
        cp.suffix = '_cross_platform'
    }
    else{
        cp.suffix = ''
    }

    if (series){
        pattern = '*'
    }
    else{
        pattern = 'GDS'
    }
    files = Sys.glob(paste(MARQ.GEO.platform.path(platform, datadir), pattern, paste('*',cp.suffix,'.orders',sep=""),sep="/"))

    return(sapply(files, function(path){ sub(".*((?:GDS|GSE)\\d+).*", '\\1', path, perl=TRUE)}, USE.NAMES = FALSE));
}

MARQ.GEO.load <- function(dataset, cross_platform = FALSE, orders = TRUE, logratios = TRUE, t = TRUE, p.values = TRUE){
    return(MA.load(MARQ.GEO.path(dataset, cross_platform), orders, logratios, t, p.values));
}

