require(RColorBrewer)
require(gplots)
require(hwriter)
require(parallel)

roastHtmlTables <- function(x,y=NULL,mat,intvar=NULL,org='human',geneid='symbol',mc.cores.x=1,mc.cores.y=1,outdir='./',stripdir='./',stripprefix='',stripid='entrez',
                            label='MaxMean',maxgs=50,padj.th=0.10,indhtml=TRUE,css='',returnData=TRUE,verbose=TRUE)
{    if (is.null(x) | is.null(mat) | is.null(intvar)) stop('Required arguments: x, mat, intvar')
    ans <- mclapply(names(x[[1]]),function(gs) ## Roast results, iterate by geneset
    {
        ans <- mclapply(names(y), function(i) ## DESeq2 results, iterate by contrast
        {
            ans <- roastHtmlTable(x[[i]][[gs]],y[[i]],mat,intvar,org,geneid,gs,i,filename=sprintf('roastGSA_%s_%s_%s.html',label,gs,i),dirname=file.path(outdir,sprintf('roastGSA/html/%s/',gs)),
                                  stripdir,stripprefix,stripid,indhtml=TRUE,maxgs=maxgs,padj.th=padj.th,verbose=verbose,css=css)
            ans
        },mc.cores=mc.cores.x)
    },mc.cores=mc.cores.y)
    if (returnData) return(ans)
}

roastHtmlTable <- function(x,y=NULL,mat,intvar,org,geneid,gs,i,filename,dirname,stripdir,stripprefix,stripid,maxgs=50,padj.th=0.10,verbose=TRUE,indhtml=TRUE,css)
{
    detable <- y
    if (indhtml & is.null(detable)) stop('indhtml is TRUE but not DE table provided')
    mygs <- x
    ## Filter for significance, if more than maxgs, maxgs
    mygs$res <- mygs$res[mygs$res$adj.pval<padj.th,]
    mygs$res <- mygs$res[order(abs(mygs$res$est),decreasing=TRUE),]
    mygs$res <- mygs$res[1:min(nrow(mygs$res),maxgs),]
    mygs$index <- mygs$index[rownames(mygs$res)]
    ## pass sortable and dragtable
    sorttable <- readLines('/Volumes/biostats/documents/bbcf_web/javascript/sorttable.js')
    dragtable <- readLines('/Volumes/biostats/documents/bbcf_web/javascript/dragtable.js')
    if (verbose) print(sprintf('Generating Roast HTML result tables for %s / %s (%s genesets selected for maxgs=%d & padj.th=%.2e)',gs,i,nrow(mygs$res),maxgs,padj.th))
    ## Write down main pathway table file
    if (indhtml) geneDEhtmlfiles <- sprintf('%s_indhtml/%s_genes.html',i,rownames(mygs$res))
    dir.create(dirname,recursive=TRUE)
    htmlrgsa2(mygs,htmlname=filename,htmlpath=dirname,
              plotpath=sprintf('%s_images/',i),indheatmap=TRUE,y=mat,intvar=intvar,ploteffsize=FALSE,mycol=redblue(100),
              geneDEhtmlfiles=geneDEhtmlfiles,sorttable=sorttable,dragtable=dragtable,whplots=rownames(mygs$res)[1:nrow(mygs$res)],
              title=sprintf('<center><h4>%s | %s</h4></center>',gs,i),css=css)
    if (indhtml) roastHtmlDETable(mygs,gs,i,org,geneid,detable,filename,dirname,stripdir,stripprefix,stripid,sorttable,dragtable,maxgs,verbose,css=css)
}


roastHtmlDETable <- function(mygs,gs,i,org,geneid,detable,filename,dirname,stripdir,stripprefix,stripid,sorttable,dragtable,maxgs,verbose=TRUE,css)
{
        outdir <- file.path(dirname,sprintf('%s_indhtml',i))
        dir.create(outdir,recursive=TRUE)
        ##st <- system.file("javascript", "sorttable.js", package = "phenoTest")
        ##dt <- system.file("javascript", "dragtable.js", package = "phenoTest")
        ##system(sprintf('cp %s %s/sorttable.js',st,outdir))
        ##system(sprintf('cp %s %s/dragtable.js',dt,outdir))
        writeLines(sorttable,sprintf('%s/sorttable.js',outdir))
        writeLines(dragtable,sprintf('%s/dragtable.js',outdir))                                
        for (j in names(mygs$index))
            {
                if ((nrow(mygs$res))>0)
                    {
                        ##intgenes <- intersect(as.character(detable$symbol),as.character(mygs$index[[j]]))
                        intgenes <- intersect(as.character(detable[,geneid]),as.character(mygs$index[[j]]))
                        ##xout <- detable[detable$symbol %in% intgenes,]
                        xout <- detable[detable[,geneid] %in% intgenes,]  
                        selcols <- c(colnames(xout)[1:4],'baseMean','log2FoldChange','stat','pvalue','padj','rej')
                        xout <- xout[,selcols]
                        xout <- xout[order(abs(xout$stat),decreasing=TRUE),]
                        isnum <- which(sapply(xout,is.numeric))
                        for (coln in isnum) xout[,coln] <- round(xout[,coln],3)
                        xout[is.na(xout)] <- 'NA'
                        ## stripcharts, first look for png, then for pdf
                        imgsrc <- sprintf('%s/%s.png',stripdir,paste0(stripprefix,xout[,stripid]))
                        ##print(sum(file.exists(file.path(outdir,imgsrc))))
                        xout$Plot <- ifelse(file.exists(file.path(outdir,imgsrc)),sprintf('<a href="%s"><img src="%s" height=50 width=50 alt="[png]"></img></a>',imgsrc,imgsrc),'-')
                        ## If stripchart in PDF is available, use it instead
                        imgsrc <- gsub('png','pdf',imgsrc)
                        xout$Plot <- ifelse(file.exists(file.path(outdir,imgsrc)),sprintf('<a href="%s"><img src="%s" height=50 width=50 alt="[png]"></img></a>',imgsrc,imgsrc),xout$Plot)
                        ## if entrez column
                        if ('entrez' %in% colnames(xout)) xout$entrez <- sprintf('<a href="https://www.ncbi.nlm.nih.gov/gene/?term=%s">%s</a>',xout$entrez,xout$entrez)
                        ## gene id link by species
                        if (org=='mouse')
                        {
                            xout[,geneid] <- sprintf('<a href="http://www.informatics.jax.org/quicksearch/summary?queryType=exactPhrase&query=%s&submit=Quick+Search">%s</a>',xout[,geneid],xout[,geneid])
                        }
                        else if (org=='human')
                            {
                                xout[,geneid] <- sprintf('<a href="http://www.informatics.jax.org/quicksearch/summary?queryType=exactPhrase&query=%s&submit=Quick+Search">%s</a>',xout[,geneid],xout[,geneid])
                            }
                        else if (org=='fly')
                        {
                            xout[,geneid] <- sprintf('<a href="http://www.informatics.jax.org/quicksearch/summary?queryType=exactPhrase&query=%s&submit=Quick+Search">%s</a>',xout[,geneid],xout[,geneid])
                        }
                        fout <- file.path(outdir, sprintf('%s_genes.html',j))
                        title <- sprintf('<center><h4>%s | %s: Genes in pathway %s\n Differential expression results</h4></center></p>',gs,i,j)
                        p <- openPage(basename(fout),dirname(fout),title=title,link.javascript=c('sorttable.js','dragtable.js') ,link.css=file.path('..',css))
                        pp <- hwrite(title,i)
                        pp <- paste(pp,hwrite(xout,table.class=list('sortable draggable'),table.align='center',row.names=FALSE),sep='\n')
                        hwrite(pp,p)
                        closePage(p)
                    }
            }
    }

write.html.mod2 <- function(x, file = paste0(htmlpath, htmlname), links = links, 
                            tiny.pic = plots, sorttable = sorttable, 
                            dragtable = dragtable,css=css,title=title)
    {
        ## Apply links and plots
        sel.links <- which(!sapply(links,is.null))
        sel.plots <- which(!sapply(tiny.pic,is.null))
        ans <- x
        for (j in sel.links) ans[,j] <- sprintf('<a href="%s">%s</a>',links[[j]],ans[,j])
        for (j in sel.plots) ans[,j] <- sprintf('<a href="%s"><img src="%s" height="100"></a>',tiny.pic[[j]],tiny.pic[[j]])
        ## Write .js files
        writeLines(sorttable,file.path(dirname(file),'sorttable.js'))
        writeLines(dragtable,file.path(dirname(file),'dragtable.js'))
        ## Write main output table
        p <- openPage(basename(file),dirname(file),link.javascript=c('sorttable.js','dragtable.js'),link.css=css)
        pp <- hwrite(title)
        pp <- paste(pp,hwrite(ans,table.class=list('sortable draggable'),table.align='center',row.names=FALSE),sep='\n')
        hwrite(pp,p)
        closePage(p)
    }

htmlrgsa2 <- function (obj, htmlpath = "", htmlname = "file.html", plotpath = "", 
    plotstats = TRUE, plotgsea = TRUE, indheatmap = TRUE, ploteffsize = TRUE, 
    links_plots = list(stats = NULL, gsea = NULL, heatmap = NULL, 
        effsize = NULL), y, whplots = NULL, geneDEhtmlfiles = NULL, 
    title = "", margins = c(15, 12), sizesHeatmap = c(12, 8), 
    typeheatmap = c("heatmap.2", "ggplot2"), intvar, adj.var = NULL, 
    mycol, varrot, psel = NULL, sorttable, dragtable, css, ...) 
{
    if (!inherits(obj, "roastgsa")) 
        stop("not a roastgsa object")
    if (ploteffsize) 
        if (missing(varrot)) 
            stop("varrot is missing")
    x <- data.frame(geneset = rownames(obj$res), obj$res)
    index <- obj$index[rownames(x)]
    psel2 <- psel
    if (plotstats | plotgsea | indheatmap) {
        dir.create(paste0(htmlpath, plotpath))
        if (is.null(whplots)) 
            whplots <- names(index)
        if (!is.na(whplots[1])) {
            stats <- sort(obj$stats)
            index <- sapply(obj$index, function(x) which(names(stats) %in% 
                x))
            for (k in whplots) {
                if (plotstats) {
                  pdf(paste0(htmlpath, plotpath, gsub("[[:punct:]]", 
                    " ", k), "_stats.pdf"))
                  plotStats(obj, whplot = k, ...)
                  dev.off()
                }
                if (plotgsea) {
                  pdf(paste0(htmlpath, plotpath, gsub("[[:punct:]]", 
                    " ", k), "_gsea.pdf"))
                  plotGSEA(obj, whplot = k, ...)
                  dev.off()
                }
                if (indheatmap) {
                  pdf(paste0(htmlpath, plotpath, gsub("[[:punct:]]", 
                    " ", k), "_heatmap.pdf"), width = sizesHeatmap[2], 
                    height = sizesHeatmap[1])
                  if (typeheatmap[1] == "ggplot2") 
                    heatmaprgsa_hm(obj, y, whplot = k, mycol = mycol, 
                      intvar = intvar, adj.var = adj.var, psel = psel2, 
                      ...)
                  else heatmaprgsa(obj, y, whplot = k, mycol = mycol, 
                    intvar = intvar, adj.var = adj.var, psel = psel2, 
                    ...)
                  dev.off()
                }
                if (ploteffsize) {
                  pdf(paste0(htmlpath, plotpath, gsub("[[:punct:]]", 
                    " ", k), "_effsize.pdf"))
                  ploteffsignaturesize(obj, varrot, whplot = k)
                  dev.off()
                }
            }
        }
    }
    if (!is.null(geneDEhtmlfiles)) 
        x$geneDEinfo <- rep("view", dim(x)[1])
    if (plotstats) 
        x$plot_stats <- NA
    if (plotgsea) 
        x$plot_gsea <- NA
    if (indheatmap) 
        x$heatmap <- NA
    if (ploteffsize) 
        x$plot_effsize <- NA
    links <- vector("list", length = ncol(x))
    names(links) <- colnames(x)
    plots <- links
    if (!is.null(links_plots$stats)) 
        links$plot_stats <- plots$plot_stats <- links_plots$stats
    else {
        if (plotstats) 
            links$plot_stats <- plots$plot_stats <- paste0(plotpath, 
                gsub("[[:punct:]]", " ", rownames(x)), "_stats.pdf")
    }
    if (!is.null(links_plots$gsea)) 
        links$plot_gsea <- plots$plot_gsea <- links_plots$gsea
    else {
        if (plotgsea) 
            links$plot_gsea <- plots$plot_gsea <- paste0(plotpath, 
                gsub("[[:punct:]]", " ", rownames(x)), "_gsea.pdf")
    }
    if (!is.null(links_plots$heatmap)) 
        links$heatmap <- plots$heatmap <- links_plots$heatmap
    else {
        if (indheatmap) 
            links$heatmap <- plots$heatmap <- paste0(plotpath, 
                gsub("[[:punct:]]", " ", rownames(x)), "_heatmap.pdf")
    }
    if (!is.null(links_plots$heatmap)) 
        links$effsize <- plots$effsize <- links_plots$effsize
    else {
        if (ploteffsize) 
            links$plot_effsize <- plots$plot_effsize <- paste0(plotpath, 
                gsub("[[:punct:]]", " ", rownames(x)), "_effsize.pdf")
    }
    if (!is.null(geneDEhtmlfiles)) {
        links$geneDEinfo <- rep(NA, dim(x)[1])
        links$geneDEinfo[1:length(geneDEhtmlfiles)] <- geneDEhtmlfiles
    }
    write.html.mod2(x, file = paste0(htmlpath, htmlname), links = links, 
        tiny.pic = plots, title = title, sorttable = sorttable, 
        dragtable = dragtable, css=css, ...)
}

