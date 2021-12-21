sgCMAP_score <-
function(gene.data, gene.up.list, gene.dn.list, perm=100, ref=NULL, scaling=c("by.geneset", "by.sample", "none")){

gene.data <- as.matrix(gene.data)
if( is.null(colnames(gene.data)) ) {
	colnames(gene.data) <- paste("Sample", 1:ncol(gene.data), sep="_")
}

if ( is.null(ref) ) {
	logfc.gene.data <- t(apply(gene.data, 1, function(x) c(x - median(x)) ) )
	
} else if (length(ref) == 1 & ref[1] == "intensity") {
	logfc.gene.data <- gene.data
	
} else if ( length(setdiff(ref, 1:ncol(gene.data))) == 0 ) {
	logfc.gene.data <- t(apply(gene.data, 1, function(x, ref=ref) c(x - mean(x[ref])) ) )
	
} else if ( !is.null(colnames(gene.data)) & length(setdiff(ref, colnames(gene.data) )) == 0 ) {
	logfc.gene.data <- t(apply(gene.data, 1, function(x, ref=ref) c(x - mean(x[ref])) ) )
	
} else {
	stop("Parameter \"ref\" is not set properly. The parameter can only be NULL, \"intensity\", column indices or column names of the reference samples.")
}

sig.genes <- c()
for (i in 1:length(gene.up.list) ){
	gene.up <- as.character(gene.up.list[[i]])
	gene.dn <- as.character(gene.dn.list[[i]])
	
	sig.genes <- union(sig.genes, union(gene.up, gene.dn))
}
sig.genes <- intersect(sig.genes, rownames(gene.data) )


sig.matrix <- matrix(0, nrow=length(sig.genes), ncol=length(gene.up.list))
rownames(sig.matrix) <- as.character(sig.genes)

nSet <- matrix(0, nrow=length(gene.up.list), ncol=4)
colnames(nSet) <- c("NumUpGene", "NumDnGene", "NumUpMatched", "NumDnMatched")


for (i in 1:length(gene.up.list) ){

gene.up <- as.character(gene.up.list[[i]])
gene.dn <- as.character(gene.dn.list[[i]])

gene.updn <- intersect(gene.up, gene.dn)
gene.up <- setdiff(gene.up, gene.updn)
gene.dn <- setdiff(gene.dn, gene.updn)

gene.up.index         <- match(as.character(gene.up), as.character(sig.genes) )
gene.up.matched.index <- gene.up.index[!is.na(gene.up.index)]

gene.dn.index         <- match(as.character(gene.dn), as.character(sig.genes) )
gene.dn.matched.index <- gene.dn.index[!is.na(gene.dn.index)]

gene.index <- rep(0, length(sig.genes))
gene.index[gene.up.matched.index] <- 1
gene.index[gene.dn.matched.index] <- -1

sig.matrix[, i] <- gene.index

nSet[i, "NumUpGene"] <- length(gene.up.list[[i]])
nSet[i, "NumDnGene"] <- length(gene.dn.list[[i]])
nSet[i, "NumUpMatched"] <- length(gene.up.matched.index)
nSet[i, "NumDnMatched"] <- length(gene.dn.matched.index)
}



if( !is.null(names(gene.up.list)) ) {
	colnames(sig.matrix) <- names(gene.up.list)
	rownames(nSet) <- names(gene.up.list)
} else {
	colnames(sig.matrix) <- paste("GeneSet_", seq(length(gene.up.list)), sep="")
	rownames(nSet) <- paste("GeneSet_", seq(length(gene.up.list)), sep="")
}



tmp <- .SIgN_connectivity_score(experiment=logfc.gene.data, query=sig.matrix, perm = perm)


score.matrix <- tmp$Score;
pval.matrix  <- tmp$pVal;



	rownames(score.matrix) <- colnames(sig.matrix)
	rownames(pval.matrix)  <- colnames(sig.matrix)


	if ( !is.null(colnames(gene.data)) ) {
		colnames(score.matrix) <- colnames(gene.data) 
		colnames(pval.matrix)  <- colnames(gene.data) 

	}
	
	
	if (scaling[1] == "by.geneset") {
		score.matrix <- apply(score.matrix, 1, function(x) {p <- ifelse(max(x) != 0, max(x), 1); q <- ifelse(min(x) != 0, abs(min(x)), 1); sapply(x, function(y) ifelse( y>=0, y/p, y/q ) ) } )
		pval.matrix  <- t(pval.matrix)


	} else if (scaling[1] == "by.sample") {
		score.matrix <- apply(score.matrix, 2, function(x) {p <- ifelse(max(x) != 0, max(x), 1); q <- ifelse(min(x) != 0, abs(min(x)), 1); sapply(x, function(y) ifelse( y>=0, y/p, y/q ) ) } ) 

	} else {
		score.matrix <- t(score.matrix)
		pval.matrix  <- t(pval.matrix)

	}

	return(list(Score=score.matrix, pValue = pval.matrix, nSet=nSet) )
}
