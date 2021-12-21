

.SIgN_connectivity_score <-
function( experiment, query, perm=100) {

	data.matrix <- as.matrix(experiment)
	query       <- as.matrix(query)


	matched.sets <- query

	sets.up     <- lapply( seq(ncol(matched.sets)),
                              function( x ) match( rownames(matched.sets)[which(matched.sets[ ,x ] == 1 )], rownames(data.matrix)) )
    
    sets.down     <- lapply( seq(ncol(matched.sets)),
                              function( x ) match( rownames(matched.sets)[which(matched.sets[ ,x ] == -1)], rownames(data.matrix)) )


    rank.matrix <- apply(data.matrix, 2, function(x) { length(x) - rank(x) + 1 } )


    scorelist <- lapply( 1:ncol(rank.matrix), function( j, rankMat, sigs.up, sigs.down )
                                                  sapply(seq_along( sigs.up ), function( n, rMat, gs.up, gs.down ) 
                                                                                   .s( rMat[gs.up[[n]], j], rMat[gs.down[[n]], j], nrow( rMat ) ), 
                                                                                   rMat=rankMat, gs.up=sigs.up, gs.down=sigs.down), 
                                                  rankMat=rank.matrix, sigs.up=sets.up, sigs.down=sets.down )
	
	
    raw.score   <- sapply(scorelist, function(x) x[3, ])
	raw.score   <- matrix(raw.score,   ncol=ncol(experiment) )




   
	score <- raw.score
	
	if ( perm > 0 ) {
	raw.score0 <- apply( rank.matrix, 2, function( r ) {
      sapply(seq_along( sets.up ), function( n ) {
        .ks( c(r[sets.up[[n]]], r[sets.down[[n]]]), length( r ) )
      })
    })
    raw.score0 <- matrix(raw.score0, ncol=ncol( experiment ))

	
	permMat <- matrix(0, nrow=length(sets.up), ncol=ncol(rank.matrix) )
	for (iter in 1:perm) { 
	
	rand.raw.score <- apply( rank.matrix, 2, function( r ) {
      sapply(seq_along( sets.up ), function( n ) {
	    rand.set <- sample( 1:nrow(rank.matrix), (length(sets.up[[n]]) + length(sets.down[[n]])), replace=FALSE)
   
		.ks( r[ rand.set ], length( r ) )
      })
    })
    rand.raw.score <- matrix(rand.raw.score, ncol=ncol( experiment ))
	
	for (i in 1:nrow(permMat) ) {
		for (j in 1:ncol(permMat)){

			if (abs(rand.raw.score[i,j]) >= abs(raw.score0[i, j])) {
				permMat[i,j] <- permMat[i,j] + 1
			}
				
		}
	}
	
	} 
	
	pvalMat <- permMat / perm

	results <- list (Score = score, pValue = pvalMat)
	
   } else {
   
	results <- list (Score = score, pValue = NULL)
	
	}
	
	results <- list (Score = score, pVal = pvalMat)
	return (results)
  }








.ks <-
function( V, n ) {
  t <- length( V )
  
  if( t == 0 )  {
    return( 0 )
  } else {
    
    if ( is.unsorted( V ) )
      V <- sort( V )
    d <- (1:t) / t - V / n
    a <- max( c(0, d))
	b <- max(c(0, (-min(d) + 1.0/t) ))
	ifelse( a >= b, a, -b )
  }
}

.s <-
function( V_up, V_down, n ) {
  ks_up   <- .ks( V_up, n )
  ks_down <- .ks( V_down, n )
  ks <- ifelse( sign( ks_up ) == sign( ks_down ), 0, ks_up - ks_down )

  return( c(ks_up, ks_down, ks) )
}


.S <-
function( scores ) {
  p <- max( scores )
  q <- min( scores )
  ifelse(
         scores == 0,
         0,
         ifelse( scores > 0, scores / p, -scores / q )
         )
}
