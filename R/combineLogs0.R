invisible('
Combine posterior parameter samples from beast
') 

burnpc = .25 # burn-in percentage?

#X <- do.call( rbind, lapply( 1:8, function(i){
## The 1:8 is here because Erik ran 8 MCMC runs in parallel. I ran only 1.
#	d = read.table( paste0( i,'/seir0.log'), header=TRUE, stringsAsFactors=FALSE)
#	i <- floor( burnpc * nrow(d))
#	tail( d, nrow(d) - i )
#}))

data_dir <- "~/Dropbox/SARS-CoV-2/data"

s <- read.table(file.path(data_dir,'alignment.log'), header = TRUE, stringsAsFactors = FALSE)
i <- floor(burnpc * nrow(s))
X <- tail(s, nrow(s) - i)

saveRDS(X , file = '~/Dropbox/SARS-CoV-2/Washington-sarscov2/beast/combinedLog.rds')




