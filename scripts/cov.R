args <- commandArgs(trailingOnly = TRUE)
keskosakeskmine <- as.numeric(args[1])
al <- as.numeric(args[2])
yl <- as.numeric(args[3])
sagedused <- as.numeric(strsplit(args[4], "_")[[1]])
osakaalud <- as.numeric(strsplit(args[5], "_")[[1]])



Leiakeskmine2=function(keskosakeskmine, al, yl, sagedused, osakaalud){
	loikekeskmine=function(mu, al, yl, sagedused, osakaalud) {
		mitu=length(sagedused)
		tulem1=rep(NA, mitu)

		musuurem = sagedused*mu*(1+dpois(yl, lambda=sagedused*mu)/(1-ppois(yl, lambda=sagedused*mu)))
		musuurem[musuurem==Inf]=0
		musuuremal = sagedused*mu*(1+dpois(al-1, lambda=sagedused*mu)/(1-ppois(al-1, lambda=sagedused*mu)))
		if (al>0) muvaiksemal = (sagedused*mu- musuuremal *(1-ppois(al-1, lambda=sagedused*mu)))/ppois(al-1, lambda=sagedused*mu) else muvaiksemal=rep(0, mitu)
		muvahepeal = (sagedused*mu-muvaiksemal*ppois(al-1, lambda=sagedused*mu)-musuurem*(1-ppois(yl, lambda=sagedused*mu)))/ (ppois(yl, lambda=sagedused*mu)-ppois(al-1, lambda=sagedused*mu))

		p = (ppois(yl, lambda=sagedused*mu)-ppois(al-1, lambda=sagedused*mu))
		osakaalud2=prop.table(osakaalud*p)

		return(sum(muvahepeal*osakaalud2))
	}

	f=function(oletustegelikule, vaadeldud, al, yl, sagedused, osakaalud) {
		h1=loikekeskmine(oletustegelikule, al, yl, sagedused, osakaalud)
		return( (h1-vaadeldud)**2 )
	}

	tulemus=optimize(f, interval=c(0, keskosakeskmine*2+1), vaadeldud=keskosakeskmine, al=al, yl=yl, sagedused=sagedused, osakaalud=osakaalud)$minimum
	return(tulemus)
}

Leiakeskmine2(keskosakeskmine, al, yl, sagedused, osakaalud)