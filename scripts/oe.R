args <- commandArgs(trailingOnly = TRUE)
n_liik1 <- as.numeric(args[1])
n_liik2 <- as.numeric(args[2])
n_koos <- as.numeric(args[3])
l_liik1 <- as.numeric(args[4])
l_liik2 <- as.numeric(args[5])
l_koos <- as.numeric(args[6])
readi_pikkus <- as.numeric(args[7])
kmeeri_pikkus <- as.numeric(args[8]) 
# ----------------------------------------------------------------
# --      Funktsioon leiab usaldusintervalli O/E suhtele.       --
# ----------------------------------------------------------------
# Sisendid:
#  kontrollitavate unikaalsete k-meeride arvud
#    n_liik1 - 1. liigi (või evolutsioonilise üksuse) unikaalsete k-meeride arv
#    n_liik2 - 2. liigi (või evolutsioonilise üksuse) unikaalsete k-meeride arv
#    n_koos - unikaalsete k-meeride arv mis esinevad mõlemas harus 
#             (kuid mida ei esine teistes harudes paiknevatel liikidel/tüvedel)
#  leitud unikaalsete k-meeride arvud
#    l_liik1 - mitu 1. liigile (või evolutsioonilisele üksusele) 
#                unikaalset k-meeri leiti üles (mitut erinevat unikaalset 
#                k-meeri nähti, kordused ei lähe arvesse).
#    l_liik2 - proovist leitud 2. liigile 
#                (või evolutsioonilisele üksusele) unikaalsete k-meeride arv
#    l_koos - proovist leitud selliste unikaalsete k-meeride arv mis peaksid
#             esinema mõlemas liigis (kontrolltiavas harus)
#    readi_pikkus - kui pikki reade kasutati
#    kmeeri_pikkus - kui pikki k-meere vaadati...
#    CI_level - usaldusintervalli soovitud katvus, CI_level=0.95 vastab 95%-usaldusintervallile
#
# ------------------------------------------------------------------
# Väljund:
#    OE - O/E suhe (kui palju ühiseid unikaalseid k-meere nägime
#              vs mida oleksime oodanud üksikute liikide andmete pealt)
#    OE_UI - ligikaudne (konservatiivne) 95%-usaldusintervall O/E suhtele
#    pvalue - ligikaudne (konservatiivne) p-väärtus testi H0: O/E=1 jaoks.
# ------------------------------------------------------------------

OE = function(n_liik1, n_liik2, n_koos, l_liik1, l_liik2, l_koos, readi_pikkus, kmeeri_pikkus, CI_level=0.95){

kmeere_readis=readi_pikkus-kmeeri_pikkus+1

 # n_liik1=100000; n_liik2=90000; n_koos=100000; l_liik1=10000; l_liik2=50000; l_koos=60000; kmeere_readis=80
  # Unikaalset k-meeri nägime/ei näinud andmete põhjal 
  #   (kui paljusid erinevaid k-meere nägime)
  #   hindame selle, mitu liigi jaoks unikaalset k-meeri kokku 
  #   nägime (ühte k-meeri võisime kohata mitu korda):
  arv1_hinnang = -log(1-l_liik1/n_liik1)*n_liik1
  arv2_hinnang = -log(1-l_liik2/n_liik2)*n_liik2
  arv_koos_hinnang = -log(1-l_koos/n_koos)*n_koos

  O=l_koos/n_koos
  E=l_liik1/n_liik1+l_liik2/n_liik2-l_liik1/n_liik1*l_liik2/n_liik2
  OE=O/E

# Funktsioon mis suudab leida dispersioone suurustele nagu
#    l_liik1  või l_koos;

disp=function(O,N, readi_pikkus, k){

  p=1-O/N

  lambda_algus=-log(p)/(readi_pikkus-k+1)
  i=1:(readi_pikkus-k)

  D_O=O*(1-O/N)+sum(((2*N-2*i)>0)*(2*N-2*i)*(  exp(-lambda_algus*(readi_pikkus-k+1+i))-exp(-2*lambda_algus*(readi_pikkus-k+1))  )  )
  return(D_O)
}

#disp=function(O,N, readi_pikkus, k){
#  p=1-O/N
#
#  lambda_algus=-log(p)/(readi_pikkus-k+1)
# i=1:(readi_pikkus-k)
#
#  D_O=O*(1-O/N)+sum((2*N-2*i+1)*(  exp(-lambda_algus*(readi_pikkus-k+1+i))-exp(-2*lambda_algus*(readi_pikkus-k+1))  )  )
#  return(D_O)
#}

# Leiame suhete dispersioonid
#  D_1 = disp(arv1_hinnang, n_liik1, kmeere_readis, kmeeri_pikkus)/n_liik1**2
#  D_2 = disp(arv2_hinnang, n_liik2, kmeere_readis, kmeeri_pikkus)/n_liik2**2
#  D_koos = disp(arv_koos_hinnang, n_koos, kmeere_readis, kmeeri_pikkus)/n_koos**2

  D_1 = disp(l_liik1, n_liik1, readi_pikkus, kmeeri_pikkus)/n_liik1**2
  D_2 = disp(l_liik2, n_liik2, readi_pikkus, kmeeri_pikkus)/n_liik2**2
  D_koos = disp(l_koos, n_koos, readi_pikkus, kmeeri_pikkus)/n_koos**2

# Edasi liigume järgmiste tähelepanekute abil:

# D(XY) = D(X)D(Y)+D(X)E(Y)^2+D(Y)*E(X)^2
# D( (p1+p2)+(-p1*p2) ) =
#   = D(p1)+D(p2) + D(p1)*D(p2) + D(p1)E(p2)^2 + D(p2)E(p1)^2
#   = D(p1)+D(p2) + D(p1)*D(p2) + D(p1)E(p2)^2 + D(p2)E(p1)^2

# Delta meetod:
# Kui
#    hinnang ~ N( par; sigma^2 )
# siis
#    f(hinnang) ~ N( f(par); sigma^2*(df/dx(par))^2 )

pp1=l_liik1/n_liik1
pp2=l_liik2/n_liik2

D_OE=D_koos/O^2+(D_1+D_2+D_1*D_2+D_1*pp2^2+D_2*pp1^2-2*(D_1*pp2+D_2*pp1))/E^2

alpha=(1-CI_level)/2
z=qnorm(alpha)

al=exp(log(O/E)+z*sqrt(D_OE))
yl=exp(log(O/E)-z*sqrt(D_OE))

#print('D_OE v22rtus:')
#print(D_OE)

#2-tailed p-value
#pvalue2 = pnorm(-abs(log(O/E)),mean=0,sd=sqrt(D_OE))*2

#1-tailed p-value
pvalue = pnorm(-(log(O/E)),mean=0,sd=sqrt(D_OE))

return(list(OE=OE, pvalue=pvalue, D_OE=D_OE))
#return(list(OE=OE, D_OE=D_OE, al=al, yl=yl, pvalue=pvalue, pvalue2=pvalue2, D_koos=D_koos, D_1=D_1, D_2=D_2))
}


OE2 = function(n_liik1, n_liik2, n_koos, l_liik1, l_liik2, l_koos, readi_pikkus, kmeeri_pikkus, CI_level=0.95){
tulem1=OE(n_liik1+4, n_liik2+4, n_koos+4, l_liik1+2, l_liik2+2, l_koos+2, readi_pikkus, kmeeri_pikkus, CI_level=0.95)
tulem2=OE(n_liik1, n_liik2, n_koos, l_liik1, l_liik2, l_koos, readi_pikkus, kmeeri_pikkus, CI_level=0.95)
#return(list(OE=tulem2$OE, D_OE=tulem1$D_OE, al=tulem1$al, yl=tulem1$yl, pvalue=tulem1$pvalue, pvalue2=tulem1$pvalue2, D_koos=tulem1$D_koos, D_1=tulem1$D_1, D_2=tulem1$D_2))
return(list(OE=tulem2$OE, pvalue=tulem1$pvalue, D_OE1=tulem1$D_OE, D_OE2=tulem2$D_OE))
}

OE2(n_liik1, n_liik2, n_koos, l_liik1, l_liik2, l_koos, readi_pikkus, kmeeri_pikkus, CI_level=0.95)
#OE(100000, 30000, 100000, 40278, 5302, 50541, 101,32)
#OE2(100000, 30000, 100000, 40278, 5302, 50541, 101,32)