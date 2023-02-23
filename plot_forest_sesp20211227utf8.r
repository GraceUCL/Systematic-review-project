plotsesp=function(sens=sensitivity_all,spec=specificity_all,digit=3,savepdf=FALSE,pdfw=16,pdfh=14,boxsize=0.025,overall=FALSE,headline=TRUE,noRawdata=TRUE,gtestgray=50){
	if(sens$level!=spec$level){stop("Sensitivity and Specificity has different levels")}
	if(sens$level.comb!=spec$level.comb){stop("Sensitivity and Specificity has different levels for pooled estimates")}	
	if(sens$bylab!=spec$bylab){stop("Sensitivity and Specificity has different group variables")}
	GN=sens$bylab
	GL=sens$bylevs
	GV=sens$byvar
	if(any(GN!=spec$bylab,GL!=spec$bylevs,GV!=spec$byvar)){stop("not all identical of group")}
	rV=function(x)format(round(x,digit), nsmall = digit)
	as.num=function(x)as.numeric(strsplit(gsub("\\]",replacement="",gsub("\\[",replacement="",x)),split=",",fixed=T)[[1]])
	sensCI.U=0;sensCI.L=1;specCI.U=0;specCI.L=1
	plotM=array(" ",c(overall*3+(length(GL)+1)*4+length(sens$n),13))
	r=1
	plotM[r,]=c("Study","TP","FP","FN","TN",
			"Sensitivity","Proportion",paste(round(sens$level*100,1),"%-CI",sep=""),
			"Specificity","Proportion",paste(round(spec$level*100,1),"%-CI",sep=""),2,1)
	r=r+1+headline
	for(g in GL){
		plotM[r,c(1,12,13)]=c(paste(GN,"=",g),2,2);r=r+1
		gint=which(GV==g);gl=length(gint)
		plotM[r:(r+gl-1),]=cbind(
			sens$studlab[gint],
			sens$event[gint],spec$n[gint]-spec$event[gint],#TP,FP
			sens$n[gint]-sens$event[gint],spec$event[gint],#FN,TN
			rep(1,gl),
			rV(sens$event[gint]/sens$n[gint]),
			paste("[",rV(sens$lower[gint]),",",rV(sens$upper[gint]),"]",sep=""),
			rep(1,gl),
			rV(spec$event[gint]/spec$n[gint]),
			paste("[",rV(spec$lower[gint]),",",rV(spec$upper[gint]),"]",sep=""),
			array(1,c(gl,2))
			);r=r+gl
		sensCI.L=min(sensCI.L,sens$lower[gint]);sensCI.U=max(sensCI.U,sens$upper[gint])
		specCI.L=min(specCI.L,spec$lower[gint]);specCI.U=max(specCI.U,spec$upper[gint])
		varint <- GL==g
		sens.propCI=rV(meta:::backtransf(
					c(sens$TE.random.w[varint],sens$lower.random.w[varint],sens$upper.random.w[varint]),
					sm=sens$sm))
		spec.propCI=rV(meta:::backtransf(
					c(spec$TE.random.w[varint],spec$lower.random.w[varint],spec$upper.random.w[varint]),
					sm=spec$sm))
		plotM[r,]=c(
			paste("Random effects model (",round(sens$level.comb*100,1),"%-CI)",sep=""),
			sens$event.w[varint],spec$n.w[varint]-spec$event.w[varint],
			sens$n.w[varint]-sens$event.w[varint],spec$event.w[varint],
			"2",sens.propCI[1],paste("[",sens.propCI[2],",",sens.propCI[3],"]",sep=""),
			"2",spec.propCI[1],paste("[",spec.propCI[2],",",spec.propCI[3],"]",sep=""),
			"2","2"
			);r=r+1
		plotM[r,c(1,2,12,13)]=c("Heterogeneity:",g,"1","2");r=r+2
		}
	sens.N.propCI=rV(meta:::backtransf(
			c(sens$TE.random,sens$lower.random,sens$upper.random),
			sm=sens$sm))
	spec.N.propCI=rV(meta:::backtransf(
			c(spec$TE.random,spec$lower.random,spec$upper.random),
			sm=spec$sm))
	if(overall){
		plotM[r,]=c(
				paste("Random effects model (",round(sens$level.comb*100,1),"%-CI)",sep=""),
				sum(sens$event),sum(spec$n)-sum(spec$event),#TP,FP
				sum(sens$n)-sum(sens$event),sum(spec$event),#FN,TN
				"2",sens.N.propCI[1],paste("[",sens.N.propCI[2],",",sens.N.propCI[3],"]",sep=""),
				"2",spec.N.propCI[1],paste("[",spec.N.propCI[2],",",spec.N.propCI[3],"]",sep=""),
				"2","1"
				);r=r+3
		plotM[r,c(1,2,12,13)]=c("Heterogeneity:","T","1","1");r=r+1
		plotM[r,c(1,2,12,13)]=c("Test for subgroup differences:","T","1","1")
		}
	
	Lscal=6
	xps=cumsum(c(0,5,0,0,0,0,Lscal,1,2,Lscal,1,2)+c(0,0,1,1,1,1,0,0,0,0,0,0)*!noRawdata )
	adjX=c(0,1,1,1,1,0.5,1,0.5,0.5,1,0.5)
	xpos<-xps
	xpos[which(adjX==1)]=xps[which(adjX==1)+1]
	xpos[which(adjX==0.5)]=(xps[which(adjX==0.5)]+xps[which(adjX==0.5)+1])/2
	ypos=0:-(dim(plotM)[1]+1)
	default.x=range(xpos);default.y=range(ypos)
	xyscal=abs(diff(default.y)/diff(default.x))
	plotM[,13][plotM[,13]!="2"]="black"
	plotM[,13][plotM[,13]=="2"]=paste("gray",gtestgray,sep="")
	polyT=function(x,y,ployS=boxsize,scal=c(Lscal,xyscal)){X=rbind(x+rep(ployS*scal[1]*c(-1,1),each=2),y+rep(ployS*prod(scal)*(pdfw/pdfh)*c(-1,1),2));return(X[,c(1,2,4,3,1)])}
	polysum=function(x,y,ci,ployH=0.3){X=rbind(c(x,ci[1],x,ci[2]),y+c(ployH,0,-ployH,0));return(X[,c(1,2,3,4,1)])}
	Llow=ypos[rev(which(plotM[,1]=="Heterogeneity:"))[1]+1-3*overall]
	rulerM=function(ra,xLP,y=Llow){
		xpp=xLP+Lscal*ra;xpps=seq(xpp[1],xpp[2],by=0.1*Lscal)
		lines(xpp,rep(y,2))
		apply(t(xpps),2,function(x)lines(rep(x,2),c(y,y-0.5)))
		apply(rbind(xpps,seq(ra[1],ra[2],by=0.1)),2,function(x)text(x = x[1], y = y-1, label = x[2],adj=c(0.5,1)))
		}
	pvalS=function(x){if(x<0.01){return(" < 0.01")}else{return(paste("=",round(x,2)))}}
	
	if(class(savepdf)!="logical")pdf(file=savepdf, width =pdfw, height=pdfh)
	##############
	par(mar = c(2, 1, 0.5, 0.5) + 0.1)
	plot(x = NULL, y = NULL, type = "n", xlim = default.x, 
        ylim = default.y, xlab = "", ylab = "", 
        xaxs = "i", yaxs = "i", axes = FALSE, frame.plot = FALSE)
	if(overall){
		lines(rep(as.num(sens.N.propCI[1]),2)*Lscal+xps[6],c(ypos[3],Llow),lty=3)
		lines(rep(as.num(spec.N.propCI[1]),2)*Lscal+xps[9],c(ypos[3],Llow),lty=3)
		}
	for(i in 1:dim(plotM)[1]){
		roltemp=plotM[i,]
		if(roltemp[1]=="Heterogeneity:"){
			if(roltemp[2]=="T" & overall){
				I2=paste(round(c(sens$I2,spec$I2)*100,0),"%",sep="")
				tau2=round(c(sens$tau2,spec$tau2),4)
				pvals=c(pvalS(sens$pval.Q),pvalS(spec$pval.Q))
				sens.t=substitute(paste(italic(l)^2,"=",A,", ",italic(tau)^2,"=",B,", ",italic(p),C) , list(A=I2[1], B=tau2[1], C=pvals[1]))
				spec.t=substitute(paste(italic(l)^2,"=",A,", ",italic(tau)^2,"=",B,", ",italic(p),C) , list(A=I2[2], B=tau2[2], C=pvals[2]))
				text(x = xps[1], y = ypos[i], label = roltemp[1],adj=c(0,1),col=roltemp[13])
				text(x = xps[9], y = ypos[i], label = sens.t,adj=c(1,1),col=roltemp[13])
				text(x = xps[12], y = ypos[i], label = spec.t,adj=c(1,1),col=roltemp[13])
				}else{
				gint=which(paste(GL)==roltemp[2])
				I2=paste(round(c(sens$I2.w[gint],spec$I2.w[gint])*100,0),"%",sep="")
				tau2=round(c(sens$tau2.w[gint],spec$tau2.w[gint]),4)
				pvals=c(pvalS(sens$pval.Q.w[gint]),pvalS(spec$pval.Q.w[gint]))
				sens.t=substitute(paste(italic(l)^2,"=",A,", ",italic(tau)^2,"=",B,", ",italic(p),C) , list(A=I2[1], B=tau2[1], C=pvals[1]))
				spec.t=substitute(paste(italic(l)^2,"=",A,", ",italic(tau)^2,"=",B,", ",italic(p),C) , list(A=I2[2], B=tau2[2], C=pvals[2]))
				text(x = xps[1], y = ypos[i], label = roltemp[1],adj=c(0,1),col=roltemp[13])
				text(x = xps[9], y = ypos[i], label = sens.t,adj=c(1,1),col=roltemp[13])
				text(x = xps[12], y = ypos[i], label = spec.t,adj=c(1,1),col=roltemp[13])
				}
			next}
		if(roltemp[1]=="Test for subgroup differences:"){
			chi=round(c(sens$Q.b.random,spec$Q.b.random),2)
			dfv=c(sens$df.Q.b,spec$df.Q.b)
			pvals=c(pvalS(sens$pval.Q.b.random),pvalS(spec$pval.Q.b.random))
			sens.t=substitute(paste(italic(chi)[B]^2,"=",A,", ",df,"=",B,", (",italic(p),C,")") , list(A=chi[1], B=dfv[1], C=pvals[1]))
			spec.t=substitute(paste(italic(chi)[B]^2,"=",A,", ",df,"=",B,", (",italic(p),C,")") , list(A=chi[2], B=dfv[2], C=pvals[2]))
			text(x = xpos[1], y = ypos[i], label = roltemp[1],adj=c(0,1),col=roltemp[13])
			text(x = xps[9], y = ypos[i], label = sens.t,adj=c(1,1),col=roltemp[13])
			text(x = xps[12], y = ypos[i], label = spec.t,adj=c(1,1),col=roltemp[13])
			next}
		for(j in 1:11){
			if(  i==1 | !(j %in% c(6,9))  ){
				if( !(j %in% 2:5)*noRawdata  ){
					text(x = xpos[j], y = ypos[i], label = roltemp[j],
						adj=c(adjX[j],1),	
						font=as.numeric(roltemp[12]),
						col=roltemp[13])
					}
				}else if(roltemp[j]=="1"){
				as.num(roltemp[j+1])*Lscal+xps[j]->xp
				XY=polyT(xp,ypos[i]-0.5)
				polygon(XY[1,],XY[2,],col="gray",border=0)
				segments(xp,ypos[i]-0.4, xp, ypos[i]-0.6,lwd=1.2)
				as.num(roltemp[j+2])*Lscal+xps[j]->xp
				segments(xp[1],ypos[i]-0.5, xp[2], ypos[i]-0.5,lwd=1.1)
				}else if(roltemp[j]=="2"){
				as.num(roltemp[j+1])*Lscal+xps[j]->xp
				as.num(roltemp[j+2])*Lscal+xps[j]->cip
				XY=polysum(xp,ypos[i]-0.5,cip)
				polygon(XY[1,],XY[2,],col="gray")
				}
			}
		}
	rulerM(round(c(sensCI.L,sensCI.U),1),xps[6])
	rulerM(round(c(specCI.L,specCI.U),1),xps[9])
	#################
	if(class(savepdf)!="logical")dev.off()
	}

