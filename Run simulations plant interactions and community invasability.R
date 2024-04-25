#########################################################################################################
###########	simulation code for plant interactions invasibility model
#########################################################################################################

### Copyright 2022 Elizabeth (Libby) J. Trevenen & Dr Michael Renton. Adapted from model described in Teste et al 2017.
### If you wish to modify or share this code, you must obtain permission from the original authors. 


############################################################################################################################
############################### parameters
############################################################################################################################

nruns= 1                  	  # the number of runs 
nyrs = 20              	        # the number of years  
worldsize=100                   # how big the world is
nsp=100                   	  # number of species  before invasion 
minsize = 0.01                  # the minimum size of a plant 
pmortnormal= 0.1                # establishment mortality rate 
pimmigration = 0.001            # immigration from species not in the world
recruitment='global'            # can be recrutited from anywhere (global) or from closest living plant (local) 
effect.type='previous'          # can be previous or neighbour 
comp.type='notsimple'           # simple = no competition, not simple (or anything else) = there is competition between neighbours
seedbankdeclinerate= 1          # the percentage of seedbank decline per year 1 = 100% (so no seedbank) 
initial.empty.prop= 1           # how much empty space will there be

percertage.invaded = 0.01       # the percentage of the landscape to be populated with an invader 
time.to.invasion = 10        # number of runs before invasion occurs 

##############################
## initial soil conditions  
##############################

init.soiltype.unconditioned = TRUE   ## set new soil to be unconditioned 
#unconditioned.soil="min"  
#unconditioned.soil="max"
unconditioned.soil="mean"            ## growth of species in unconditioned soil is set to mean (standard growth)
#unconditioned.soil="sterile" 

# Number of species at the start of the simulation 
no.sp=0  # if = 0 then it is an empty world to start with and plants immigrate in naturally
all.sp=1:nsp
initialspeciespool = no.sp; initialspeciespooln = "no.sp"

##############################
###### Plant interaction matricies 
##############################

source('create PSFI matricies.R')


###############################
ani = FALSE		# TRUE is to make a plot 		

### set up neighbourhood
pos=expand.grid(x=1:worldsize,y=1:worldsize)
neighbourlist = matrix(0,nrow=worldsize^2,ncol=4)
rr = pos$x
cc = pos$y
rru = rr-1 ### uppper row. row 0 -1
rru[rru==0]=worldsize ## if i am in the top row the next row up will be on the bottom. the world wraps around 
ccu = cc-1
ccu[ccu==0]=worldsize 
rrd = rr+1
rrd[rrd==worldsize+1]=1 
ccd = cc+1
ccd[ccd==worldsize+1]=1
#cbind(rr,cc,rru,ccu,rrd,ccd)
neighbourlist[,1] = rru+(cc-1)* worldsize
neighbourlist[,2] = rrd+(cc-1)* worldsize
neighbourlist[,3] = rr+(ccu-1)* worldsize
neighbourlist[,4] = rr+(ccd-1)* worldsize


########## define functions
rgr.func = function(i){
	if (planttype[i]==0) return(0)
	if (soiltype[i]==0) {
		if (unconditioned.soil=="min") return(min(rgrm))
		if (unconditioned.soil=="max") return(max(rgrm))
		if (unconditioned.soil=="mean") return(rowMeans(rgrm)[planttype[i]])
		if (unconditioned.soil=="sterile") return(no.soil.biota[planttype[i],1])
	}
	if (effect.type=='self') return(rgrm[planttype[i],planttype[i]])
	if (effect.type=='previous') return(rgrm[planttype[i],soiltype[i]]) 
	if (effect.type=='neighbour') {
		grs = rgrm[planttype[i],planttype[neighbourlist[i,1:4]]]
		grs = c(grs, rep(exp(4), 4-length(grs)))
		return(prod(grs)^(1/4)) # effect of neighbour 
	}
}

## plant grows until it reaches max size 
## there is a curve in here used to relfect a more realistic growh curve, however when tested it was made no difference. 
constrainfunc = function(plantsize ,potgrowth, maxsize){ 
	if (comp.type=='simple') {					
		plantsize = plantsize+potgrowth
		plantsize[plantsize>maxsize]=maxsize
		return (plantsize)
	} else {
		psc = plantsize 
		psc[1:(worldsize-1),] = psc[1:(worldsize-1),] + plantsize[2:worldsize,] 
		psc[2:(worldsize),] = psc[2:(worldsize),] + plantsize[1:(worldsize-1),] 
		psc[,1:(worldsize-1)] = psc[,1:(worldsize-1)] + plantsize[,2:worldsize] 
		psc[,2:(worldsize)] = psc[,2:(worldsize)] + plantsize[,1:(worldsize-1)]  
		maxsz = matrix(maxsize*5, nrow=worldsize,ncol=worldsize)
		maxsz[1,]=maxsz[1,]-maxsize
		maxsz[worldsize,]=maxsz[worldsize,]-maxsize
		maxsz[,1]=maxsz[,1]-maxsize
		maxsz[,worldsize]=maxsz[,worldsize]-maxsize
		actgrowth = potgrowth*((maxsz-psc)/maxsz)
		actgrowth[actgrowth>potgrowth]=potgrowth[actgrowth>potgrowth]
		actgrowth[actgrowth<0]=0	
		return(plantsize+actgrowth)
	}
}

#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################

### loop options flags
iiii=0  # counter set to 0 
ptm <- proc.time()	

for (hypn in 1:39){ 

iiii=iiii+1
print(paste('hyp:',hypn))

rgrm = hyplist[[hypn]]
head(rgrm)
rgrm = exp(rgrm) 
head(rgrm)


### set up overall records and start runs loop 
allrichrec=NULL 
allshrec=NULL 
allsimprec=NULL
allinvsimprec=NULL
allallrec=list()
alltotalbiomassrec = NULL

############################################# 
############################################# 
# beginning of loop 
############################################# 
############################################# 

for (run in 1:nruns){
print(paste('run',run,'elapsed time:'))
#print(proc.time()-ptm) 


############################
### pre-invasion conditions 
############################

planttype=matrix(sample(initialspeciespool,worldsize^2,replace=TRUE),nrow=worldsize)
planttype[runif(worldsize^2)<initial.empty.prop] = 0 
plantsize=matrix(minsize ,nrow=worldsize,ncol=worldsize) 
plantsize[planttype==0]=0 # adding in spaces 
if (init.soiltype.unconditioned) {soiltype=matrix(rep(0,worldsize^2,replace=TRUE),nrow=worldsize)} else {soiltype=matrix(sample(1:nsp,worldsize^2,replace=TRUE),nrow=worldsize)} # soil type is all unconditioned or  soil type is random across the world 
soilstrength=matrix(0.1,nrow=worldsize,ncol=worldsize)
pos=expand.grid(x=1:worldsize,y=1:worldsize)


## set up recorders for each simulation
richrec = NULL
shrec = NULL
simprec = NULL
invsimprec=NULL
biomassbyspeciesrec = NULL
totalbiomassrec = NULL

if (recruitment=='global'){
	seedbank=rep(0,nsp+1)
} else {
	seedbank = array(0,dim=c(worldsize*worldsize,nsp+1))
}

############################################# 
############################################# 
# beginning of simulation
############################################# 
#############################################

for (t in 1:nyrs){ 				
	#print(t)   

	

	########################################################################
	## add in invader 
	########################################################################

	if(t==time.to.invasion){

	# find empty cells and place the invader into a percentage of them 
	planttype2=sample(which(planttype==0),I(percertage.invaded*worldsize^2))
	planttype[planttype2]=101
	plantsize[planttype2]=minsize
	}


	########################################################################
	## immigration - invader cannot re-immigrate in 
	########################################################################
	spaces=which(planttype==0)
	plantsize[spaces]=minsize  
	ifimmigrant=rbinom(length(spaces), 1,pimmigration)
	immigrant = spaces[ifimmigrant==1] 
	planttype[immigrant]= sample(1:nsp,length(immigrant),replace=TRUE)

	########################################################################
	## recruitment and seedbank dynamics - invader can recruit 
	########################################################################
	spaces=which(planttype==0)

	if (t>1){


	if (recruitment=='global' & seedbankdeclinerate==1) {    ## ie no seedbank, and global recruitment 
		newrecruits = unlist(sapply(spaces, function(i){
			poss=neighbourlist[i,]
			if (all(planttype[poss]==0)) return(0)
			return (sample(planttype[poss],1,prob=plantsize[poss])) 
		}))
		planttype[spaces]=newrecruits


	} else {    ## local AND persistent seedbank
		#print ('yes')
		seedbank = seedbank * (1-seedbankdeclinerate)
		for (e in 1:worldsize^2) {		## replenish seedbank locally
			poss=neighbourlist[e,]
			if (any(planttype[poss]!=0)) 	{
				pts = planttype[poss]
				pss = plantsize[poss][pts>0]
				seedbank[e,planttype[poss][pts>0]] = seedbank[e,planttype[poss][pts>0]] + pss
			}
		}
		newrecruits = unlist(sapply(spaces, function(i){
			if (sum(seedbank[i,])==0) return (0) else return (sample(0:nsp,1,prob=c(minsize,seedbank[i,]))) 
		}))
		planttype[spaces]=newrecruits

	}
	plantsize[spaces]=minsize 
	plantsize[planttype==0]=0
	}

	########################################################################
	## mortality
	########################################################################
	pmort = pmortnormal
	chance.death = unlist(sapply(1:worldsize^2 , function(i) pmort ))	
	chance.death [plantsize>minsize ] = pmort 
	chance.death [planttype == 0 ] = 0
	chance.death[is.na(chance.death)]= pmort 
	ifdie = chance.death>runif(worldsize^2)
	soiltype[ifdie] = planttype[ifdie]
	planttype[ifdie]=0
	plantsize[planttype==0]=0

	########################################################################
	## growth
	########################################################################
	rgr = unlist(sapply(1:worldsize^2 , rgr.func ))	
	rgr[is.na(rgr)]=0
	potgrowth = plantsize*rgr - plantsize
	plantsize = constrainfunc(plantsize ,potgrowth , 100)

	## plot and record
	if (ani) plot(pos$x,pos$y,col=planttype,pch=(planttype>8)*2+15,xlab=t,ylab="",xaxt='n',yaxt='n',cex=log(plantsize)/4)
	richrec = c(richrec,length(unique(as.vector(planttype)))-1)
	biomassbyspecies = tapply(plantsize ,planttype, sum)
	biomassbyspeciesrec = rbind(biomassbyspeciesrec , data.frame(biomass=biomassbyspecies,species=names(biomassbyspecies),year=t))
	totalbiomassrec = c(totalbiomassrec , sum(biomassbyspecies ))  
	## inverse simpsions index
	invsimp_diversity = function(x){
	ps=x/sum(x)
	1/sum(ps^2)
	}
	invsimprec= c(invsimprec, invsimp_diversity(biomassbyspecies)) 


} ## end of each year
	print(run)

	write.csv(biomassbyspeciesrec , paste('biomassbyspecies',run,hypnames[hypn],effect.type,seedbankdeclinerate,recruitment,pmortnormal,pimmigration,comp.type,'.csv',sep=''))
	allrichrec=rbind(allrichrec,richrec)
	alltotalbiomassrec = rbind(alltotalbiomassrec ,totalbiomassrec )
	allinvsimprec=rbind(allinvsimprec, invsimprec)

}   ## end of each run

	write.csv(allrichrec , paste('allrichrec2',hypnames[hypn],effect.type,seedbankdeclinerate,recruitment,pmortnormal,pimmigration,comp.type,'.csv',sep=''))
	write.csv(alltotalbiomassrec , paste('alltotalbiomassrec2',hypnames[hypn],effect.type,seedbankdeclinerate,recruitment,pmortnormal,pimmigration,comp.type,'.csv',sep=''))
	write.csv(allinvsimprec, paste('allinvsimprec2',hypnames[hypn],effect.type,seedbankdeclinerate,recruitment,pmortnormal,pimmigration,comp.type,'.csv',sep=''))

}  ## end of all loops
	biomassbyspeciesrec = NULL
	totalbiomassrec = NULL

#########################################################################################################
#########################################################################################################
#########################################################################################################


### Copyright 2022 Elizabeth (Libby) J. Trevenen & Dr Michael Renton. Adapted from model described in Teste et al 2017.
### If you wish to modify or share this code, you must obtain permission from the original authors. 


