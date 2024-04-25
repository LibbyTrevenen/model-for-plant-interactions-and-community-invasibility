### Copyright 2019 Elizabeth (Libby) J. Trevenen & Dr Michael Renton. Adapted from model described in Teste et al 2017.
### If you wish to modify or share this code, you must obtain permission from the original authors. 


########################################################################################
###########		making the matricies
########################################################################################

# growth rates

p = 0.6           # fast growth -average increase in growth rate due to plant soil interaction
n = 0.4		# slow growth -average decrease in growth rate due to plant soil interaction 
m = 0.5		# standard growth


########################################################################################
### list of matricies

########################################################################################
########################################################################################
# Null scenarios 
########################################################################################
########################################################################################

#######################################################################################
##   Null 
#######  
Null = matrix(m,nrow=101,ncol=101)
#windows()
###image(Null[seq(101,1,by=-1),])
###image(Null,zlim=c(n,p))
diag(Null)[c(101,101)] = p
Null

########################################################################################
########################################################################################
# conspecifics
########################################################################################
########################################################################################

##############################
##  negative conspecific 
####### 

Null.NC = matrix(m,nrow=101,ncol=101)
diag(Null.NC)=n
diag(Null.NC)[c(101,101)] = p
###image(Null.NC[seq(101,1,by=-1),])
###image(Null.NC,zlim=c(n,p))

#############################
## 3 Null.PC
####### positive conspecific 

Null.PC = matrix(m, nrow=101,ncol=101)
diag(Null.PC)=p
diag(Null.PC)[c(101,101)] = p
###image(Null.PC[seq(101,1,by=-1),])
###image(Null.PC,zlim=c(n,p))
Null.PC

########################################################################################
########################################################################################
### Ring PSI NC
########################################################################################
########################################################################################

########################################################################################
##### positive ring systems 

##############################
##  Pos.Ring.10.NC
####### general ring help 10 species on etiher side but dosent help itself. 
neff = 10
ring.10 = matrix(m,nrow=101,ncol=101)
diag(ring.10)=n
for (i in 1:101) {
	ring.10[i,seq(i,i+neff-1)%%101+1] = p
}
###image(ring.10[seq(101,1,by=-1),])
###image(ring.10,zlim=c(n,p))
ring.10[101,1:25] = m 
ring.10[1:101,101] = m 
diag(ring.10)[c(101,101)] = p
Pos.Ring.10.NC=ring.10



##############################
##  Pos.Ring.5.NC
####### general ring help 5 species on etiher side but dosent help itself. 
neff = 5
ring.5 = matrix(m,nrow=101,ncol=101)
diag(ring.5)=n
for (i in 1:101) {
	ring.5[i,seq(i,i+neff-1)%%101+1] = p
}
###image(ring.5[seq(101,1,by=-1),])
###image(ring.5,zlim=c(n,p))
ring.5[1:101,101] = m 
ring.5[101,1:25] = m 
diag(ring.5)[c(101,101)] = p
Pos.Ring.5.NC=ring.5


##############################
##  Pos.Ring.20.NC
####### general ring help 20 species on etiher side but dosent help itself. 
neff = 20
ring.20 = matrix(m,nrow=101,ncol=101)
diag(ring.20)=n
for (i in 1:101) {
	ring.20[i,seq(i,i+neff-1)%%101+1] = p
}
###image(ring.20,zlim=c(n,p))
ring.20[101,1:25] = m  
ring.20[1:101,101] = m 
diag(ring.20)[c(101,101)] = p
Pos.Ring.20.NC=ring.20


########################################################################################
##### negative ring systems 

##############################
##  Neg.Ring.10
####### general ring hinder 10 species on etiher side but dosent help itself. 
neff = 10
Neg.Ring.10 = matrix(m,nrow=101,ncol=101)
diag(Neg.Ring.10)=n
for (i in 1:101) {
	Neg.Ring.10[i,seq(i,i+neff-1)%%101+1] = n
}
###image(Neg.Ring.10[seq(101,1,by=-1),])
###image(Neg.Ring.10,zlim=c(n,p))

Neg.Ring.10[1:101,101] = m
Neg.Ring.10[101,1:25] = m 
diag(Neg.Ring.10)[c(101,101)] = p
Neg.Ring.10.NC=Neg.Ring.10

##############################
##  Neg.Ring.5
####### general ring hinder 10 species on etiher side but dosent help itself. 
neff = 5
Neg.Ring.5 = matrix(m,nrow=101,ncol=101)
diag(Neg.Ring.5)=n
for (i in 1:101) {
	Neg.Ring.5[i,seq(i,i+neff-1)%%101+1] = n
}
###image(Neg.Ring.5[seq(101,1,by=-1),])
###image(Neg.Ring.5,zlim=c(n,p))
Neg.Ring.5[101,1:25] = m 
Neg.Ring.5[1:101,101] = m  
diag(Neg.Ring.5)[c(101,101)] = p

Neg.Ring.5.NC=Neg.Ring.5

##############################
##  Neg.Ring.20.NC
####### general ring hinder 20 species on etiher side but dosent help itself. 
neff = 20
Neg.Ring.20 = matrix(m,nrow=101,ncol=101)
diag(Neg.Ring.20)=n
for (i in 1:101) {
	Neg.Ring.20[i,seq(i,i+neff-1)%%101+1] = n
}
###image(Neg.Ring.20[seq(101,1,by=-1),])
###image(Neg.Ring.20,zlim=c(n,p))
Neg.Ring.20[1:101,101] = m
Neg.Ring.20[101,1:25] = m 
diag(Neg.Ring.20)[c(101,101)] = p
Neg.Ring.20.NC=Neg.Ring.20

########################################################################################
########################################################################################
### Modular PSI NC
########################################################################################
########################################################################################

########################################################################################
##### modular systems with positive groups 

##############################
##  Pos.Mod.10.NC
####### 10 groups of 10 species. 
Pos.Mod.10 = matrix(m,nrow=101,ncol=101)
for (gr in 1:10) {
	Pos.Mod.10[((gr-1)*10+1):(gr*10),  ((gr-1)*10+1):(gr*10)] = p
}
diag(Pos.Mod.10)=n
diag(Pos.Mod.10)[c(101,101)] = p
###image(Pos.Mod.10[seq(101,1,by=-1),][80:101,1:20])
#windows()
###image(Pos.Mod.10,zlim=c(n,p))
Pos.Mod.10.NC=Pos.Mod.10

##############################
##  Pos.Mod.5.NC
####### 20 groups of 5 species
Pos.Mod.5 = matrix(m,nrow=101,ncol=101)
for (gr in 1:20) {
	Pos.Mod.5[((gr-1)*5+1):(gr*5),  ((gr-1)*5+1):(gr*5)] = p
}
diag(Pos.Mod.5)=n
diag(Pos.Mod.5)[c(101,101)] = p
###image(Pos.Mod.5[seq(101,1,by=-1),])
###image(Pos.Mod.5,zlim=c(n,p))
Pos.Mod.5.NC=Pos.Mod.5


##############################
##  Pos.Mod.20.NC
####### 5 groups of 20 species
Pos.Mod.20 = matrix(m,nrow=101,ncol=101)
for (gr in 1:5) {
	Pos.Mod.20[((gr-1)*20+1):(gr*20),  ((gr-1)*20+1):(gr*20)] = p
}
diag(Pos.Mod.20)=n
diag(Pos.Mod.20)[c(101,101)] = p
###image(Pos.Mod.20[seq(101,1,by=-1),])
###image(Pos.Mod.20,zlim=c(n,p))
Pos.Mod.20.NC=Pos.Mod.20


########################################################################################
##### negative modular group systems 

##############################
##  Neg.Mod.10.NC
####### 10 groups of 10 species. Species dont help own species and negatively affecting 10 species "around them" 
Neg.Mod.10 = matrix(m,nrow=101,ncol=101)
for (gr in 1:10) {
	Neg.Mod.10[((gr-1)*10+1):(gr*10),  ((gr-1)*10+1):(gr*10)] = n
}
diag(Neg.Mod.10)=n
diag(Neg.Mod.10)[c(101,101)] = p
###image(Neg.Mod.10[seq(101,1,by=-1),])
#windows()
###image(Neg.Mod.10,zlim=c(n,p))
Neg.Mod.10.NC=Neg.Mod.10

##############################
##  Neg.Mod.5.NC
#######20 groups of 5 species. Species dont help own species and negatively affecting 20 species "around them" 
Neg.Mod.5 = matrix(m,nrow=101,ncol=101)
for (gr in 1:20) {
	Neg.Mod.5[((gr-1)*5+1):(gr*5),  ((gr-1)*5+1):(gr*5)] = n
}
diag(Neg.Mod.5)=n
diag(Neg.Mod.5)[c(101,101)] = p
###image(Neg.Mod.5[seq(101,1,by=-1),])
###image(Neg.Mod.5,zlim=c(n,p))
Neg.Mod.5.NC=Neg.Mod.5

##############################
##  Neg.Mod.20.NC
####### 5 groups of 20 species. Species dont help own species and negatively affecting 20 species "around them" 
Neg.Mod.20 = matrix(m,nrow=101,ncol=101)
for (gr in 1:5) {
	Neg.Mod.20[((gr-1)*20+1):(gr*20),  ((gr-1)*20+1):(gr*20)] = n
}
diag(Neg.Mod.20)=n
diag(Neg.Mod.20)[c(101,101)] = p
###image(Neg.Mod.20[seq(101,1,by=-1),])
###image(Neg.Mod.20,zlim=c(n,p))
Neg.Mod.20.NC=Neg.Mod.20

########################################################################################
########################################################################################
### Nested PSI NC
########################################################################################
########################################################################################

########################################################################################
##### positive nested  PSIs 

##############################
##  Pos.Nest.10.NC
####### every species grows well in 10 species soil. 
a= matrix(p, nrow=101, ncol=11)
b= matrix(m, nrow=101, ncol=90)
nested.10.NC= cbind(b,a)
diag(nested.10.NC)=n
diag(nested.10.NC)[c(101,101)] = p
nested.10.NC[101,90:100] = m
nested.10.NC[1:100,101] = m
###image(nested.10.NC,zlim=c(n,p))
Pos.Nest.10.NC=nested.10.NC

##############################
##  Pos.Nest.20.NC
####### every species grows well in 20 species soil NC. 
a= matrix(p, nrow=101, ncol=21)
b= matrix(m, nrow=101, ncol=80)
nested.20.NC= cbind(b,a)
diag(nested.20.NC)=n 
diag(nested.20.NC)[c(101,101)] = p
nested.20.NC[101,80:100] = m
nested.20.NC[1:100,101] = m
###image(nested.20.NC,zlim=c(n,p))
Pos.Nest.20.NC=nested.20.NC


##############################
##  Pos.Nest.5.NC
####### every species grows well in 5 species soil NC. 
a= matrix(p, nrow=101, ncol=6)
b= matrix(m, nrow=101, ncol=95)
nested.5.NC= cbind(b,a)
diag(nested.5.NC)=n
diag(nested.5.NC)[c(101,101)] = p
nested.5.NC[101,80:100] = m
nested.5.NC[1:100,101] = m
###image(nested.5.NC,zlim=c(n,p))
Pos.Nest.5.NC=nested.5.NC


########################################################################################
##### Negative Nested  PSIs

##############################
##  Neg.Nest.10.NC
####### every species grows well in 10 species soil. 
a= matrix(n, nrow=101, ncol=11)
b= matrix(m, nrow=101, ncol=90)
neg.nested.10.NC= cbind(b,a)
diag(neg.nested.10.NC)=n
diag(neg.nested.10.NC)[c(101,101)] = p
neg.nested.10.NC[101,90:100] = m
neg.nested.10.NC[1:100,101] = m
##image(neg.nested.10.NC,zlim=c(n,p))
Neg.Nest.10.NC=neg.nested.10.NC

##############################
##  Neg.Nest.20.NC
####### every species grows well in 20 species soil NC. 
a= matrix(n, nrow=101, ncol=21)
b= matrix(m, nrow=101, ncol=80)
neg.nested.20.NC= cbind(b,a)
diag(neg.nested.20.NC)=n 
diag(neg.nested.20.NC)[c(101,101)] = p
neg.nested.20.NC[101,80:100] = m
neg.nested.20.NC[1:100,101] = m
###image(neg.nested.20.NC,zlim=c(n,p))
Neg.Nest.20.NC=neg.nested.20.NC


##############################
##  Neg.Nest.5.NC
####### every species grows well in 5 species soil NC. 
a= matrix(n, nrow=101, ncol=6)
b= matrix(m, nrow=101, ncol=95)
neg.nested.5.NC= cbind(b,a)
diag(neg.nested.5.NC)=n 
diag(neg.nested.5.NC)[c(101,101)] = p
neg.nested.5.NC[101,80:100] = m
neg.nested.5.NC[1:100,101] = m
##image(neg.nested.5.NC,zlim=c(n,p))
Neg.Nest.5.NC=neg.nested.5.NC


########################################################################################
########################################################################################
### Ring PSI MC
########################################################################################
########################################################################################

########################################################################################
##### positive ring systems 

##############################
##  Pos.Ring.10
####### general ring help 10 species but dosent help itself. 
neff = 10
ring.10 = matrix(m,nrow=101,ncol=101)
diag(ring.10)=m


for (i in 1:101) {
	ring.10[i,seq(i,i+neff-1)%%101+1] = p
}
ring.10[1:101,101] = m
ring.10[101,1:25] = m   
diag(ring.10)[c(101,101)] = p
###image(ring.10,zlim=c(n,p))
Pos.Ring.10=ring.10

##############################
##  Pos.Ring.5
####### general ring help 5 species but dosent help itself. 
neff = 5
ring.5 = matrix(m,nrow=101,ncol=101)
diag(ring.5)=m
for (i in 1:101) {
	ring.5[i,seq(i,i+neff-1)%%101+1] = p
}

ring.5[1:101,101] = m
ring.5[101,1:25] = m   
diag(ring.5)[c(101,101)] = p
###image(ring.5,zlim=c(n,p))
Pos.Ring.5=ring.5


##############################
##  Pos.Ring.20
####### general ring help 20 species but dosent help itself. 
neff = 20
ring.20 = matrix(m,nrow=101,ncol=101)
diag(ring.20)=m
for (i in 1:101) {
	ring.20[i,seq(i,i+neff-1)%%101+1] = p
}
ring.20[1:101,101] = m
ring.20[101,1:25] = m 
diag(ring.20)[c(101,101)] = p
  
###image(ring.20,zlim=c(n,p))
Pos.Ring.20=ring.20

########################################################################################
##### negative ring systems 

##############################
##  Neg.Ring.10
####### general ring hinders 10 species but dosent help itself. 
neff = 10
neg.ring.10 = matrix(m,nrow=101,ncol=101)
diag(neg.ring.10)=m
for (i in 1:101) {
	neg.ring.10[i,seq(i,i+neff-1)%%101+1] = n
}
neg.ring.10[1:101,101] = m
neg.ring.10[101,1:25] = m 
diag(neg.ring.10)[c(101,101)] = p

###image(neg.ring.10,zlim=c(n,p))
Neg.Ring.10=neg.ring.10

##############################
##  Neg.Ring.5
####### general ring hinders 5 species but dosent help itself. 
neff = 5
neg.ring.5 = matrix(m,nrow=101,ncol=101)
diag(neg.ring.5)=m
for (i in 1:101) {
	neg.ring.5[i,seq(i,i+neff-1)%%101+1] = n
}
neg.ring.5[1:101,101] = m
neg.ring.5[101,1:25] = m 
diag(neg.ring.5)[c(101,101)] = p

###image(neg.ring.5,zlim=c(n,p))
Neg.Ring.5=neg.ring.5

##############################
##  Neg.Ring.20
####### general ring hinders 20 species but dosent help itself. 
neff = 20
neg.ring.20 = matrix(m,nrow=101,ncol=101)
diag(neg.ring.20)=m
for (i in 1:101) {
	neg.ring.20[i,seq(i,i+neff-1)%%101+1] = n
}
neg.ring.20[1:101,101] = m
neg.ring.20[101,1:25] = m 
diag(neg.ring.20)[c(101,101)] = p

###image(neg.ring.20,zlim=c(n,p))
Neg.Ring.20=neg.ring.20

########################################################################################
########################################################################################
### Modular PSI MC
########################################################################################
########################################################################################

########################################################################################
##### modular systems with positive groups 

##############################
##  Pos.Mod.5
####### 20 groups of 5 species
friends.5 = matrix(m,nrow=101,ncol=101)
for (gr in 1:20) {
	friends.5[((gr-1)*5+1):(gr*5),  ((gr-1)*5+1):(gr*5)] = p
}
diag(friends.5)=m
diag(friends.5)[c(101,101)] = p
###image(friends.5,zlim=c(n,p))
Pos.Mod.5=friends.5

##############################
##  Pos.Mod.10
####### 10 groups of 10 species. Species dont help own species but are positively helping 10 species 'around them' 
friends.10 = matrix(m,nrow=101,ncol=101)
for (gr in 1:10) {
	friends.10[((gr-1)*10+1):(gr*10),  ((gr-1)*10+1):(gr*10)] = p
}
diag(friends.10)=m
diag(friends.10)[c(101,101)] = p
###image(friends.10,zlim=c(n,p))
Pos.Mod.10=friends.10

##############################
##  Pos.Mod.20
####### 5 groups of 20 species
friends.20 = matrix(m,nrow=101,ncol=101)
for (gr in 1:5) {
	friends.20[((gr-1)*20+1):(gr*20),  ((gr-1)*20+1):(gr*20)] = p
}
diag(friends.20)=m
diag(friends.20)[c(101,101)] = p
###image(friends.20,zlim=c(n,p))
Pos.Mod.20=friends.20


########################################################################################
##### negative modular group systems 

##############################
##  Neg.Mod.10
####### 10 groups of 10 species. 
enemies.10 = matrix(m,nrow=101,ncol=101)
for (gr in 1:10) {
	enemies.10[((gr-1)*10+1):(gr*10),  ((gr-1)*10+1):(gr*10)] = n
}
diag(enemies.10)=m
diag(enemies.10)[c(101,101)] = p
###image(enemies.10,zlim=c(n,p))
Neg.Mod.10=enemies.10

##############################
##  Neg.Mod.5
#######20 groups of 5 species. 
enemies.5 = matrix(m,nrow=101,ncol=101)
for (gr in 1:20) {
	enemies.5[((gr-1)*5+1):(gr*5),  ((gr-1)*5+1):(gr*5)] = n
}
diag(enemies.5)=m
diag(enemies.5)[c(101,101)] = p
###image(enemies.5,zlim=c(n,p))
Neg.Mod.5=enemies.5

##############################
##  Neg.Mod.20
####### 5 groups of 20 species. 
enemies.20 = matrix(m,nrow=101,ncol=101)
for (gr in 1:5) {
	enemies.20[((gr-1)*20+1):(gr*20),  ((gr-1)*20+1):(gr*20)] = n
}
diag(enemies.20)=m
diag(enemies.20)[c(101,101)] = p
###image(enemies.20,zlim=c(n,p))
Neg.Mod.20=enemies.20


########################################################################################
########################################################################################
### Nested PSI MC
########################################################################################
########################################################################################

########################################################################################
##### positive Pos.Nest  PSIs 

##############################
##  Pos.Nest.10
####### every species grows well under 10 species soil. 
a= matrix(p, nrow=101, ncol=11)
b= matrix(m, nrow=101, ncol=90)
nested.10.NC= cbind(b,a)
diag(nested.10.NC)=m
diag(nested.10.NC)[c(101,101)] = p
nested.10.NC[101,90:100] = m
nested.10.NC[1:100,101] = m
###image(nested.10.NC,zlim=c(n,p))
Pos.Nest.10=nested.10.NC



##############################
##  Pos.Nest.20
####### every species grows well under 20 species soil MC. 
a= matrix(p, nrow=101, ncol=21)
b= matrix(m, nrow=101, ncol=80)
nested.20.NC= cbind(b,a)
diag(nested.20.NC)=m 
diag(nested.20.NC)[c(101,101)] = p
nested.20.NC[101,80:100] = m
nested.20.NC[1:100,101] = m
###image(nested.20.NC,zlim=c(n,p))
Pos.Nest.20=nested.20.NC


##############################
##  Pos.Nest.5
####### every species grows well under 5 species soil MC. 
a= matrix(p, nrow=101, ncol=6)
b= matrix(m, nrow=101, ncol=95)
nested.5.NC= cbind(b,a)
diag(nested.5.NC)=m
diag(nested.5.NC)[c(101,101)] = p
nested.5.NC[101,80:100] = m
nested.5.NC[1:100,101] = m
###image(nested.5.NC,zlim=c(n,p))
Pos.Nest.5=nested.5.NC


########################################################################################
##### Negative Nested  PSIs

##############################
##  Neg.Nest.10
####### every species grows well under 10 species soil. 
a= matrix(n, nrow=101, ncol=11)
b= matrix(m, nrow=101, ncol=90)
neg.nested.10.NC= cbind(b,a)
diag(neg.nested.10.NC)=m
diag(neg.nested.10.NC)[c(101,101)] = p
neg.nested.10.NC[101,90:100] = m
neg.nested.10.NC[1:100,101] = m
###image(neg.nested.10.NC,zlim=c(n,p))
Neg.Nest.10=neg.nested.10.NC


##############################
##  Neg.Nest.20
####### every species grows well under 20 species soil MC. 
a= matrix(n, nrow=101, ncol=21)
b= matrix(m, nrow=101, ncol=80)
neg.nested.20.NC= cbind(b,a)
diag(neg.nested.20.NC)=m 
diag(neg.nested.20.NC)[c(101,101)] = p
neg.nested.20.NC[101,80:100] = m
neg.nested.20.NC[1:100,101] = m
###image(neg.nested.20.NC,zlim=c(n,p))
Neg.Nest.20=neg.nested.20.NC

##############################
##  Neg.Nest.5
####### every species grows well under 5 species soil MC. 
a= matrix(n, nrow=101, ncol=6)
b= matrix(m, nrow=101, ncol=95)
neg.nested.5.NC= cbind(b,a)
diag(neg.nested.5.NC)=m 
diag(neg.nested.5.NC)[c(101,101)] = p
neg.nested.5.NC[101,80:100] = m
neg.nested.5.NC[1:100,101] = m
###image(neg.nested.5.NC,zlim=c(n,p))
Neg.Nest.5=neg.nested.5.NC


########################################################################################
########################################################################################
########################################################################################
########################################################################################

hyplist=list()
hyplist[[1]]=Null 
hyplist[[2]]=Null.NC
hyplist[[3]]=Null.PC
hyplist[[4]]=Pos.Ring.10.NC
hyplist[[5]]=Pos.Ring.5.NC    
hyplist[[6]]=Pos.Ring.20.NC   
hyplist[[7]]=Neg.Ring.10.NC
hyplist[[8]]=Neg.Ring.5.NC
hyplist[[9]]=Neg.Ring.20.NC
hyplist[[10]]=Pos.Mod.10.NC
hyplist[[11]]=Pos.Mod.5.NC
hyplist[[12]]=Pos.Mod.20.NC
hyplist[[13]]=Neg.Mod.10.NC
hyplist[[14]]=Neg.Mod.5.NC
hyplist[[15]]=Neg.Mod.20.NC
hyplist[[16]]=Pos.Nest.10.NC
hyplist[[17]]=Pos.Nest.20.NC 
hyplist[[18]]=Pos.Nest.5.NC 
hyplist[[19]]=Neg.Nest.10.NC  
hyplist[[20]]=Neg.Nest.20.NC
hyplist[[21]]=Neg.Nest.5.NC
hyplist[[22]]=Pos.Ring.10
hyplist[[23]]=Pos.Ring.5
hyplist[[24]]=Pos.Ring.20
hyplist[[25]]=Neg.Ring.10
hyplist[[26]]=Neg.Ring.5
hyplist[[27]]=Neg.Ring.20   
hyplist[[28]]=Pos.Mod.10
hyplist[[29]]=Pos.Mod.5
hyplist[[30]]=Pos.Mod.20
hyplist[[31]]=Neg.Mod.10
hyplist[[32]]=Neg.Mod.5
hyplist[[33]]=Neg.Mod.20
hyplist[[34]]=Pos.Nest.10
hyplist[[35]]=Pos.Nest.20
hyplist[[36]]=Pos.Nest.5
hyplist[[37]]=Neg.Nest.10
hyplist[[38]]=Neg.Nest.20
hyplist[[39]]=Neg.Nest.5

hypnames = c( 
"Null", 
"Null.NC",
"Null.PC",
"Pos.Ring.10.NC",
"Pos.Ring.5.NC",    
"Pos.Ring.20.NC",   
"Neg.Ring.10.NC",
"Neg.Ring.5.NC",
"Neg.Ring.20.NC",
"Pos.Mod.10.NC",
"Pos.Mod.5.NC",
"Pos.Mod.20.NC",
"Neg.Mod.10.NC",
"Neg.Mod.5.NC",
"Neg.Mod.20.NC",
"Pos.Nest.10.NC",
"Pos.Nest.20.NC",
"Pos.Nest.5.NC", 
"Neg.Nest.10.NC",  
"Neg.Nest.20.NC",
"Neg.Nest.5.NC",
"Pos.Ring.10",
"Pos.Ring.5",
"Pos.Ring.20",
"Neg.Ring.10",
"Neg.Ring.5",
"Neg.Ring.20",  
"Pos.Mod.10",
"Pos.Mod.5",
"Pos.Mod.20",
"Neg.Mod.10",
"Neg.Mod.5",
"Neg.Mod.20",
"Pos.Nest.10",
"Pos.Nest.20",
"Pos.Nest.5",
"Neg.Nest.10",
"Neg.Nest.20",
"Neg.Nest.5")






########################################################################################
########################################################################################
########################################################################################
########################################################################################

### Copyright 2019 Elizabeth (Libby) J. Trevenen & Dr Michael Renton. Adapted from model described in Teste et al 2017.
### If you wish to modify or share this code, you must obtain permission from the original authors. 




