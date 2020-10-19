from collections import defaultdict 
d=defaultdict()
   
# Initializing the value  

d[2]=1
d[3]=4
d[1]=2

d[7]=3
d[5]=3
d[4]=3
d[9]=3
d[8]=3


returnList=[i for i,j in sorted(d.items(), key = lambda kv:(kv[1], kv[0]))]
#returnList=[i for i,j in sorted(d.items(), key = lambda kv: kv[1])]
print(returnList)    

