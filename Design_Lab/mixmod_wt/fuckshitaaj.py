modctr = 0

def __modularity(commu, status, graph):
    global modctr
    modctr += 1
    #print("modularity called", modctr, "edgewt: ", [graph[1][nbr].get('weight',1) for nbr in graph[1]])
    #print("From modularity, node_c: ", status.node_c)
    #print("From modularity, node_l: ", status.node_l)

    layer=status.layer
    node_l=status.node_l
    node_c=status.node_c       
    top=status.top
    bot=status.bot
    edge_l=status.edge_l
    edge_c=status.edge_c
    couple=status.couple
    mu = status.mu
    
    #Construct a dict = {node1: layer, node2: layer......}
    nodelayer = {}
    for l in layer:
        for node in layer[l]:
            nodelayer[node] = l

    #calculate Intra_inter ----------------------------------------------------------
    f=0    
    intra_inter={}
    for c in commu:
        intra_inter[c]=set()
        
        for n in commu[c]:
            for l in layer:
                if n in layer[l]:
                    intra_inter[c].add(l)

    #--------------------------------------------------------------------------------

    #calculate |E1| , |E2| , |E12|---------------------------------------------------
    E={}
    E12=0
    for l in layer:
        E[l]=0
        for n in layer[l]:
            for nei in node_l.get(n,set()):
                E[l]+= graph[n][nei].get('weight',1)
        E[l] = E[l]/2

    for n in node_c:
        for nei in node_c[n]:
            E12+=graph[n][nei].get('weight',1)
    E12 = E12/2
    #--------------------------------------------------------------------------------

   
    modularity=0    
    x1={}
    x2={}    
    
    for c in commu:
        x1[c]=0
        x2[c]=0
        modc_layer=0
        if len(intra_inter[c])>1:
            for l in layer:
                Aij=0.0
                hihj=0.0
                
                #compute summation Aij------------------------
                for n in commu[c]:
                    if n in layer[l]:
                        for nei in node_l.get(n,set()):
                            if nei in commu[c]:
                                Aij+=graph[n][nei].get('weight',1)
                #Aij = Aij/2

                #compute summation hihj--------------------------
                for n1 in commu[c]:
                    if n1 in layer[l]:
                        for n2 in commu[c]:
                            if n2 in layer[l]:
                                if(n1==n2):
                                   if(n1 not in node_l[n1]): 
                                       continue
                                hi = sum([graph[n1][nbr].get('weight',1) for nbr in node_l.get(n1,set())])
                                hj =sum([graph[n2][nbr].get('weight',1) for nbr in node_l.get(n2,set())])
                                hihj+= (hi*hj)
                #hihj = hihj/2
                #-------------------------------------------------
                try:
                    mod = (1.0/(2*E[l]))*(Aij - (hihj*1.0/(2*E[l])))
                except:
                    mod=0
                
                if f==1:# f is always 0
                    modc_layer+=mu*mod
                else:    
                    modc_layer+=mod
                    x1[c]+=mod
        #----------------------------------------------------------------------
            modc_couple=0
            Aij=0
            cicj=0
            
            #compute Aij------------------------------
            for n in commu[c]:
                if n in node_c:
                    for nei in node_c[n]:
                        if nei in commu[c]:
                            Aij+= graph[n][nei].get('weight',1)
            Aij = Aij/2

            #compute cicj-----------------------------
            tempmod=0.0
            for n1 in commu[c]:
                for n2 in commu[c]:
                    if(n1==n2 or nodelayer[n1]==nodelayer[n2] ) :
                        #If nodes in same layer, continue
                        continue
                    norm = 0;
                    note =0
                    if n1 in node_c:
                        ci = sum([graph[n1][nbr].get('weight',1) for nbr in node_c[n1]])

                        #extra penalty -> if a node is having coupling degree , then we penalize it for 
                        #            having link in its own layer with nodes that are not in the community
                        extra_penalty = sum(graph[n1][nbr].get('weight',1) for nbr in node_l[n1] \
                                        if nbr not in commu[c])

                        if extra_penalty>0:
                            lay = nodelayer[n1]
                            norm+=2*E[lay]
                            ci+=extra_penalty

                        norm+=E12
                        note =1
                    else:
                        ci = sum([graph[n1][nbr].get('weight',1) for nbr in node_l.get(n1,set())])
                        lay = nodelayer[n1]
                        norm+=2*E[lay]
                    if n2 in node_c:
                        cj = sum([graph[n2][nbr].get('weight',1) for nbr in node_c[n2]])
                        
                        #extra penalty -> if a node is having coupling degree , then we penalize it for 
                        #            having link in its own layer with nodes that are not in the community
                        extra_penalty = sum(graph[n2][nbr].get('weight',1) for nbr in node_l[n2] \
                                        if nbr not in commu[c])
                        if (extra_penalty>0):
                            lay = nodelayer[n2]
                            norm+=2*E[lay]
                            cj+=extra_penalty

                        if(note==0):
                            norm+=E12
                    else:
                        cj = sum([graph[n2][nbr].get('weight',1) for nbr in node_l.get(n2,set())])
                        lay = nodelayer[n2]
                        norm+= 2*E[lay]
                    if(n1 in node_c[n2]):
                        aij=1
                    else:
                        aij=0
                    try:
                        tempmod+= (aij-(ci*cj*1.0)/(norm*1.0))/(2*norm)
                        #cicj += (ci*cj*1.0)/(norm*norm*1.0)
                        #print("ci = ",ci," cj = ",cj," norm = ",norm)
                    except:
                        tempmod+=0
                        #cicj +=0
            '''
            cicj = cicj/2

            try:
                mod = (Aij*1.0/E12 - (cicj))
            except:
                mod=0
            '''
            mod = tempmod
            if f==1:
                modc_couple+=2*(1-mu)*mod
            else:
                modc_couple+=mod
                x2[c]+=mod
            #print modc_couple            
            ##print "ha hh"

            #-----------------------------------------------------------------------
            modularity+=modc_layer+modc_couple
            
        else:
            l=list(intra_inter[c])
            l=l[0];
            Aij=0.0
            hihj=0.0
            
            #compute summation AIJ------------------------
            for n in commu[c]:
                for nei in node_l.get(n,set()):
                    if nei in commu[c]:
                        Aij+=graph[n][nei].get('weight',1)
            #Aij = Aij/2

            #compute summation hihj--------------------------
            for n1 in commu[c]:
                for n2 in commu[c]:
                    if(n1==n2):
                       if(n1 not in node_l[n1]):   #implies not a loop
                           continue
                    hi = sum([graph[n1][nbr].get('weight',1) for nbr in node_l.get(n1,set())])
                    ci = sum([graph[n1][nbr].get('weight',1) for nbr in node_c.get(n1,set())])
                    hj =sum([graph[n2][nbr].get('weight',1) for nbr in node_l.get(n2,set())])
                    cj = sum([graph[n2][nbr].get('weight',1) for nbr in node_c.get(n2,set())])
                    hihj+= ((hi+ci)*(hj+cj))

            #hihj = hihj/2
            
            #-------------------------------------------------
            try:
                norm = 1.0/(2*E[l] + E12)
                mod = norm*(Aij*1.0 - (norm*hihj))
            except:
                mod=0
            
            if f==1: #f is always 0
                modc_layer+=mu*mod
            else:    
                modc_layer+=mod
                x1[c]+=mod

            modularity+=modc_layer        
                                    
    ##print x1,x2    
    return 0.333*modularity    

