import numpy as np
import matplotlib.pyplot as plt
from numpy import mean, absolute
from math import *
import operator
from collections import OrderedDict
from dijkstra import *
import copy

node_x = {1: 0, 2: 100, 3: 0, 4: 80, 5: 0, 6: 30, 7: 0}
node_y = {1: 100, 2: 100, 3: 80, 4: 80, 5: 30, 6: 30, 7: 0}
edge_st = {1: 7, 2: 7, 3: 5, 4: 5, 5: 6, 6: 3, 7: 3, 8: 4, 9: 1}
edge_end = {1: 5, 2: 6, 3: 3, 4: 6, 5: 4, 6: 1, 7: 4, 8: 2, 9: 2}


a_x = [2, -6, 8, -5, 20, 35, 43, 67, 90]
a_y = [10, 20, 40, 65, 78, 90, 87, 75, 95]

#a_x = [10,  18, 43, 35]
#a_y = [12,  20,  87, 90]

def length(x1, x2, y1, y2):
	d = sqrt(pow(x1-x2,2)+pow(y1-y2,2))
	return d

G = {1:{2:length(node_x[1],node_x[2],node_y[1],node_y[2]),3:length(node_x[1],node_x[3],node_y[1],node_y[3])}
,2:{1:length(node_x[1],node_x[2],node_y[1],node_y[2]),4:length(node_x[4],node_x[2],node_y[4],node_y[2])}
,3:{5:length(node_x[3],node_x[5],node_y[3],node_y[5]),4:length(node_x[3],node_x[4],node_y[3],node_y[4]),1:length(node_x[1],node_x[3],node_y[1],node_y[3])}
,4:{2:length(node_x[4],node_x[2],node_y[4],node_y[2]),3:length(node_x[3],node_x[4],node_y[3],node_y[4]),6:length(node_x[4],node_x[6],node_y[4],node_y[6])}
,5:{7:length(node_x[5],node_x[7],node_y[5],node_y[7]),6:length(node_x[5],node_x[6],node_y[5],node_y[6]),3:length(node_x[3],node_x[5],node_y[3],node_y[5])}
,6:{5:length(node_x[5],node_x[6],node_y[5],node_y[6]),4:length(node_x[4],node_x[6],node_y[4],node_y[6]),7:length(node_x[7],node_x[6],node_y[7],node_y[6])}
,7:{5:length(node_x[5],node_x[7],node_y[5],node_y[7]),6:length(node_x[7],node_x[6],node_y[7],node_y[6])}}



def mad(data, axis=None):
    return mean(absolute(data - mean(data, axis)), axis)


def Shortest_Distance(a_x, a_y):
	queue = {}
	q = []
	Q = {}
	d = 0
	min_R = []
	min_D = []
	b = []
	proj = {0:{},1:{}}
	for key in range(len(a_x)):
		#last = d
		queue = {}
		proj[0].update({key:[]})
		proj[1].update({key:[]})
		for i in edge_st:
			distance_a = (pow((node_x[edge_st[i]] - node_x[edge_end[i]]), 2) + pow((node_y[edge_st[i]] - node_y[edge_end[i]]), 2))**0.5
			distance_b = ((a_x[key] - node_x[edge_st[i]])** 2 + (a_y[key] - node_y[edge_st[i]])** 2)**0.5
			distance_c = ((a_x[key] - node_x[edge_end[i]])*(a_x[key] - node_x[edge_end[i]]) + (a_y[key] - node_y[edge_end[i]])*(a_y[key] - node_y[edge_end[i]]))**0.5
			p0 = (distance_a + distance_b + distance_c) / 2
			if (distance_b + distance_c == distance_a):
				d = 0
				proj[0][key].append(a_x[key])
				proj[1][key].append(a_y[key])
			elif(pow(distance_c,2) >= pow(distance_a,2) + pow(distance_b, 2)):
				d = distance_b
				proj[0][key].append(node_x[edge_st[i]])
				proj[1][key].append(node_y[edge_st[i]])
			elif(pow(distance_b,2) >= pow(distance_a,2) + pow(distance_c, 2)):
				d = distance_c
				proj[0][key].append(node_x[edge_end[i]])
				proj[1][key].append(node_y[edge_end[i]])
			else:
				d = 2 * ((p0*(p0 - distance_a)*(p0 - distance_b)*(p0 - distance_c))**0.5)/distance_a
				tmp_x = (a_x[key]*pow((node_x[edge_st[i]] - node_x[edge_end[i]]), 2)+a_y[key]*(node_y[edge_st[i]] - node_y[edge_end[i]])*(node_x[edge_st[i]] - node_x[edge_end[i]])+(node_x[edge_end[i]]*node_y[edge_st[i]]-node_x[edge_st[i]]*node_y[edge_end[i]])*(node_y[edge_st[i]]-node_y[edge_end[i]]))/(pow((node_x[edge_st[i]] - node_x[edge_end[i]]), 2) + pow((node_y[edge_st[i]] - node_y[edge_end[i]]), 2))
				tmp_y = (a_x[key]*(node_x[edge_st[i]] - node_x[edge_end[i]])*(node_y[edge_st[i]]-node_y[edge_end[i]])+a_y[key]*((node_y[edge_st[i]] - node_y[edge_end[i]])**2)+(node_x[edge_st[i]]*node_y[edge_end[i]]-node_x[edge_end[i]]*node_y[edge_st[i]])*(node_x[edge_st[i]]-node_x[edge_end[i]]))/(pow((node_x[edge_st[i]] - node_x[edge_end[i]]), 2) + pow((node_y[edge_st[i]] - node_y[edge_end[i]]), 2))
				proj[0][key].append(tmp_x)
				proj[1][key].append(tmp_y)
			queue.update({i:d})
		Q.update({key:queue})
		sorted_queue = sorted(queue.items(), key = operator.itemgetter(1))
		q.append(sorted_queue[0])
		#print q[key][0]
        #(min_road,min_d) = sorted_queue.pop(0)# k->road segment num; v->the min distance
        #(min_road, min_d) = sorted_queue[0]
		b_x = (a_x[key]*pow((node_x[edge_st[q[key][0]]] - node_x[edge_end[q[key][0]]]), 2)+a_y[key]*(node_y[edge_st[q[key][0]]] - node_y[edge_end[q[key][0]]])*(node_x[edge_st[q[key][0]]] - node_x[edge_end[q[key][0]]])+(node_x[edge_end[q[key][0]]]*node_y[edge_st[q[key][0]]]-node_x[edge_st[q[key][0]]]*node_y[edge_end[q[key][0]]])*(node_y[edge_st[q[key][0]]]-node_y[edge_end[q[key][0]]]))/(pow((node_x[edge_st[q[key][0]]] - node_x[edge_end[q[key][0]]]), 2) + pow((node_y[edge_st[q[key][0]]] - node_y[edge_end[q[key][0]]]), 2))
		b_y = (a_x[key]*(node_x[edge_st[q[key][0]]] - node_x[edge_end[q[key][0]]])*(node_y[edge_st[q[key][0]]]-node_y[edge_end[q[key][0]]])+a_y[key]*((node_y[edge_st[q[key][0]]] - node_y[edge_end[q[key][0]]])**2)+(node_x[edge_st[q[key][0]]]*node_y[edge_end[q[key][0]]]-node_x[edge_end[q[key][0]]]*node_y[edge_st[q[key][0]]])*(node_x[edge_st[q[key][0]]]-node_x[edge_end[q[key][0]]]))/(pow((node_x[edge_st[q[key][0]]] - node_x[edge_end[q[key][0]]]), 2) + pow((node_y[edge_st[q[key][0]]] - node_y[edge_end[q[key][0]]]), 2))
		b.append((b_x, b_y))
	# b: the projection coordinate on ideal road segment
	# q: # of road segment and shortest distance
	# Q: every distance to road segment
	# proj: every coordinate on each road segment
	return b,q,Q,proj
#print Shortest_Distance(a_x, a_y)[2]
'''
def cal_sigma():
	l = []
	for i in range (8):
		a_xx = a_x[i]
		a_yy = a_y[i]
		(k,v,b_xx, b_yy,queue) = Shortest_Distance(a_xx, a_yy)
		l.append(v)
	sigma = mad(l)
	return sigma
'''
def cal_Emission_Pr():
	Pr_E = {}
	dic = {}
	sigma = 4.07
	b, q, Q, proj = Shortest_Distance(a_x, a_y)
	for t in range (len(a_x)):
		Pr_E.update({t+1:{}})#t->time, i->road seg, p->emission probability
		#print queue
		for i in edge_st.keys():
			p = (1/(sqrt(2*pi)*sigma))*exp(-0.5*pow(Q[t][i]/sigma, 2))
			Pr_E[t+1].update({i:log(p)})
		#Pr_E[t+1].update({t+1:dic})#t->time, i->road seg, p->emission probability
	return Pr_E
#print cal_Emission_Pr()
#===================================================================================
#           				 Shortest Path
#===================================================================================
def Dijkstra(G,start,end=None):
	D = {}	# dictionary of final distances
	P = {}	# dictionary of predecessors
	Q = priorityDictionary()	# estimated distances of non-final vertices
	Q[start] = 0
	
	for v in Q:
		D[v] = Q[v]
		if v == end: break
		
		for w in G[v]:
			vwLength = D[v] + G[v][w]
			if w in D:
				if vwLength < D[w]:
					raise ValueError, "Dijkstra: found better path to already-final vertex"
			elif w not in Q or vwLength < Q[w]:
				Q[w] = vwLength
				P[w] = v
	
	return (D,P)
			
def shortestPath(G,start,end):
	"""
	Find a single shortest path from the given start vertex to the given end vertex.
	The input has the same conventions as Dijkstra().
	The output is a list of the vertices in order along the shortest path.
	"""

	D,P = Dijkstra(G,start,end)
	Path = []
	l = 0
	while 1:
		Path.append(end)
		if end == start: break
		end = P[end]
	Path.reverse()
	for i in range(len(Path)-1):
		l += G[Path[i]][Path[i+1]]
	return Path,l

#print Dijkstra(G_0,1)
#print shortestPath(G_0,1,5)

def cal_beta():
	l = []
	b, q, Q, proj = Shortest_Distance(a_x, a_y)
	d_x = [10, 20, 25, 35, 15, 8, 24, 13]
	#d_x = [40, 65, 8]
	for i in range(0,len(a_x)-1):
		d_z = sqrt(pow(a_x[i]-a_x[i+1],2)+pow(a_y[i]-a_y[i+1],2))
		'''
		G0 = copy.deepcopy(G)
		d_x = sqrt(pow(b[i][0]-b[i-1][0],2)+pow(b[i][1]-b[i-1][1],2))
		G0.update({'x_t':{edge_st[q[i][0]]:length(node_x[edge_st[q[i][0]]],b[i][0],node_y[edge_st[q[i][0]]],b[i][1]),edge_end[q[i][0]]:length(node_x[edge_end[q[i][0]]],b[i][0],node_y[edge_end[q[i][0]]],b[i][1])}})
		G0.update({'x_t+1':{edge_st[q[i+1][0]]:length(node_x[edge_st[q[i+1][0]]],b[i+1][0],node_y[edge_st[q[i+1][0]]],b[i+1][1]),edge_end[q[i+1][0]]:length(node_x[edge_end[q[i+1][0]]],b[i+1][0],node_y[edge_end[q[i+1][0]]],b[i+1][1])}})
		G0[edge_st[q[i][0]]].update({'x_t':length(node_x[edge_st[q[i][0]]],b[i][0],node_y[edge_st[q[i][0]]],b[i][1])})
		G0[edge_end[q[i][0]]].update({'x_t':length(node_x[edge_end[q[i][0]]],b[i][0],node_y[edge_end[q[i][0]]],b[i][1])})
		del G0[edge_st[q[i][0]]][edge_end[q[i][0]]]
		del G0[edge_end[q[i][0]]][edge_st[q[i][0]]]
		G0[edge_st[q[i+1][0]]].update({'x_t+1':length(node_x[edge_st[q[i+1][0]]],b[i+1][0],node_y[edge_st[q[i+1][0]]],b[i+1][1])})
		G0[edge_end[q[i+1][0]]].update({'x_t+1':length(node_x[edge_end[q[i+1][0]]],b[i+1][0],node_y[edge_end[q[i+1][0]]],b[i+1][1])})
		if q[i][0] != q[i+1][0]:
			del G0[edge_st[q[i+1][0]]][edge_end[q[i+1][0]]] 
			del G0[edge_end[q[i+1][0]]][edge_st[q[i+1][0]]] 
		else:
			G0['x_t'].update({'x_t+1':length(b[i][0], b[i+1][0], b[i][1], b[i+1][1])})
			G0['x_t+1'].update({'x_t':length(b[i][0], b[i+1][0], b[i][1], b[i+1][1])})
		lo = shortestPath(G0, 'x_t', 'x_t+1')
		#print lo
		'''
		delta = absolute(d_z-d_x[i])
		l.append(delta)
	beta = (1/log(2,e))*mad(l)
	return beta
#print cal_beta()

def cal_Transition_Pr():
	Pr_T = {}
	dic1 = {}
	dic2 = {}
	d_z = []
	b_x = []
	b_y = []
	b_x1 = []
	b_y1 = []
	beta = cal_beta()
	proj = Shortest_Distance(a_x, a_y)[3]
	G0 = copy.deepcopy(G)
	#b,q,Q = Shortest_Distance(a_x, a_y)
	for i in range(len(a_x)-1):
		#d_z = sqrt(pow(a_x[i]-a_x[i+1],2)+pow(a_y[i]-a_y[i+1],2))
		#d_x = sqrt(pow(b[i][0]-b[i+1][0],2)+pow(b[i][1]-b[i+1][1],2))
		#d_t = abs(d_z - d_x)
		#Pr_T.append((i,(1/beta)*exp(-d_t/beta)))
		d_z.append(sqrt(pow(a_x[i]-a_x[i+1],2)+pow(a_y[i]-a_y[i+1],2)))
	i = 0
	while(i<len(a_x)-1):	
		b_x = proj[0][i]
		b_y = proj[1][i]
		b_x1= proj[0][i+1]
		b_y1= proj[1][i+1]
		for m in range(len(edge_st)):
			for n in range (len(edge_st)):
				G0.update({'x_t':{edge_st[m+1]:length(node_x[edge_st[m+1]],b_x[m],node_y[edge_st[m+1]],b_y[m]),edge_end[m+1]:length(node_x[edge_end[m+1]],b_x[m],node_y[edge_end[m+1]],b_y[m])}})
				G0.update({'x_t+1':{edge_st[n+1]:length(node_x[edge_st[n+1]],b_x1[n],node_y[edge_st[n+1]],b_y1[n]),edge_end[n+1]:length(node_x[edge_end[n+1]],b_x1[n],node_y[edge_end[n+1]],b_y1[n])}})
				G0[edge_st[m+1]].update({'x_t':length(node_x[edge_st[m+1]],b_x[m],node_y[edge_st[m+1]],b_y[m])})
				G0[edge_end[m+1]].update({'x_t':length(node_x[edge_end[m+1]],b_x[m],node_y[edge_end[m+1]],b_y[m])})
				del G0[edge_st[m+1]][edge_end[m+1]]
				del G0[edge_end[m+1]][edge_st[m+1]]
				G0[edge_st[n+1]].update({'x_t+1':length(node_x[edge_st[n+1]],b_x1[n],node_y[edge_st[n+1]],b_y1[n])})
				G0[edge_end[n+1]].update({'x_t+1':length(node_x[edge_end[n+1]],b_x1[n],node_y[edge_end[n+1]],b_y1[n])})
				if n != m:
					del G0[edge_st[n+1]][edge_end[n+1]] 
					del G0[edge_end[n+1]][edge_st[n+1]] 
				d_x = shortestPath(G0, 'x_t', 'x_t+1')[1]
				#print shortestPath(G0, 'x_t', 'x_t+1')
				G0 = copy.deepcopy(G)
				#d_x = sqrt(pow(b_x[m]-b_x1[n],2)+pow(b_y[m]-b_y1[n],2))
				delta = absolute(d_z[i]-d_x)
				dic1.update({n+1: log((1/beta)*exp(-delta/beta))}) #i-time, m+1->the first road segment, n+1->the next road seg
			dic2.update({m+1:dic1})
			dic1 = {}
		Pr_T.update({i+1:dic2})
		dic2 = {}
		i+=1	
	return Pr_T
#print cal_Transition_Pr()


def find_shortest_path():
	p_T = cal_Transition_Pr()
	p_E_tmp = cal_Emission_Pr()
	p_E = copy.deepcopy(p_E_tmp)
	max_p = -exp(200)
	path_rec = {} 
	tmp = 0
	p_path = {}
	dic1 = []
	t = 1
	while(t<len(a_x)):
		path_rec.update({t:{}})
		for i in range(1,1+len(edge_st)):
			for j in range(1, 1+len(edge_st)):
				prob = p_E[t][j]+p_T[t][j][i]
				#print prob
				if prob >= max_p:
					max_p = prob
					tmp = j
					path_rec[t].update({i:j})	#i->end
				else:
					max_p = max_p
			p_E[t+1][i] = p_E[t+1][i]+max_p
			max_p = -exp(200)
		#print t+1,p_E[t+1]
		if t == len(a_x)-1:
			#print p_E[t+1]
			end = max(p_E[t+1].iteritems(), key=operator.itemgetter(1))[0]
		t+=1
	i = len(a_x)-1
	dic1.append(end)
	while i > 0:
		dic1.append(path_rec[i][dic1[-1]])
		i-=1
	dic1.reverse()
	return dic1
print find_shortest_path()


x = []
y = []
b_xx = []
b_yy = []
path = find_shortest_path()
p_x = []
p_y = []
new_x = []
new_y = []
b,q,new = Shortest_Distance(a_x, a_y)[0], Shortest_Distance(a_x, a_y)[1], Shortest_Distance(a_x, a_y)[3]
# plot road map
for i in range(1,len(edge_st)+1):
    x.append(node_x[edge_st[i]])
    x.append(node_x[edge_end[i]])
    y.append(node_y[edge_st[i]])
    y.append(node_y[edge_end[i]])
for i in range(len(b)):
    b_xx.append(b[i][0])
    b_yy.append(b[i][1])
# plot real path
for i in range(len(path)):
	p_x.append(node_x[edge_st[path[i]]])
	p_x.append(node_x[edge_end[path[i]]])
	p_y.append(node_y[edge_st[path[i]]])
	p_y.append(node_y[edge_end[path[i]]])
for i in range(len(a_x)):
	new_x.append(new[0][i][path[i]-1])
	new_y.append(new[1][i][path[i]-1])
# plot road map, trajectory and projection coordinate on road segments 
#plt.xlim([-20, 150])
#plt.ylim([-20, 120])
plt.plot(x, y,color = 'b', linewidth = 2, label = 'road map')
plt.plot(b_xx,b_yy, 'o', color = 'r', label = 'projection on roads')
plt.plot(a_x, a_y, 'o', color = 'k', label = 'trajectory')
plt.plot(p_x, p_y, color = 'g', linewidth = 3, label = 'Map_Matching')
plt.plot(new_x, new_y, '*', color = 'y', linewidth = 5, label = 'candidate')
plt.title('Least distance Map_Matching Algorithm')
plt.legend(bbox_to_anchor=(1,0.05), loc=4, borderaxespad=0.)
plt.show()

