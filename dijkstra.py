from __future__ import generators

class priorityDictionary(dict):
	def __init__(self):
		'''Initialize priorityDictionary by creating binary heap of pairs (value,key).
Note that changing or removing a dict entry will not remove the old pair from the heap
until it is found by smallest() or until the heap is rebuilt.'''
		self.__heap = []
		dict.__init__(self)

	def smallest(self):
		'''Find smallest item after removing deleted items from front of heap.'''
		if len(self) == 0:
			raise IndexError, "smallest of empty priorityDictionary"
		heap = self.__heap
		while heap[0][1] not in self or self[heap[0][1]] != heap[0][0]:
			lastItem = heap.pop()
			insertionPoint = 0
			while 1:
				smallChild = 2*insertionPoint+1
				if smallChild+1 < len(heap) and heap[smallChild] > heap[smallChild+1] :
					smallChild += 1
				if smallChild >= len(heap) or lastItem <= heap[smallChild]:
					heap[insertionPoint] = lastItem
					break
				heap[insertionPoint] = heap[smallChild]
				insertionPoint = smallChild
		return heap[0][1]
	
	def __iter__(self):
		'''Create destructive sorted iterator of priorityDictionary.'''
		def iterfn():
			while len(self) > 0:
				x = self.smallest()
				yield x
				del self[x]
		return iterfn()
	
	def __setitem__(self,key,val):
		'''Change value stored in dictionary and add corresponding pair to heap.
Rebuilds the heap if the number of deleted items gets large, to avoid memory leakage.'''
		dict.__setitem__(self,key,val)
		heap = self.__heap
		if len(heap) > 2 * len(self):
			self.__heap = [(v,k) for k,v in self.iteritems()]
			self.__heap.sort()  # builtin sort probably faster than O(n)-time heapify
		else:
			newPair = (val,key)
			insertionPoint = len(heap)
			heap.append(None)
			while insertionPoint > 0 and newPair < heap[(insertionPoint-1)//2]:
				heap[insertionPoint] = heap[(insertionPoint-1)//2]
				insertionPoint = (insertionPoint-1)//2
			heap[insertionPoint] = newPair
	
	def setdefault(self,key,val):
		'''Reimplement setdefault to pass through our customized __setitem__.'''
		if key not in self:
			self[key] = val
		return self[key]
'''

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

G = {'s': {'u':10, 'x':5},
	'u': {'v':1, 'x':2},
	'v': {'y':4},
	'x':{'u':3,'v':9,'y':2},
	'y':{'s':7,'v':6}}

print Dijkstra(G,'s')
print shortestPath(G,'s','v')
'''
