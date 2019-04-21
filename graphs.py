import numpy as np
import random

class Vertex:
    def __init__(self, vertex):
        self.name = vertex
        self.neighbors = []
        
    def add_neighbor(self, neighbor):
        if isinstance(neighbor, Vertex):
            if neighbor.name not in self.neighbors:
                self.neighbors.append(neighbor.name)
                neighbor.neighbors.append(self.name)
                self.neighbors = sorted(self.neighbors)
                neighbor.neighbors = sorted(neighbor.neighbors)
        else:
            return False
        
    def __repr__(self):
        return str(self.neighbors)

class Graph:
    def __init__(self,num):
        self.vertices = {}
        self.adjacency_matrix = np.zeros(shape=(num,num))
        self.man_adj_mat = np.zeros(shape=(num,num))
    
    def add_vertex(self, vertex):
        if isinstance(vertex, Vertex):
            self.vertices[vertex.name] = vertex.neighbors

            
    def add_vertices(self, vertices):
        for vertex in vertices:
            if isinstance(vertex, Vertex):
                self.vertices[vertex.name] = vertex.neighbors
            
    def add_edge(self, vertex_from, vertex_to, i, j):
        wt = random.randint(1,20)
    	self.adjacency_matrix[i-1,j-1] = self.adjacency_matrix[j-1,i-1] = wt
        if isinstance(vertex_from, Vertex) and isinstance(vertex_to, Vertex):
            vertex_from.add_neighbor(vertex_to)
            if isinstance(vertex_from, Vertex) and isinstance(vertex_to, Vertex):
                self.vertices[vertex_from.name] = vertex_from.neighbors
                self.vertices[vertex_to.name] = vertex_to.neighbors
                #self.adjacency_matrix[i-1,j-1] = self.adjacency_matrix[j-1,i-1] = random.randint(1,20)        
        return wt
    
    def adjacencyList(self):
        if len(self.vertices) >= 1:
                #return [str(key) + ":" + str(self.vertices[key]) for key in self.vertices.keys()]  
                return [ (self.vertices[key]) for key in self.vertices.keys()]  
        else:
            return dict()

    def add_edge_mst(self, vertex_from, vertex_to):
        if isinstance(vertex_from, Vertex) and isinstance(vertex_to, Vertex):
            vertex_from.add_neighbor(vertex_to)
            if isinstance(vertex_from, Vertex) and isinstance(vertex_to, Vertex):
                self.vertices[vertex_from.name] = vertex_from.neighbors
                self.vertices[vertex_to.name] = vertex_to.neighbors

