# Python program to find root which gives minimum
# height to tree

# This class represents a undirected graph using
# adjacency list representation
from queue import Queue


class MHGraph:

    # Constructor of mhgraph, initialize adjacency list
    # and degree vector
    def __init__(self, bags=None):
        if bags is None:
            self.bags = []
        else:
            self.bags = list(bags)
        self.V = len(bags)
        self.adj = dict((i, []) for i in range(self.V))
        self.degree = list()
        for i in range(self.V):
            self.degree.append(0)

    def addVertex(self, bag):
        if bag not in self.bags:
            self.V += 1
            self.bags.append(bag)
            self.degree.append(0)
            self.adj.update({(self.V, [])})
        else:
            print("Bag already exists.")

    # addEdge method adds vertex to adjacency list and
    # increases degree by 1
    def addEdge(self, v, w):
        v_ind = self.bags.index(v)
        w_ind = self.bags.index(w)
        self.adj[v_ind].append(w_ind)  # Adds w to v's list
        self.adj[w_ind].append(v_ind)  # Adds v to w's list
        self.degree[v_ind] += 1  # increment degree of v by 1
        self.degree[w_ind] += 1  # increment degree of w by 1

    # Method to return roots which gives minimum height to tree
    def rootForMinimumHeight(self):
        q = Queue()

        # First enqueue all leaf nodes in queue
        for i in range(self.V):
            if self.degree[i] == 1:
                q.put(i)

        # loop until total vertex remains less than 2
        while (self.V > 2):
            p = q.qsize()
            self.V -= p
            for i in range(p):
                t = q.get()

                # for each neighbour, decrease its degree and
                # if it become leaf, insert into queue
                for j in self.adj[t]:
                    self.degree[j] -= 1
                    if self.degree[j] == 1:
                        q.put(j)

        #  Copying the result from queue to result vector
        res = list()
        while (q.qsize() > 0):
            res.append(q.get())

        return res

    def get_bag(self, i):
        return self.bags[i]

# # Driver code
# g = MHGraph(6)
# g.addEdge(0, 3)
# g.addEdge(1, 3)
# g.addEdge(2, 3)
# g.addEdge(4, 3)
# g.addEdge(5, 4)
#
#
# # Function call
# res = g.rootForMinimumHeight()
# for i in res:
#     print(i, end=" ")
