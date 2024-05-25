import math, sys, numpy as np, pandas as pd
import symnmfmodule as sm

np.random.seed(0)

class Point:
    def __init__(self, c:list[float], d:int, clust: int):
        self.coords = c
        self.dim = d
        self.cluster = clust

DEF_MAX_ITER = 300
def main():
    if(len(sys.argv) != 3):
        print("An error has occurred")
        exit(0)

    try:
        k = int(sys.argv[1])
        fileName = str(sys.argv[2])
    except:
        print("An error has occurred")
        exit(0)

    N,d, data = parseFile(fileName)

    # Run SymNMF and find clustering for it
    W = sm.norm(data)
    Hinit = initH(W, N,k)
    
    final_mat = sm.symnmf(Hinit, W)
    
    for i in range(N):
        data[i].cluster = np.array(final_mat[i]).argmax()
    
    SYMNMF_clustering = calcClustering(data, k)
    
    # Reset clustering
    for i in range(N):
        data[i].cluster = -1
    
    # Run K-means and find clustering for it
    final_cents = KMeans(k,N,d,DEF_MAX_ITER, data)
    for i in range(N):
        data[i].cluster = FindClosestCentroid(data[i], final_cents, d)

    KMEANS_clustering = calcClustering(data, k)
    
    print("nmf:", calcSilCoef(SYMNMF_clustering, N,k))
    print("kmeans:", calcSilCoef(KMEANS_clustering, N,k))

# Recieves clustered data and calculates the silhouette coefficient
def calcSilCoef(clustering:list[list[Point]], N:int,k:int) -> float:
    sill_score = 0

    for i in range(k):
        for j in range(len(clustering[i])):
            a = calcMeanDist(clustering[i][j], clustering[i])
            b = calcMinOuter(clustering[i][j], i, clustering)
            sill_score += (b-a)/(max(a,b))

    sill_score /= float(N)

    return sill_score


def calcMinOuter(x:Point, cluster:int, clustering:list[list[Point]]):
    b = 0
    bmin = 9999999 # Arbitrary "INFINITE" value

    for i in range(len(clustering)):
        if(i != cluster):
            b = calcMeanDist(x, clustering[i])

            if(b < bmin):
                bmin = b
    return bmin

# Calculates the mean distance from a point within a cluster
def calcMeanDist(x:Point, cluster:list[Point]):
    meanDist:float = 0;
    
    for i in range(len(cluster)):
        meanDist += dist(x, cluster[i], x.dim)
    meanDist /= float(len(cluster))

    return meanDist
    

# Returns a 2-dimensional array where each row represents a cluster
def calcClustering(data: list[Point], k) -> list[list[Point]]:
    N = len(data)

    clustering = []
    for i in range(k):
        row = []
        for j in range(N):
            if(data[j].cluster == i):
                temp = Point(data[j].coords.copy(), data[j].dim, i)
                row.append(temp)
        clustering.append(row)
    
    return clustering



# Initialize H randomly as said in the assignment
def initH(W:list[list[float]], n:int, k:int) -> list[list[float]]:
    avg:float = 0
    for i in range(n):
        for j in range(n):
            avg += W[i][j]
    avg /= math.pow(n,2) #The average of all entries of W

    bound = 2.0*math.sqrt(avg/float(k))

    # Initialize a zero nxk matrix
    H = []
    for i in range(n):
        row = []
        for j in range(k):
            row.append(0)
        H.append(row)
    
    for i in range(n):
        for j in range(k):
            H[i][j] = np.random.uniform(0, bound)
    
    return H

# Parses a csv file and returns N: number of observations, d: dimension, data: Point array of data
def parseFile(filename: str) -> tuple[int,int,list[Point]]:
    df = pd.read_csv(filename, header=None)
    N,d = df.shape

    # Convert pandas df to python matrix
    df = list(df.to_numpy())
    for i in df:
        i = list(i)

    # Initialize Point array
    data = [Point([0]*d,d,-1)]*N
    for i in range(len(df)):
        data[i] = Point(list(map(float,df[i])), d, -1)

    
    return tuple([N,d,data])

def printMat(mat :list[list[any]]):
    n = len(mat)
    m = len(mat[0])

    for i in range(n):
        for j in range(m):
            print("%.4f" % (mat[i][j]) ,end='')
            if(j != m-1):
                print(',',end='')
        print()


#Main logic of kmeans algorithm (Returns new centroids)
def KMeans(K, N, d, iter, data: list[Point]):
    Epsilon = 0.0001 #Convergence threshold
    iteration = 0 #Counts iterations
    DeltaVector = [0]*K #Saves distance difference between previous and current centroids

    #Init centroids to be first K dataPoints
    centroids = data[: K].copy()
    
    while True:
        #Reset means and counters for next iteration
        PtCtr = [0]*K #Counts points in each cluster (each index represents a cluster)
        KMEANS = [Point([0]*d,d,-1)]*K #The K Means (each index represents a cluster)

        #Assign each point to cluster
        for i in range(N):
            data[i].cluster = FindClosestCentroid(data[i], centroids, d)
            PtCtr[data[i].cluster] += 1 #Point assigned to cluster - increase count
            KMEANS[data[i].cluster] = ADD(data[i], KMEANS[data[i].cluster], d) #Add point to according mean's sum

        #Normalize mean for each cluster and compute difference vector
        for i in range(K):
            KMEANS[i] = MULT(KMEANS[i], 1.0/float(PtCtr[i]), d)
            DeltaVector[i] = dist(centroids[i], KMEANS[i], d)

        #Update the centroids
        centroids = KMEANS.copy()

        #Check if converged
        EXIT_FLAG = True
        for i in range(K):
            if(DeltaVector[i] >= Epsilon):
                EXIT_FLAG = False

        iteration += 1
        if(iteration == iter or EXIT_FLAG == True):
            break

    return centroids

#Finds and return closest centroid to given point x (denoted by cluster number)
def FindClosestCentroid(x, centroids, dim):
    assigned = 0
    minDist = dist(x, centroids[0], dim)
    for i in range(len(centroids)):
        curDist = dist(x, centroids[i], dim)
        if(curDist < minDist):
            minDist = curDist
            assigned = i
    
    return assigned

#Returns distance between two points
def dist(x,y, dim):
    dist = 0
    for i in range(dim):
        dist += pow(x.coords[i] - y.coords[i], 2)

    dist = math.sqrt(dist)
    return dist

#Returns a new point that is the sum of x and y
def ADD(x,y,dim):
    res = Point([0]*dim, dim, -1)
    for i in range(dim):
        res.coords[i] = x.coords[i] + y.coords[i]

    return res

#Returns a new point that is x multiplied by the scalar z
def MULT(x,z, dim):
    res = Point([0]*dim, dim, -1)
    for i in range(dim):
        res.coords[i] = x.coords[i] * z

    return res


if(__name__ =="__main__"):
    main()