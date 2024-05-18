import math, sys, numpy as np, pandas as pd,time
import symnmfmodule as sm

np.random.seed(0)

class Point:
    def __init__(self, c:list[float], d:int, clust: int):
        self.coords = c
        self.dim = d
        self.cluster = clust


def main():
    argc = len(sys.argv)
    if(argc != 4):
        print("ERR")
        exit(0)

    k = int(sys.argv[1])
    goal = str(sys.argv[2])
    fileName = str(sys.argv[3])

    N, d, data = parseFile(fileName)
    
    W = sm.norm(data)
    Hinitial = initH(W, N, k)

    OPT = sm.symnmf(Hinitial, W)



def initH(W:list[list[float]], n:int, k:int) -> list[list[float]]:
    avg:float = 0
    for i in range(n):
        for j in range(n):
            avg += W[i][j]
    avg /= math.pow(n,2) #The average of all entries of W

    bound = 2.0*math.sqrt(avg/float(k))
    H = [[0]*k]*n
    
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


        
if __name__ == '__main__':
    main()