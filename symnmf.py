import math, sys, numpy as np, pandas as pd
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

    try:
        k = int(sys.argv[1])
        goal = str(sys.argv[2])
        fileName = str(sys.argv[3])
    except:
        print("ERR")
        exit(0)

    N,d, data = parseFile(fileName)

    if(goal == "sym"):
        output = sm.sym(data)
    elif(goal == "ddg"):
        output = sm.ddg(data)
    elif(goal == "norm"):
        output = sm.norm(data)
    elif(goal == "symnmf"):
        W = sm.norm(data)
        Hinit = initH(W, N,k)

        output = sm.symnmf(Hinit,W)
    else:
        print("ERR")
        exit(0)
    
    printMat(output)



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

        
if __name__ == '__main__':
    main()