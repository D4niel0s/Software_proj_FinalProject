import math, sys, numpy as np
np.random.seed(0)

class Point:
    def __init__(self, c:list[float], d:int):
        self.coords = c
        self.dim = d


def main():
    argc = len(sys.argv)
    if(argc != 4):
        print("ERR")
        exit(0)
    
    k = float(sys.argv[1])
    goal = str(sys.argv[2])
    fileName = str(sys.argv[3])

    



def initH(W, n,k):
    avg = 0
    for i in range(n):
        for j in range(n):
            avg += W[i][j]
    avg /= math.pow(n,2) #The average of all entries of W

    bound = 2.0*math.sqrt(avg/float(k))
    H = [(0 for j in range(n)) for i in range(n)]
    
    for i in range(n):
        for j in range(n):
            H[i][j] = np.random.uniform(0, bound)

    return H


if __name__ == '__main__':
    main()