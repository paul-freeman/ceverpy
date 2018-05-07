"""Python implementation of the Dynamic Time Warping algorithm"""
import numpy as np

def dtw(list1, list2):
    """Basic dynamic time warping algorithm"""
    costs = np.empty((len(list1)+1, len(list2)+1), dtype='float64')
    costs[0, 0] = -(abs(list1[0] - list2[0]))
    costs[0, 1:] = np.inf
    costs[1:, 0] = np.inf
    path = np.empty((len(list1)+1, len(list2)+1), dtype='uint8')
    path[0, 0] = 0
    path[0, 1:] = 0
    path[1:, 0] = 0
    for i in range(1, costs.shape[0]):
        for j in range(1, costs.shape[1]):
            distance = abs(list1[i-1] - list2[j-1])
            if costs[i-1, j] < costs[i, j-1] and  costs[i-1, j] - costs[i-1, j-1] < distance:
                costs[i, j], path[i, j] = costs[i-1, j] + distance, 1
                continue
            if costs[i, j-1] - costs[i-1, j-1] < distance:
                costs[i, j], path[i, j] = costs[i, j-1] + distance, 2
                continue
            costs[i, j], path[i, j] = costs[i-1, j-1] + distance * 2, 3
    i = costs.shape[0] - 1
    j = costs.shape[1] - 1
    index1 = [i]
    index2 = [j]
    while i > 1 or j > 1:
        if path[i, j] == 1:
            i, j = i-1, j
        elif path[i, j] == 2:
            i, j = i, j-1
        else: # path[i, j] == 3:
            i, j = i-1, j-1
        index1.append(i)
        index2.append(j)
    index1.reverse()
    index2.reverse()
    return (float(costs[-1, -1]),
            np.asarray(index1, dtype=float),
            np.asarray(index2, dtype=float),
            costs[1:, 1:])
