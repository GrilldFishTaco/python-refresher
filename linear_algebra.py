import numpy as np

# Problem 1
a = np.array([1, 2, 3])
b = np.array([4, 5, 6])

print(a + b)
print(b - a)

# Problem 2
A = np.array([[1, 2], [3, 4]])
B = np.array([[5, 6], [7, 8]])

print(A + B)
print(B - A)

# Problem 3
a = np.array([1, 2, 3])
b = np.array([4, 5, 6])

print(np.dot(a, b))

# Problem 4
A = np.array([[1, 2, 3], [4, 5, 6]])
B = np.array([[7, 8, 9, 10], [11, 12, 13, 14], [15, 16, 17, 18]])

print(np.dot(a, b))

# Problem 5
a = np.array([[1], [2], [3]])

print(np.linalg.norm(a))

# Problem 6
A = np.array([[1, 2], [3, 4]])

print(A.T)
