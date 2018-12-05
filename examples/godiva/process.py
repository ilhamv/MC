import matplotlib.pyplot as plt

with open("k.txt") as f:
    content = f.readlines()
k = [float(x) for x in content]
print(k)

plt.plot(k)
plt.xlabel("iteration #")
plt.ylabel("k")
plt.savefig("k.png",dpi=1000)
