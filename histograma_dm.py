import math

datos = []
with open('data/modulos_dm.dat', 'r') as f:
    for line in f:
        datos.append(float(line.strip()))

T = 1.5
n = len(datos)
v_max = max(datos) * 1.2
nbins = 20
bin_width = v_max / nbins

bins = []
counts = []

for i in range(nbins):
    bins.append((i + 0.5) * bin_width)
    counts.append(0)

for v in datos:
    ibin = int(v / bin_width)
    if ibin >= nbins:
        ibin = nbins - 1
    counts[ibin] += 1

for i in range(nbins):
    counts[i] = counts[i] / (n * bin_width)

with open('data/histograma_dm.dat', 'w') as f:
    for i in range(nbins):
        f.write(f"{bins[i]:.4f} {counts[i]:.6f}\n")

with open('data/maxwell_teorica.dat', 'w') as f:
    for i in range(101):
        v = i * v_max / 100
        mb = 4 * math.pi * (1 / (2 * math.pi * T))**1.5 * v**2 * math.exp(-v**2 / (2 * T))
        f.write(f"{v:.4f} {mb:.6f}\n")

print(f"N={n} particulas")
print(f"Media: {sum(datos)/n:.4f}")
print(f"Teorico: {math.sqrt(8*T/math.pi):.4f}")
