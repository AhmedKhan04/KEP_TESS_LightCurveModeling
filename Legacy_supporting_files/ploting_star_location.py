from astroquery.vizier import Vizier
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

kic_values = np.loadtxt(r'vizier_KIC_cross.txt', dtype=str)
print(kic_values)
kic_ids = kic_values

"""
kic_ids = [
    6231538, 6606229, 6963490, 7300184, 7900367, 7915515, 9851822, 10350769, 10355055, 10415087,
    11027806, 11769929, 12216817, 1849235, 3424493, 3648131, 3965879, 4072890, 8264404, 8330102,
    8453431, 8493159, 8585472, 8845312, 9050337, 4544967, 4569150, 4577647, 4940217, 5108514,
    5358323, 5900260, 9051991, 9306095, 9368220, 9368524, 9649801, 9775887, 9532219, 9594189,
    9897710, 9717148, 2581626, 2972514, 3119295, 3953144, 4036687, 4066203, 4243668, 4374279,
    4466691, 4547067, 4651526, 4995588, 5027750, 5284701, 5286485, 5353653, 5357882, 5534340,
    5707205, 5788165, 6304420, 6442207, 6444630, 6778487, 6836820, 6955650, 7048016, 7124161,
    7347529, 7521682, 7601767, 7617649, 7668283, 7750215, 7905603, 7937097, 7948091, 7984934,
    8052082, 8087649, 8090059, 8144212, 8150307, 8245366, 8248967, 8249829, 8315263, 8393922,
    8516900, 8648251, 8649814, 8960514, 8963394, 9075949, 9077483, 9137819, 9202969, 9214444,
    9364179, 9594857, 9614153, 9700322, 9706609, 9724292, 9942562, 10451090, 11143576, 11704101, 11852985
]
"""

Vizier.ROW_LIMIT = -1  
catalog = "V/133/kic"  

result = Vizier.query_constraints(catalog=catalog, KICID=",".join(map(str, kic_ids)))


table = result[0].to_pandas()


df = table[['KIC', 'RAJ2000', 'DEJ2000']]
df.columns = ['KIC', 'RA', 'Dec']

# Plot sky positions - Cartesian
plt.figure(figsize=(12, 7))
plt.scatter(df['RA'], df['Dec'], s=15, label='Analyzed Stars', alpha=0.7)
plt.xlabel('Right Ascension (degrees)')
plt.ylabel('Declination (degrees)')
plt.title('Sky Position of Analyzed Stars')
plt.legend()
plt.grid(True)
plt.gca().invert_xaxis()
plt.show()

# Plot sky positions - Aitoff projection
plt.figure(figsize=(12, 7))
plt.subplot(111, projection="aitoff")
plt.scatter(np.radians(df['RA'] - 180), np.radians(df['Dec']), s=15)
plt.title('Sky Map in Aitoff Projection')
plt.grid(True)
plt.show()
