import pandas as pd
import os
direct = '/'.join(os.getcwd().split('/')[:-1]) + "/src"
file = direct +  "/data.json"
df = pd.read_json(file)
df.to_csv(direct + "/data.txt", index = False)