import json
import pandas as pd
df = pd.read_json('data.json')
print(df)

short = df[df['arm'] == 'short']
# print(short)
print(short.shape)

longarm = df[df['arm'] == 'long']
# print(longarm)
print(longarm.shape)
