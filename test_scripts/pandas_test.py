import pandas as pd

df1 = pd.DataFrame()
df2 = pd.DataFrame()

df_dict = {1:[df1], 2:[df2]}

df3 = pd.DataFrame.from_dict(df_dict)
print(df3)
